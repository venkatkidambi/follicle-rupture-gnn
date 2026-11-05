"""Prefix-only inverse fitting and frozen-parameter free forecasts."""
from copy import deepcopy
import numpy as np
import torch

from .mechanics import rollout
from .model import RegionalModel, decode, regularization, NAMES


def trajectory_loss(predicted, observed):
    valid = torch.isfinite(observed).all(-1)
    if not valid.any():
        raise ValueError('No valid observations')
    residual = predicted[valid] - observed[valid]
    return residual.square().mean()


def fit(graph, observations, times, train_frames, dynamics, boundary,
        kind='gnn', epochs=80, seed=0, max_step=.02, learning_rate=.03, penalty=1e-5):
    if not 3 <= train_frames < len(times) or epochs < 1:
        raise ValueError('Need >=3 training frames, future holdout, and positive epochs')
    torch.manual_seed(seed)
    model = RegionalModel(graph, kind=kind, regions=int(graph.region.max())+1)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    # Slice before any objective, normalization or encoder operation.
    observed = torch.as_tensor(observations[:train_frames], dtype=graph.x0.dtype)
    train_time = torch.as_tensor(times[:train_frames], dtype=graph.x0.dtype)
    scale = trajectory_loss(graph.x0.expand_as(observed)[1:], observed[1:]).detach().clamp_min(1e-6)
    history, best, best_weights = [], float('inf'), None
    for epoch in range(epochs):
        optimizer.zero_grad()
        predicted = rollout(graph, model(graph), train_time, dynamics, boundary, max_step)
        loss = trajectory_loss(predicted.x[1:], observed[1:]) / scale
        objective = loss + penalty*regularization(model.regional_values(graph))
        if not torch.isfinite(objective):
            raise FloatingPointError('Nonfinite training objective')
        score = float(objective.detach())
        if score < best:
            best, best_weights = score, deepcopy(model.state_dict())
        objective.backward()
        if any(p.grad is not None and not torch.isfinite(p.grad).all() for p in model.parameters()):
            raise FloatingPointError('Nonfinite mechanical gradients')
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.)
        optimizer.step()
        history.append(float(loss.detach()))
    model.load_state_dict(best_weights)
    return model, dict(seed=seed, epochs=epochs, normalized_training_loss=best,
                      loss_history=history, selection='lowest training objective only')


def metrics(predicted, observed, origin, start, direction_floor=.1):
    """Frame-zero displacement errors; direction excludes tiny observed motion."""
    p = np.asarray(predicted)[start:] - origin
    y = np.asarray(observed)[start:] - origin
    valid = np.isfinite(y).all(-1) & np.isfinite(p).all(-1)
    if not valid.any():
        return dict(valid_samples=0)
    p, y = p[valid], y[valid]
    pn, yn = np.linalg.norm(p, axis=-1), np.linalg.norm(y, axis=-1)
    directional = (pn > direction_floor) & (yn > direction_floor)
    cosine = (p[directional]*y[directional]).sum(-1) / (pn[directional]*yn[directional])
    angles = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    return dict(valid_samples=int(valid.sum()),
                vector_rmse=float(np.sqrt(np.mean(np.sum((p-y)**2, axis=-1)))),
                magnitude_mae=float(np.mean(abs(pn-yn))),
                direction_deg=float(angles.mean()) if len(angles) else None,
                direction_samples=len(angles), direction_floor=direction_floor,
                stationary_predictions=int(np.sum(pn <= direction_floor)))


def motion_baselines(observations, times, train_frames):
    """Stationary and per-cell linear motion from the training prefix only."""
    x, times = np.asarray(observations), np.asarray(times)
    velocity = np.zeros_like(x[0])
    for node in range(x.shape[1]):
        valid = np.isfinite(x[:train_frames, node]).all(1)
        t = times[:train_frames][valid] - times[0]
        if (t*t).sum() > 0:
            velocity[node] = (t[:, None]*(x[:train_frames, node][valid]-x[0, node])).sum(0)/(t*t).sum()
    return dict(stationary=np.broadcast_to(x[0], x.shape),
                linear_prefix=x[0] + (times-times[0])[:, None, None]*velocity)


def sensitivity(graph, values, observations, times, train_frames, dynamics, boundary, max_step=.02):
    """Training-only finite differences in log regional parameters.

    Singular values expose nonidentifiability. Linearized standard errors are
    conditional on model/boundaries and iid residuals, not calibrated intervals.
    """
    values = values.detach()
    observed = np.asarray(observations[:train_frames])
    mask = np.isfinite(observed[1:]).all(-1)
    columns = []
    with torch.no_grad():
        base = rollout(graph, decode(values, graph), times[:train_frames], dynamics, boundary, max_step).x.numpy()
        for index in range(values.numel()):
            delta = torch.zeros_like(values).flatten()
            delta[index] = .01
            plus, minus = values*torch.exp(delta.reshape_as(values)), values*torch.exp(-delta.reshape_as(values))
            p = rollout(graph, decode(plus, graph), times[:train_frames], dynamics, boundary, max_step).x.numpy()
            m = rollout(graph, decode(minus, graph), times[:train_frames], dynamics, boundary, max_step).x.numpy()
            columns.append(((p[1:]-m[1:])[mask]/.02).ravel())
    jacobian = np.stack(columns, axis=1)
    _, singular, vt = np.linalg.svd(jacobian, full_matrices=False)
    cutoff = max(singular[0]*1e-5, 1e-10)
    rank = int(np.sum(singular > cutoff))
    residual = (base[1:]-observed[1:])[mask].ravel()
    variance = float(residual@residual/max(1, len(residual)-rank))
    se = None
    if rank == values.numel():
        covariance = (vt.T * (variance/singular**2)) @ vt
        se = np.sqrt(np.diag(covariance)).reshape(values.shape).tolist()
    return dict(parameter_order=[f'region_{r}/{name}' for r in range(len(values)) for name in NAMES],
                log_parameter_singular_values=singular.tolist(), effective_rank=rank,
                parameter_count=values.numel(), relative_rank_cutoff=1e-5,
                log_parameter_standard_errors=se,
                weakest_log_parameter_direction=vt[-1].tolist(),
                uncertainty='local iid-residual approximation; null directions unbounded; no topology/boundary uncertainty')
