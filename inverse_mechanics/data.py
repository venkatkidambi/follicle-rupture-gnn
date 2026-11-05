"""Initial-frame MATLAB geometry and calibrated persistent trajectories."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from scipy.io import loadmat


KNOWN_FILE = 'trackedNodes_nuc_48_2_OA7_1.mat'
SPACING_UM = np.array([0.4649160043116627, 0.4649160043116627, 1.0])


@dataclass
class Tracks:
    positions: np.ndarray  # T,N,3, micrometers; missing observations are NaN
    times: np.ndarray      # seconds, acquisition clock retained after slicing
    ids: np.ndarray        # persistent original 1-based MATLAB row identities
    metadata: dict


def load_tracks(path, frame_limit=35, timestamps=None, spacing_um=None, duration_s=None):
    """Select nodes only at frame zero, exactly as MATLAB does.

    Future track availability never changes model membership. Missing training
    or scoring observations are masked, not interpolated. Row identity is an
    acquisition assumption; there is no identity matching from future motion.
    """
    path = Path(path)
    raw = np.asarray(loadmat(path)['trackedNodes'], dtype=float)
    if raw.ndim != 3 or raw.shape[1] != 3 or raw.shape[2] < 3:
        raise ValueError('Expected trackedNodes N x 3 x T, T >= 3')
    known = path.name == KNOWN_FILE
    spacing = np.asarray(SPACING_UM if spacing_um is None and known else spacing_um, dtype=float)
    if spacing.shape != (3,) or not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError('Supply three positive finite xyz spacings in micrometers')
    time_source = 'explicit timestamps'
    if timestamps is None:
        if duration_s is None and known and raw.shape[2] == 45:
            duration_s = 1800.
            time_source = 'repository acquisition convention: 1800 s / 44 intervals; not MAT metadata'
        elif duration_s is None:
            raise ValueError('Unknown acquisition: supply timestamps or duration_s')
        else:
            time_source = 'explicit duration, assumed uniform acquisition'
        timestamps = np.linspace(0, duration_s, raw.shape[2])
    times = np.asarray(timestamps, dtype=float)
    if times.shape != (raw.shape[2],) or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError('Need strictly increasing finite timestamps for the full acquisition')
    if not isinstance(frame_limit, int) or isinstance(frame_limit, bool) or not 3 <= frame_limit <= len(times):
        raise ValueError('frame_limit must be an integer between 3 and the acquisition frame count')
    valid = np.isfinite(raw).all(1) & (raw != 0).any(1)
    ids = np.flatnonzero(valid[:, 0])
    # MATLAB unique(..., rows, stable) on the initial calibrated geometry.
    _, first = np.unique(raw[ids, :, 0], axis=0, return_index=True)
    ids = ids[np.sort(first)]
    if len(ids) < 4:
        raise ValueError('Need four distinct valid reference nodes')
    xyz = raw[ids, :, :frame_limit].transpose(2, 0, 1) * spacing
    xyz[~valid[ids, :frame_limit].T] = np.nan
    origin = np.median(xyz[0], axis=0)
    _, _, vt = np.linalg.svd(xyz[0] - origin, full_matrices=False)
    axes = vt.T
    for col in range(3):
        if axes[np.argmax(abs(axes[:, col])), col] < 0:
            axes[:, col] *= -1
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1
    axes[:, 0] *= -1  # experimentalLongAxisSign = -1 in MATLAB
    xyz = (xyz - origin) @ axes
    return Tracks(xyz, times[:frame_limit] - times[0], ids + 1, dict(
        source=str(path), spacing_um=spacing.tolist(), origin_um=origin.tolist(),
        initial_axes=axes.tolist(), time_source=time_source,
        acquisition_frames=len(times), frame_indices=list(range(frame_limit)),
        selection='valid unique frame-zero rows only; later missing observations masked',
        identity='persistent MATLAB row IDs assumed, not independently verified',
        registration='one fixed translation and orientation from frame zero; no future registration'))


def candidate_edges(x, neighbors=6):
    """MATLAB buildExperimentalEdgeList: unique sorted symmetric kNN pairs."""
    x = np.asarray(x)
    if x.ndim != 2 or x.shape[1] != 3 or len(x) < 2 or not np.isfinite(x).all():
        raise ValueError('Need finite N x 3 reference geometry')
    if not isinstance(neighbors, int) or neighbors < 1:
        raise ValueError('neighbors must be a positive integer')
    k = min(neighbors, len(x) - 1)
    d2 = ((x[:, None] - x[None]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    if (d2 == 0).any():
        raise ValueError('Duplicate reference nodes')
    order = np.argsort(d2, axis=1, kind='stable')[:, :k]
    return np.unique(np.sort(np.column_stack((np.repeat(np.arange(len(x)), k), order.ravel())), axis=1), axis=0)


@dataclass
class Graph:
    x0: torch.Tensor
    edges: torch.Tensor
    lengths: torch.Tensor
    region: torch.Tensor
    features: torch.Tensor
    anchor: torch.Tensor
    active: torch.Tensor
    weak: torch.Tensor


def make_graph(x, edges=None, regions=3):
    """Reference-only masks: MATLAB axial contraction, tip patch and P90 anchor.

    x is already in model coordinates. Three axial bins bound the number of
    learned mechanical values. No per-edge free parameters or ID features.
    """
    x = torch.as_tensor(x, dtype=torch.float64)
    if regions < 1 or regions > 4:
        raise ValueError('Use one to four regions')
    edges = torch.as_tensor(candidate_edges(x.numpy()) if edges is None else edges, dtype=torch.long)
    if edges.ndim != 2 or edges.shape[1] != 2 or len(edges) == 0:
        raise ValueError('Need E x 2 edges')
    if edges.min() < 0 or edges.max() >= len(x) or (edges[:, 0] == edges[:, 1]).any():
        raise ValueError('Invalid edge endpoints')
    if len(torch.unique(torch.sort(edges, dim=1).values, dim=0)) != len(edges):
        raise ValueError('Duplicate undirected edges')
    lengths = torch.linalg.vector_norm(x[edges[:, 1]] - x[edges[:, 0]], dim=1)
    if not torch.isfinite(x).all() or (lengths <= 1e-12).any():
        raise ValueError('Invalid or coincident geometry')
    # Already oriented x is the long axis for experimental data.
    score = x[:, 0]
    ap = (score - score.min()) / (score.max() - score.min()).clamp_min(1e-12)
    midpoint = ap[edges].mean(1)
    region = (midpoint * regions).long().clamp_max(regions - 1)
    anchor = score >= torch.quantile(score, .9)
    weak_candidates = (score[edges] <= torch.quantile(score, .1)).all(1)
    weak = torch.zeros(len(edges), dtype=torch.bool)
    eligible = torch.where(weak_candidates)[0]
    weak[eligible[torch.argsort(midpoint[eligible])[:22]]] = True
    active = ((score[edges].mean(1) <= torch.quantile(score, .9)) &
              ~anchor[edges].all(1) & ~weak &
              ((score[edges[:, 1]] - score[edges[:, 0]]).abs() / lengths >= .35))
    features = torch.cat((x - x.mean(0), ap[:, None], anchor[:, None].to(x.dtype)), dim=1)
    return Graph(x, edges, lengths, region, features, anchor, active, weak)
