"""Explicit overdamped mechanics, ported from the original MATLAB RHS.

Lengths/time are dimensionless in fitting; node drag is fixed to one. All
force scales are therefore relative to drag, not absolute material moduli.
"""
from dataclasses import dataclass
from typing import NamedTuple
import math

import torch


FAMILIES = ('elastic', 'maxwell', 'kelvin-voigt', 'zener', 'burgers')


def smoothstep(x):
    x = x.clamp(0, 1)
    return x.square() * (3 - 2 * x)


def natural_lengths(t, rest, rate, floor=.5, start=0., ramp=.1):
    """Exact MATLAB piecewise integrated smoothstep (including assignment order).

    Use duration >= ramp in fitting; the MATLAB overlapping-ramp case can be
    nonphysical. Inactive rates remain exactly zero without NaN gradients.
    """
    ramp = max(float(ramp), 1e-8)
    rate = torch.as_tensor(rate, dtype=rest.dtype, device=rest.device).expand_as(rest)
    start = torch.as_tensor(start, dtype=rest.dtype, device=rest.device)
    floor = torch.as_tensor(floor, dtype=rest.dtype, device=rest.device)
    elapsed = (torch.as_tensor(t, dtype=rest.dtype, device=rest.device) - start).clamp_min(0)
    duration = (1 - floor) / rate.clamp_min(1e-30)
    active = rate > 0
    u = (elapsed / ramp).clamp(0, 1)
    integral = u ** 3 - .5 * u ** 4
    effective = torch.where(active & (elapsed < ramp), ramp * integral, 0.)
    activation = torch.where(active & (elapsed < ramp), smoothstep(u), 0.)
    plateau = active & (elapsed >= ramp) & (elapsed <= duration)
    effective = torch.where(plateau, elapsed - .5 * ramp, effective)
    activation = torch.where(plateau, 1., activation)
    down = active & (elapsed > duration) & (elapsed < duration + ramp)
    u = ((elapsed - duration) / ramp).clamp(0, 1)
    effective = torch.where(down, duration - .5*ramp + ramp*(u - u**3 + .5*u**4), effective)
    activation = torch.where(down, 1 - smoothstep(u), activation)
    effective = torch.where(active & (elapsed >= duration + ramp), duration, effective)
    return rest * torch.maximum(1 - rate * effective, floor), -rest * rate * activation


def constitutive(family, extension, extension_rate, sigma, z, k1, k2, c_m, c_s):
    """Signed force and internal-state rates; sigma has force units, z length."""
    zero = torch.zeros_like(extension)
    if family == 'elastic':
        return k1 * extension, zero, zero
    if family == 'maxwell':
        return sigma, k1 * extension_rate - k1 / c_m * sigma, zero
    if family == 'kelvin-voigt':
        return k1 * extension + c_m * extension_rate, zero, zero
    if family == 'zener':
        return k1 * extension + sigma, k2 * extension_rate - k2 / c_m * sigma, zero
    if family == 'burgers':
        zdot = (sigma - k1 * z) / c_s
        return sigma, k2 * (extension_rate - zdot) - k2 / c_m * sigma, zdot
    raise ValueError(f'Unknown constitutive family: {family}')


class State(NamedTuple):
    x: torch.Tensor
    sigma: torch.Tensor
    z: torch.Tensor
    damage: torch.Tensor
    release: torch.Tensor


@dataclass(frozen=True)
class Boundary:
    drag: float = 1.
    anchor_k: float = 1.
    anchor_drag: float = 0.
    adhesion_k: float = .3
    stuck_drag: float = .05
    sliding_drag: float = .01
    load_threshold: float = .02
    release_tau: float = .15
    gate_damage: float = .5
    # Minimal model intentionally uses free 3D geometry and explicit reference
    # tethers. It does not prescribe a fitted/moving future surface or traction.


@dataclass(frozen=True)
class Dynamics:
    family: str = 'zener'
    floor: float = .5
    ramp: float = .2
    start: float = 0.
    damage_tau: float = .08
    damage_width: float = .01  # 0 = hard MATLAB threshold, >0 = smooth surrogate
    cascade: bool = False
    enable_time: float = 0.
    damage_enabled: bool = True
    stiffness_floor: float = 0.
    damping_floor: float = 1e-6


def initial_state(graph):
    e, n = len(graph.edges), len(graph.x0)
    return State(graph.x0, graph.x0.new_zeros(e), graph.x0.new_zeros(e),
                 graph.x0.new_zeros(e), graph.x0.new_zeros(n))


def edge_sum(graph, vectors):
    i, j = graph.edges.T
    return torch.zeros_like(graph.x0).index_add(0, i, vectors).index_add(0, j, -vectors)


def coupled_velocity(graph, direction, damping, drag, force):
    """Solve SPD Kelvin–Voigt drag balance; no lagged neighbor velocities.

    Dense solve for fixtures; matrix-free preconditioned CG for larger graphs.
    Iterations retain autograd. A failed residual check raises rather than
    silently accepting an inaccurate force balance.
    """
    i, j = graph.edges.T
    def apply(v):
        rate = ((v[i] - v[j]) * direction).sum(1)
        return drag[:, None] * v + edge_sum(graph, (damping * rate)[:, None] * direction)
    if len(graph.x0) <= 32:
        b = graph.x0.new_zeros((len(i), 3 * len(graph.x0)))
        rows = torch.arange(len(i))[:, None]
        b[rows, 3 * i[:, None] + torch.arange(3)] = direction
        b[rows, 3 * j[:, None] + torch.arange(3)] = -direction
        matrix = torch.diag(drag.repeat_interleave(3)) + b.T @ (damping[:, None] * b)
        return torch.linalg.solve(matrix, force.flatten()).reshape_as(force)
    diagonal = drag[:, None].expand_as(force).clone()
    local = damping[:, None] * direction.square()
    diagonal = diagonal.index_add(0, i, local).index_add(0, j, local)
    x = torch.zeros_like(force)
    r = force
    z = r / diagonal
    p = z
    rz = (r*z).sum()
    tolerance = 1e-10 * max(float(torch.linalg.vector_norm(force).detach()), 1e-10)
    for _ in range(200):
        if float(torch.linalg.vector_norm(r).detach()) <= tolerance:
            break
        ap = apply(p)
        alpha = rz / (p*ap).sum().clamp_min(1e-300)
        x = x + alpha*p
        r = r - alpha*ap
        z = r / diagonal
        new_rz = (r*z).sum()
        p = z + (new_rz / rz.clamp_min(1e-300))*p
        rz = new_rz
    if float(torch.linalg.vector_norm(apply(x)-force).detach()) > 10*tolerance:
        raise RuntimeError('Coupled drag solve did not converge')
    return x


def rhs(t, state, graph, parameters, dynamics=Dynamics(), boundary=Boundary()):
    x, sigma, z, damage, release = state
    i, j = graph.edges.T
    separation = x[j] - x[i]
    lengths = torch.linalg.vector_norm(separation, dim=1)
    direction = separation / lengths.clamp_min(1e-18)[:, None]
    natural, natural_dot = natural_lengths(t, graph.lengths, parameters['rate'] * graph.active,
                                          dynamics.floor, dynamics.start, dynamics.ramp)
    extension = lengths - natural
    integrity = 1 - damage.clamp(0, 1)
    g = integrity.clamp_min(dynamics.stiffness_floor)
    gd = integrity.clamp_min(dynamics.damping_floor)
    k1, k2 = g * parameters['k1'], g * parameters['k2']
    cm = (gd * parameters['cm']).clamp_min(1e-12)
    cs = (gd * parameters['cs']).clamp_min(1e-12)
    passive, _, _ = constitutive(dynamics.family, extension, -natural_dot, sigma, z, k1, k2, cm, cs)
    internal = edge_sum(graph, (integrity * passive)[:, None] * direction)
    a = release.clamp(0, 1)
    free = (~graph.anchor).to(x.dtype)
    displacement = x - graph.x0
    tether = (1-a) * boundary.adhesion_k * free + boundary.anchor_k * graph.anchor
    drag = (boundary.drag + boundary.anchor_drag * graph.anchor +
            ((1-a)*boundary.stuck_drag + a*boundary.sliding_drag)*free)
    force = internal - tether[:, None] * displacement
    velocity = (coupled_velocity(graph, direction, integrity*cm, drag, force)
                if dynamics.family == 'kelvin-voigt' else force / drag[:, None])
    extension_dot = ((velocity[j]-velocity[i])*direction).sum(1) - natural_dot
    passive, sigma_dot, z_dot = constitutive(dynamics.family, extension, extension_dot, sigma, z, k1, k2, cm, cs)
    # Hard validation mode matches MATLAB threshold/cascade RHS. Smooth mode
    # has a compact transition above threshold: exactly zero subthreshold damage.
    margin = extension / graph.lengths - parameters['threshold']
    if dynamics.damage_width > 0:
        activity = smoothstep(margin / dynamics.damage_width)
    else:
        activity = (margin >= 0).to(x.dtype)
    activity = activity * graph.weak
    if dynamics.cascade:
        if dynamics.damage_width == 0:
            trigger = ((activity.max() > 0) | ((damage * graph.weak).max() > 0)).to(x.dtype)
        else:
            trigger = torch.maximum(activity.max(), smoothstep((damage*graph.weak).max()/.01))
        activity = trigger * graph.weak
    if not dynamics.damage_enabled or float(t) < dynamics.enable_time:
        activity = activity * 0
    damage_dot = integrity * activity / dynamics.damage_tau
    selected_damage = (damage.clamp(0, 1) * graph.weak).max()
    gate = smoothstep((selected_damage-boundary.gate_damage)/(1-boundary.gate_damage))
    trial = -boundary.adhesion_k * displacement - boundary.stuck_drag * velocity
    load = torch.sqrt(trial.square().sum(1) + 1e-24) - 1e-12
    overload = .5 * (1 + torch.tanh((load-boundary.load_threshold) / (.25*boundary.load_threshold)))
    release_dot = gate * overload * (1-a) / boundary.release_tau * free
    return State(velocity, sigma_dot, z_dot, damage_dot, release_dot)


def rollout(graph, parameters, times, dynamics=Dynamics(), boundary=Boundary(), max_step=.02, state=None):
    """Differentiable explicit midpoint integration with fixed bounded substeps.

    Damage/release are monotone and bounded when h <= 0.25*tau. Hard threshold
    event timing is step-resolved, not MATLAB's event-localized ode15s solve.
    Internal stresses start relaxed unless an explicit state is supplied.
    """
    if dynamics.family not in FAMILIES or max_step <= 0 or boundary.drag <= 0:
        raise ValueError('Invalid family, step or drag')
    if dynamics.damage_tau <= 0 or boundary.release_tau <= 0 or not 0 <= boundary.gate_damage < 1:
        raise ValueError('Invalid damage/release configuration')
    if boundary.load_threshold <= 0 or not 0 < dynamics.floor <= 1 or dynamics.damage_width < 0:
        raise ValueError('Invalid contraction floor or damage thresholds')
    times = torch.as_tensor(times, dtype=graph.x0.dtype)
    if times.ndim != 1 or len(times) < 2 or not torch.isfinite(times).all() or (torch.diff(times) <= 0).any():
        raise ValueError('Need finite increasing times')
    state = initial_state(graph) if state is None else state
    history = [state]
    step_limit = min(max_step, .25*dynamics.damage_tau, .25*boundary.release_tau)
    for left, right in zip(times[:-1], times[1:]):
        steps = math.ceil(float(right-left) / step_limit)
        h = (right-left) / steps
        for substep in range(steps):
            t = left + substep*h
            first = rhs(t, state, graph, parameters, dynamics, boundary)
            midpoint = State(*(v + .5*h*dv for v, dv in zip(state, first)))
            second = rhs(t+.5*h, midpoint, graph, parameters, dynamics, boundary)
            state = State(*(v + h*dv for v, dv in zip(state, second)))
        if not all(torch.isfinite(v).all() for v in state):
            raise FloatingPointError('Unstable rollout: reduce max_step or mechanical parameter bounds')
        if (state.damage < -1e-12).any() or (state.damage > 1+1e-12).any():
            raise FloatingPointError('Damage left [0,1]; reduce step')
        history.append(state)
    return State(*(torch.stack([s[index] for s in history]) for index in range(5)))
