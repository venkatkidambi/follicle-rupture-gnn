"""Small known 3D network with local failure and reproducible observations."""
from dataclasses import replace
import numpy as np
import torch

from .data import make_graph
from .mechanics import Boundary, Dynamics, rollout
from .model import decode


def example(family='zener', noise=0., seed=0, max_step=.02):
    axis = np.linspace(0, 3, 9)
    x = np.column_stack((axis, .18*np.sin(2*axis), .15*np.cos(3*axis)))
    edges = np.array([(i, i+d) for i in range(len(x)) for d in (1, 2) if i+d < len(x)])
    graph = make_graph(x, edges, regions=3)
    graph.weak = graph.edges[:, 0] == 0
    graph.active = ~graph.weak
    graph.anchor = torch.arange(len(x)) == len(x)-1
    graph.features[:, -1] = graph.anchor
    values = torch.tensor([[.45, .5, .20], [.65, .7, .16], [.35, .4, .23]], dtype=torch.float64)
    parameters = decode(values, graph, threshold=.05)
    times = torch.linspace(0, 3., 25, dtype=torch.float64)
    dynamics = Dynamics(family=family, ramp=.15, damage_width=.01)
    boundary = Boundary(anchor_k=2., adhesion_k=.7, load_threshold=.008)
    with torch.no_grad():
        truth = rollout(graph, parameters, times, dynamics, boundary, max_step=max_step)
    observed = truth.x.numpy().copy()
    observed[1:] += np.random.default_rng(seed).normal(0, noise, observed[1:].shape)
    return graph, values, times, dynamics, boundary, truth, observed


def rupture_summary(damage, times, threshold=.5):
    damage, times = np.asarray(damage), np.asarray(times)
    crossed = damage >= threshold
    onset = [float(times[np.flatnonzero(row)[0]]) if row.any() else None for row in crossed.T]
    finite = [v for v in onset if v is not None]
    return dict(damage_definition=f'first saved sample D >= {threshold}',
                first_rupture_time=min(finite) if finite else None,
                edge_onset=onset, failed_count_per_frame=crossed.sum(1).tolist(),
                max_damage=float(damage.max()))
