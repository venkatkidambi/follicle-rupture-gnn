"""Shared message passing, symmetric edge features, few regional coefficients."""
import torch
from torch import nn


BOUNDS = ((.04, 1.2), (.12, 1.8), (.015, .35))  # k/drag, relaxation, contraction
NAMES = ('stiffness', 'relaxation', 'contraction')


def decode(values, graph, threshold=.05):
    """Fixed branch ratios remove unsupported constitutive degrees of freedom."""
    k, tau, rate = values[graph.region].unbind(1)
    return dict(k1=k, k2=k, cm=k*tau, cs=k*tau, rate=rate,
                threshold=torch.full_like(k, threshold))


class RegionalModel(nn.Module):
    """Three scalars per region, from a GNN or a direct fitting baseline.

    Topology, branch ratios, initial internal states, boundary mechanics and
    damage susceptibility are fixed. Geometry is the only encoder input.
    """
    def __init__(self, graph, kind='gnn', regions=3, hidden=12, threshold=.05):
        super().__init__()
        if kind not in ('gnn', 'regional', 'uniform'):
            raise ValueError('kind must be gnn, regional or uniform')
        if int(graph.region.max()) >= regions:
            raise ValueError('Graph region exceeds model region count')
        self.kind, self.regions, self.threshold = kind, regions, threshold
        self.register_buffer('lower', torch.tensor([v[0] for v in BOUNDS], dtype=torch.float64))
        self.register_buffer('upper', torch.tensor([v[1] for v in BOUNDS], dtype=torch.float64))
        initial = torch.tensor([.3, .6, .12], dtype=torch.float64)
        logits = torch.logit((initial-self.lower)/(self.upper-self.lower))
        if kind == 'gnn':
            self.encoder = nn.Sequential(nn.Linear(5, hidden), nn.Tanh())
            self.message = nn.Sequential(nn.Linear(2*hidden+1, hidden), nn.Tanh())
            self.update = nn.Sequential(nn.Linear(2*hidden, hidden), nn.Tanh())
            self.head = nn.Sequential(nn.Linear(2*hidden+1, hidden), nn.Tanh(), nn.Linear(hidden, 3))
            self.double()
            nn.init.normal_(self.head[-1].weight, std=.01)
            with torch.no_grad():
                self.head[-1].bias.copy_(logits)
        else:
            count = 1 if kind == 'uniform' else regions
            self.logits = nn.Parameter(logits.repeat(count, 1))

    def regional_values(self, graph):
        if self.kind == 'gnn':
            i, j = graph.edges.T
            h = self.encoder(graph.features)
            length = graph.lengths[:, None]
            # Same message and update weights at both endpoints, two rounds.
            for _ in range(2):
                symmetric = torch.cat((h[i]+h[j], (h[i]-h[j]).abs(), length), dim=1)
                messages = self.message(symmetric)
                total = torch.zeros_like(h).index_add(0, i, messages).index_add(0, j, messages)
                degree = h.new_zeros(len(h)).index_add(0, i, torch.ones_like(length[:, 0])).index_add(0, j, torch.ones_like(length[:, 0]))
                h = self.update(torch.cat((h, total/degree.clamp_min(1)[:, None]), dim=1))
            edge = torch.cat((h[i]+h[j], (h[i]-h[j]).abs(), length), dim=1)
            pooled = edge.new_zeros((self.regions, edge.shape[1])).index_add(0, graph.region, edge)
            count = torch.bincount(graph.region, minlength=self.regions).clamp_min(1)
            logits = self.head(pooled / count[:, None])
        else:
            logits = self.logits.expand(self.regions, -1)
        return self.lower + (self.upper-self.lower)*torch.sigmoid(logits)

    def forward(self, graph):
        return decode(self.regional_values(graph), graph, self.threshold)


def regularization(values):
    # Weak relative spatial smoothness; no penalty favouring damage or rupture.
    return (torch.diff(torch.log(values), dim=0).square().mean()
            if len(values) > 1 else values.sum()*0)
