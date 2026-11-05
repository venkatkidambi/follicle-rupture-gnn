import json
from dataclasses import replace
from pathlib import Path
import unittest

import numpy as np
import torch

from inverse_mechanics.data import make_graph
from inverse_mechanics.mechanics import (Boundary, Dynamics, State, FAMILIES,
    natural_lengths, constitutive, initial_state, rhs, rollout, edge_sum, coupled_velocity)
from inverse_mechanics.model import decode, RegionalModel
from inverse_mechanics.synthetic import example


torch.set_num_threads(1)


def fixture_setup():
    graph = make_graph([[0,0,0],[1,.2,.1],[2,-.1,.3]], [[0,1],[1,2]], regions=1)
    graph.anchor = torch.tensor([False,False,True])
    graph.active = torch.ones(2, dtype=torch.bool)
    graph.weak = torch.tensor([True,False])
    p = {key: torch.tensor(value, dtype=torch.float64) for key,value in dict(
        k1=[.7,.9], k2=[.8,1.1], cm=[.3,.4], cs=[.25,.35], rate=[.2,.15], threshold=[.02,.02]).items()}
    state = State(graph.x0+torch.tensor([[0,0,0],[.08,-.01,.02],[.05,.01,-.01]],dtype=torch.float64),
                  torch.tensor([.12,-.03],dtype=torch.float64), torch.tensor([.02,.01],dtype=torch.float64),
                  torch.tensor([.7,0.],dtype=torch.float64), torch.tensor([.1,.2,0.],dtype=torch.float64))
    state = State(*(v.double() for v in state))
    boundary = Boundary(anchor_drag=.2, adhesion_k=.7, load_threshold=.008)
    return graph, p, state, boundary


class MechanicsTests(unittest.TestCase):
    def test_constitutive_hand_values(self):
        a = [torch.tensor(v, dtype=torch.float64) for v in (.2,.3,.4,.1,2.,3.,4.,5.)]
        expected = {'elastic':(.4,0,0), 'maxwell':(.4,.4,0),
                    'kelvin-voigt':(1.6,0,0), 'zener':(.8,.6,0), 'burgers':(.4,.48,.04)}
        for family in FAMILIES:
            np.testing.assert_allclose([float(v) for v in constitutive(family,*a)], expected[family], atol=1e-12)

    def test_matlab_rhs_and_contraction(self):
        path = Path('tests/python/fixtures/matlab_reference.json')
        self.assertTrue(path.exists(), 'Generate MATLAB reference fixture first')
        data = json.loads(path.read_text())
        graph,p,state,boundary = fixture_setup()
        for case in data['cases']:
            derivative = rhs(.4,state,graph,p,Dynamics(family=case['family'],ramp=.1,damage_width=0),boundary)
            actual = torch.cat([v.flatten() for v in derivative]).numpy()
            np.testing.assert_allclose(actual, case['derivative'], atol=2e-7, rtol=2e-6)
            self.assertLess(np.abs(case['net_force']).max(), 1e-12)
        for case in data['contraction']:
            length, rate = natural_lengths(case['time'], graph.lengths,p['rate'],ramp=.1)
            np.testing.assert_allclose(length,case['length'],atol=1e-12)
            np.testing.assert_allclose(rate,case['rate'],atol=1e-12)

    def test_matlab_integrated_states(self):
        data = json.loads(Path('tests/python/fixtures/matlab_reference.json').read_text())
        graph,p,state,boundary = fixture_setup()
        p['threshold'] = torch.full((2,), float('inf'))
        for case in data['cases']:
            result = rollout(graph,p,case['times'],Dynamics(family=case['family'],ramp=.1,damage_width=0),boundary,.0005,state)
            actual = torch.cat([v.reshape(len(case['times']),-1) for v in result],dim=1).numpy()
            np.testing.assert_allclose(actual,case['history'],atol=3e-7,rtol=3e-6)

    def test_force_conservation_and_failure(self):
        graph,p,state,boundary = fixture_setup()
        free = Boundary(anchor_k=0,adhesion_k=0,stuck_drag=0,sliding_drag=0)
        for family in FAMILIES:
            d = rhs(.4,state,graph,p,Dynamics(family=family),free)
            np.testing.assert_allclose(d.x.sum(0),0,atol=1e-12)
            failed = state._replace(damage=torch.ones_like(state.damage))
            d = rhs(.4,failed,graph,p,Dynamics(family=family),free)
            np.testing.assert_allclose(d.x,0,atol=1e-12)

    def test_damage_release_irreversibility(self):
        graph,p,state,boundary = fixture_setup()
        state = state._replace(damage=torch.tensor([.7,0.],dtype=torch.float64))
        result = rollout(graph,p,torch.linspace(.4,.8,15),Dynamics(damage_width=.01),boundary,.002,state)
        self.assertTrue((torch.diff(result.damage,dim=0)>=-1e-12).all())
        self.assertTrue((torch.diff(result.release,dim=0)>=-1e-12).all())
        self.assertGreater(float(result.release[-1,0]-result.release[0,0]),0)
        self.assertEqual(float(result.release[-1,-1]),0.)
        self.assertTrue((result.damage<=1).all())

    def test_subthreshold_and_no_damage(self):
        graph,p,state,boundary = fixture_setup()
        state = state._replace(damage=torch.zeros_like(state.damage),release=torch.zeros_like(state.release))
        p['threshold'] = torch.full((2,),10.)
        d = rhs(.4,state,graph,p,Dynamics(),boundary)
        self.assertEqual(float(d.damage.sum()+d.release.sum()),0.)

    def test_gradient_finite_difference(self):
        graph,values,t,d,b,_,_ = example()
        values = values.clone().requires_grad_()
        def objective(v):
            result = rollout(graph,decode(v,graph),t[:7],d,b,.015)
            return result.x[-1].square().sum()
        gradient = torch.autograd.grad(objective(values),values)[0]
        for index in ((0,0),(1,1),(2,2)):
            delta = torch.zeros_like(values); delta[index] = 1e-5
            finite = (objective(values+delta)-objective(values-delta))/(2e-5)
            self.assertAlmostEqual(float(gradient[index]),float(finite.detach()),places=6)
        self.assertGreater(float(gradient.abs().min()),1e-8)

    def test_damage_threshold_gradient(self):
        graph,p,state,boundary = fixture_setup()
        natural,_ = natural_lengths(.4,graph.lengths,p['rate'],ramp=.2)
        extension = torch.linalg.vector_norm(state.x[graph.edges[:,1]]-state.x[graph.edges[:,0]],dim=1)-natural
        p['threshold'] = (extension/graph.lengths-.005).detach().requires_grad_()
        d = rhs(.4,state,graph,p,Dynamics(damage_width=.01),boundary)
        grad = torch.autograd.grad(d.damage.sum(),p['threshold'])[0]
        self.assertLess(float(grad[0]),0.)
        self.assertEqual(float(grad[1]),0.)

    def test_midpoint_convergence(self):
        graph,values,t,d,b,_,_ = example()
        p = decode(values,graph)
        with torch.no_grad():
            solutions = [rollout(graph,p,t,d,b,h) for h in (.02,.01,.005)]
        coarse = torch.linalg.vector_norm(solutions[0].x-solutions[2].x)
        fine = torch.linalg.vector_norm(solutions[1].x-solutions[2].x)
        self.assertLess(float(fine),float(coarse)*.5)

    def test_gnn_symmetry_and_gradients(self):
        graph,_,times,d,b,_,_ = example()
        torch.manual_seed(2)
        model = RegionalModel(graph)
        expected = model.regional_values(graph)
        reversed_graph = replace(graph,edges=graph.edges.flip(1))
        torch.testing.assert_close(expected,model.regional_values(reversed_graph))
        result = rollout(graph,model(graph),times[:4],d,b)
        result.x[-1].square().sum().backward()
        self.assertGreater(sum(float(p.grad.abs().sum()) for p in model.parameters()),0)
        self.assertTrue(all(torch.isfinite(p.grad).all() for p in model.parameters()))

    def test_drag_scale_nonidentifiability(self):
        graph,values,t,d,b,_,_ = example()
        p = decode(values,graph)
        factor = 3.
        scaled = {k:(v*factor if k in ('k1','k2','cm','cs') else v) for k,v in p.items()}
        b2 = replace(b, **{key:getattr(b,key)*factor for key in ('drag','anchor_k','anchor_drag','adhesion_k','stuck_drag','sliding_drag','load_threshold')})
        with torch.no_grad():
            a = rollout(graph,p,t,d,b)
            c = rollout(graph,scaled,t,d,b2)
        torch.testing.assert_close(a.x,c.x,atol=1e-11,rtol=1e-10)


if __name__ == '__main__':
    unittest.main()
