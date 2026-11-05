from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
import numpy as np
import torch
from scipy.io import savemat

from inverse_mechanics.data import load_tracks, candidate_edges, make_graph
from inverse_mechanics.fit import fit, motion_baselines, metrics
from inverse_mechanics.mechanics import coupled_velocity, edge_sum
from inverse_mechanics.synthetic import example
from inverse_mechanics.model import RegionalModel


torch.set_num_threads(1)


class DataAndFitTests(unittest.TestCase):
    def test_initial_geometry_and_clock_ignore_future(self):
        rng = np.random.default_rng(0)
        raw = rng.uniform(1,10,(10,3,45))
        raw[0] = 0
        raw[1,:,0] = raw[2,:,0]  # duplicate reference row
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'trackedNodes_nuc_48_2_OA7_1.mat'
            savemat(path,{'trackedNodes':raw})
            first = load_tracks(path)
            raw[3:,:,24:] = np.nan
            savemat(path,{'trackedNodes':raw})
            second = load_tracks(path)
        np.testing.assert_array_equal(first.ids,second.ids)
        np.testing.assert_allclose(first.positions[:24],second.positions[:24])
        np.testing.assert_array_equal(candidate_edges(first.positions[0]),candidate_edges(second.positions[0]))
        self.assertEqual(len(first.times),35)
        self.assertAlmostEqual(first.times[-1],34*1800/44)
        self.assertTrue(np.isnan(second.positions[24:]).any())
        self.assertNotIn(1,first.ids)
        self.assertNotIn(3,first.ids)

    def test_unknown_acquisition_needs_calibration(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'unknown.mat'
            savemat(path,{'trackedNodes':np.ones((5,3,4))})
            with self.assertRaises(ValueError):
                load_tracks(path,frame_limit=4)

    def test_future_never_changes_fit_or_baselines(self):
        graph,_,time,dynamics,boundary,_,observed=example()
        changed=observed.copy(); changed[8:]=1e6
        a,_=fit(graph,observed,time,8,dynamics,boundary,epochs=2,seed=7)
        b,_=fit(graph,changed,time,8,dynamics,boundary,epochs=2,seed=7)
        torch.testing.assert_close(a.regional_values(graph),b.regional_values(graph),rtol=0,atol=0)
        for name,array in motion_baselines(observed,time,8).items():
            np.testing.assert_array_equal(array,motion_baselines(changed,time,8)[name])

    def test_permutation_equivariance(self):
        graph,*_=example()
        model=RegionalModel(graph)
        permutation=torch.randperm(len(graph.x0))
        inverse=torch.argsort(permutation)
        other=replace(graph,x0=graph.x0[permutation],features=graph.features[permutation],
                      anchor=graph.anchor[permutation],edges=inverse[graph.edges])
        torch.testing.assert_close(model.regional_values(graph),model.regional_values(other))

    def test_large_coupled_force_balance(self):
        rng=np.random.default_rng(1)
        graph=make_graph(rng.normal(size=(40,3)))
        i,j=graph.edges.T
        direction=(graph.x0[j]-graph.x0[i])/graph.lengths[:,None]
        damping=torch.full_like(graph.lengths,.3,requires_grad=True)
        drag=torch.ones(40,dtype=torch.float64)
        force=torch.tensor(rng.normal(size=(40,3)),dtype=torch.float64)
        velocity=coupled_velocity(graph,direction,damping,drag,force)
        residual=drag[:,None]*velocity+edge_sum(graph,(damping*((velocity[i]-velocity[j])*direction).sum(1))[:,None]*direction)-force
        self.assertLess(float(residual.detach().norm()),1e-8)
        velocity.square().sum().backward()
        self.assertTrue(torch.isfinite(damping.grad).all())

    def test_missing_metrics_and_stationary_direction(self):
        y=np.ones((3,4,3)); y[2,0]=np.nan
        p=np.zeros((3,4,3))
        report=metrics(p,y,np.zeros((4,3)),1)
        self.assertEqual(report['valid_samples'],7)
        self.assertIsNone(report['direction_deg'])
        self.assertEqual(report['stationary_predictions'],7)


if __name__=='__main__':
    unittest.main()
