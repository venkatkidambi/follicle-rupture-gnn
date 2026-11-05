"""Reproducible synthetic recovery, discrimination and misspecification audit."""
import argparse
from dataclasses import replace, asdict
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
import torch

from .fit import fit, metrics, sensitivity
from .mechanics import FAMILIES, rollout
from .model import decode
from .synthetic import example, rupture_summary


def source_digest():
    paths = sorted(Path('inverse_mechanics').glob('*.py')) + sorted(Path('tests/python').rglob('*.py'))
    paths += [Path('tests/python/fixtures/matlab_reference.json')]
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path).encode()); digest.update(path.read_bytes())
    return digest.hexdigest()


def run_validation(out, epochs=60):
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use an empty output directory')
    out.mkdir(parents=True,exist_ok=True)
    torch.set_num_threads(1)
    suite=unittest.defaultTestLoader.discover('tests/python')
    checks=unittest.TextTestRunner(verbosity=1).run(suite)
    if not checks.wasSuccessful():
        raise RuntimeError('Mechanical tests failed; experimental inference is disabled')
    report=dict(source_digest=source_digest(), tests_run=checks.testsRun, tests_passed=True,
                dtype='float64', epochs=epochs, train_frames=19, recovery={}, discrimination={}, controls={})
    zener_case=None
    for family in FAMILIES:
        graph,values,t,d,b,truth,observed=example(family)
        model,history=fit(graph,observed,t,19,d,b,kind='regional',epochs=epochs,penalty=0)
        with torch.no_grad():
            predicted=rollout(graph,model(graph),t,d,b)
            recovered=model.regional_values(graph)
            refined=rollout(graph,decode(values,graph),t,d,b,.005)
        error=(recovered-values).abs()/values
        true_rupture=rupture_summary(truth.damage,t)
        fitted_rupture=rupture_summary(predicted.damage,t)
        actual_failed=truth.damage[-1].numpy()>=.5
        predicted_failed=predicted.damage[-1].numpy()>=.5
        union=np.logical_or(actual_failed,predicted_failed).sum()
        record=dict(true_parameters=values.tolist(),fitted_parameters=recovered.tolist(),
                    mean_relative_parameter_error=float(error.mean()),max_relative_parameter_error=float(error.max()),
                    training=history,holdout=metrics(predicted.x,observed,observed[0],19,direction_floor=1e-4),
                    true_rupture=true_rupture,fitted_rupture=fitted_rupture,
                    failed_edge_iou=float(np.logical_and(actual_failed,predicted_failed).sum()/max(union,1)),
                    max_dt_refinement_position_difference=float((truth.x-refined.x).abs().max()),
                    refined_rupture=rupture_summary(refined.damage,t))
        report['recovery'][family]=record
        print(f'Synthetic {family}: parameter error {float(error.mean()):.3g}, holdout {record["holdout"]["vector_rmse"]:.4g}',flush=True)
        if family=='zener':
            zener_case=(graph,values,t,d,b,truth,observed,model)
    graph,values,t,d,b,truth,observed,model=zener_case
    # Discrimination uses the same known Zener trajectories for every family.
    for family in FAMILIES:
        hypothesis=replace(d,family=family)
        if family=='zener':
            candidate=model
        else:
            candidate,_=fit(graph,observed,t,19,hypothesis,b,kind='regional',epochs=epochs,penalty=0)
        with torch.no_grad():
            prediction=rollout(graph,candidate(graph),t,hypothesis,b)
        report['discrimination'][family]=metrics(prediction.x,observed,observed[0],19,direction_floor=1e-4)
        print('Discrimination',family,report['discrimination'][family]['vector_rmse'],flush=True)
    for name,kind,noise,seed,hypothesis,boundary in [
        ('gnn','gnn',0.,0,d,b), ('no_damage','regional',0.,0,replace(d,damage_enabled=False),b),
        ('wrong_adhesion','regional',0.,0,d,replace(b,adhesion_k=0.,stuck_drag=0.)),
        ('noise_seed_1','regional',.002,1,d,b), ('noise_seed_2','regional',.002,2,d,b)]:
        noisy=observed.copy()
        noisy[1:]+=np.random.default_rng(seed).normal(0,noise,noisy[1:].shape)
        candidate,history=fit(graph,noisy,t,19,hypothesis,boundary,kind=kind,epochs=epochs,seed=seed,penalty=0)
        with torch.no_grad():
            prediction=rollout(graph,candidate(graph),t,hypothesis,boundary)
        report['controls'][name]=dict(parameters=candidate.regional_values(graph).detach().tolist(),training=history,
            clean_holdout=metrics(prediction.x,observed,observed[0],19,direction_floor=1e-4),
            rupture=rupture_summary(prediction.damage,t))
        print('Control',name,report['controls'][name]['clean_holdout']['vector_rmse'],flush=True)
    # Incorrect topology is fixed before fitting, not optimized using holdout.
    keep=torch.arange(len(graph.edges))!=5
    wrong=replace(graph,edges=graph.edges[keep],lengths=graph.lengths[keep],region=graph.region[keep],
                  active=graph.active[keep],weak=graph.weak[keep])
    candidate,history=fit(wrong,observed,t,19,d,b,kind='regional',epochs=epochs,penalty=0)
    with torch.no_grad():
        prediction=rollout(wrong,candidate(wrong),t,d,b)
    report['controls']['wrong_topology']=dict(training=history,
        clean_holdout=metrics(prediction.x,observed,observed[0],19,direction_floor=1e-4))
    report['sensitivity']=sensitivity(graph,model.regional_values(graph),observed,t,19,d,b)
    report['scope']='Fixed known topology and damage mask; branch ratios fixed; no claim of joint topology or arbitrary constitutive recovery.'
    report['identifiability']='Absolute force scale is unidentifiable (tested drag-rescaling invariance). Recovery metrics do not prove unique parameters; elastic viscosity is structurally unidentifiable.'
    report['synthetic_audit_complete']=True
    (out/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(out/'known_zener.npz',observed=observed,times=t.numpy(),edges=graph.edges.numpy(),
                        parameters=values.numpy(),damage=truth.damage.numpy(),release=truth.release.numpy())
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('results/gnn_synthetic'))
    parser.add_argument('--epochs',type=int,default=60)
    args=parser.parse_args()
    run_validation(args.out,args.epochs)


if __name__=='__main__':
    main()
