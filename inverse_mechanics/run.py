"""Synthetic-first experimental comparison, with free future-frame forecasts."""
import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path

import numpy as np
import torch
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .data import load_tracks, make_graph, candidate_edges
from .fit import fit, metrics, motion_baselines, sensitivity
from .mechanics import FAMILIES, Dynamics, Boundary, rollout
from .synthetic import rupture_summary
from .validate import source_digest


def plot_result(out, observed, predicted, times, graph, parameters, damage, train_frames):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    fig=plt.figure(figsize=(13,9),layout='constrained')
    ax=fig.add_subplot(221,projection='3d')
    selected=np.linspace(0,len(graph.x0)-1,min(18,len(graph.x0)),dtype=int)
    for node in selected:
        ax.plot(*observed[:,node].T,color='black',alpha=.5,lw=.8)
        ax.plot(*predicted[:,node].T,color='tab:orange',alpha=.8,lw=.8)
    ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)',title='Observed (black) and simulated (orange)')
    ax=fig.add_subplot(222)
    error=np.linalg.norm(predicted-observed,axis=-1)
    ax.plot(times,np.sqrt(np.nanmean(error**2,axis=1)))
    ax.axvline(times[train_frames-1],color='black',ls='--')
    ax.set(xlabel='Time (s)',ylabel='Position RMSE (µm)',title='Rollout error; dashed line ends training')
    ax=fig.add_subplot(223,projection='3d')
    lines=graph.x0.numpy()[graph.edges.numpy()]
    collection=Line3DCollection(lines,cmap='viridis',linewidths=.7)
    collection.set_array(parameters['k1'].numpy())
    # Fixed physical range avoids turning roundoff-level regional differences
    # into an apparent heterogeneous mechanical map.
    from .model import BOUNDS
    collection.set_clim(*BOUNDS[0])
    ax.add_collection3d(collection)
    ax.auto_scale_xyz(*graph.x0.numpy().T)
    ax.set(title='Candidate graph: effective stiffness',xlabel='x / L',ylabel='y / L',zlabel='z / L')
    fig.colorbar(collection,ax=ax,shrink=.5,label='k / fixed drag (model units)')
    ax=fig.add_subplot(224)
    im=ax.imshow(damage.T,aspect='auto',origin='lower',extent=[times[0],times[-1],0,damage.shape[1]],vmin=0,vmax=1,cmap='magma')
    ax.axvline(times[train_frames-1],color='cyan',ls='--')
    ax.set(xlabel='Time (s)',ylabel='Candidate edge index',title='Conditional damage; fixed susceptibility mask')
    fig.colorbar(im,ax=ax,label='Damage D')
    fig.savefig(out/'comparison.png',dpi=160)
    plt.close(fig)


def run(args):
    validation=json.loads(args.validation.read_text())
    if not validation.get('tests_passed') or not validation.get('synthetic_audit_complete'):
        raise ValueError('Complete the synthetic audit before experimental fitting')
    if validation.get('source_digest')!=source_digest():
        raise ValueError('Source changed since validation; rerun synthetic audit')
    if args.out.exists() and any(args.out.iterdir()):
        raise ValueError('Choose an empty output directory; previous results are preserved')
    torch.set_num_threads(1)
    tracks=load_tracks(args.data,frame_limit=args.frame_limit)
    if not 3<=args.train_frames<len(tracks.times):
        raise ValueError('Training prefix must leave at least one future frame')
    edges=candidate_edges(tracks.positions[0])
    lengths=np.linalg.norm(tracks.positions[0,edges[:,1]]-tracks.positions[0,edges[:,0]],axis=1)
    length_scale=float(np.median(lengths))
    time_scale=args.time_scale
    if time_scale<=0:
        raise ValueError('Positive fixed time scale required')
    observed=tracks.positions/length_scale
    times=tracks.times/time_scale
    graph=make_graph(observed[0],edges)
    adjacency=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(graph.x0),)*2)
    components=int(connected_components(adjacency,directed=False)[0])
    boundary=Boundary(anchor_k=args.anchor_k,adhesion_k=args.adhesion_k)
    dynamics=Dynamics(ramp=.2,damage_width=.01)
    args.out.mkdir(parents=True,exist_ok=True)
    report=dict(metadata=tracks.metadata, nodes=len(graph.x0), candidate_edges=len(edges),connected_components=components,
        train_frame_indices=list(range(args.train_frames)),heldout_frame_indices=list(range(args.train_frames,len(times))),
        length_scale_um=length_scale,time_scale_s=time_scale,force_scale='node drag fixed to 1; no absolute force calibration',
        graph='exact MATLAB initial symmetric six-neighbor candidate graph, not verified junctions',
        boundary=asdict(boundary),dynamics=asdict(dynamics),source_digest=source_digest(),
        validation=str(args.validation),models={},motion_baselines={},
        evaluation='frame-zero displacements in one fixed initial coordinate transform; no future registration or teacher forcing',
        limitations=['Fixed topology and fixed tip susceptibility; neither was discovered.',
                     'Constitutive family is global; equal spring/branch ratios are fixed assumptions.',
                     'Three axial regions; only stiffness, relaxation and contraction learned.',
                     'Free 3D reference tethers; no MATLAB ellipsoid projection, moving surface, external traction or inertia.',
                     'Smooth damage surrogate; no exact event-time gradients or experimental rupture ground truth.',
                     'One trajectory and local conditional uncertainty do not establish unique mechanics.'])
    for name,predicted in motion_baselines(tracks.positions,tracks.times,args.train_frames).items():
        report['motion_baselines'][name]=metrics(predicted,tracks.positions,tracks.positions[0],args.train_frames)
    hypotheses=[(f'gnn_{family}','gnn',replace(dynamics,family=family)) for family in args.families]
    hypotheses += [('uniform_zener','uniform',dynamics),('regional_zener','regional',dynamics),
                   ('gnn_zener_no_damage','gnn',replace(dynamics,damage_enabled=False))]
    for name,kind,hypothesis in hypotheses:
        model,training=fit(graph,observed,times,args.train_frames,hypothesis,boundary,
                           kind=kind,epochs=args.epochs,seed=args.seed,max_step=args.max_step)
        with torch.no_grad():
            parameters=model(graph)
            predicted=rollout(graph,parameters,times,hypothesis,boundary,args.max_step)
            refined=rollout(graph,parameters,times,hypothesis,boundary,args.max_step/2)
        physical=predicted.x.numpy()*length_scale
        record=dict(kind=kind,family=hypothesis.family,damage_enabled=hypothesis.damage_enabled,training=training,
            regional_parameters=model.regional_values(graph).detach().tolist(),
            learned_mechanical_values=3*(1 if kind=='uniform' else 3),
            heldout=metrics(physical,tracks.positions,tracks.positions[0],args.train_frames),
            rupture=rupture_summary(predicted.damage,tracks.times),
            max_step_refinement_difference_um=float((predicted.x-refined.x).abs().max())*length_scale,
            max_adhesion_release=float(predicted.release.max()),
            mean_final_integrity=float((1-predicted.damage[-1]).mean()))
        report['models'][name]=record
        model_out=args.out/name; model_out.mkdir()
        np.savez_compressed(model_out/'forecast.npz',times_s=tracks.times,ids=tracks.ids,edges=edges,
            edge_cell_ids=tracks.ids[edges],reference_lengths_um=lengths,region=graph.region.numpy(),
            observed_um=tracks.positions,predicted_um=physical,damage=predicted.damage.numpy(),
            integrity=1-predicted.damage.numpy(),release=predicted.release.numpy(),sigma=predicted.sigma.numpy(),
            z=predicted.z.numpy(),active_mask=graph.active.numpy(),susceptibility_mask=graph.weak.numpy(),
            **{key:value.numpy() for key,value in parameters.items()})
        torch.save(model.state_dict(),model_out/'weights.pt')
        if name=='gnn_zener':
            record['sensitivity']=sensitivity(graph,model.regional_values(graph),observed,times,args.train_frames,hypothesis,boundary,args.max_step)
            plot_result(model_out,tracks.positions,physical,tracks.times,graph,parameters,predicted.damage.numpy(),args.train_frames)
        (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        print(f'{name}: heldout RMSE {record["heldout"]["vector_rmse"]:.3f} µm; magnitude MAE {record["heldout"]["magnitude_mae"]:.3f} µm; rupture {record["rupture"]["first_rupture_time"]}',flush=True)
    report['mechanical_network_discovery_established']=False
    report['conclusion']='Comparative effective-model fit only. Fixed graph/damage prior, parameter ambiguity, and baseline comparisons preclude a claim of discovered mechanical topology or experimental rupture.'
    (args.out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path('input_node_data/tracked_nodes/trackedNodes_nuc_48_2_OA7_1.mat'))
    parser.add_argument('--out',type=Path,default=Path('results/gnn_experiment'))
    parser.add_argument('--validation',type=Path,default=Path('results/gnn_synthetic/validation.json'))
    parser.add_argument('--frame-limit',type=int,default=35)
    parser.add_argument('--train-frames',type=int,default=24)
    parser.add_argument('--epochs',type=int,default=60)
    parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--max-step',type=float,default=.02)
    parser.add_argument('--time-scale',type=float,default=300.)
    parser.add_argument('--anchor-k',type=float,default=1.)
    parser.add_argument('--adhesion-k',type=float,default=.3)
    parser.add_argument('--families',nargs='+',choices=FAMILIES,default=list(FAMILIES))
    args=parser.parse_args()
    run(args)


if __name__=='__main__':
    main()
