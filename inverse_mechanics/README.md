# Minimal physics-informed network inference

The pipeline is **reference graph → shared GNN → bounded regional mechanical
parameters → explicit differentiable simulator → multiframe trajectory loss**.
This first stage tests effective mechanics on a fixed candidate graph. It does
not claim to identify actual cell junctions, absolute moduli, or unique rupture
mechanisms. See [RESULTS.md](RESULTS.md) for the completed measurements.

## Running

From the repository root, with `requirements-inverse.txt` installed:

```sh
rtk proxy .venv/bin/python -m unittest discover -s tests/python -v
rtk proxy .venv/bin/python -m inverse_mechanics.validate --epochs 60 --out results/gnn_synthetic
rtk proxy .venv/bin/python -m inverse_mechanics.run --epochs 60 --out results/gnn_experiment
```

Use new output directories for subsequent runs. Experimental fitting requires
a completed synthetic report matching the Python source and numerical tests.
This gate establishes that the audit was run; it does **not** certify unique
mechanical recovery. Every synthetic failure or ambiguity remains in the report.

The experimental defaults retain frames 0–34, fit frames 0–23, and forecast
frames 24–34. `--families zener` runs only that GNN family plus the three
Zener controls. `--epochs`, `--seed`, `--max-step`, `--anchor-k`, `--adhesion-k`,
`--train-frames` and `--time-scale` expose the few useful numerical choices.
These choices must be fixed without optimizing on the forecast frames.

To regenerate independent fixtures with MATLAB R2024a or a compatible release:

```sh
rtk proxy /Applications/MATLAB_R2024a.app/bin/matlab -batch "addpath('tests'); matlab_gnn_reference('tests/python/fixtures/matlab_reference.json')"
```

The checked-in JSON allows routine Python tests without MATLAB. The generator
calls the original, unchanged MATLAB RHS and `ode15s`; it is not a second
implementation of the Python equations.

## Small module layout

| Module | Responsibility |
| --- | --- |
| `data.py` | Calibration, persistent row IDs, MATLAB initial orientation, kNN graph and reference masks |
| `mechanics.py` | Contraction, five constitutive families, force balance, internal states, damage, adhesion, midpoint rollout |
| `model.py` | Two message-passing rounds, symmetric edge pooling, three regional coefficients; direct and uniform controls |
| `fit.py` | Prefix-only optimization, displacement metrics, motion baselines, local sensitivity |
| `synthetic.py` | Known 3D network, localized failure, noisy observations and rupture summaries |
| `validate.py` | Numerical test gate, recovery/discrimination/misspecification audit |
| `run.py` | Experimental comparisons and portable arrays, reports and one static figure |

NumPy and PyTorch implement the model. SciPy reads existing MATLAB data and
checks graph connectivity. Matplotlib writes the scientific comparison figure.
No PyTorch Geometric, dashboards or configuration framework are needed.

## Data, geometry and leakage

`trackedNodes` is N × 3 × T. Original one-based MATLAB row IDs remain attached
to every node and edge. Their persistence is an acquisition assumption, not a
newly verified tracking result. As in MATLAB, frame-zero invalid/zero-placeholder
rows and duplicate coordinates are excluded. Later missing observations become
NaN and are masked from losses and scoring; future availability never selects
model nodes. Tracks are not interpolated or reidentified.

For `trackedNodes_nuc_48_2_OA7_1.mat`, xyz calibration is
`[0.4649160043116627, 0.4649160043116627, 1]` µm per input unit, from
`src/experimental_coordinate_units.m`. The existing acquisition convention is
1,800 seconds over 45 frames, i.e. 44 intervals. Frame 34 is therefore at
1,390.909 seconds; retained frames are never stretched to 1,800 seconds.
`load_tracks` accepts explicit full-acquisition timestamps and spacing for
other datasets. Timing is not available from this MAT array itself.

The initial median translation and principal-axis orientation reproduce
`follicle_rupture_msd_fixed.m`, including the long-axis sign. They are applied
once to every frame. No future registration transform, inferred moving surface,
observed velocity, or observed future position enters a mechanical rollout.

Each node chooses its nearest six reference neighbors; undirected pairs are
sorted and deduplicated. Lengths are measured from that same initial geometry.
There is no silent spanning-tree completion or added edge. Connectivity is
reported. The original MATLAB tuned graph is a different object and is not used
as the candidate graph.

Three equal-width axial regions bound spatial complexity. MATLAB-inspired
reference masks select the P90 anchor cap, axial contraction links with alignment
at least 0.35, and up to 22 links whose endpoints lie in the P0–10 tip cap.
The tip mask is an explicit anatomical/model prior, not inferred damage location.

## Mechanics

State is `(x, sigma, z, D, A)`: absolute 3D position, internal edge force,
Burgers internal extension, irreversible edge damage, and node adhesion release.
Initial sigma, z, D and A are zero unless a state is explicitly supplied.
Reference lengths are immutable. `natural_lengths` ports the MATLAB integrated
smoothstep ramp, plateau, ramp-down and rest-length floor, including `Ln_dot`.
Training bounds keep contraction duration longer than the ramp.

For extension `e = length − Ln` and rate `edot = length_dot − Ln_dot`:

| Family | Passive axial force | Internal-state evolution |
| --- | --- | --- |
| Elastic | `k1 e` | none |
| Maxwell | `sigma` | `sigma_dot = k1 edot − (k1/cm) sigma` |
| Kelvin–Voigt | `k1 e + cm edot` | none |
| Zener | `k1 e + sigma` | `sigma_dot = k2 edot − (k2/cm) sigma` |
| Burgers | `sigma` | `z_dot = (sigma − k1 z)/cs`; `sigma_dot = k2(edot − z_dot) − (k2/cm) sigma` |

Damage scales both springs by `max(1−D, stiffness_floor)` and both dashpots by
`max(1−D, damping_floor)`, with MATLAB's small viscosity floor. A final `1−D`
multiplier on transmitted force guarantees no residual load from internal stress
after complete failure. This deliberately preserves MATLAB's double integrity
scaling for instantaneous elastic terms; it is not a single generic edge gate.

Damage is driven by extension relative to **initial** reference length. Hard mode
matches the MATLAB extension threshold and optional patch cascade RHS. Default
training uses a compact smoothstep over 0.01 strain **above** threshold: damage
is exactly zero below threshold, and `D_dot = (1−D) activity / tau_damage >= 0`.
This is a documented differentiable surrogate, not an exact event-time gradient.
No-damage controls disable this evolution. The fixed threshold is 0.05 strain.
Per-edge failure continues only under overload unless cascade mode is enabled;
damage never heals. The smooth cascade is an approximation to MATLAB's latch.

Reference-position adhesion springs preload the network. The maximum selected
edge damage gates irreversible, load-sensitive adhesion release. Adhesion drag
interpolates between stuck and sliding values, and stored network stresses then
drive retraction. Anchored nodes never release. A zero rupture gate prevents
adhesion detachment before damage. Boundaries are explicit reference tethers and
isotropic drag; no external geometry is prescribed.

At every RHS evaluation the model solves `0 = F_internal + F_boundary − C v`.
Most families admit a local division by positive drag. Kelvin–Voigt additionally
couples endpoint velocities: its global positive-definite drag system is solved
directly for small graphs and by residual-checked differentiable preconditioned
conjugate gradients for large graphs. Node velocity is not learned or lagged.

Fixed-step explicit midpoint integration retains autograd through every substep.
Damage/release times limit the maximum step; high stiffness may require further
refinement. Nonfinite rollouts and failed velocity solves raise errors. Reports
compare each forecast with half-step integration. Hard events are step-resolved,
whereas MATLAB's cell-scale `ode15s` implementation localizes its cascade event.

The following MATLAB features are intentionally outside this minimal version:
ellipsoid penalty/projection and normal damping, externally imposed traction,
moving surfaces, inertia, and exact event restarts. Comparisons use MATLAB's
`surfaceConstraintMode='none'` with equivalent explicit tethers and drag.

## GNN and inverse problem

Five node features contain initial centered xyz, normalized axial position and
anchor membership. Shared message/update networks perform two rounds. Edge
features use endpoint sums, absolute differences and reference length; reversing
an edge or permuting cell order leaves predictions equivalent. Regional pooling
produces at most nine physical values, with sigmoid bounds:

| Parameter | Bounds in model units |
| --- | --- |
| Stiffness / fixed node drag | 0.04–1.2 |
| Relaxation time | 0.12–1.8 |
| Fractional natural-length contraction rate | 0.015–0.35 |

`k2=k1`, `cm=k1*tau`, and `cs=k1*tau` fix otherwise poorly constrained branch
ratios. Elastic relaxation has no effect and must not be described as recovered.
One constitutive family is selected globally per fit. GNN weights are numerous,
but the mechanical simulator receives only three values per region, not free
parameters on every edge. On a single graph a direct regional fit is an especially
important baseline; the GNN has no established cross-follicle generalization.

Length unit is median initial edge length; default time unit is 300 s. Node drag
is fixed to one. Absolute force, stiffness and viscosity cannot be recovered
without an independent drag/force calibration. The tests explicitly demonstrate
the invariance under simultaneous force/drag rescaling.

Optimization uses multistep free rollouts from frame zero and masked position MSE
(equivalently displacement MSE with fixed initial positions). Normalization uses
training displacement only. A weak log-parameter regional smoothness penalty is
added. Adam selection uses the lowest **training** objective. Forecasts continue
from the same initial state with frozen learned parameters and no teacher forcing.
Reported future scores do not choose a family or tune settings automatically.

Uniform and directly fitted regional Zener models, a GNN Zener no-damage model,
stationary motion, and training-prefix linear motion provide controls. Held-out
errors include 3D vector RMSE, displacement-magnitude MAE and angular error.
Angles require both observed and predicted displacement above 0.1 µm; the report
includes valid direction counts and stationary predictions rather than assigning
an arbitrary direction to zero motion.

## Validation, outputs and uncertainty

Tests compare all five families against original MATLAB derivatives and short
`ode15s` trajectories, including natural contraction, damping, partially damaged
states, adhesion and anchors. They also check failure forces, irreversibility,
finite-difference gradients, symmetric/permutation-equivariant inference,
large-graph force balance, timestep convergence and future-data isolation.

The synthetic audit fits each known family, evaluates family discrimination on
the same Zener data, and runs GNN, no-damage, incorrect adhesion/topology and two
noise controls. It reports parameter error, sampled rupture time, failed-edge IoU
and timestep sensitivity. The synthetic susceptibility mask and topology are
known; localization scores test mechanical state evolution conditional on that
mask and do not establish discovery of an unknown rupture site.

`validation.json` records the audit. Experimental `report.json` records units,
splits, mechanics, all model scores, rupture progression and numerical refinement.
Each model directory contains `weights.pt` and `forecast.npz`, with observed and
simulated xyz, persistent IDs, reference lengths, edge cell-ID pairs, regional
parameters, effective integrity, damage, internal states and adhesion release.
The GNN Zener directory also includes `comparison.png` with 3D trajectories,
errors, candidate connectivity/stiffness and damage progression.

Training-only finite differences in log regional parameters provide singular
values, weak directions and a numerical rank. Local linearized standard errors
are reported only for full rank; otherwise unconstrained directions are left
unbounded. These are conditional iid-residual approximations, not calibrated
confidence intervals. Noise controls expose sensitivity to observations, but
neither accounts fully for identity, boundary, topology or model discrepancy.

## Recovery and future extensions

Before replacement, old Python code, tests and generated Python outputs were
saved to Git commit `34a7c414706a88d2c91b451c79e77615b277773a`, retained by branch
`checkpoint/python-before-gnn-rebuild-20261104`. All 348 selected files were
verified byte-for-byte against that commit. `checkpoint.json` lists the paths.
Use `git archive` to recover them into a separate directory without overwriting
the rebuild. The normal HEAD and staging index were not changed by checkpointing.

`preservation.json` hashes the 190 pre-existing MATLAB, dataset and non-Python
result files. They remain unchanged. Original MATLAB reference tests are retained.
Only the old Python implementation/tests, its obsolete generated output paths,
and obsolete README links were replaced.

Topology refinement, learned damage susceptibility, unknown rupture localization,
edge-specific constitutive mixtures, free branch ratios, learned initial stress,
and full uncertainty propagation are **future extensions**. They should be added
one at a time only when synthetic identifiability and experimental predictive
performance justify the extra freedom.
