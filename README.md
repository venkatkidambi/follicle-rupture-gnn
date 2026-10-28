# Follicle Rupture Shell Simulation

This repository contains the core MATLAB simulation of follicle rupture on an
ellipsoidal epithelial shell. It is intentionally minimal so the mechanics can
be rebuilt from first principles.

The model currently includes:

- generated ellipsoidal or spreadsheet-based surface nodes
- passive viscoelastic connections between neighboring nodes
- localized natural-length contraction and weakest-edge rupture
- optional soft anchoring and surface constraints
- adaptive ODE integration with simulation data saved to `results/`

Post-processing includes shell snapshots, signed long-axis percentile-ring
displacement and velocity plots, and an MP4 animation of the rupture motion.
An additional force figure shows rupture-patch tension and signed retraction
forces in the same percentile bands.
Mechanics tests in `tests/` cover units, force scaling, damage, contraction,
and an analytical overdamped retraction solution.

## Run

Open MATLAB in the repository root and run:

```matlab
run_simulation
```

The simulation writes `results/follicleShellSim_test2.mat`, then exports the
post-processing figures as `results/node_position_snapshots.png` and
`results/percentile_ring_velocity.png`, plus
`results/rupture_retraction_force.png`. Force measurements and the estimated
rupture-onset interval are saved in `results/rupture_retraction_force.mat`
as `forceInfo` (forces in newtons). It also renders the rupture motion to
`results/follicle_rupture_motion.mp4` through the final saved simulation time.
Generated results are ignored by git; `results/.gitkeep` preserves the output
directory.

The default `mechanicsProfile = 'legacy'` retains the original mechanical
coefficients and historical artificial inertia. Both profiles now use the
source image calibration by default, so results differ from the earlier
uncalibrated runs. To try the physical estimates,
set `mechanicsProfile = 'cell-scale'` in `run_simulation.m`. It writes
`follicleShellSim_cell_scale.mat` and adds `_cell_scale` to figure/video names.
Alternatively, run only the mechanics without rendering:

```matlab
addpath('src');
simFile = follicle_rupture_msd_fixed('cell-scale');
% Optional second argument selects a different output MAT-file.
% Third argument selects 'source' (default) or 'input-units' calibration.
```

Each node represents one Drosophila follicle cell. Spatial calibration now
comes from `Step14_findingAxis_v2.m`; mechanical coefficients remain estimates.
To reproduce the old coordinate assumption, set
`coordinateCalibration = 'input-units'` in `run_simulation.m`, or call:

```matlab
follicle_rupture_msd_fixed('legacy', [], 'input-units');
```

This adds `_input_units` to output names. Earlier local results are archived
in `results/before_source_calibration/`. Coordinate calibration and mechanics
profile are independent: source calibration applies to both mechanics profiles.

## Repository Layout

- `run_simulation.m`: adds `src/` to the MATLAB path and starts the model.
- `src/follicle_rupture_msd_fixed.m`: setup, geometry, integration, and output.
- `src/follicle_shell_rhs.m`: shared mechanical equations and force evaluation.
- `src/contracting_natural_lengths.m`: smooth natural-length contraction.
- `src/cell_scale_mechanics.m`: common force scaling and estimated cell mass.
- `src/experimental_coordinate_units.m`: dataset-specific raw xyz calibration.
- `src/shell_jacobian_pattern.m`: sparse coupling pattern for implicit integration.
- `src/integrate_cell_scale.m`: first-rupture event and overdamped integration.
- `src/overdamped_shell_rhs.m`: solves node velocities from instantaneous force balance.
- `plot_node_position_snapshots.m`: initial, middle, and final shell positions.
- `plot_percentile_ring_velocity.m`: ring-average displacement and speed.
- `plot_rupture_retraction_force.m`: rupture-patch tension and ring-average axial tissue pull.
- `export_follicle_rupture_motion_video.m`: MP4 animation of shell rupture.
- `input_node_data/`: optional experimental node-coordinate spreadsheets.
- `results/`: local generated output.

## Main Parameters

Edit the user-parameter section near the top of
`src/follicle_rupture_msd_fixed.m`.

```matlab
cellDiameter_um = 10;
coordinateInputMode = 'microns'; % generated/other inputs; known tracked file uses source calibration
displayMode = 'microns';       % 'cells' or 'microns'
generatedEllipsoidSemiAxes = [50 15 15];

generatedNodeCount = 300;
tEnd = 1800;
numSave = 300;
ruptureSaveEnd = 900;          % dense post-processing window, seconds
ruptureSaveInterval = 1;      % output spacing in that window, seconds
odeSolver = 'ode45';          % legacy; cell-scale selects ode15s
```

Mechanics coefficients and saved data use SI units. Coordinate inputs are
converted to meters once when geometry is loaded. SI labels alone do not
make the legacy inertia physical; see the audit below.

## Tracked Node Data

The default run uses only frame 1 of the ignored local file:

```text
input_node_data/tracked_nodes/trackedNodes_nuc_48_2_OA7_1.mat
```

It loads `trackedNodes(:,:,1)` as the initial node cloud and converts raw
image coordinates to meters **before centering, PCA, ellipsoid fitting, and
six-nearest-neighbor mesh construction**. Later tracked frames are not used
by the mechanics. The snapshot comparison recovers the full affine transform
from the calibrated initial nodes and applies it to all experimental frames.

`Step14_findingAxis_v2.m` loads this exact filename at line 15, defines the
spacings at lines 33–37, and multiplies raw xyz coordinates at lines 99–101:

| Raw coordinate | Spacing in meters | Spacing in µm |
| --- | ---: | ---: |
| x | `2.3245800215583136e-7 / 0.5` | 0.464916004312 |
| y | `2.3245800215583136e-7 / 0.5` | 0.464916004312 |
| z | `9.9999999999999995e-7` | 1 |

The source's `ScalingX*zscale` equals `ScalingZ`; no further z correction
is applied. This acquisition-specific conversion is selected only for
`trackedNodes_nuc_48_2_OA7_1.mat`. Other filenames, spreadsheets, and generated
geometry retain `coordinateInputMode`. To add a different acquisition, supply
its own calibration in `src/experimental_coordinate_units.m`.

The calibrated local data contain 357 valid cells and 1,235 links. Fitted
full ellipsoid axes are approximately 492.42 × 177.69 × 177.83 µm, and median
link length is 26.35 µm. No 500 µm length constraint from the analysis script
is imposed on the simulation fit. The old uniform-micron assumption produced
1,050 × 383 × 180 µm and 1,249 links; anisotropic scaling changes neighbors
and may also change rupture-patch membership and dynamics.

Saved `param.units.coordinateInputMode` is `voxels`,
`param.units.coordinateToMeters` and `param.experimentalCoordinateScale`
contain the row vector of xyz spacings, and `param.units.calibrationSource`
records the script/lines. The units schema version is 2 for these calibrated
results; physical outputs remain in meters, seconds, and newtons.

Time calibration is separate: `Step15_3Danalysis_v3.m:308–312` uses supplied
`tvec` values or defaults to one time unit per frame. It does not establish
seconds per frame. The acquisition duration was subsequently confirmed as
**30 minutes across 45 frames**. Tuning assumes uniform sampling from 0 to
1,800 s: **1,800 / 44 = 40.9090909 s between frames**, not 1,800 / 45.
The simulation's existing default time parameters are unchanged.
The `input_node_data/tracked_nodes/` folder remains ignored by Git.

## Parameter fitting against tracked cells

To also fit **static cell connections and the weak rupture patch**, run:

```matlab
fit = tune_network_simulation(160); % coarse candidate budget; resumable
fitFile = fullfile('results','network_tuning','best_fit.mat');
summarize_tuning(fitFile);
verify_tuning(fitFile);
simFile = run_tuned_simulation(fitFile);
% Extend a completed/interrupted search with a larger total budget:
% fit = tune_network_simulation(240, true);
```

To regenerate the complete output set from the latest saved fit, run:

```matlab
outputs = regenerate_tuned_results;
% Or select a specific saved fit:
% outputs = regenerate_tuned_results('results/network_tuning/best_fit.mat');
```

This reruns the forward simulation with all fitted coefficients, connections,
and weak links, then exports snapshots, ring displacement/velocity, rupture
and retraction forces, and the motion MP4. Network-fit outputs use the
`_network_tuned` suffix in `results/`; earlier parameter-only outputs retain
`_tuned`. It also refreshes the fit comparison plots/report, checks tighter
solver tolerances and force balance, and decodes the video to verify its
frame count and coverage through 1,800 s. The corresponding
`visualization_info_network_tuned.mat` records source paths and plot/video
metadata.

This search starts from the existing parameter-only fit (or the checked-in
preset) and alternates connection density/orientation, rupture-patch
location/size, the seven mechanical coefficients, and explicit local edge
swaps. Candidate connections come from the symmetric 12-neighbor graph of
the **initial measured positions**, plus original links. A minimum spanning
tree preserves connectivity in generated graphs; edge swaps preserve a
minimum node degree of three. Global neighbor counts span 4, 6, 8, and 10;
axial distance weights span 0.5, 1, and 2. Patch proposals cover the
nonanchored surface and then explore nearby centers, with 1–32 connected
weak links. These are search assumptions, not measured cell contacts.

The graph and weak patch stay fixed throughout each forward simulation.
Rest lengths come from initial positions. Weak links fail through the
existing extension-threshold event and damage law; the rupture time is not
a fitted clock. No node position, velocity, motion schedule, or trajectory
feedback is introduced. Seven global mechanical scalars remain the only
fitted actuation/material coefficients. The original training/validation
split and rigid registration stay fixed; **only training 3D position error**
selects graphs, patches, and coefficients. The five best coarse candidates
are checked at fine solver accuracy against the previous parameter-only
fit, which remains eligible as the incumbent.

Outputs live separately in `results/network_tuning/`: the full fit and
checkpoint, `connections.csv` with original cell IDs and weak-link labels,
mechanical parameters, cell errors, trajectory and network comparison
figures, and a report. The independent forward result is
`results/follicleShellSim_network_tuned.mat`. Use `best_fit.mat` for network
forward runs: the scalar parameter CSV alone does not encode the fitted
connections or weak patch. The tuned-run, summary, comparison, and verification
commands prefer the network fit when it exists; pass the parameter-only
`best_fit.mat` explicitly to use the earlier fit. Existing parameter-only
results and `run_simulation` remain unchanged. The inferred topology is a
bounded, heuristic model fit; improved trajectory agreement does not
uniquely identify biological contacts or the biological rupture site.

Run the bounded parameter search from the repository root in MATLAB:

```matlab
fit = tune_simulation(48);      % maximum number of coarse candidate solves
fitFile = fullfile('results','parameter_tuning','best_fit.mat');
summarize_tuning(fitFile);      % trajectory comparison PNG and numerical report
simFile = run_tuned_simulation(fitFile); % independent forward run
verify_tuning(fitFile);         % tighter integration and force-balance checks
```

The fitter uses the calibrated `cell-scale` model. It reads an existing
`results/follicleShellSim_cell_scale.mat` baseline, generating one if absent.
It preserves the initial geometry, cell identities, mesh, spring constants,
mass assumption, and shared force equations. `run_simulation.m` and its
default coefficients remain unchanged; fitted runs are saved separately as
`results/follicleShellSim_tuned.mat`.
The checked-in `tuned_parameters.csv` records the first fitted preset for
this acquisition. `run_tuned_simulation` prefers the local network fit,
then the parameter-only fit, and falls back to that preset otherwise. To run the recorded
preset explicitly, use `run_tuned_simulation('tuned_parameters.csv')`.

**No measured displacement is prescribed in the simulation.** The forward
solver takes only model parameters, initial geometry, and output times.
Node velocities emerge from mechanical force balance in the overdamped
mass-spring-damper system. The existing global natural-length contraction
law supplies active loading. Its onset and rate are fitted scalars; there
are no fitted per-cell motion schedules, target-position forces, or imposed
rupture times. Rupture still occurs at an extension-threshold crossing.

| Fitted coefficient | Search interval | Units |
| --- | ---: | --- |
| Global natural-length shortening rate | 2e-5–4e-4 | 1/s |
| Global contraction onset | 0–1200 | s |
| Maxwell relaxation time `c/k2` | 5–300 | s |
| Adhesion release-load multiplier | 0.15–8 | dimensionless |
| Adhesion release time | 2–180 | s |
| Weak-patch rupture strain | 0.01–0.15 | dimensionless |
| Global contraction rise time | 60–900 | s |

These are bounded modeling assumptions, not experimentally measured ranges.
The equilibrium and transient springs stay at 0.002 N/m and cell mass stays
at 5e-13 kg. The fitted relaxation time changes the edge dashpot. Absolute
force scale cannot be identified from these trajectories alone: multiplying
all forces and damping by a common factor preserves overdamped motion.
Consequently, fitted trajectories do not independently validate predicted
forces, spring constants, or cell mass.

The comparison includes the 357 cells present in frame 1; the other 27 rows
have no initial position and are excluded. A rigid rotation and translation
align each experimental frame's far-body cap with the initial cap. This
uses training cells only and never changes simulated positions, rescales the
specimen, or fits deformation. It assumes the cap's coherent motion is
specimen drift; real cap movement could also be removed. Both registered
and raw observations are retained in the fit artifact.

Every fifth cell in initial axial order is held out, giving 286 training
cells and 71 validation cells. Even-numbered frames supply the training
targets; odd frames after frame 1 supply 22 held-out frames. The objective
is mean squared **3D position error**, so the large late retraction is not
downweighted as an outlier. Only training error selects candidates. The
joint holdout score evaluates excluded cells at excluded times. Pose
correction at those times still uses observed training-cap cells, so this is
interpolation validation within one specimen, not forecasting or independent
biological validation. An unmoving initial configuration is scored alongside
the original model to expose fits that merely suppress all motion.

Search uses deterministic initial guesses, bounded coordinate steps, cached
evaluations, sparse stiff integration, and coarse solver tolerances before
rerunning the best three candidates at the standard fine tolerances. It uses
no Optimization Toolbox, parallel pool, or videos during fitting. A trial
time budget guards expensive candidates; timed-out trials are computationally
censored rather than evidence of bad physics. Checkpoints preserve completed
solves. To extend an unchanged search to a larger cumulative budget:

```matlab
fit = tune_simulation(72, true);
summarize_tuning;
```

The local `results/parameter_tuning/` folder contains `best_fit.mat`
(parameters, observations, predictions, split masks, errors, and history),
`parameters.csv`, `search_history.mat`, `trajectory_comparison.png`, and
`tuning_report.txt`. `cell_errors.csv` records errors by original tracked-cell
ID, and `verification.mat` records numerical convergence and force balance.
These generated outputs and the raw experimental MAT
file remain ignored by Git. The fit is a local search; it does not establish
parameter uniqueness, confidence intervals, or the adequacy of all model
mechanisms.

### Initial fit result

The recorded preset improves the original model substantially but is **not
a close match to the complete experimental motion**. After the documented
rigid pose correction, training / joint holdout 3D RMSE are approximately
15.40 / 19.75 µm, compared with 72.30 / 78.74 µm for the original cell-scale
model and 19.90 / 21.09 µm for an unmoving initial configuration. The holdout
improvement over no motion is only about 6.3%; the large improvement over
the original model alone would overstate the calibration's success.

| Final mean axial motion | Experiment, pose corrected | Fitted model |
| --- | ---: | ---: |
| Initial 0–10% band | 61.27 µm | 63.93 µm |
| Initial 25–35% band | 59.80 µm | 27.99 µm |
| Initial 45–55% band | 52.66 µm | 18.39 µm |

The tip's final displacement is close, while the later middle-band motion
is substantially underestimated. A targeted search over the existing global
activation rise time did not improve training error. Damping and release
coefficients were also weakly distinguished by the tested trajectories.
The preset is a provisional, reproducible starting point; it does not show
that parameter adjustment alone can reproduce all three bands. Its predicted
rupture event is about 694.17 s, which is an internal model event rather than
an experimentally established rupture timestamp.

The 79 recorded coarse evaluations include nine early timeouts that were
subsequently rerun successfully; the top three candidates were checked with
fine integration. All 15 self-contained tests pass. Tenfold tighter solver
tolerances and a halved maximum step change node positions by at most
0.000491 µm and the rupture time by 0.000075 s. The maximum checked net-force
residual is 2.03e-24 N. Numerical error is therefore negligible compared with
the experimental mismatch.

## Approximately concentric rupture variant

To reduce the one-sided clustering while retaining cell-to-cell variation:

```matlab
simFile = run_concentric_simulation;
```

This starts from the saved tuned parameters and initial coordinates, runs a
new forward solve, and writes `results/follicleShellSim_concentric.mat`,
`concentric_comparison.png`, `concentric_comparison.csv`, and
`follicle_rupture_motion_concentric.mp4`. The comparison includes an end-on
view, angular drift, angular concentration, and axial retraction. Original
and tuned runs are preserved. If no tuned result exists, the runner first
creates one using `run_tuned_simulation`.

The original 22 weak links were selected by extreme axial position. In this
specimen their midpoint azimuths occupy only about 91 degrees, leaving most
of the circumference intact. The concentric variant instead defines a cap
by ellipsoid-normalized polar angle (the closest 10% of nodes to the pole)
and weakens predominantly circumferential links throughout that cap. It
retains meridional links to transmit retraction and rebuilds the contraction
mask so weak links cannot shorten. The same strain-threshold rupture event
and damage cascade still determine when the patch opens.

A second, explicit modeling assumption supplies **soft circumferential
substrate support**. Each cell's initial radial direction defines a
half-plane through the long axis. A finite restoring force penalizes motion
away from that half-plane, projected onto the local shell tangent. Axial
and radial motion within the half-plane remain free; no equal cell speeds,
prescribed ring radii, or measured target positions are imposed. Its default
stiffness is 0.5 times mean equilibrium-link stiffness. This additional
support is not an experimentally established material property.

```matlab
% Arguments: source simulation, support/equilibrium stiffness ratio, render video
run_concentric_simulation('results/follicleShellSim_tuned.mat', 0.5, false);
% Use ratio 0 to isolate the effect of the new weak-region geometry.
```

The soft support is needed because broadening rupture alone still permits
circumferential collapse. Breaking *all* cap links was rejected because it
leaves some tip nodes mechanically disconnected. The selected variant
checks that every node retains an intact connection. Natural geometric
irregularity and local speed differences remain, including some angular
outliers; this is not an enforced circular boundary.

This variant changes the mechanics and has **not been refitted** to the
experimental trajectories. Reduced lateral clustering also reduces axial
retraction, so the earlier fit scores and rupture time do not apply to it.
The original runners omit the new support and retain their prior equations.

## Spreadsheet Node Data

To use spreadsheet coordinates, place the file in `input_node_data/` and set:

```matlab
useExperimentalNodeData = true;
experimentalNodeDataFileName = 'nodeDataSample1.xlsx';
```

The first four spreadsheet columns must contain:

```text
Cell ID, x, y, z
```

The coordinate columns follow `coordinateInputMode` and
`cellDiameter_um`.

## Model State

The ODE state is:

- `q`: node displacement
- `v`: node velocity
- `sigma`, `zeta`: viscoelastic edge states
- `D`: edge damage from `0` (intact) to `1` (failed)
- `A`: irreversible node-adhesion release from `0` (attached) to `1` (detached)

The selected weak edge patch lies at the extreme of the principal-axis 0-10%
tip. Predominantly axial, intact edges from the tip toward the 90th-percentile
far-body cap begin contracting together with one smooth global ramp; there is no
spatial activation schedule. Zener linkages retain elastic stress during this
slow pull, while elastic reference-position adhesion tethers suppress bulk
motion. Active traction is disabled. When one weak edge reaches its extension
threshold, the selected local patch opens and unlocks irreversible,
load-dependent adhesion release. The lower tip threshold and load redistribution
then produce bodyward release through the mechanically coupled network. The
far-body cap remains softly anchored. The rupture edges do not contract. Both
endpoints of a
rupture edge must belong to the initial 0-10% node band; among those candidates,
the 22 edges nearest the tip are selected. Other edges are effectively
unbreakable.

The editable controls are grouped in the `USER PARAMETERS` block at the top
of `src/follicle_rupture_msd_fixed.m`. Edge- and node-specific arrays are
derived from these scalar settings later in the function:

```matlab
param.prestressStartTime = 3;         % seconds
param.anchorPercentile = 90;          % softly anchored far-body cap
param.contractionMinAxialAlignment = 0.35;
param.activationRiseTime = 900;      % global slow-pull onset, seconds
param.contractionRateFrac = 4e-4;     % plateau rate, original lengths/second
param.minRestLengthFrac = 0.50;       % maximum 50% natural-length shortening
param.activeTraction_on = false;      % recoil uses stored elastic energy
param.adhesionStiffness = 1;          % reference-position tether, N/m
param.adhesionLoadThresholdTip = 5e-7;    % lower tip release load, N
param.adhesionLoadThresholdBody = 3e-6;   % higher body release load, N
param.adhesionTransitionWidthFrac = 0.25; % smooth force-transition width
param.adhesionReleaseTau = 5;              % irreversible release time, seconds
param.adhesionDragStuck = 0.05;            % attached damping, N*s/m
param.adhesionDragSliding = 0.01;          % detached drag, N*s/m
param.extensionCritFrac = 0.05;      % rupture-patch extension threshold
param.rupturePatchCascade_on = true; % first break opens the selected patch
```

These are the legacy coefficients before the optional common force scaling.
The 900 s rise, 1,250 s equivalent full-rate shortening duration, and smooth
ramp down reach the 50% natural-length floor at 2,153 s, after the default
1,800 s run ends. Natural-length shortening is not the same as tissue strain.

Saved output includes `Lnatural_time`, the natural length of every edge at
every saved time, in addition to displacement `Y`, velocity `V`, edge damage
`D_time`, and node adhesion release `A_time`. After simulation,
`run_simulation` opens a node-position snapshot
figure and a second figure showing signed long-axis displacement and velocity
for the initial 0-10%, 25-35%, and 45-55% principal-axis rings. Positive values
indicate motion away from the rupture tip.

## Rupture and Retraction Forces

The separate force figure reports two model-derived measurements:

- **Rupture-patch tensile load:** the sum of positive axial link forces over
  the selected weak edges, including their current damage. The marked rupture
  force estimate is the load at the last saved intact sample before any patch
  edge exceeds damage `1e-6`. The title reports the interval between that sample
  and the first damaged sample; finer saved output improves time resolution.
  This is a sum of link tensions, not a vector resultant across a tissue
  section. Rupture is still triggered by extension, not a force threshold.
- **Retraction force:** the mean internal link force per node projected onto
  the initial long axis for each initial 0–10%, 25–35%, and 45–55% band.
  Positive means away from the rupture tip. Link forces include the model's
  surface projection; external adhesion, drag, anchors, surface constraints,
  and active traction are excluded. This measures tissue pull, not the net
  force after resistance. The entire loading/retraction history is shown.

Both panels display µN. No rupture force is reported if the patch never
ruptures or is already damaged in the first saved sample. Empty percentile
bands return `NaN` rather than a zero force.

The simulation saves signed `edgeForce_time` (`Nt` by `Ne`),
`nodeInternalForce_time`, and `nodeNetForce_time` (both `Nt` by `3N`, in
node-wise x/y/z order), all in newtons. These are evaluated with the same
mechanical equations and full viscoelastic states used by the solver; the net
force includes all external loads and resistance. Older MAT-files must be
regenerated because they do not contain the force histories or the internal
viscoelastic states needed to reconstruct them.

To inspect a newly generated result without rerunning the simulation:

```matlab
[forceFigure, forceInfo] = plot_rupture_retraction_force( ...
    'results/follicleShellSim_test2.mat');
```

## Physical scale and inertia audit

The original wrapper multiplies physical acceleration `F/m` by
`epsilon = m*(k1+k2)/c^2`, in addition to the ordinary time/length unit
conversion. Thus the simulated inertial mass is `m/epsilon = c^2/(k1+k2)`:
**2,500 kg per node** at the original settings. Editing `mVal` alone cannot
correct this, because mass cancels. This behavior remains only in `legacy`
for exact reproducibility; it is not a physical cell mass.

The `cell-scale` profile estimates wet mass from density × volume:
`1000 kg/m^3 × 500 µm^3 = 5e-13 kg = 0.5 ng`. Volume and density are explicit
assumptions, independent of `cellDiameter_um`, which remains a coordinate and
display conversion. A tracked nuclear coordinate does not measure cell volume.

| Quantity | Legacy | Cell-scale | Meaning |
| --- | ---: | ---: | --- |
| Wet mass per node | nominal 1 kg | 5e-13 kg | one cell; estimated |
| Equilibrium / transient spring | 2 / 2 N/m | 0.002 / 0.002 N/m | 2 nN/µm per branch |
| Maxwell / Burgers dashpot | 100 / 100 N·s/m | 0.1 / 0.1 N·s/m | effective tissue dissipation |
| Maxwell relaxation `c/k2` | 50 s | 50 s | unchanged |
| Adhesion stiffness | 1 N/m | 0.001 N/m | effective tether |
| Adhesion release load, tip/body | 500 / 3000 nN | 0.5 / 3 nN | assumed release law |
| Background drag | 0.1 N·s/m | 1e-4 N·s/m | tissue resistance, not water drag |
| Attached / detached adhesion drag | 0.05 / 0.01 N·s/m | 5e-5 / 1e-5 N·s/m | effective friction |
| Anchor stiffness / damping | 10 N/m / 1 N·s/m | 0.01 N/m / 0.001 N·s/m | boundary tether |
| Surface penalty coefficient | 5e-6 N | 5e-9 N | multiplies dimensionless ellipsoid error |
| Normal damping | 1 N·s/m | 0.001 N·s/m | effective surface resistance |

All dimensional forces, stiffnesses, dashpots, load thresholds and load
regularization are multiplied by the same `cellForceScale = 1e-3`, including
the disabled traction amplitude. This preserves relative force strengths,
the contraction program, rupture strain, and all `c/k` times. It does **not**
preserve inertial transients when replacing the legacy inertia.

At the estimated mass, `m/c_drag = 5e-9 s`. The physical profile therefore
uses the overdamped limit `F_net = 0`: node velocities are solved from
instantaneous force balance, then `ode15s` integrates positions and material
states. The default Zener model needs only local normal/tangential drag solves;
pure dampers and Kelvin-Voigt links use a coupled sparse solve. Saved velocities
are the force-balance velocities, not finite differences. Actual mass and the inertial relaxation
estimate remain in `param`; `effectiveInertialMass = 0` explicitly records
the approximation. Solver position states are absolute dimensionless
coordinates to avoid roundoff during initial loading; saved `Y` remains
displacement in meters. Incomplete physical runs raise an error before saving.
With the patch cascade enabled, the first extension crossing is located as a
terminal event, then integration restarts with irreversible patch failure.
This keeps finite-difference Jacobian probes from activating the intact patch;
the threshold and 0.5 s damage time remain unchanged. `param.ruptureTime`
records that event (infinity if none); the force plot still uses its documented
saved-sample bracket, for consistency with legacy results.

### Evidence and limits of the estimates

Actomyosin-driven contraction and posterior wall breakdown are supported for
Drosophila ovulation by [Cho et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11441566/).
That study supports the mechanism, not a numerical spring constant or a
calibrated rupture force for this mesh. The current extension-triggered patch
cascade is a phenomenological substitute for local biochemical weakening.

As a cross-system order-of-magnitude comparison, the junction model in
[Huebner et al. (2021), Appendix 1 §6](https://elifesciences.org/articles/65390)
discusses stiffnesses of 0.1–1 nN/µm and effective damping around
100 nN·s/µm. Our 2 nN/µm spring and 100 nN·s/µm dashpot are close to those
scales, but that paper studies Xenopus, not mature Drosophila follicle cells.
These are plausible effective coefficients, **not measured follicle values**.

For a dimensional check only, `k = E*A/L` with assumed `L = 10 µm` and
`A = 10 µm^2` gives `E ≈ 2 kPa` per spring branch, while
`c = eta*A/L` gives `eta ≈ 1e5 Pa·s`. Neither the cross-section nor these
inferred moduli have been measured here. In particular, a link between cell
centers is not a single cadherin bond. [Crest et al. (2017)](https://elifesciences.org/articles/24958)
measured roughly 30–70 kPa in the **basement membrane** of developing
Drosophila follicles; those measurements should not be substituted directly
for mature cell-link stiffness.

At the calibrated median link length, 5% extension is about 1.32 µm and the
cell-scale equilibrium spring contributes about 2.63 nN before damage. The transient
branch depends on loading history. The saved patch sum is not a tissue
section resultant or an experimentally established rupture threshold.

The model still lacks explicit oocyte pressure/contact, cell volumes and
shape changes, junction remodeling, and a biochemical degradation law.
Surface projection and the far-body anchor represent external reactions.
The generated `[50 15 15]` geometry is a legacy demonstration, not a
calibrated mature follicle. Confirm stage and acquisition timing, then fit
forces and recoil times to independent measurements before claiming
quantitative biological prediction. The spatial conversion now follows the
supplied Step14 analysis script.

## Verification

Run the self-contained tests (no tracked data required):

```matlab
assertSuccess(runtests('tests'));
```

The cleanup was compared against a fresh execution of the original source
on the local 357-cell, 1,249-link input: all 12 saved numerical arrays were
bit-for-bit identical. Alternate material models, rupture modes, input modes,
and saved metadata were retained because they are supported functionality or
provenance, rather than dead parameters. Physical-profile results must be
compared separately because correcting inertia changes the dynamics.

### Source-calibration verification

All 12 tests pass (nine mechanics tests plus three calibration tests). Both
mechanics profiles were rerun through 1,800 s using the source xyz spacings.
Geometry checks verify that conversion occurs before rotation and mesh
construction, and that the experimental snapshot affine transform includes
the same unequal xyz spacings. The old `input-units` option reproduces all
12 pre-calibration numerical arrays bit-for-bit. Original outputs are archived
under `results/before_source_calibration/`.

Calibrated runs produced the following model outputs, without retuning
mechanics or timing parameters:

| Measurement | Legacy mechanics | Cell-scale mechanics |
| --- | ---: | ---: |
| First sampled rupture interval | 499–500 s | 497–498 s |
| Last intact patch tensile load | 19.500 µN | 19.377 nN |
| Final mean axial displacement, 0–10% band | 266.508 µm | 386.433 µm |
| Final mean axial displacement, 25–35% band | 151.818 µm | 217.006 µm |
| Final mean axial displacement, 45–55% band | 103.343 µm | 149.365 µm |

The cell-scale rupture event is at approximately 497.264 s and its maximum
saved net-force residual is 4.65e-24 N. Both profiles' figures, the comparison
figure, and both 240-frame videos were regenerated; video read-back checks
passed. These values are model predictions, not experimental measurements.

### Historical audit before spatial calibration

The earlier local audit (1,800 s, original coordinate assumption) gave the
following values. These are retained as historical evidence, not as the
expected results of the newly calibrated default:

| Measurement | Legacy | Cell-scale |
| --- | ---: | ---: |
| First sampled rupture interval | 461–462 s | 461–462 s |
| Last intact patch tensile load | 27.659 µN | 27.747 nN |
| Final mean axial displacement, 0–10% band | 137.434 µm | 196.886 µm |
| Final mean axial displacement, 25–35% band | 87.790 µm | 120.189 µm |
| Final mean axial displacement, 45–55% band | 54.514 µm | 78.287 µm |

The physical event occurred at 461.284539 s. The ordering and direction of
ring retraction are retained, but removing the artificial inertia increases
final ring displacement by approximately 37–44% and changes recoil velocities.
The new profile should therefore not be treated as numerically interchangeable
with the old one. That historical comparison is archived at
`results/before_source_calibration/mechanics_profile_comparison.png`.
The current `results/mechanics_profile_comparison.png` compares both mechanics
profiles on calibrated geometry; physical figures and the full 240-frame video
have `_cell_scale` suffixes.

The nine mechanics tests in that audit passed, including force balance for all six material
models, an analytical retraction solution, and the exponential patch-failure
law after the rupture event. All three legacy plot information structures
matched the original code, and an exported MP4 was decoded successfully.
Tightening `RelTol`/`AbsTol` from `1e-5`/`1e-8` to `1e-6`/`1e-9` and halving
`MaxStep` changed physical rupture time by 0.000084 s and saved displacement
by at most 0.000246 µm. The maximum default-run net-force residual was
1.32e-23 N. These checks establish numerical consistency, not biological
calibration.
