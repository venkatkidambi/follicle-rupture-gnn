# Follicle Rupture Shell Simulation

This repository contains a MATLAB simulation of follicle rupture driven by
actomyosin-like contraction on an ellipsoidal epithelial shell. It is written as
a small research codebase: most settings live near the top of one main file, and
the simulation output is saved as a `.mat` file that can be plotted or turned
into a movie.

The current preferred model uses **active rest-length shortening**: selected
cell-cell junctions contract because their preferred length becomes shorter.
The older **active tension** model is still available for comparison, but it is
off by default.

## Quick Start

1. Open MATLAB in the repository root folder.
2. Run:

   ```matlab
   run_simulation
   ```

3. Look in `results/` for the outputs:

   - `follicleShellSim_test2.mat`: saved simulation data
   - `follicleShellSim_test2_dashboard.png`: diagnostic dashboard, if
     postprocessing is enabled
   - `follicleRupture_motion.mp4`: motion movie, if video export is enabled

The wrapper `run_simulation.m` adds `src/` to the MATLAB path and then runs the
main simulation file:

```matlab
src/follicle_rupture_msd_fixed.m
```

## Files

- `run_simulation.m`: easiest entry point.
- `src/follicle_rupture_msd_fixed.m`: main model, parameters, ODE solve, and
  save step.
- `src/confirm_actomyosin_working.m`: postprocessing dashboard for checking
  motion, recruitment timing, and active shortening.
- `src/export_follicle_rupture_motion_video.m`: standalone MP4 renderer for a
  saved `.mat` simulation.
- `results/`: generated outputs. This folder can be deleted and regenerated.

## What The Simulation Represents

The follicle is represented as points on an ellipsoid surface. Nearby points are
connected by edges, which act like cell-cell junctions. Each edge can stretch,
relax, shorten, and fail.

The state variables are:

- `q`: node displacement
- `v`: node velocity
- `sigma`, `zeta`: internal viscoelastic edge states
- `D`: edge damage, where `0` means intact and `1` means failed
- `A_edge`: active recruitment level, where `0` means inactive and `1` means
  fully recruited

The code builds a surface mesh, chooses posterior/tip regions, defines active
and damage masks, integrates the equations of motion, and saves the node motion
over time.

## Main Parameters To Edit

Most user-facing settings are near the top of
`src/follicle_rupture_msd_fixed.m` in the `USER-PARAMETERS` section.

Common settings:

```matlab
N = 300;                         % number of surface nodes
tEnd = 900;                      % simulated time, seconds
numSave = 300;                   % number of saved frames
runPostprocess = true;           % make dashboard after simulation
makeRuptureMotionVideo = true;   % export MP4 after simulation
odeSolver = 'ode45';             % 'ode45', 'ode15s', or 'ode23t'
```

For a fast test run, use:

```matlab
N = 100;
tEnd = 60;
numSave = 60;
runPostprocess = false;
makeRuptureMotionVideo = false;
```

For a nicer final run, increase `N`, `tEnd`, and `numSave`, then turn the
dashboard or movie back on.

## Preferred Mode: Rest-Length Shortening

This is the recommended mode.

Instead of adding an extra pulling force to an edge, the model shortens the
edge's preferred length:

```matlab
param.act_on = false;
param.activeRestLength_on = true;
```

Conceptually, this means:

1. An edge has an original rest length, `L0_edge0`.
2. Active actomyosin recruitment increases `A_edge`.
3. The active rest length becomes shorter:

   ```matlab
   L0_active = L0_edge0 * (1 - activeFrac)
   ```

4. The passive viscoelastic edge model reacts to the mismatch between the
   current length and this shorter preferred length.

Important parameters:

```matlab
param.activeShorteningMax = 0.32;       % max fractional shortening
param.activeShorteningTau = 1;          % shortening timescale, seconds
param.activeShorteningActivation = 'exponential';
param.activeShorteningMask = param.motorMask;
```

`activeShorteningMax = 0.32` means an active edge can try to become up to 32%
shorter than its original length. The mask limits shortening to the motor
region, so not every edge contracts.

The shortening activation can be:

- `'exponential'`: turns on quickly and asymptotically approaches full
  shortening.
- `'sigmoid'`: smooth switch-like onset.
- `'accelerating'`: starts slowly and speeds up during the ramp.

### Optional Prescribed Rest-Length Wave

The code also contains a prescribed wave for rest-length shortening:

```matlab
param.activeRestLengthWave_on = false;
```

When enabled, different edges begin shortening at different prescribed times.
This can be useful for controlled experiments where you want a known posterior
to anterior activation pattern. In the current default configuration, this is
off because local recruitment controls the activation spread.

## Recruitment Adaptation

Recruitment is the adaptive part of the preferred model. Instead of every motor
edge turning on at the same time, each edge has its own activation state,
`A_edge`.

Default:

```matlab
param.activeRecruitment_on = true;
param.seedMask = param.tipMask;
param.recruitmentMask = param.motorMask;
param.recruitmentSignalMode = 'neighborA_plus_strainRate';
```

The posterior tip seeds activation. Recruitment then spreads through the edge
network, but only inside the recruitable motor domain.

The ODE for `A_edge` has two main contributions:

- Seed activation: edges in `seedMask` autonomously turn on.
- Neighbor recruitment: inactive edges turn on when nearby edges are already
  active.

In the default signal mode, strain rate also helps recruitment:

```matlab
param.recruitmentSignalMode = 'neighborA_plus_strainRate';
```

This means the model asks two questions:

1. Are neighboring edges recruited?
2. Is this local region mechanically changing fast enough to strengthen the
   signal?

Useful recruitment parameters:

```matlab
param.k_seed = 1/60;                         % autonomous seed rate
param.k_neighbor = 1/90;                     % neighbor-only rate
param.k_strainRateRecruit = 1/120;           % neighbor + strain-rate rate
param.recruitmentThreshold = 0.03;           % suppress tiny neighbor leakage
param.recruitmentStrainRateThreshold = 1e-4; % ignore small strain-rate noise
param.recruitmentStrainRateScale = 1e-3;     % strain-rate saturation scale
param.recruitmentStrainRateGain = 2;         % strain-rate amplification
```

To make recruitment spread faster, increase `k_seed`,
`k_neighbor`, or `k_strainRateRecruit`. To make it harder for recruitment to
spread, increase `recruitmentThreshold` or
`recruitmentStrainRateThreshold`.

The saved file includes `A_time`, so you can inspect recruitment over time.
The dashboard plots recruitment phase delays between posterior-to-anterior
bands.

## Legacy Mode: Active Tension

The older force-based actomyosin model is still in the code for comparison.
It adds active tension or active traction directly instead of changing the
edge's preferred length.

To use legacy active tension, switch:

```matlab
param.act_on = true;
param.activeRestLength_on = false;
param.activeRecruitment_on = false;   % optional, if you want a purely legacy run
```

Main legacy parameters:

```matlab
param.Tact_main = 1e-7;                % main-body active tension/traction
param.Tact_tip = 1e-8;                 % weaker tip cohesion
param.activeForceMode = 'posterior-traction';
param.balanceActiveTraction = true;
param.activeStrokeLimit = 100e-6;
param.activeStrokeTau = 2e-6;
```

There are two force modes:

- `'posterior-traction'`: applies nodal traction along the posterior axis.
- `'edge-tension'`: adds tension along active edges.

The default legacy setting is `'posterior-traction'`. With
`balanceActiveTraction = true`, the code subtracts the mean active force so the
whole object is less likely to drift as a rigid body.

### Optional Active Tension Wave

Legacy active tension can also be applied as a wave:

```matlab
param.activeWave_on = true;
param.activeWaveMode = 'continuous';   % or 'rings'
param.activeWaveDelayMax = 240;
param.activeWavePulse_on = true;
param.activeWavePulseDuration = 10;
```

Use this when you want force-based actomyosin to turn on first near the
posterior side and later toward the anterior side.

## Passive Edge Models

The passive edge model controls how junctions respond to stretch and stretch
rate. The current default is:

```matlab
param.connectionModel = 'burgers';
```

Available models:

- `'spring'`: elastic force only.
- `'damper'`: viscous force only.
- `'maxwell'`: spring and damper in series.
- `'kelvin-voigt'`: spring and damper in parallel.
- `'zener'`: standard linear solid.
- `'burgers'`: Kelvin-Voigt branch in series with Maxwell branch.

You can also assign one model per edge using `param.connectionModel_edge` as a
cell array, string array, or categorical array.

## Damage And Rupture

Damage is stored in `D`. An edge with `D = 0` is intact; an edge with `D = 1`
is failed. Damage reduces passive stiffness and damping on that edge.

The current setup selects a small posterior damage region:

```matlab
param.damageMode = 'tipmost-edge';   % 'seam', 'tip', or 'tipmost-edge'
param.tipDamage_on = true;
param.sig_fail_tip = 1e-9;
param.tau_dmg = 0.005;
param.eps_gate_on = true;
param.eps_gate = 0.001;
```

`tipmost-edge` damages only the most posterior edge. Use `'tip'` to allow the
whole posterior cap to damage, or `'seam'` to damage a band near the posterior
selection cutoff.

Damage only accumulates when the edge is overloaded, and with `eps_gate_on`
enabled it must also be opening enough to pass the strain gate.

The saved file includes `D_time`, which is used by the movie renderer to show
broken edges.

## Performance Optimizations

Several choices in the code are there to make the simulation faster and more
stable.

### 1. Vectorized Edge Calculations

The right-hand side computes stretch, stretch rate, active shortening, damage,
and passive stress for all edges as vectors. This avoids a MATLAB loop over
every edge at every ODE step.

The passive connection laws are handled by:

```matlab
connectionElement1DVectorized(...)
```

Internally, each model type is selected with logical masks, so mixed per-edge
models can still be evaluated in batches.

### 2. Sparse Edge-To-Node Accumulation

The code precomputes sparse matrices:

```matlab
param.edgeAccumI
param.edgeAccumJ
```

These matrices accumulate edge forces onto nodes efficiently. This is faster
and cleaner than repeatedly looping through every edge and manually adding
forces to the two endpoint nodes.

The recruitment neighbor graph is also sparse:

```matlab
param.edgeNeighbor
param.edgeNeighborCount
```

This lets the code compute average neighboring activation with one sparse
matrix multiply.

### 3. Nondimensionalized ODE State

The solver integrates a nondimensional state rather than raw SI units. This
keeps state variables closer in numerical scale, which helps adaptive ODE
solvers choose more reasonable steps.

The physical right-hand side is still written in physical units, then converted
back in:

```matlab
nondimRHS_shell(...)
```

### 4. Adaptive ODE Solvers

The solver still adapts internally even though the saved output uses
`numSave` time points.

Available solvers:

```matlab
odeSolver = 'ode45';
odeSolver = 'ode15s';
odeSolver = 'ode23t';
```

Use `ode45` for general testing. Try `ode15s` or `ode23t` if the run becomes
stiff, especially when damage is very fast or viscoelastic parameters create
short timescales.

### 5. Saved Output Is Thinned

`numSave` controls how many frames are saved, not how many ODE steps the solver
takes:

```matlab
numSave = 300;
```

Lowering `numSave` reduces file size and speeds up postprocessing and movie
export, but it does not necessarily make the ODE solve much faster.

### 6. Optional Postprocessing And Video

Dashboard generation and MP4 export can take noticeable time. Turn them off
while tuning parameters:

```matlab
runPostprocess = false;
makeRuptureMotionVideo = false;
```

After you have a simulation you like, turn them back on or call the renderer
directly:

```matlab
export_follicle_rupture_motion_video( ...
    'results/follicleShellSim_test2.mat', ...
    'results/follicleRupture_motion.mp4');
```

The movie renderer also batches edge drawing into single `plot3` calls per edge
class, which is much faster than drawing each edge as a separate graphics
object.

## Troubleshooting

If the run is too slow:

- Lower `N`.
- Lower `tEnd`.
- Turn off `runPostprocess` and `makeRuptureMotionVideo`.
- Try fewer saved frames with `numSave`.
- Use `ode15s` or `ode23t` if the solver reports stiffness or takes tiny steps.

If the follicle flies away or drifts:

- Keep `param.balanceActiveTraction = true` in legacy active-tension mode.
- Increase `param.c_drag`.
- Keep `param.anchor_on = true`.

If nothing ruptures:

- Lower `param.sig_fail_tip`.
- Increase `param.activeShorteningMax`.
- Use a larger damage region such as `param.damageMode = 'tip'`.
- Check the dashboard to confirm recruitment and shortening are actually
  turning on.

If everything ruptures immediately:

- Increase `param.sig_fail_tip`.
- Increase `param.tau_dmg`.
- Lower `param.activeShorteningMax` or active tension.
- Use a smaller damage region such as `param.damageMode = 'tipmost-edge'`.

## Recommended Workflow

1. Start with the preferred rest-length shortening mode.
2. Use a small, fast run to check that the code runs.
3. Inspect the dashboard to confirm that recruitment and shortening look
   reasonable.
4. Increase resolution and duration for the final run.
5. Use the legacy active-tension mode only when you specifically want to compare
   against the older force-driven mechanism.

