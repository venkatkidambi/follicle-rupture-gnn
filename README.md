# Follicle Rupture Shell Simulation

This repository contains a stripped-down MATLAB simulation of follicle rupture
on an ellipsoidal epithelial shell.

The current model is intentionally bare bones:

- surface nodes are connected by passive viscoelastic edges
- selected rupture edges accumulate damage
- damage weakens edge stiffness and damping until an edge can fail
- there is no active actomyosin force
- there is no active rest-length shortening
- there is no adaptive neighbor-based recruitment or shortening spread

This makes the codebase a baseline rupture-and-weakening model that can be
extended later with a new shortening mechanism.

## Quick Start

Open MATLAB in the repository root and run:

```matlab
run_simulation
```

Outputs are written to `results/`:

- `follicleShellSim_test2.mat`: saved simulation data
- `follicleShellSim_test2_dashboard.png`: rupture/motion diagnostic dashboard
- `follicleRupture_motion.mp4`: optional motion movie

## Files

- `run_simulation.m`: entry point that adds `src/` to the MATLAB path.
- `src/follicle_rupture_msd_fixed.m`: main model, parameters, ODE solve, and
  save step.
- `src/confirm_actomyosin_working.m`: dashboard generator. The filename is
  retained for compatibility, but the diagnostics now focus on rupture damage
  and motion rather than actomyosin recruitment.
- `src/export_follicle_rupture_motion_video.m`: MP4 renderer for a saved `.mat`
  simulation.
- `input_node_data/`: optional experimental node coordinate spreadsheets.
- `results/`: generated outputs.

## Model State

The ODE state is:

- `q`: node displacement
- `v`: node velocity
- `sigma`, `zeta`: viscoelastic edge states
- `D`: edge damage, where `0` means intact and `1` means failed

There is no `A_edge` recruitment state in this stripped version.

Each edge keeps its original geometric rest length:

```matlab
L0_ref = param.L0_edge0;
stretch = L - L0_ref;
stretch_dot = geometricStretchDot;
```

Damage weakens the passive connection element:

```matlab
g  = max(1 - D, param.g_floor_stiff);
gd = max(1 - D, param.g_floor_damp);
```

`g` scales stiffness. `gd` scales damping, with a small floor to avoid
numerical singularities.

## Main Parameters To Edit

Most settings live near the top of `src/follicle_rupture_msd_fixed.m`.

```matlab
N = 300;                         % number of surface nodes
tEnd = 1800;                     % simulated time, seconds
numSave = 300;                   % number of saved frames
runPostprocess = true;           % make dashboard after simulation
makeRuptureMotionVideo = true;   % export MP4 after simulation
odeSolver = 'ode45';             % 'ode45', 'ode15s', or 'ode23t'
```

For a fast test run:

```matlab
N = 100;
tEnd = 60;
numSave = 60;
runPostprocess = false;
makeRuptureMotionVideo = false;
```

## Rupture And Weakening

The rupture region is selected in the damage section:

```matlab
param.damageMode = 'tipmost-edge'; % 'seam', 'tip', or 'tipmost-edge'
```

The default model cuts the selected damage edge immediately:

```matlab
param.instantDamage_on   = true;
param.instantDamageTime  = 0;
param.instantDamageTau   = 1e-4;
```

Tension-gated damage remains available:

```matlab
param.tipDamage_on = true;
param.sig_fail_tip = 1e-9;
param.tau_dmg      = 0.005;
param.dmg_pow      = 2;
```

Only edges in `param.damageMask` are eligible for localized rupture damage.

## Experimental Node Data

Experimental node data support is still diagnostic. Set:

```matlab
useExperimentalNodeData = true;
```

Expected spreadsheet layout:

```text
Cell ID, x, y, z
```

The loader centers and orients coordinates, fits an ellipsoid for projection,
and builds a local nearest-neighbor edge graph when `experimentalEdgeMode` is
`'knn'`.

## Notes For Future Shortening Work

This version deliberately removes the old adaptive neighbor-based shortening
path. A future shortening mechanism can be added back against a clean baseline
by introducing a new rest-length law and, if needed, a new state variable.
