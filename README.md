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
cellDiameter_um = 10;          % positive cell diameter in microns
coordinateInputMode = 'cells'; % 'cells', 'microns', or 'meters'
displayMode = 'microns';       % 'cells' or 'microns'
generatedEllipsoidSemiAxes = [50 15 15];

N = 300;                         % number of surface nodes
tEnd = 1800;                     % simulated time, seconds
numSave = 300;                   % number of saved frames
runPostprocess = true;           % make dashboard after simulation
makeRuptureMotionVideo = true;   % export MP4 after simulation
odeSolver = 'ode45';             % 'ode45', 'ode15s', or 'ode23t'
```

### Cell-scale coordinate input

`coordinateInputMode` controls how both `generatedEllipsoidSemiAxes` and the
`x`, `y`, and `z` columns of an experimental spreadsheet are interpreted:

- `'cells'`: one coordinate unit is one cell diameter
- `'microns'`: one coordinate unit is 1 µm
- `'meters'`: coordinates are already in meters

Every coordinate array is converted to meters once, when it enters the
simulation. The generated default `[50 15 15]` in `'cells'` mode with a
10 µm cell diameter therefore produces the existing physical semi-axes of
500 µm × 150 µm × 150 µm.

For example, with 10 µm cells:

```matlab
cellDiameter_um = 10;
coordinateInputMode = 'cells';
generatedEllipsoidSemiAxes = [50 15 15]; % 500 x 150 x 150 µm
```

With 12 µm cells, 50 input cell diameters represent 600 µm:

```matlab
cellDiameter_um = 12;
coordinateInputMode = 'cells';
generatedEllipsoidSemiAxes = [50 15 15]; % 600 x 180 x 180 µm
```

To express the original 500 µm × 150 µm × 150 µm geometry directly in
microns instead:

```matlab
coordinateInputMode = 'microns';
generatedEllipsoidSemiAxes = [500 150 150];
```

`displayMode = 'cells'` plots position and displacement in cell diameters and
velocity in cells/s. `displayMode = 'microns'` plots them in µm and µm/s.
The selection affects only logging and plots, never the saved physical data.

### Unit convention

All mechanics use SI units. `nodeCoords`, displacement matrix `Y`, rest
lengths `L0_edge0`, and time vector `tSol` are saved in meters or seconds as
appropriate. Each MAT file stores the input, display, cell-diameter, and
internal-unit metadata in `param.units`.

| Parameter/data category | Internal and saved unit |
| --- | --- |
| Coordinates, positions, displacements, ellipsoid axes, edge/rest lengths, `zeta` | m |
| Time, damage timescales, saved `tSol` | s |
| Velocity | m/s |
| Force, edge force state `sigma`, failure threshold, surface-penalty force | N |
| Edge, anchor, and other linear stiffness | N/m |
| Edge, surface-normal, drag, and anchor damping | N·s/m |
| Node mass `mVal` | kg |
| Strain, damage, masks, percentiles, nondimensional solver state | dimensionless |
| Damage/priming rates | 1/s |

> **Warning:** changing the coordinate input unit is only a coordinate
> conversion. It is not valid to reinterpret lengths without also preserving
> the SI meanings of stiffness, damping, force, and mass. These mechanical
> parameters remain N/m, N·s/m, N, and kg and are never automatically rescaled.

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

The spreadsheet coordinate columns use the same `coordinateInputMode` and
`cellDiameter_um` configuration as generated coordinates. For example, a
spreadsheet value of `50` becomes 500 µm with 10 µm cells in `'cells'` mode,
while a value of `500` becomes 500 µm in `'microns'` mode.

The loader centers and orients coordinates, fits an ellipsoid for projection,
and builds a local nearest-neighbor edge graph when `experimentalEdgeMode` is
`'knn'`.

## Unit Validation

Run the lightweight conversion checks in MATLAB:

```matlab
addpath('src');
run('tests/validate_unit_conversions.m');
```

The checks cover 10 µm and 12 µm cell diameters, equivalent micron/meter
inputs, display conversion in both modes, invalid configurations, and the
legacy MAT-file plotting fallback.

## Notes For Future Shortening Work

This version deliberately removes the old adaptive neighbor-based shortening
path. A future shortening mechanism can be added back against a clean baseline
by introducing a new rest-length law and, if needed, a new state variable.
