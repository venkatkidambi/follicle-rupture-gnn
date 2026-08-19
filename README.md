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

Post-processing includes shell snapshots and percentile-ring velocity plots.
Video generation and standalone tests are not included at this stage.

## Run

Open MATLAB in the repository root and run:

```matlab
run_simulation
```

The simulation writes `results/follicleShellSim_test2.mat`. Generated results
are ignored by git; `results/.gitkeep` preserves the output directory.

## Repository Layout

- `run_simulation.m`: adds `src/` to the MATLAB path and starts the model.
- `src/follicle_rupture_msd_fixed.m`: the complete simulation, including its
  local geometry, mechanics, unit, input, and progress helpers.
- `plot_node_position_snapshots.m`: initial, middle, and final shell positions.
- `plot_percentile_ring_velocity.m`: ring-average speed and rupture damage.
- `input_node_data/`: optional experimental node-coordinate spreadsheets.
- `results/`: local generated output.

## Main Parameters

Edit the user-parameter section near the top of
`src/follicle_rupture_msd_fixed.m`.

```matlab
cellDiameter_um = 10;
coordinateInputMode = 'cells'; % 'cells', 'microns', or 'meters'
displayMode = 'microns';       % 'cells' or 'microns'
generatedEllipsoidSemiAxes = [50 15 15];

N = 300;
tEnd = 1800;
numSave = 300;
ruptureSaveEnd = 120;          % dense post-processing window, seconds
ruptureSaveInterval = 0.25;    % output spacing in that window, seconds
odeSolver = 'ode45';           % 'ode45', 'ode15s', or 'ode23t'
```

All mechanics and saved physical data use SI units. Coordinate inputs are
converted to meters once when geometry is loaded.

## Tracked Node Data

The default run uses only frame 1 of the ignored local file:

```text
input_node_data/tracked_nodes/trackedNodes_nuc_48_2_OA7_1.mat
```

It loads `trackedNodes(:,:,1)` as the initial node cloud, interprets its
coordinates as microns, builds a six-nearest-neighbor edge mesh, and then
runs the same contraction-driven tip-patch rupture. Later tracked frames are
not used by the mechanics.

The `input_node_data/tracked_nodes/` folder is ignored by Git.

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

The selected weak edge patch lies at the extreme of the principal-axis 0-10%
tip. Its intact shoulder edges contract, loading the patch in tension and
pulling its rim nodes back along the shell after rupture. The rupture edges do
not contract, which avoids collapsing all of their endpoints into a cluster
before they fail. Both endpoints of a rupture edge must belong to the initial
0-10% node band; among those candidates, the eight edges nearest the tip are
selected. Other edges keep their geometric natural lengths and are effectively
unbreakable.

The main controls are in `src/follicle_rupture_msd_fixed.m`:

```matlab
param.contractionStartTime = 3;       % seconds
param.contractionRateFrac_edge(...) = 0.20; % original edge lengths/second
param.minRestLengthFrac_edge(...) = 0.30;
param.extensionCritFrac_edge(...) = 0.15; % rupture-patch extension threshold
```

Saved output includes `Lnatural_time`, the natural length of every edge at
every saved time, in addition to displacement `Y`, velocity `V`, and damage
`D_time`. After simulation, `run_simulation` opens a node-position snapshot
figure and a second figure showing mean velocity magnitude for the initial
0-10%, 25-35%, and 45-55% principal-axis rings.
