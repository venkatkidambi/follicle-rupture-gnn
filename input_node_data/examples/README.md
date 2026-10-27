# Input Node Data Examples

Place experimental node spreadsheets in `input_node_data/`.

The simulation expects the first four columns to be:

1. `Cell ID`
2. `x`
3. `y`
4. `z`

Real spreadsheets in `input_node_data/` are ignored by git.

Tracked MATLAB data may also be placed in `input_node_data/tracked_nodes/`.
The current model expects a numeric `trackedNodes` array shaped `N x 3 x T`
and uses frame 1 only. All-zero frame rows are treated as absent-track
placeholders and omitted.

The supplied `trackedNodes_nuc_48_2_OA7_1.mat` dataset uses raw image
coordinates. By default, its x/y spacing is 0.464916004312 µm and z spacing
is 1 µm, as specified in `Step14_findingAxis_v2.m`. The conversion occurs
before centering and orientation. Other files use `coordinateInputMode`;
this acquisition's calibration is not automatically applied to them.
