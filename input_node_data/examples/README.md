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
