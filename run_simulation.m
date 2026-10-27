repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
mechanicsProfile = 'legacy'; % original mechanics; 'cell-scale' uses cell mechanics
coordinateCalibration = 'source'; % Step14 xyz spacing; 'input-units' restores old scaling
simFile = follicle_rupture_msd_fixed(mechanicsProfile, [], coordinateCalibration);
outputSuffix = '';
if strcmp(mechanicsProfile, 'cell-scale')
    outputSuffix = '_cell_scale';
end
if strcmp(coordinateCalibration, 'input-units')
    outputSuffix = [outputSuffix '_input_units'];
end

[snapshotFigure, ~] = plot_node_position_snapshots(simFile);
[velocityFigure, ~] = plot_percentile_ring_velocity(simFile);
[forceFigure, forceInfo] = plot_rupture_retraction_force(simFile);

snapshotFigure.Units = 'inches';
snapshotFigure.Position(3:4) = [16, 6.5];
exportgraphics(snapshotFigure, ...
    fullfile(repoRoot, 'results', ['node_position_snapshots' outputSuffix '.png']), ...
    'Resolution', 300);
exportgraphics(velocityFigure, ...
    fullfile(repoRoot, 'results', ['percentile_ring_velocity' outputSuffix '.png']), ...
    'Resolution', 300);
exportgraphics(forceFigure, ...
    fullfile(repoRoot, 'results', ['rupture_retraction_force' outputSuffix '.png']), ...
    'Resolution', 300);
save(fullfile(repoRoot, 'results', ['rupture_retraction_force' outputSuffix '.mat']), 'forceInfo');

export_follicle_rupture_motion_video(simFile, ...
    fullfile(repoRoot, 'results', ['follicle_rupture_motion' outputSuffix '.mp4']));
