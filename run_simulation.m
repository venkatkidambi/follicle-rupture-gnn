repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
follicle_rupture_msd_fixed;

simFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
[snapshotFigure, ~] = plot_node_position_snapshots(simFile);
[velocityFigure, ~] = plot_percentile_ring_velocity(simFile);
[forceFigure, forceInfo] = plot_rupture_retraction_force(simFile);

snapshotFigure.Units = 'inches';
snapshotFigure.Position(3:4) = [16, 6.5];
exportgraphics(snapshotFigure, ...
    fullfile(repoRoot, 'results', 'node_position_snapshots.png'), ...
    'Resolution', 300);
exportgraphics(velocityFigure, ...
    fullfile(repoRoot, 'results', 'percentile_ring_velocity.png'), ...
    'Resolution', 300);
exportgraphics(forceFigure, ...
    fullfile(repoRoot, 'results', 'rupture_retraction_force.png'), ...
    'Resolution', 300);
save(fullfile(repoRoot, 'results', 'rupture_retraction_force.mat'), 'forceInfo');

export_follicle_rupture_motion_video(simFile, ...
    fullfile(repoRoot, 'results', 'follicle_rupture_motion.mp4'));
