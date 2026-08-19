repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'src'));
follicle_rupture_msd_fixed;

simFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
plot_node_position_snapshots(simFile);
plot_percentile_ring_velocity(simFile);
