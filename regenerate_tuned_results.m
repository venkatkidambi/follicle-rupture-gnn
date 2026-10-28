function outputs = regenerate_tuned_results(fitFile)
%REGENERATE_TUNED_RESULTS Simulate, render, and verify the latest saved fit.
% Uses the fitted graph and rupture patch as well as all mechanical scalars.
repo = fileparts(mfilename('fullpath'));
if nargin < 1
    fitFile = fullfile(repo, 'results', 'network_tuning', 'best_fit.mat');
    if ~isfile(fitFile)
        fitFile = fullfile(repo, 'results', 'parameter_tuning', 'best_fit.mat');
    end
end
loaded = load(fitFile, 'fit');
fit = loaded.fit;
outputs.fitFile = fitFile;
outputs.simFile = run_tuned_simulation(fitFile);
simulation = load(outputs.simFile, 'param', 'tSol', 'Y', 'nodeNetForce_time');
% Integration adds event metadata; every input field must still match the fit.
fields = fieldnames(fit.param);
for idx = 1:numel(fields)
    name = fields{idx};
    if strcmp(name, 'ruptureTime'), continue; end
    assert(isequaln(simulation.param.(name), fit.param.(name)), ...
        'Forward simulation changed fitted parameter: %s', name);
end
assert(simulation.tSol(1) == 0 && simulation.tSol(end) == 1800);
assert(all(isfinite(simulation.Y), 'all'));
outputs.maxNetForceResidual_N = max(vecnorm( ...
    reshape(simulation.nodeNetForce_time.', 3, []), 2, 1));
assert(outputs.maxNetForceResidual_N < 1e-18, 'Force balance failed.');
outputs.ruptureTime_s = simulation.param.ruptureTime;

[~, stem] = fileparts(outputs.simFile);
suffix = erase(stem, 'follicleShellSim');
outDir = fullfile(repo, 'results');
outputs.snapshot = fullfile(outDir, ['node_position_snapshots' suffix '.png']);
outputs.velocity = fullfile(outDir, ['percentile_ring_velocity' suffix '.png']);
outputs.force = fullfile(outDir, ['rupture_retraction_force' suffix '.png']);
outputs.video = fullfile(outDir, ['follicle_rupture_motion' suffix '.mp4']);

[fig, snapshotInfo] = plot_node_position_snapshots(outputs.simFile);
fig.Units = 'inches';
fig.Position(3:4) = [16, 6.5];
exportgraphics(fig, outputs.snapshot, 'Resolution', 300);
close(fig);
[fig, velocityInfo] = plot_percentile_ring_velocity(outputs.simFile);
exportgraphics(fig, outputs.velocity, 'Resolution', 300);
close(fig);
[fig, forceInfo] = plot_rupture_retraction_force(outputs.simFile);
exportgraphics(fig, outputs.force, 'Resolution', 300);
close(fig);
save(fullfile(outDir, ['rupture_retraction_force' suffix '.mat']), 'forceInfo');
fprintf('Saved snapshot, velocity, and force plots from %s\n', outputs.simFile);

summarize_tuning(fitFile);
outputs.verification = verify_tuning(fitFile);
[~, videoInfo] = export_follicle_rupture_motion_video(outputs.simFile, outputs.video);
reader = VideoReader(outputs.video);
decodedFrames = 0;
while hasFrame(reader)
    frame = readFrame(reader);
    assert(~isempty(frame), 'Decoded an empty video frame.');
    decodedFrames = decodedFrames + 1;
end
assert(decodedFrames == numel(videoInfo.frameIndices), 'Video frame count mismatch.');
assert(videoInfo.frameTimes(1) == simulation.tSol(1) && ...
    videoInfo.frameTimes(end) == simulation.tSol(end), 'Video omits simulation endpoints.');
assert(isequal(videoInfo.ruptureEdges(:), find(fit.param.damageMask(:) > 0)));
outputs.decodedVideoFrames = decodedFrames;
outputs.generatedAt = char(datetime('now', 'TimeZone', 'UTC'));
save(fullfile(outDir, ['visualization_info' suffix '.mat']), ...
    'outputs', 'snapshotInfo', 'velocityInfo', 'forceInfo', 'videoInfo');
fprintf('Verified %d video frames spanning 0 to 1800 s.\n', decodedFrames);
disp(outputs);
end
