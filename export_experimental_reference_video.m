function outputs = export_experimental_reference_video()
%EXPORT_EXPERIMENTAL_REFERENCE_VIDEO Render measured tracks in simulation style.
% Raw motion uses only the fixed initial calibration/orientation. A companion
% removes rigid cap motion using the same registration as the fit comparison.
% Missing tracks stay hidden; measured frames are held, never interpolated.
repo = fileparts(mfilename('fullpath'));
addpath(fullfile(repo, 'src'));
outDir = fullfile(repo, 'results');
referenceFile = fullfile(outDir, 'visualization_info_network_tuned.mat');
reference = load(referenceFile, 'videoInfo');
simulationFile = fullfile(outDir, 'follicleShellSim_network_tuned.mat');
simulation = load(simulationFile, 'param');
base = simulation.param;
data = tuning_prepare_data(base, 1800);
source = load(base.experimentalNodeDataFile, base.experimentalNodeDataVariable);
raw = double(source.(base.experimentalNodeDataVariable));
raw = raw(base.experimentalCellIds,:,:);
affine = [raw(:,:,1), ones(base.N,1)] \ base.nodeCoords;
assert(max(abs([raw(:,:,1),ones(base.N,1)]*affine-base.nodeCoords), [], 'all') < 1e-12);

% Match the existing simulation video's frame times and playback rate.
tSol = reference.videoInfo.frameTimes(:);
sourceFrameIndices = zeros(size(tSol));
for idx = 1:numel(tSol)
    sourceFrameIndices(idx) = find(data.time <= tSol(idx) + 1e-9, 1, 'last');
end
assert(isequal(unique(sourceFrameIndices), (1:numel(data.time)).'), ...
    'Reference playback must include every measured frame.');
nodeCoords = base.nodeCoords;
edgeList = zeros(0,2); % Tracks do not measure cell contacts or link failures.
param = struct('postAxis',base.postAxis,'postSign',base.postSign, ...
    'units',base.units,'damageMask',zeros(0,1));
opts = struct('fps',reference.videoInfo.frameRate,'targetFrames',numel(tSol));
renderFile = [tempname '.mat'];
cleanup = onCleanup(@() delete(renderFile));
outputs = struct();
modes = {'raw','registered'};
for modeIdx = 1:numel(modes)
    mode = modes{modeIdx};
    if strcmp(mode,'raw')
        position = data.rawPosition;
        label = 'Experiment (raw tracks)';
    else
        position = data.position;
        label = 'Experiment (pose corrected)';
    end
    Y = zeros(numel(tSol),3*base.N);
    opts.frameTitles = cell(numel(tSol),1);
    for idx = 1:numel(tSol)
        frame = sourceFrameIndices(idx);
        displacement = position(:,:,frame)-nodeCoords;
        displacement(~data.valid(:,frame),:) = NaN;
        Y(idx,:) = reshape(displacement.',1,[]);
        opts.frameTitles{idx} = sprintf('%s | measured t = %.1f s | frame %d/%d', ...
            label, data.time(frame), frame, numel(data.time));
    end
    save(renderFile, 'nodeCoords','edgeList','Y','tSol','param');
    outFile = fullfile(outDir, ['experimental_tracked_motion_' mode '.mp4']);
    [~, info] = export_follicle_rupture_motion_video(renderFile, outFile, opts);
    assert(isequal(info.bandMasks, reference.videoInfo.bandMasks), 'Ring membership changed.');
    assert(isempty(info.ruptureEdges));
    reader = VideoReader(outFile);
    decoded = 0;
    while hasFrame(reader)
        frameImage = readFrame(reader);
        decoded = decoded+1;
        if decoded == 1
            imwrite(frameImage, fullfile(outDir, ['experimental_tracked_motion_' mode '_preview.png']));
        end
    end
    assert(decoded == numel(tSol), 'Video frame count mismatch.');
    assert(abs(reader.Duration-reference.videoInfo.duration) < 1/opts.fps);
    info.matFile = base.experimentalNodeDataFile;
    info.sourceVariable = base.experimentalNodeDataVariable;
    info.sourceCellIds = base.experimentalCellIds;
    info.sourceFrameIndices = sourceFrameIndices;
    info.measuredFrameTimes = data.time(sourceFrameIndices);
    info.motionMode = mode;
    info.missingTracksHidden = true;
    info.temporalInterpolation = 'None; hold last measured frame';
    info.decodedFrames = decoded;
    outputs.(mode) = info;
end
outputs.initialAffineTransform = affine;
outputs.registrationDescription = data.registrationDescription;
outputs.referenceSimulationFile = simulationFile;
outputs.acquisitionTiming = '45 uniformly spaced frames across 1800 s (1800/44 s intervals)';
save(fullfile(outDir, 'experimental_reference_info.mat'), 'outputs', 'data');
fprintf('Verified both experimental videos: %d measured frames, %d playback frames each.\n', ...
    numel(data.time), numel(tSol));
end
