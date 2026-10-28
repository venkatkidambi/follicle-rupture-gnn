function [outFile, info] = export_follicle_rupture_motion_video(matFile, outFile, opts)
%EXPORT_FOLLICLE_RUPTURE_MOTION_VIDEO Render the saved shell motion as MP4.
%   export_follicle_rupture_motion_video() uses the default simulation file
%   and renders through its final saved time.
%   export_follicle_rupture_motion_video(matFile, outFile) selects the input
%   MAT-file and output MP4. An optional opts structure can override fps,
%   quality, targetFrames, endTime, viewAz, viewEl, and damageThreshold.
%   frameTitles optionally supplies one title per saved time for tracked data.

repoRoot = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(matFile)
    matFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end
if nargin < 2 || isempty(outFile)
    outFile = fullfile(repoRoot, 'results', 'follicle_rupture_motion.mp4');
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

matFile = char(matFile);
outFile = char(outFile);
if ~isfile(matFile)
    error('export_follicle_rupture_motion_video:fileNotFound', ...
        'Simulation result not found: %s', matFile);
end

loadFields = intersect({'nodeCoords', 'edgeList', 'Y', 'tSol', ...
    'D_time', 'param'}, who('-file', matFile), 'stable');
S = load(matFile, loadFields{:});
requiredFields = {'nodeCoords', 'edgeList', 'Y', 'tSol'};
for fieldIdx = 1:numel(requiredFields)
    if ~isfield(S, requiredFields{fieldIdx})
        error('export_follicle_rupture_motion_video:missingData', ...
            'The result file is missing variable "%s".', requiredFields{fieldIdx});
    end
end

nodeCoords = S.nodeCoords;
edgeList = S.edgeList;
Y = S.Y;
tSol = S.tSol(:);
N = size(nodeCoords, 1);
Nt = numel(tSol);
Ne = size(edgeList, 1);
if size(nodeCoords,2) ~= 3 || size(Y,1) ~= Nt || size(Y,2) ~= 3*N || ...
        size(edgeList,2) ~= 2
    error('export_follicle_rupture_motion_video:inconsistentDimensions', ...
        'Expected nodeCoords [N x 3], edgeList [Ne x 2], and Y [Nt x 3N].');
end

opts = videoDefaults(opts, S);
if isfield(opts, 'frameTitles')
    assert(iscellstr(opts.frameTitles) && numel(opts.frameTitles) == Nt, ...
        'frameTitles must contain one character vector per saved time.');
end
frameCandidates = find(tSol <= opts.endTime);
if isempty(frameCandidates)
    frameCandidates = 1;
end
frameCount = min(numel(frameCandidates), max(2, round(opts.targetFrames)));
framePositions = unique(round(linspace(1, numel(frameCandidates), frameCount)));
frameIndices = frameCandidates(framePositions);

[pcaBasis, center0, bandMasks, percentileRanges] = ...
    initialPrincipalCoordinates(nodeCoords, S);
if isfield(S, 'param') && isfield(S.param, 'damageMask')
    ruptureEdges = find(S.param.damageMask(:) > 0);
elseif isfield(S, 'D_time')
    ruptureEdges = find(max(S.D_time, [], 1).' >= opts.damageThreshold);
else
    ruptureEdges = zeros(0,1);
end

displayScale = 1e6;
displayUnit = '\mum';
if isfield(S, 'param') && isfield(S.param, 'units') && ...
        isfield(S.param.units, 'displayFromMeters') && ...
        isfield(S.param.units, 'displayLengthUnit')
    displayScale = S.param.units.displayFromMeters;
    displayUnit = S.param.units.displayLengthUnit;
end
plotLimits = fixedVideoLimits(nodeCoords, Y, frameIndices, center0, ...
    pcaBasis, displayScale);

outputDir = fileparts(outFile);
if ~isempty(outputDir) && ~isfolder(outputDir)
    mkdir(outputDir);
end

frameSize = [960 720];
fig = figure('Color', 'w', 'Units', 'pixels', ...
    'Position', [100 100 frameSize], 'Resize', 'off');
ax = axes(fig, 'Position', [0.03 0.03 0.94 0.91]);
writer = VideoWriter(outFile, 'MPEG-4');
writer.FrameRate = opts.fps;
writer.Quality = opts.quality;
writerIsOpen = false;

colors.edges = [0.68 0.70 0.74];
colors.rupture = [0.75 0.00 0.75];
colors.otherNodes = [0.55 0.57 0.61];
colors.rings = [0.88 0.20 0.18; 0.20 0.62 0.30; 0.20 0.42 0.88];

hold(ax, 'on');
intactEdgePlot = plot3(ax, nan, nan, nan, '-', 'Color', colors.edges, ...
    'LineWidth', 0.55, 'HandleVisibility', 'off');
ruptureEdgePlot = plot3(ax, nan, nan, nan, '-', 'Color', colors.rupture, ...
    'LineWidth', 2.6, 'HandleVisibility', 'off');
otherNodePlot = scatter3(ax, nan, nan, nan, 18, colors.otherNodes, ...
    'filled', 'MarkerEdgeColor', 'none', 'HandleVisibility', 'off');
ringNodePlots = gobjects(size(bandMasks,2), 1);
for bandIdx = 1:size(bandMasks,2)
    ringNodePlots(bandIdx) = scatter3(ax, nan, nan, nan, 34, ...
        colors.rings(bandIdx,:), 'filled', ...
        'MarkerEdgeColor', [0.15 0.15 0.15], 'LineWidth', 0.25, ...
        'DisplayName', sprintf('%g-%g percentile', ...
        percentileRanges(bandIdx,1), percentileRanges(bandIdx,2)));
end
failedNodePlot = scatter3(ax, nan, nan, nan, 90, colors.rupture, ...
    'o', 'LineWidth', 1.8, 'HandleVisibility', 'off');
failedMidpointPlot = scatter3(ax, nan, nan, nan, 120, colors.rupture, ...
    'x', 'LineWidth', 2.4, 'HandleVisibility', 'off');

xlim(ax, plotLimits(1,:));
ylim(ax, plotLimits(2,:));
zlim(ax, plotLimits(3,:));
axis(ax, 'equal');
axis(ax, 'vis3d');
view(ax, opts.viewAz, opts.viewEl);
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, sprintf('Principal long axis (%s)', displayUnit));
ylabel(ax, sprintf('PC2 (%s)', displayUnit));
zlabel(ax, sprintf('PC3 (%s)', displayUnit));
titleHandle = title(ax, '', 'Interpreter', 'none');
legend(ax, ringNodePlots, 'Location', 'northeast');

otherNodeMask = ~any(bandMasks, 2);
try
    open(writer);
    writerIsOpen = true;
    for outputFrameIdx = 1:numel(frameIndices)
        k = frameIndices(outputFrameIdx);
        displacement = reshape(Y(k,:), 3, []).';
        absolutePosition = nodeCoords + displacement;
        P = (absolutePosition - center0) * pcaBasis * displayScale;

        if isfield(S, 'D_time')
            failedEdges = S.D_time(k,:).' >= opts.damageThreshold;
        else
            failedEdges = false(Ne,1);
        end
        intactEdges = ~failedEdges;
        selectedIntact = ruptureEdges(~failedEdges(ruptureEdges));
        selectedFailed = ruptureEdges(failedEdges(ruptureEdges));

        [edgeX, edgeY, edgeZ] = videoEdgeCoordinates(P, edgeList, find(intactEdges));
        set(intactEdgePlot, 'XData', edgeX, 'YData', edgeY, 'ZData', edgeZ);
        [edgeX, edgeY, edgeZ] = videoEdgeCoordinates(P, edgeList, selectedIntact);
        set(ruptureEdgePlot, 'XData', edgeX, 'YData', edgeY, 'ZData', edgeZ);
        set(otherNodePlot, 'XData', P(otherNodeMask,1), ...
            'YData', P(otherNodeMask,2), 'ZData', P(otherNodeMask,3));
        for bandIdx = 1:size(bandMasks,2)
            bandMask = bandMasks(:,bandIdx);
            set(ringNodePlots(bandIdx), 'XData', P(bandMask,1), ...
                'YData', P(bandMask,2), 'ZData', P(bandMask,3));
        end

        if ~isempty(selectedFailed)
            failedNodes = unique(edgeList(selectedFailed,:));
            failedMidpoints = 0.5 * (P(edgeList(selectedFailed,1),:) + ...
                P(edgeList(selectedFailed,2),:));
            set(failedNodePlot, 'XData', P(failedNodes,1), ...
                'YData', P(failedNodes,2), 'ZData', P(failedNodes,3));
            set(failedMidpointPlot, 'XData', failedMidpoints(:,1), ...
                'YData', failedMidpoints(:,2), 'ZData', failedMidpoints(:,3));
        else
            set(failedNodePlot, 'XData', nan, 'YData', nan, 'ZData', nan);
            set(failedMidpointPlot, 'XData', nan, 'YData', nan, 'ZData', nan);
        end

        titleHandle.String = sprintf(['Follicle rupture | t = %.2f s | ', ...
            'rupture edges failed = %d/%d'], tSol(k), numel(selectedFailed), ...
            numel(ruptureEdges));
        if isfield(opts, 'frameTitles')
            titleHandle.String = opts.frameTitles{k};
        end
        drawnow;
        writeVideo(writer, getframe(fig));
    end
    close(writer);
    writerIsOpen = false;
    close(fig);
catch exportError
    if writerIsOpen
        close(writer);
    end
    if isgraphics(fig)
        close(fig);
    end
    rethrow(exportError);
end

info.matFile = matFile;
info.outFile = outFile;
info.frameIndices = frameIndices;
info.frameTimes = tSol(frameIndices);
info.frameRate = opts.fps;
info.duration = numel(frameIndices) / opts.fps;
info.ruptureEdges = ruptureEdges;
info.percentileRanges = percentileRanges;
info.bandMasks = bandMasks;
info.tipNodeMask = bandMasks(:,1);
fprintf('Saved rupture video: %s (%d frames, %.2f s)\n', ...
    outFile, numel(frameIndices), info.duration);
end


function opts = videoDefaults(opts, S)
defaults.fps = 30;
defaults.quality = 95;
defaults.targetFrames = 240;
defaults.viewAz = -35;
defaults.viewEl = 20;
defaults.damageThreshold = 0.5;
defaults.endTime = S.tSol(end);

names = fieldnames(defaults);
for fieldIdx = 1:numel(names)
    name = names{fieldIdx};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
if opts.fps <= 0 || opts.quality < 0 || opts.quality > 100 || ...
        opts.targetFrames < 1 || opts.endTime < S.tSol(1)
    error('export_follicle_rupture_motion_video:invalidOptions', ...
        'Invalid fps, quality, targetFrames, or endTime option.');
end
end


function [pcaBasis, center0, bandMasks, percentileRanges] = ...
        initialPrincipalCoordinates(nodeCoords, S)
N = size(nodeCoords, 1);
center0 = mean(nodeCoords, 1);
centered0 = nodeCoords - center0;
covariance = (centered0.' * centered0) / max(N - 1, 1);
[eigenvectors, eigenvalues] = eig(covariance, 'vector');
[~, eigenOrder] = sort(eigenvalues, 'descend');
longAxis = eigenvectors(:, eigenOrder(1));
secondAxis = eigenvectors(:, eigenOrder(2));

if isfield(S, 'param') && isfield(S.param, 'postAxis') && isfield(S.param, 'postSign')
    posteriorDirection = zeros(3,1);
    posteriorDirection(S.param.postAxis) = S.param.postSign;
    if dot(longAxis, posteriorDirection) < 0
        longAxis = -longAxis;
    end
end
thirdAxis = cross(longAxis, secondAxis);
thirdAxis = thirdAxis / norm(thirdAxis);
secondAxis = cross(thirdAxis, longAxis);
secondAxis = secondAxis / norm(secondAxis);
pcaBasis = [longAxis, secondAxis, thirdAxis];

longScore = centered0 * longAxis;
percentileRanges = [0 10; 25 35; 45 55];
percentileEdges = prctile(longScore, reshape(percentileRanges.', 1, []));
bandMasks = false(N, size(percentileRanges,1));
for bandIdx = 1:size(percentileRanges,1)
    lowerEdge = percentileEdges(2*bandIdx - 1);
    upperEdge = percentileEdges(2*bandIdx);
    bandMasks(:,bandIdx) = longScore >= lowerEdge & longScore <= upperEdge;
end
end


function plotLimits = fixedVideoLimits(nodeCoords, Y, frameIndices, center0, ...
        pcaBasis, displayScale)
dataMin = [inf inf inf];
dataMax = [-inf -inf -inf];
for frameIdx = frameIndices(:).'
    displacement = reshape(Y(frameIdx,:), 3, []).';
    P = (nodeCoords + displacement - center0) * pcaBasis * displayScale;
    dataMin = min(dataMin, min(P, [], 1));
    dataMax = max(dataMax, max(P, [], 1));
end
dataRange = dataMax - dataMin;
dataRange(~isfinite(dataRange) | dataRange <= 0) = 1;
plotLimits = [dataMin; dataMax].';
plotLimits(:,1) = plotLimits(:,1) - 0.08 * dataRange.';
plotLimits(:,2) = plotLimits(:,2) + 0.08 * dataRange.';
end


function [X, Y, Z] = videoEdgeCoordinates(P, edgeList, edgeIndices)
if isempty(edgeIndices)
    X = nan;
    Y = nan;
    Z = nan;
    return;
end
i = edgeList(edgeIndices,1);
j = edgeList(edgeIndices,2);
X = [P(i,1), P(j,1), nan(numel(edgeIndices),1)].';
Y = [P(i,2), P(j,2), nan(numel(edgeIndices),1)].';
Z = [P(i,3), P(j,3), nan(numel(edgeIndices),1)].';
X = X(:);
Y = Y(:);
Z = Z(:);
end
