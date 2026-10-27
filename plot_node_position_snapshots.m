function [hFig, info] = plot_node_position_snapshots(simFile)
%PLOT_NODE_POSITION_SNAPSHOTS Compare tracked and simulated shell positions.
%   plot_node_position_snapshots() reads the default simulation result.
%   plot_node_position_snapshots(simFile) reads the specified MAT-file.
%
% The top row shows the beginning, middle, and end frames from the tracked
% experimental input. The bottom row shows the corresponding beginning,
% middle, and end snapshots from the simulation. The simulation's initial
% node cloud defines the common principal-axis coordinate system.

if nargin < 1 || isempty(simFile)
    repoRoot = fileparts(mfilename('fullpath'));
    simFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end

if ~(ischar(simFile) || (isstring(simFile) && isscalar(simFile)))
    error('plot_node_position_snapshots:invalidFile', ...
        'simFile must be a character vector or scalar string.');
end
simFile = char(simFile);
if ~isfile(simFile)
    error('plot_node_position_snapshots:fileNotFound', ...
        'Simulation result not found: %s', simFile);
end

S = load(simFile, 'nodeCoords', 'edgeList', 'Y', 'tSol', 'D_time', 'param');
requiredFields = {'nodeCoords', 'edgeList', 'Y', 'tSol', 'D_time', 'param'};
for fieldIdx = 1:numel(requiredFields)
    if ~isfield(S, requiredFields{fieldIdx})
        error('plot_node_position_snapshots:missingData', ...
            'The result file is missing variable "%s".', requiredFields{fieldIdx});
    end
end

nodeCoords = S.nodeCoords;
tSol = S.tSol(:);
Y = S.Y;
edgeList = S.edgeList;
D_time = S.D_time;
N = size(nodeCoords, 1);
if size(nodeCoords, 2) ~= 3 || size(Y, 2) ~= 3*N || size(Y, 1) ~= numel(tSol)
    error('plot_node_position_snapshots:inconsistentDimensions', ...
        'Expected nodeCoords [N x 3], Y [Nt x 3N], and tSol [Nt x 1].');
end
if size(edgeList, 2) ~= 2 || size(D_time, 1) ~= numel(tSol) || ...
        size(D_time, 2) ~= size(edgeList, 1)
    error('plot_node_position_snapshots:inconsistentEdgeDimensions', ...
        'Expected edgeList [Ne x 2] and D_time [Nt x Ne].');
end
if numel(tSol) < 3
    error('plot_node_position_snapshots:tooFewFrames', ...
        'At least three saved time points are required.');
end

frameIdx = [1, round((numel(tSol) + 1)/2), numel(tSol)];
frameNames = {'Beginning', 'Middle', 'End'};
panelLabels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

% Establish a right-handed PCA basis from the undeformed node cloud. This
% uses EIG rather than PCA, so the Statistics Toolbox is not required.
center0 = mean(nodeCoords, 1);
centered0 = nodeCoords - center0;
covariance = (centered0.' * centered0) / max(N - 1, 1);
[eigenvectors, eigenvalues] = eig(covariance, 'vector');
[~, eigenOrder] = sort(eigenvalues, 'descend');
longAxis = eigenvectors(:, eigenOrder(1));
secondAxis = eigenvectors(:, eigenOrder(2));

% Resolve the arbitrary PCA sign using the model's posterior convention.
% Increasing aligned x therefore points toward the posterior when metadata
% are available.
if isfield(S, 'param') && isfield(S.param, 'postAxis') && isfield(S.param, 'postSign')
    posteriorDirection = zeros(3,1);
    posteriorDirection(S.param.postAxis) = S.param.postSign;
    if dot(longAxis, posteriorDirection) < 0
        longAxis = -longAxis;
    end
else
    [~, dominantComponent] = max(abs(longAxis));
    if longAxis(dominantComponent) < 0
        longAxis = -longAxis;
    end
end

thirdAxis = cross(longAxis, secondAxis);
thirdAxis = thirdAxis / norm(thirdAxis);
secondAxis = cross(thirdAxis, longAxis);
secondAxis = secondAxis / norm(secondAxis);
pcaBasis = [longAxis, secondAxis, thirdAxis];

initialAligned = centered0 * pcaBasis;
initialLongScore = initialAligned(:,1);
percentileEdges = prctile(initialLongScore, [0 10 25 35 45 55]);

bandMasks = false(N,3);
bandMasks(:,1) = initialLongScore >= percentileEdges(1) & ...
                 initialLongScore <= percentileEdges(2);
bandMasks(:,2) = initialLongScore >= percentileEdges(3) & ...
                 initialLongScore <= percentileEdges(4);
bandMasks(:,3) = initialLongScore >= percentileEdges(5) & ...
                 initialLongScore <= percentileEdges(6);
otherMask = ~any(bandMasks, 2);

if isfield(S.param, 'damageMask')
    ruptureEdges = find(S.param.damageMask(:) > 0);
else
    ruptureEdges = find(max(D_time, [], 1).' > 0.5);
end

displayScale = 1e6;
displayUnit = '\mum';
if isfield(S, 'param') && isfield(S.param, 'units') && ...
        isfield(S.param.units, 'displayFromMeters') && ...
        isfield(S.param.units, 'displayLengthUnit')
    displayScale = S.param.units.displayFromMeters;
    displayUnit = S.param.units.displayLengthUnit;
end

alignedFrames = cell(1,3);
for panelIdx = 1:3
    displacement = reshape(Y(frameIdx(panelIdx),:), 3, []).';
    absolutePosition = nodeCoords + displacement;
    alignedFrames{panelIdx} = (absolutePosition - center0) * pcaBasis * displayScale;
end

[experimentalFrames, experimentalFrameIdx, experimentalValidMasks, ...
    experimentalBandMasks, experimentalMetadata] = loadExperimentalFrames( ...
    S.param, nodeCoords, bandMasks, center0, pcaBasis, displayScale);

allAligned = vertcat(alignedFrames{:});
for panelIdx = 1:3
    allAligned = [allAligned; ...
        experimentalFrames{panelIdx}(experimentalValidMasks{panelIdx},:)]; %#ok<AGROW>
end

dataMin = min(allAligned, [], 1);
dataMax = max(allAligned, [], 1);
dataRange = max(dataMax - dataMin, eps);
plotMin = dataMin - 0.08*dataRange;
plotMax = dataMax + 0.08*dataRange;

colors.other = [0.72 0.74 0.78];
colors.bands = [0.88 0.20 0.18; ...
                0.20 0.62 0.30; ...
                0.20 0.42 0.88];
bandLabels = {'0-10 percentile', '25-35 percentile', '45-55 percentile'};

hFig = figure('Color', 'w', 'Name', 'Experimental and simulated node position snapshots');
layout = tiledlayout(hFig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, ...
    'Tracked experiment (top) and simulation (bottom) in initial principal-axis coordinates');

axesHandles = gobjects(2,3);
for panelIdx = 1:3
    ax = nexttile(layout, panelIdx);
    axesHandles(1,panelIdx) = ax;
    hold(ax, 'on');
    P = experimentalFrames{panelIdx};
    validMask = experimentalValidMasks{panelIdx};
    currentBandMasks = experimentalBandMasks & validMask;
    otherExperimentalMask = validMask & ~any(experimentalBandMasks, 2);

    scatter3(ax, P(otherExperimentalMask,1), P(otherExperimentalMask,2), ...
        P(otherExperimentalMask,3), 18, colors.other, 'filled', ...
        'MarkerFaceAlpha', 0.50, 'DisplayName', 'Other nodes');
    for bandIdx = 1:3
        mask = currentBandMasks(:,bandIdx);
        scatter3(ax, P(mask,1), P(mask,2), P(mask,3), ...
            30, colors.bands(bandIdx,:), 'filled', ...
            'MarkerEdgeColor', [0.15 0.15 0.15], ...
            'LineWidth', 0.35, 'DisplayName', bandLabels{bandIdx});
    end

    drawLongAxis(ax, P(validMask,:));
    formatSnapshotAxes(ax, plotMin, plotMax, displayUnit);
    title(ax, sprintf('%s Experimental — %s\nframe %d of %d', ...
        panelLabels{panelIdx}, frameNames{panelIdx}, experimentalFrameIdx(panelIdx), ...
        experimentalMetadata.frameCount));
end

for panelIdx = 1:3
    ax = nexttile(layout, panelIdx + 3);
    axesHandles(2,panelIdx) = ax;
    hold(ax, 'on');
    P = alignedFrames{panelIdx};

    scatter3(ax, P(otherMask,1), P(otherMask,2), P(otherMask,3), ...
        18, colors.other, 'filled', 'MarkerFaceAlpha', 0.50, ...
        'DisplayName', 'Other nodes');
    for bandIdx = 1:3
        mask = bandMasks(:,bandIdx);
        scatter3(ax, P(mask,1), P(mask,2), P(mask,3), ...
            30, colors.bands(bandIdx,:), 'filled', ...
            'MarkerEdgeColor', [0.15 0.15 0.15], ...
            'LineWidth', 0.35, 'DisplayName', bandLabels{bandIdx});
    end

    drawLongAxis(ax, P);
    title(ax, sprintf('%s Simulation — %s: t = %.4g s', ...
        panelLabels{panelIdx + 3}, frameNames{panelIdx}, tSol(frameIdx(panelIdx))));
    formatSnapshotAxes(ax, plotMin, plotMax, displayUnit);
end

legend(axesHandles(2,1), 'Location', 'southoutside', 'NumColumns', 3);

info.simFile = simFile;
info.frameIndices = frameIdx;
info.frameTimes = tSol(frameIdx);
info.initialCenter = center0;
info.pcaBasis = pcaBasis;
info.longAxis = longAxis;
info.percentileEdges = percentileEdges;
info.bandMasks = bandMasks;
info.bandNodeIndices = {find(bandMasks(:,1)), find(bandMasks(:,2)), find(bandMasks(:,3))};
info.ruptureEdges = ruptureEdges;
info.ruptureDamageAtFrames = D_time(frameIdx, ruptureEdges);
info.experimentalFile = experimentalMetadata.file;
info.experimentalVariable = experimentalMetadata.variable;
info.experimentalFrameIndices = experimentalFrameIdx;
info.experimentalFrameCount = experimentalMetadata.frameCount;
info.experimentalValidMasks = experimentalValidMasks;
info.experimentalAffineTransform = experimentalMetadata.affineTransform;
end


function [alignedFrames, frameIdx, validMasks, experimentalBandMasks, metadata] = ...
        loadExperimentalFrames(param, nodeCoords, bandMasks, center0, pcaBasis, ...
        displayScale)
requiredParamFields = {'useExperimentalNodeData', 'experimentalNodeDataFile', ...
    'experimentalNodeDataVariable', 'experimentalNodeDataFrame'};
for fieldIdx = 1:numel(requiredParamFields)
    if ~isfield(param, requiredParamFields{fieldIdx})
        error('plot_node_position_snapshots:missingExperimentalMetadata', ...
            'Simulation metadata is missing param.%s.', requiredParamFields{fieldIdx});
    end
end
if ~param.useExperimentalNodeData
    error('plot_node_position_snapshots:noExperimentalInput', ...
        'This six-panel figure requires a simulation initialized from tracked experimental data.');
end

experimentalFile = char(param.experimentalNodeDataFile);
if ~isfile(experimentalFile)
    [~, sourceName, sourceExtension] = fileparts(experimentalFile);
    repoRoot = fileparts(mfilename('fullpath'));
    fallbackFile = fullfile(repoRoot, 'input_node_data', 'tracked_nodes', ...
        [sourceName sourceExtension]);
    if isfile(fallbackFile)
        experimentalFile = fallbackFile;
    else
        error('plot_node_position_snapshots:experimentalFileNotFound', ...
            'Tracked experimental data not found: %s', experimentalFile);
    end
end

experimentalVariable = char(param.experimentalNodeDataVariable);
source = load(experimentalFile, experimentalVariable);
if ~isfield(source, experimentalVariable)
    error('plot_node_position_snapshots:missingExperimentalVariable', ...
        'Experimental file does not contain variable "%s".', experimentalVariable);
end
tracked = source.(experimentalVariable);
if ~isnumeric(tracked) || ndims(tracked) > 3 || size(tracked,2) ~= 3
    error('plot_node_position_snapshots:invalidExperimentalData', ...
        'Experimental variable "%s" must be numeric [N x 3 x T].', experimentalVariable);
end

frameCount = size(tracked,3);
if frameCount < 3
    error('plot_node_position_snapshots:tooFewExperimentalFrames', ...
        'At least three tracked experimental frames are required.');
end
frameIdx = [1, round((frameCount + 1)/2), frameCount];

if isfield(param, 'experimentalCellIds') && ...
        numel(param.experimentalCellIds) == size(nodeCoords,1)
    trackIds = param.experimentalCellIds(:);
else
    trackIds = (1:size(nodeCoords,1)).';
end
if any(trackIds < 1 | trackIds > size(tracked,1) | trackIds ~= round(trackIds))
    error('plot_node_position_snapshots:invalidTrackIds', ...
        'Saved experimental track IDs do not index the tracked node array.');
end

initialFrame = param.experimentalNodeDataFrame;
if initialFrame < 1 || initialFrame > frameCount
    error('plot_node_position_snapshots:invalidInitialExperimentalFrame', ...
        'The saved experimental initialization frame is unavailable.');
end
calibrationRaw = double(tracked(trackIds,:,initialFrame));
if any(~isfinite(calibrationRaw), 'all') || any(all(calibrationRaw == 0, 2))
    error('plot_node_position_snapshots:invalidCalibrationTracks', ...
        'Tracks used to initialize the simulation are invalid in the source data.');
end

% Recover the exact affine preprocessing that mapped raw tracked coordinates
% to the saved initial simulation nodes (unit conversion, centering, and
% principal-axis orientation). Applying it to all tracked frames preserves
% the measured experimental motion in the simulation's display coordinates.
affineTransform = [calibrationRaw, ones(size(calibrationRaw,1),1)] \ nodeCoords;

trackCount = size(tracked,1);
experimentalBandMasks = false(trackCount,3);
for bandIdx = 1:3
    experimentalBandMasks(trackIds(bandMasks(:,bandIdx)),bandIdx) = true;
end

alignedFrames = cell(1,3);
validMasks = cell(1,3);
for panelIdx = 1:3
    rawFrame = double(tracked(:,:,frameIdx(panelIdx)));
    validMasks{panelIdx} = all(isfinite(rawFrame), 2) & ~all(rawFrame == 0, 2);
    processedFrame = [rawFrame, ones(trackCount,1)] * affineTransform;
    alignedFrames{panelIdx} = (processedFrame - center0) * pcaBasis * displayScale;
end

metadata.file = experimentalFile;
metadata.variable = experimentalVariable;
metadata.frameCount = frameCount;
metadata.affineTransform = affineTransform;
end


function drawLongAxis(ax, P)
% PC1 is the x axis; draw it through the current valid cloud centroid.
currentCenter = mean(P, 1);
axisX = [min(P(:,1)), max(P(:,1))];
plot3(ax, axisX, currentCenter(2)*[1 1], currentCenter(3)*[1 1], ...
    'k-', 'LineWidth', 2.2, 'DisplayName', 'Principal long axis');
end


function formatSnapshotAxes(ax, plotMin, plotMax, displayUnit)
xlabel(ax, sprintf('Principal long axis (%s)', displayUnit));
ylabel(ax, sprintf('PC2 (%s)', displayUnit));
zlabel(ax, sprintf('PC3 (%s)', displayUnit));
grid(ax, 'on');
box(ax, 'on');
axis(ax, 'equal');
xlim(ax, [plotMin(1), plotMax(1)]);
ylim(ax, [plotMin(2), plotMax(2)]);
zlim(ax, [plotMin(3), plotMax(3)]);
% Negative PC1 is the end closest to the viewer.
view(ax, -35, 20);
end
