function [hFig, info] = plot_percentile_ring_velocity(simFile)
%PLOT_PERCENTILE_RING_VELOCITY Plot signed long-axis ring motion.
%   plot_percentile_ring_velocity() reads the default simulation result.
%   plot_percentile_ring_velocity(simFile) reads the specified MAT-file.
%
% Ring membership is fixed from the undeformed node cloud, matching
% plot_node_position_snapshots. If the result contains the solver velocity
% state V, it is used directly. Older files fall back to finite differences
% of the saved displacement Y.

if nargin < 1 || isempty(simFile)
    repoRoot = fileparts(mfilename('fullpath'));
    simFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end

if ~(ischar(simFile) || (isstring(simFile) && isscalar(simFile)))
    error('plot_percentile_ring_velocity:invalidFile', ...
        'simFile must be a character vector or scalar string.');
end
simFile = char(simFile);
if ~isfile(simFile)
    error('plot_percentile_ring_velocity:fileNotFound', ...
        'Simulation result not found: %s', simFile);
end

fileVariables = who('-file', simFile);
loadVariables = {'nodeCoords', 'Y', 'tSol', 'param'};
if ismember('V', fileVariables)
    loadVariables{end+1} = 'V';
end
S = load(simFile, loadVariables{:});
requiredFields = {'nodeCoords', 'Y', 'tSol'};
for fieldIdx = 1:numel(requiredFields)
    if ~isfield(S, requiredFields{fieldIdx})
        error('plot_percentile_ring_velocity:missingData', ...
            'The result file is missing variable "%s".', requiredFields{fieldIdx});
    end
end

nodeCoords = S.nodeCoords;
Y = S.Y;
tSol = S.tSol(:);
N = size(nodeCoords, 1);
Nt = numel(tSol);
if size(nodeCoords,2) ~= 3 || size(Y,1) ~= Nt || size(Y,2) ~= 3*N
    error('plot_percentile_ring_velocity:inconsistentDimensions', ...
        'Expected nodeCoords [N x 3], Y [Nt x 3N], and tSol [Nt x 1].');
end
if Nt < 2 || any(diff(tSol) <= 0)
    error('plot_percentile_ring_velocity:invalidTime', ...
        'tSol must contain at least two strictly increasing times.');
end

velocitySource = 'saved solver velocity';
if isfield(S, 'V') && isequal(size(S.V), size(Y))
    V = S.V;
else
    velocitySource = 'finite difference of displacement';
    V = zeros(size(Y));
    V(1,:) = (Y(2,:) - Y(1,:)) / (tSol(2) - tSol(1));
    V(end,:) = (Y(end,:) - Y(end-1,:)) / (tSol(end) - tSol(end-1));
    V(2:end-1,:) = (Y(3:end,:) - Y(1:end-2,:)) ./ ...
        (tSol(3:end) - tSol(1:end-2));
end

% Use the same initial long-axis convention and ring boundaries as the
% node-position snapshot figure.
centeredNodes = nodeCoords - mean(nodeCoords, 1);
covariance = (centeredNodes.' * centeredNodes) / max(N - 1, 1);
[principalVectors, principalValues] = eig(covariance, 'vector');
[~, principalOrder] = sort(principalValues, 'descend');
longAxis = principalVectors(:, principalOrder(1));
if isfield(S, 'param') && isfield(S.param, 'postAxis') && isfield(S.param, 'postSign')
    posteriorDirection = zeros(3,1);
    posteriorDirection(S.param.postAxis) = S.param.postSign;
    if dot(longAxis, posteriorDirection) < 0
        longAxis = -longAxis;
    end
end

initialLongScore = centeredNodes * longAxis;
percentileRanges = [0 10; 25 35; 45 55];
percentileEdges = prctile(initialLongScore, reshape(percentileRanges.', 1, []));
bandMasks = false(N, size(percentileRanges,1));
for bandIdx = 1:size(percentileRanges,1)
    lowerEdge = percentileEdges(2*bandIdx - 1);
    upperEdge = percentileEdges(2*bandIdx);
    bandMasks(:,bandIdx) = initialLongScore >= lowerEdge & ...
        initialLongScore <= upperEdge;
end

velocityComponents = reshape(V.', 3, N, Nt);
nodeSpeed = squeeze(sqrt(sum(velocityComponents.^2, 1))).';
nodeAxialVelocity = squeeze(sum(velocityComponents .* ...
    reshape(longAxis, 3, 1, 1), 1)).';
displacementComponents = reshape(Y.', 3, N, Nt);
nodeDisplacement = squeeze(sqrt(sum(displacementComponents.^2, 1))).';
nodeAxialDisplacement = squeeze(sum(displacementComponents .* ...
    reshape(longAxis, 3, 1, 1), 1)).';
meanBandSpeed = zeros(Nt, size(bandMasks,2));
meanBandDisplacement = zeros(Nt, size(bandMasks,2));
meanBandAxialVelocity = zeros(Nt, size(bandMasks,2));
meanBandAxialDisplacement = zeros(Nt, size(bandMasks,2));
for bandIdx = 1:size(bandMasks,2)
    meanBandSpeed(:,bandIdx) = mean(nodeSpeed(:,bandMasks(:,bandIdx)), 2);
    meanBandDisplacement(:,bandIdx) = ...
        mean(nodeDisplacement(:,bandMasks(:,bandIdx)), 2);
    meanBandAxialVelocity(:,bandIdx) = ...
        mean(nodeAxialVelocity(:,bandMasks(:,bandIdx)), 2);
    meanBandAxialDisplacement(:,bandIdx) = ...
        mean(nodeAxialDisplacement(:,bandMasks(:,bandIdx)), 2);
end

displayScale = 1e6;
displayUnit = '\mum';
if isfield(S, 'param') && isfield(S.param, 'units') && ...
        isfield(S.param.units, 'displayFromMeters') && ...
        isfield(S.param.units, 'displayLengthUnit')
    displayScale = S.param.units.displayFromMeters;
    displayUnit = S.param.units.displayLengthUnit;
end
meanBandSpeedDisplay = meanBandSpeed * displayScale;
meanBandDisplacementDisplay = meanBandDisplacement * displayScale;
meanBandAxialVelocityDisplay = meanBandAxialVelocity * displayScale;
meanBandAxialDisplacementDisplay = meanBandAxialDisplacement * displayScale;

colors = [0.88 0.20 0.18; 0.20 0.62 0.30; 0.20 0.42 0.88];
bandLabels = {'0-10 percentile', '25-35 percentile', '45-55 percentile'};
[~, simFileName, simFileExtension] = fileparts(simFile);
hFig = figure('Color', 'w', 'Name', 'Signed percentile-ring axial motion');
layout = tiledlayout(hFig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
displacementAx = nexttile(layout);
velocityAx = nexttile(layout);

hold(displacementAx, 'on');
hold(velocityAx, 'on');
for bandIdx = 1:size(meanBandAxialVelocityDisplay,2)
    plot(displacementAx, tSol, meanBandAxialDisplacementDisplay(:,bandIdx), '-', ...
        'Color', colors(bandIdx,:), 'LineWidth', 2.0, ...
        'DisplayName', bandLabels{bandIdx});
    plot(velocityAx, tSol, meanBandAxialVelocityDisplay(:,bandIdx), '-', ...
        'Color', colors(bandIdx,:), 'LineWidth', 2.0, ...
        'DisplayName', bandLabels{bandIdx});
end
yline(displacementAx, 0, ':', 'Color', [0.35 0.35 0.35], ...
    'HandleVisibility', 'off');
yline(velocityAx, 0, ':', 'Color', [0.35 0.35 0.35], ...
    'HandleVisibility', 'off');
ylabel(displacementAx, sprintf('Mean axial displacement (%s)', displayUnit));
ylabel(velocityAx, sprintf('Mean axial velocity (%s/s)', displayUnit));

% Show the complete saved simulation, including the long post-rupture tail.
xLimits = [tSol(1), tSol(end)];
xlim(displacementAx, xLimits);
xlim(velocityAx, xLimits);
peakVisibleDisplacement = max(abs(meanBandAxialDisplacementDisplay), [], 'all');
if peakVisibleDisplacement > 0
    ylim(displacementAx, 1.05 * peakVisibleDisplacement * [-1 1]);
end
peakVisibleSpeed = max(abs(meanBandAxialVelocityDisplay), [], 'all');
if peakVisibleSpeed > 0
    ylim(velocityAx, 1.05 * peakVisibleSpeed * [-1 1]);
end

xlabel(velocityAx, 'Time (s)');
title(displacementAx, 'Signed displacement away from tip (+) / toward tip (-)');
title(velocityAx, 'Signed velocity away from tip (+) / toward tip (-)');
title(layout, sprintf('Long-axis motion of initial principal-axis percentile rings — %s', ...
    [simFileName simFileExtension]), 'Interpreter', 'none');
grid(displacementAx, 'on');
grid(velocityAx, 'on');
box(displacementAx, 'on');
box(velocityAx, 'on');
legend(displacementAx, 'Location', 'best');

info.simFile = simFile;
info.velocitySource = velocitySource;
info.longAxis = longAxis;
info.percentileRanges = percentileRanges;
info.percentileEdges = percentileEdges;
info.bandMasks = bandMasks;
info.bandNodeIndices = arrayfun(@(idx) find(bandMasks(:,idx)), ...
    1:size(bandMasks,2), 'UniformOutput', false);
info.nodeDisplacement = nodeDisplacement;
info.meanBandDisplacement = meanBandDisplacement;
info.meanBandDisplacementDisplay = meanBandDisplacementDisplay;
info.nodeAxialDisplacement = nodeAxialDisplacement;
info.meanBandAxialDisplacement = meanBandAxialDisplacement;
info.meanBandAxialDisplacementDisplay = meanBandAxialDisplacementDisplay;
info.nodeSpeed = nodeSpeed;
info.meanBandSpeed = meanBandSpeed;
info.meanBandSpeedDisplay = meanBandSpeedDisplay;
info.nodeAxialVelocity = nodeAxialVelocity;
info.meanBandAxialVelocity = meanBandAxialVelocity;
info.meanBandAxialVelocityDisplay = meanBandAxialVelocityDisplay;
info.plotTimeLimit = xLimits;
end
