function [hFig, info] = plot_node_position_snapshots(simFile)
%PLOT_NODE_POSITION_SNAPSHOTS Plot the shell at its initial, middle, and final times.
%   plot_node_position_snapshots() reads the default simulation result.
%   plot_node_position_snapshots(simFile) reads the specified MAT-file.
%
% The initial node cloud defines a common principal-axis coordinate system.
% PC1 (the football's long axis) becomes the horizontal x coordinate in all
% three panels. Node identities in the requested initial PC1 percentile
% bands are preserved through time, so their motion can be compared.

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
requiredFields = {'nodeCoords', 'edgeList', 'Y', 'tSol', 'D_time'};
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
allAligned = zeros(3*N,3);
for panelIdx = 1:3
    displacement = reshape(Y(frameIdx(panelIdx),:), 3, []).';
    absolutePosition = nodeCoords + displacement;
    alignedFrames{panelIdx} = (absolutePosition - center0) * pcaBasis * displayScale;
    rowIdx = (panelIdx-1)*N + (1:N);
    allAligned(rowIdx,:) = alignedFrames{panelIdx};
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
colors.edges = [0.60 0.62 0.66];
colors.rupture = [0.75 0.00 0.75];
bandLabels = {'0-10 percentile', '25-35 percentile', '45-55 percentile'};

hFig = figure('Color', 'w', 'Name', 'Node position snapshots');
layout = tiledlayout(hFig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('Node positions in initial principal-axis coordinates — %s', ...
    string(simFile)), 'Interpreter', 'none');

for panelIdx = 1:3
    ax = nexttile(layout);
    hold(ax, 'on');
    P = alignedFrames{panelIdx};

    damageNow = D_time(frameIdx(panelIdx),:).';
    failedEdges = damageNow >= 0.5;
    intactEdges = ~failedEdges;
    plotEdgeSegments(ax, P, edgeList(intactEdges,:), colors.edges, 0.45);

    selectedIntact = ruptureEdges(~failedEdges(ruptureEdges));
    if ~isempty(selectedIntact)
        plotEdgeSegments(ax, P, edgeList(selectedIntact,:), ...
            colors.rupture, 3.0, 'Edge selected to rupture');
    end

    selectedFailed = ruptureEdges(failedEdges(ruptureEdges));
    if ~isempty(selectedFailed)
        failedNodes = unique(edgeList(selectedFailed,:));
        scatter3(ax, P(failedNodes,1), P(failedNodes,2), P(failedNodes,3), ...
            95, colors.rupture, 'o', 'LineWidth', 2.2, ...
            'DisplayName', 'Rupture-patch nodes');
        failedMidpoints = 0.5*(P(edgeList(selectedFailed,1),:) + ...
                               P(edgeList(selectedFailed,2),:));
        scatter3(ax, failedMidpoints(:,1), failedMidpoints(:,2), failedMidpoints(:,3), ...
            150, colors.rupture, 'x', 'LineWidth', 3.0, ...
            'DisplayName', 'Ruptured edges (removed)');
    end

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

    % PC1 is the x axis in this coordinate system. Draw it through the
    % current cloud centroid so any overall translation remains visible.
    currentCenter = mean(P, 1);
    axisX = [min(P(:,1)), max(P(:,1))];
    plot3(ax, axisX, currentCenter(2)*[1 1], currentCenter(3)*[1 1], ...
        'k-', 'LineWidth', 2.2, 'DisplayName', 'Principal long axis');

    xlabel(ax, sprintf('Principal long axis (%s)', displayUnit));
    ylabel(ax, sprintf('PC2 (%s)', displayUnit));
    zlabel(ax, sprintf('PC3 (%s)', displayUnit));
    if isempty(ruptureEdges)
        ruptureStatus = 'no rupture edge selected';
    else
        ruptureStatus = sprintf('rupture damage D = %.3f', max(damageNow(ruptureEdges)));
    end
    title(ax, sprintf('%s: t = %.4g s\n%s', frameNames{panelIdx}, ...
        tSol(frameIdx(panelIdx)), ruptureStatus));
    grid(ax, 'on');
    box(ax, 'on');
    axis(ax, 'equal');
    xlim(ax, [plotMin(1), plotMax(1)]);
    ylim(ax, [plotMin(2), plotMax(2)]);
    zlim(ax, [plotMin(3), plotMax(3)]);
    % Keep the original principal-axis coordinates and change only the
    % viewpoint: negative PC1 is now the end closest to the viewer.
    view(ax, -35, 20);
end

legend(nexttile(layout, 1), 'Location', 'southoutside', 'NumColumns', 3);

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
end


function h = plotEdgeSegments(ax, P, edges, color, lineWidth, displayName)
% Plot many disconnected edges as one NaN-separated line object.
if isempty(edges)
    h = gobjects(0);
    return;
end

nEdges = size(edges,1);
edgeX = [P(edges(:,1),1), P(edges(:,2),1), nan(nEdges,1)].';
edgeY = [P(edges(:,1),2), P(edges(:,2),2), nan(nEdges,1)].';
edgeZ = [P(edges(:,1),3), P(edges(:,2),3), nan(nEdges,1)].';
h = plot3(ax, edgeX(:), edgeY(:), edgeZ(:), '-', ...
    'Color', color, 'LineWidth', lineWidth);

if nargin < 6 || isempty(displayName)
    h.HandleVisibility = 'off';
else
    h.DisplayName = displayName;
end
end
