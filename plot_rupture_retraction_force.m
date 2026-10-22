function [hFig, info] = plot_rupture_retraction_force(simFile)
%PLOT_RUPTURE_RETRACTION_FORCE Plot patch tension and signed tissue pull.
% Forces are evaluated by the simulation RHS and saved in SI units. The
% rupture estimate is the total tensile link load in the selected patch at
% the last saved sample before damage exceeds 1e-6. It is a sum of link
% tensions, not a resultant force across a tissue section or a failure law.
% Retraction force is the mean projected internal link force per node along
% the initial long axis, positive away from the tip. External drag, adhesion,
% anchoring, surface constraint and active traction are excluded.

if nargin < 1 || isempty(simFile)
    repoRoot = fileparts(mfilename('fullpath'));
    simFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end
if ~(ischar(simFile) || (isstring(simFile) && isscalar(simFile)))
    error('plot_rupture_retraction_force:invalidFile', ...
        'simFile must be a character vector or scalar string.');
end
if ~isfile(simFile)
    error('plot_rupture_retraction_force:fileNotFound', ...
        'Simulation result not found: %s', simFile);
end
required = {'tSol', 'nodeCoords', 'param', 'D_time', ...
    'edgeForce_time', 'nodeInternalForce_time'};
if ~all(ismember(required, who('-file', simFile)))
    error('plot_rupture_retraction_force:missingForceData', ...
        ['This result lacks saved force histories. Rerun run_simulation ', ...
         'to evaluate forces with the viscoelastic history states.']);
end
S = load(simFile, required{:});
t = S.tSol(:);
N = size(S.nodeCoords, 1);
Nt = numel(t);
if ~isfield(S.param, 'damageMask')
    error('plot_rupture_retraction_force:missingPatch', ...
        'The result must contain param.damageMask.');
end
patchMask = S.param.damageMask(:) > 0;
Ne = numel(patchMask);
if Nt < 2 || any(~isfinite(t)) || any(diff(t) <= 0) || ...
        size(S.nodeCoords,2) ~= 3 || ...
        ~isequal(size(S.edgeForce_time), [Nt Ne]) || ...
        ~isequal(size(S.D_time), [Nt Ne]) || ...
        ~isequal(size(S.nodeInternalForce_time), [Nt 3*N])
    error('plot_rupture_retraction_force:inconsistentData', ...
        'Expected increasing times and matching node, edge, damage and force arrays.');
end
if any(~isfinite(S.edgeForce_time), 'all') || ...
        any(~isfinite(S.nodeInternalForce_time), 'all') || ...
        any(~isfinite(S.D_time), 'all')
    error('plot_rupture_retraction_force:nonfiniteData', ...
        'Force and damage histories must be finite.');
end

% Match the initial principal-axis bands used in the motion figures.
centeredNodes = S.nodeCoords - mean(S.nodeCoords, 1);
[vectors, values] = eig(centeredNodes.' * centeredNodes, 'vector');
[~, idx] = max(values);
longAxis = vectors(:,idx);
if isfield(S.param, 'postAxis') && isfield(S.param, 'postSign')
    if longAxis(S.param.postAxis) * S.param.postSign < 0
        longAxis = -longAxis;
    end
end
score = centeredNodes * longAxis;
percentileRanges = [0 10; 25 35; 45 55];
bounds = prctile(score, reshape(percentileRanges.', 1, []));
bandMasks = false(N, 3);
components = reshape(S.nodeInternalForce_time.', 3, N, Nt);
nodeAxialForce = reshape(sum(components .* reshape(longAxis, 3, 1, 1), 1), N, Nt).';
meanBandAxialForce = nan(Nt, 3);
for band = 1:3
    bandMasks(:,band) = score >= bounds(2*band-1) & score <= bounds(2*band);
    if any(bandMasks(:,band))
        meanBandAxialForce(:,band) = mean(nodeAxialForce(:,bandMasks(:,band)), 2);
    end
end

patchTension = sum(max(S.edgeForce_time(:,patchMask), 0), 2);
damageTolerance = 1e-6;
firstDamage = find(any(S.D_time(:,patchMask) > damageTolerance, 2), 1);
ruptureForce = NaN;
ruptureTimeBracket = [NaN NaN];
if ~any(patchMask)
    ruptureLabel = 'No rupture patch selected';
elseif isempty(firstDamage)
    ruptureLabel = 'No rupture detected in saved samples';
elseif firstDamage == 1
    ruptureLabel = 'Patch already damaged at first sample; preload unavailable';
else
    ruptureForce = patchTension(firstDamage-1);
    ruptureTimeBracket = t([firstDamage-1 firstDamage]).';
    ruptureLabel = sprintf('Last intact sample: %.3g \\muN; onset between %.3g and %.3g s', ...
        ruptureForce * 1e6, ruptureTimeBracket);
end

hFig = figure('Color', 'w', 'Name', 'Rupture and retraction forces');
hFig.Units = 'inches';
hFig.Position(3:4) = [11 8];
layout = tiledlayout(hFig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
patchAx = nexttile(layout);
plot(patchAx, t, patchTension * 1e6, 'Color', [0.55 0.20 0.65], 'LineWidth', 2);
hold(patchAx, 'on');
if isfinite(ruptureForce)
    plot(patchAx, t(firstDamage-1), ruptureForce * 1e6, 'ko', ...
        'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
end
ylabel(patchAx, 'Sum of tensile link forces (\muN)');
title(patchAx, {'Rupture-patch tensile load', ruptureLabel});

retractionAx = nexttile(layout);
hold(retractionAx, 'on');
colors = [0.88 0.20 0.18; 0.20 0.62 0.30; 0.20 0.42 0.88];
for band = 1:3
    plot(retractionAx, t, meanBandAxialForce(:,band) * 1e6, ...
        'Color', colors(band,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('%d-%d percentile', percentileRanges(band,:)));
end
yline(retractionAx, 0, ':', 'HandleVisibility', 'off');
ylabel(retractionAx, 'Mean axial force per node (\muN)');
title(retractionAx, {'Retraction: internal tissue pull away from tip (+) / toward tip (-)', ...
    'Projected link forces; excludes external loads and drag'});
legend(retractionAx, 'Location', 'best');
xlabel(retractionAx, 'Time (s)');
for ax = [patchAx retractionAx]
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [t(1) t(end)]);
    if ~isempty(firstDamage)
        xline(ax, t(firstDamage), '--', 'First sampled damage', ...
            'HandleVisibility', 'off');
    end
end
linkaxes([patchAx retractionAx], 'x');
title(layout, 'Rupture and retraction forces');

info.simFile = simFile;
info.time = t;
info.forceUnit = 'N';
info.patchMask = patchMask;
info.patchTension = patchTension;
info.damageTolerance = damageTolerance;
info.firstDamageIndex = firstDamage;
info.ruptureForce = ruptureForce;
info.ruptureTimeBracket = ruptureTimeBracket;
info.longAxis = longAxis;
info.percentileRanges = percentileRanges;
info.bandMasks = bandMasks;
info.meanBandAxialForce = meanBandAxialForce;
end
