function confirm_actomyosin_working(matFile)
% Rupture-focused postprocessing diagnostics for follicleShellSim_*.mat.
% The function name is kept for compatibility with older run scripts.

if nargin < 1 || isempty(matFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    matFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end

close all;

velocityRingBandsPct = [0 10; 25 35; 45 55];
nVelocityRings = size(velocityRingBandsPct, 1);

S = load(matFile);
nodeCoords = S.nodeCoords;
edgeList   = S.edgeList;
tSol       = S.tSol(:);
Y          = S.Y;
param      = S.param;
units      = follicle_unit_system(param);
displayFromMeters = units.displayFromMeters;
nodeCoordsDisplay = nodeCoords * displayFromMeters;

N  = size(nodeCoords,1);
Nt = numel(tSol);
E  = edgeList;
Ne = size(E,1);

q   = reshape(Y, Nt, N, 3);
pos = q + reshape(nodeCoords, 1, N, 3);

if ~isfield(param,'postAxis') || ~isfield(param,'postSign')
    error('param.postAxis and param.postSign are required to reconstruct masks.');
end

ax = param.postAxis;
sg = param.postSign;
edgeMid0 = 0.5*(nodeCoords(E(:,1),:) + nodeCoords(E(:,2),:));
sNodes = sg * nodeCoords(:,ax);
sMid   = sg * edgeMid0(:,ax);

if isfield(param,'bodyMask')
    bodyEdge = logical(param.bodyMask(:));
elseif isfield(param,'motorMask')
    bodyEdge = logical(param.motorMask(:));
else
    mainPct = getParamOrDefault(param, 'mainPct', 80);
    sCutMain = prctile(sNodes, mainPct);
    bodyEdge = sMid < sCutMain;
end

if isfield(param,'damageMask')
    damageEdge = logical(param.damageMask(:));
else
    damageEdge = false(Ne,1);
end

if isfield(S,'D_time') && isequal(size(S.D_time), [Nt, Ne])
    damageTime = min(max(S.D_time, 0), 1);
else
    damageTime = zeros(Nt, Ne);
end

posteriorToAnteriorNodeFrac = (max(sNodes) - sNodes) ./ max(range(sNodes), eps);
posteriorToAnteriorEdgeFrac = (max(sNodes) - sMid) ./ max(range(sNodes), eps);
ringBands = velocityRingBandsPct ./ 100;
ringLabels = cell(1, nVelocityRings);
ringIdNode = zeros(N,1);
ringIdEdge = zeros(Ne,1);

for r = 1:nVelocityRings
    lo = ringBands(r,1);
    hi = ringBands(r,2);
    if r == nVelocityRings
        nodeInRing = posteriorToAnteriorNodeFrac >= lo & posteriorToAnteriorNodeFrac <= hi;
        edgeInRing = posteriorToAnteriorEdgeFrac >= lo & posteriorToAnteriorEdgeFrac <= hi;
    else
        nodeInRing = posteriorToAnteriorNodeFrac >= lo & posteriorToAnteriorNodeFrac < hi;
        edgeInRing = posteriorToAnteriorEdgeFrac >= lo & posteriorToAnteriorEdgeFrac < hi;
    end
    ringIdNode(nodeInRing) = r;
    ringIdEdge(edgeInRing) = r;
    ringLabels{r} = sprintf('%g-%g%%', velocityRingBandsPct(r,1), velocityRingBandsPct(r,2));
end

sPos = sg * squeeze(pos(:,:,ax));
vAnterior = -timeGradientRows(sPos, tSol);
ringVelocity = nan(Nt, nVelocityRings);
ringCounts = zeros(1, nVelocityRings);
meanDamageRing = nan(Nt, nVelocityRings);

for r = 1:nVelocityRings
    nodeInRing = ringIdNode == r;
    edgeInRing = ringIdEdge == r;
    ringCounts(r) = nnz(nodeInRing);
    if any(nodeInRing)
        ringVelocity(:,r) = mean(vAnterior(:,nodeInRing), 2);
    end
    if any(edgeInRing)
        meanDamageRing(:,r) = mean(damageTime(:,edgeInRing), 2);
    end
end

ringVelocityMag = abs(ringVelocity);
if any(damageEdge)
    selectedDamage = mean(damageTime(:,damageEdge), 2);
else
    selectedDamage = zeros(Nt,1);
end

figDashboard = figure('Color','w','Position',[100 100 1500 850]);
ringColors = lines(nVelocityRings);

subplot(2,3,1); hold on; grid on; box on;
velLeg = cell(1, nVelocityRings);
for r = 1:nVelocityRings
    plot(tSol, ringVelocity(:,r) * displayFromMeters, ...
        'LineWidth', 1.8, 'Color', ringColors(r,:));
    velLeg{r} = sprintf('%s, n=%d', ringLabels{r}, ringCounts(r));
end
yline(0, ':');
xlabel('Time (s)');
ylabel(sprintf('Anterior-axis velocity (%s)', units.displayVelocityUnit));
title('Band velocity: anteriorward positive');
legend(velLeg, 'Location','best');

subplot(2,3,2); hold on; grid on; box on;
otherNode = ringIdNode == 0;
scatter3(nodeCoordsDisplay(otherNode,1), nodeCoordsDisplay(otherNode,2), nodeCoordsDisplay(otherNode,3), ...
    8, [0.78 0.78 0.78], 'filled');
locLeg = {'outside bands'};
for r = 1:nVelocityRings
    inRing = ringIdNode == r;
    scatter3(nodeCoordsDisplay(inRing,1), nodeCoordsDisplay(inRing,2), nodeCoordsDisplay(inRing,3), ...
        16, ringColors(r,:), 'filled');
    locLeg{end+1} = ringLabels{r}; %#ok<AGROW>
end
axis equal; view(3);
xlabel(sprintf('x (%s)', units.displayLengthUnit));
ylabel(sprintf('y (%s)', units.displayLengthUnit));
zlabel(sprintf('z (%s)', units.displayLengthUnit));
title('Velocity bands');
legend(locLeg, 'Location','bestoutside');

subplot(2,3,3); hold on; grid on; box on;
for r = 1:nVelocityRings
    plot(tSol, meanDamageRing(:,r), 'LineWidth', 1.8, 'Color', ringColors(r,:));
end
xlabel('Time (s)');
ylabel('Mean damage D');
title('Damage by edge band');
legend(ringLabels, 'Location','best');

subplot(2,3,4); hold on; grid on; box on;
for r = 1:nVelocityRings
    plot(tSol, ringVelocityMag(:,r) * displayFromMeters, ...
        'LineWidth', 1.8, 'Color', ringColors(r,:));
end
yline(0, ':');
xlabel('Time (s)');
ylabel(sprintf('Velocity magnitude (%s)', units.displayVelocityUnit));
title('Band velocity magnitude');
legend(velLeg, 'Location','best');

subplot(2,3,5); hold on; grid on; box on;
plot(tSol, selectedDamage, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Mean damage D');
title('Selected rupture damage');
ylim([0 1]);

subplot(2,3,6); hold on; grid on; box on;
scatter3(nodeCoordsDisplay(:,1), nodeCoordsDisplay(:,2), nodeCoordsDisplay(:,3), 8, 'filled');
rng(1);
showK = min(3000, Ne);
pick = randperm(Ne, showK);
for kk = pick
    i = E(kk,1); j = E(kk,2);
    Xi = nodeCoordsDisplay(i,:); Xj = nodeCoordsDisplay(j,:);
    if bodyEdge(kk)
        plot3([Xi(1) Xj(1)], [Xi(2) Xj(2)], [Xi(3) Xj(3)], '-', ...
            'Color', [0.2 0.45 0.85]);
    end
    if damageEdge(kk)
        plot3([Xi(1) Xj(1)], [Xi(2) Xj(2)], [Xi(3) Xj(3)], '-', ...
            'Color', [0.85 0.05 0.05], 'LineWidth', 2);
    end
end
axis equal;
xlabel(sprintf('x (%s)', units.displayLengthUnit));
ylabel(sprintf('y (%s)', units.displayLengthUnit));
zlabel(sprintf('z (%s)', units.displayLengthUnit));
title('Body and damage edge masks');

[matDir, matBase] = fileparts(matFile);
pngFile = fullfile(matDir, [matBase '_dashboard.png']);
try
    exportgraphics(figDashboard, pngFile, 'Resolution', 300);
catch
    print(figDashboard, pngFile, '-dpng', '-r300');
end
fprintf('[POSTPROCESS] saved dashboard PNG: %s\n', pngFile);

end


function dXdt = timeGradientRows(X, t)
Nt = size(X,1);
dXdt = zeros(size(X));

if Nt < 2
    return;
end

dtForward = max(t(2) - t(1), eps);
dXdt(1,:) = (X(2,:) - X(1,:)) ./ dtForward;

dtBackward = max(t(end) - t(end-1), eps);
dXdt(end,:) = (X(end,:) - X(end-1,:)) ./ dtBackward;

if Nt > 2
    dtMid = max(t(3:end) - t(1:end-2), eps);
    dXdt(2:end-1,:) = (X(3:end,:) - X(1:end-2,:)) ./ dtMid;
end
end


function val = getParamOrDefault(param, fieldName, defaultVal)
if isfield(param, fieldName) && ~isempty(param.(fieldName))
    val = param.(fieldName);
else
    val = defaultVal;
end
end
