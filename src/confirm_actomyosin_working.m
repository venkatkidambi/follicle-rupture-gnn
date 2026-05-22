function confirm_actomyosin_working(matFile)
% confirm_actomyosin_working.m
% Postprocessing diagnostics for follicleShellSim_*.mat
%
% Key feature: motor masks exactly match SIM logic:
%   sNodes = sg * nodeCoords(:,ax)
%   sMid   = sg * edgeMid(:,ax)
%   sCutMain = prctile(sNodes, 80)
%   motorMask = (sMid < sCutMain), optionally filtered by edge orientation
%
% Outputs one dashboard figure with:
%   (1) posterior-band velocity vs time
%   (2) posterior-band velocity locations
%   (3) recruitment phase delay between rings
%   (4) posterior-band velocity magnitude vs time
%   (5) active rest-length shortening vs time
%   (6) motor region sanity (subset of edges with motors ON)

if nargin < 1 || isempty(matFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    matFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end

close all;

%% ---------------- USER SETTINGS ----------------
velocityRingBandsPct = [0 10; 25 35; 45 55]; % percent distance from posterior tip toward anterior
nVelocityRings = size(velocityRingBandsPct, 1);

%% ---------------- LOAD ----------------
S = load(matFile);

nodeCoords = S.nodeCoords;
edgeList   = S.edgeList;
tSol       = S.tSol(:);
Y          = S.Y;
param      = S.param;

% Motor region definition (MUST match sim)
mainPct = getParamOrDefault(param, 'mainPct', 80);  % posterior motor-off cutoff

N  = size(nodeCoords,1);
Nt = numel(tSol);
E  = edgeList;
Ne = size(E,1);

%% ---------------- POSITIONS OVER TIME ----------------
q   = reshape(Y, Nt, N, 3);                         % displacement
pos = q + reshape(nodeCoords, 1, N, 3);             % absolute positions

%% ---------------- MASKS (MATCH SIM EXACTLY) ----------------
edgeMid0 = 0.5*(nodeCoords(E(:,1),:) + nodeCoords(E(:,2),:));

% Must exist in param
if ~isfield(param,'postAxis') || ~isfield(param,'postSign')
    error('param.postAxis and/or param.postSign missing. Cannot reconstruct masks.');
end

ax = param.postAxis;          % 1/2/3
sg = param.postSign;          % +1 or -1

% Posterior score: larger = more posterior (by your definition)
sNodes = sg * nodeCoords(:,ax);
sMid   = sg * edgeMid0(:,ax);

sCutMain = prctile(sNodes, mainPct);

if isfield(param,'motorMask')
    motorEdge = logical(param.motorMask(:));
else
    motorEdge = (sMid <  sCutMain);
end

fprintf('[MASK] motorEdge fraction = %.3f\n', mean(motorEdge));

%% ---------------- ACTIVE REST-LENGTH SHORTENING OVER TIME ----------------
[activeFracTime, ~] = activeRestLengthShorteningPost(tSol, param, Ne);
if any(motorEdge)
    activeFracMotor = activeFracTime(:, motorEdge);
else
    activeFracMotor = activeFracTime;
end
meanActiveShrink = mean(activeFracMotor, 2);
maxActiveShrink = max(activeFracMotor, [], 2);

%% ---------------- POSTERIOR-BAND RING VELOCITY ----------------
% Group nodes by percent distance from posterior tip toward anterior.
% Velocity is along the anterior axis: positive means moving anteriorward.
posteriorToAnteriorNodeFrac = (max(sNodes) - sNodes) ./ max(range(sNodes), eps);
posteriorToAnteriorEdgeFrac = (max(sNodes) - sMid) ./ max(range(sNodes), eps);
ringBands = velocityRingBandsPct ./ 100;
ringLabels = cell(1, nVelocityRings);
ringIdNode = zeros(N,1);
for r = 1:nVelocityRings
    lo = ringBands(r,1);
    hi = ringBands(r,2);
    if r == nVelocityRings
        inRing = posteriorToAnteriorNodeFrac >= lo & posteriorToAnteriorNodeFrac <= hi;
    else
        inRing = posteriorToAnteriorNodeFrac >= lo & posteriorToAnteriorNodeFrac < hi;
    end
    ringIdNode(inRing) = r;
    ringLabels{r} = sprintf('%g-%g%%', velocityRingBandsPct(r,1), velocityRingBandsPct(r,2));
end

sPos = sg * squeeze(pos(:,:,ax));         % [Nt x N], larger = more posterior
vAnterior = -timeGradientRows(sPos, tSol); % [Nt x N], m/s along anterior axis
ringVelocity = nan(Nt, nVelocityRings);
ringCounts = zeros(1, nVelocityRings);
for r = 1:nVelocityRings
    inRing = ringIdNode == r;
    ringCounts(r) = nnz(inRing);
    if any(inRing)
        ringVelocity(:,r) = mean(vAnterior(:,inRing), 2);
    end
end
ringVelocityMag = abs(ringVelocity);

%% ---------------- RECRUITMENT PHASE DELAY DIAGNOSTIC ----------------
ringIdEdge = zeros(Ne,1);
for r = 1:nVelocityRings
    lo = ringBands(r,1);
    hi = ringBands(r,2);
    if r == nVelocityRings
        inRing = posteriorToAnteriorEdgeFrac >= lo & posteriorToAnteriorEdgeFrac <= hi;
    else
        inRing = posteriorToAnteriorEdgeFrac >= lo & posteriorToAnteriorEdgeFrac < hi;
    end
    ringIdEdge(inRing) = r;
end

if isfield(S,'A_time') && isequal(size(S.A_time), [Nt, Ne])
    recruitmentTime = min(max(S.A_time, 0), 1);
    recruitmentThreshold = getParamOrDefault(param, 'recruitmentThreshold', 0.3);
    recruitmentLabel = 'A_{edge}';
else
    maxShortening = max(getParamOrDefault(param, 'activeShorteningMax', 1), eps);
    recruitmentTime = min(max(activeFracTime ./ maxShortening, 0), 1);
    recruitmentThreshold = 0.3;
    recruitmentLabel = 'normalized shortening';
end

meanRecruitmentRing = nan(Nt, nVelocityRings);
for r = 1:nVelocityRings
    inRing = ringIdEdge == r;
    if any(inRing)
        meanRecruitmentRing(:,r) = mean(recruitmentTime(:,inRing), 2);
    end
end

recruitmentOnset = nan(1, nVelocityRings);
for r = 1:nVelocityRings
    idx = find(meanRecruitmentRing(:,r) >= recruitmentThreshold, 1, 'first');
    if ~isempty(idx)
        recruitmentOnset(r) = tSol(idx);
    end
end
recruitmentDelay = recruitmentOnset - recruitmentOnset(1);

%% =========================
% ===== FIGURE DASHBOARD ====
% =========================
figDashboard = figure('Color','w','Position',[100 100 1500 850]);

ringColors = lines(nVelocityRings);

% (1) posterior-band velocity
subplot(2,3,1); hold on; grid on; box on;
velLeg = cell(1, nVelocityRings);
for r = 1:nVelocityRings
    plot(tSol, ringVelocity(:,r)*1e6, 'LineWidth', 1.8, 'Color', ringColors(r,:));
    velLeg{r} = sprintf('%s, n=%d', ringLabels{r}, ringCounts(r));
end
yline(0, ':');
if isfield(param,'tActOn'), xline(param.tActOn,'--','tActOn'); end
xlabel('Time (s)');
ylabel('Anterior-axis velocity (\mum/s)');
title('Band velocity: anteriorward positive');
legend(velLeg, 'Location','best');

% (2) velocity ring locations
subplot(2,3,2); hold on; grid on; box on;
otherNode = ringIdNode == 0;
scatter3(nodeCoords(otherNode,1)*1e6, nodeCoords(otherNode,2)*1e6, nodeCoords(otherNode,3)*1e6, ...
    8, [0.78 0.78 0.78], 'filled');
locLeg = {'anterior half'};
for r = 1:nVelocityRings
    inRing = ringIdNode == r;
    scatter3(nodeCoords(inRing,1)*1e6, nodeCoords(inRing,2)*1e6, nodeCoords(inRing,3)*1e6, ...
        16, ringColors(r,:), 'filled');
    locLeg{end+1} = ringLabels{r}; %#ok<AGROW>
end
axis equal; view(3);
xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
title('Velocity bands: percent from posterior tip');
legend(locLeg, 'Location','bestoutside');

% (3) recruitment phase delay between rings
subplot(2,3,3); hold on; grid on; box on;
for r = 1:nVelocityRings
    plot(tSol, meanRecruitmentRing(:,r), 'LineWidth', 1.8, 'Color', ringColors(r,:));
end
yline(recruitmentThreshold, ':', 'threshold');
if isfield(param,'tActOn'), xline(param.tActOn,'--','tActOn'); end
xlabel('Time (s)');
ylabel(sprintf('Mean %s', recruitmentLabel));
title('Recruitment phase delay by ring');
phaseLeg = cell(1, nVelocityRings);
for r = 1:nVelocityRings
    if isfinite(recruitmentDelay(r))
        phaseLeg{r} = sprintf('%s, delay %.1fs', ringLabels{r}, recruitmentDelay(r));
    else
        phaseLeg{r} = sprintf('%s, no crossing', ringLabels{r});
    end
end
legend(phaseLeg, 'Location','best');

% (4) posterior-band velocity magnitude
subplot(2,3,4); hold on; grid on; box on;
for r = 1:nVelocityRings
    plot(tSol, ringVelocityMag(:,r)*1e6, 'LineWidth', 1.8, 'Color', ringColors(r,:));
end
yline(0, ':');
if isfield(param,'tActOn'), xline(param.tActOn,'--','tActOn'); end
xlabel('Time (s)');
ylabel('Velocity magnitude (\mum/s)');
title('Band velocity magnitude');
legend(velLeg, 'Location','best');

% (5) active rest-length shortening
subplot(2,3,5); hold on; grid on; box on;
if isfield(param,'activeRestLength_on') && param.activeRestLength_on
    plot(tSol, 100*meanActiveShrink, 'LineWidth',2);
    plot(tSol, 100*maxActiveShrink, '--', 'LineWidth',1.5);
    if isfield(param,'tActOn'), xline(param.tActOn,'--','tActOn'); end
    xlabel('Time (s)'); ylabel('Rest-length shortening (%)');
    title('Active preferred-length shortening');
    legend('mean selected edges','max selected edge','Location','best');
else
    text(0.1,0.5,'Active rest-length shortening is disabled.', ...
         'Units','normalized','FontSize',11);
    axis off;
end

% (6) initial geometry + motor edges (visual sanity)
subplot(2,3,6); hold on; grid on; box on;
scatter3(nodeCoords(:,1)*1e6, nodeCoords(:,2)*1e6, nodeCoords(:,3)*1e6, 8, 'filled');
rng(1);
showK = min(3000, Ne);
pick = randperm(Ne, showK);
for kk = pick
    i = E(kk,1); j = E(kk,2);
    Xi = nodeCoords(i,:); Xj = nodeCoords(j,:);
    if motorEdge(kk)
        plot3([Xi(1) Xj(1)]*1e6, [Xi(2) Xj(2)]*1e6, [Xi(3) Xj(3)]*1e6, '-');
    end
end
axis equal;
xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
title('Motor region sanity: edges with motors ON (subset)');

%% ---------------- EXPORT DASHBOARD PNG ----------------
[matDir, matBase] = fileparts(matFile);
pngFile = fullfile(matDir, [matBase '_dashboard.png']);
try
    exportgraphics(figDashboard, pngFile, 'Resolution', 300);
catch
    print(figDashboard, pngFile, '-dpng', '-r300');
end
fprintf('[POSTPROCESS] saved dashboard PNG: %s\n', pngFile);

end


function [activeFrac, activeFracDot] = activeRestLengthShorteningPost(tSol, param, Ne)
Nt = numel(tSol);
activeFrac = zeros(Nt, Ne);
activeFracDot = zeros(Nt, Ne);

if ~isfield(param,'activeRestLength_on') || ~param.activeRestLength_on
    return;
end

if isfield(param,'tActiveRestLengthOn_edge') && numel(param.tActiveRestLengthOn_edge) == Ne
    tOn = param.tActiveRestLengthOn_edge(:).';
else
    tOn = getParamOrDefault(param, 'tActOn', 0) * ones(1, Ne);
end

activationMode = lower(strtrim(char(getParamOrDefault(param, 'activeShorteningActivation', 'exponential'))));

switch activationMode
    case 'exponential'
        tau = max(getParamOrDefault(param, 'activeShorteningTau', 1), eps);
        dtOn = tSol - tOn;
        isOn = dtOn >= 0;
        A = zeros(Nt, Ne);
        Adot = zeros(Nt, Ne);
        A(isOn) = 1 - exp(-dtOn(isOn) ./ tau);
        Adot(isOn) = exp(-dtOn(isOn) ./ tau) ./ tau;
    case {'accelerating', 'power', 'reverse-exponential'}
        tau = max(getParamOrDefault(param, 'activeShorteningTau', 1), eps);
        p = max(getParamOrDefault(param, 'activeShorteningAccelPower', 2), 1);
        dtOn = tSol - tOn;
        u = max(dtOn ./ tau, 0);
        isRamping = u < 1;
        A = min(u .^ p, 1);
        Adot = zeros(Nt, Ne);
        Adot(isRamping & u > 0) = (p ./ tau) .* u(isRamping & u > 0).^(p - 1);
    case {'sigmoid', 'delayed-sigmoid'}
        delay = getParamOrDefault(param, 'activeShorteningSigmoidDelay', 0);
        sigmoidTau = max(getParamOrDefault(param, 'activeShorteningSigmoidTau', ...
            getParamOrDefault(param, 'activeShorteningTau', 1)), eps);
        tMid = tOn + delay;
        A = 1 ./ (1 + exp(-(tSol - tMid) ./ sigmoidTau));
        Adot = (A .* (1 - A)) ./ sigmoidTau;
        tau = sigmoidTau;
    otherwise
        error('Unknown activeShorteningActivation "%s". Use "exponential", "sigmoid", or "accelerating".', activationMode);
end

if isfield(param,'activeRestLengthWavePulse_on') && param.activeRestLengthWavePulse_on && ...
        isfield(param,'tActiveRestLengthOff_edge') && numel(param.tActiveRestLengthOff_edge) == Ne
    tOff = param.tActiveRestLengthOff_edge(:).';
    offTau = tau;
    if isfield(param,'activeRestLengthWaveOffTau') && isfinite(param.activeRestLengthWaveOffTau) && param.activeRestLengthWaveOffTau > 0
        offTau = param.activeRestLengthWaveOffTau;
    end
    Aoff = 1 ./ (1 + exp((tSol - tOff) ./ max(offTau, eps)));
    AoffDot = -(Aoff .* (1 - Aoff)) ./ max(offTau, eps);
    Adot = Adot .* Aoff + A .* AoffDot;
    A = A .* Aoff;
end

if isfield(param,'activeShorteningMask') && ~isempty(param.activeShorteningMask)
    mask = param.activeShorteningMask(:).';
else
    mask = ones(1, Ne);
end

if isscalar(mask)
    mask = mask * ones(1, Ne);
elseif numel(mask) ~= Ne
    error('param.activeShorteningMask must be scalar or have one value per edge.');
end

maxShortening = getParamOrDefault(param, 'activeShorteningMax', 0);
rawFrac = maxShortening * A .* mask;
rawFracDot = maxShortening * Adot .* mask;

activeFrac = min(max(rawFrac, 0), 0.95);
activeFracDot = rawFracDot;
activeFracDot(rawFrac <= 0 | rawFrac >= 0.95) = 0;
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
