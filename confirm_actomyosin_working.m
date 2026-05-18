function confirm_actomyosin_working(matFile)
% confirm_actomyosin_working.m
% Postprocessing + rupture localization diagnostics for follicleShellSim_*.mat
%
% Key feature: tip/motor masks exactly match SIM logic:
%   sNodes = sg * nodeCoords(:,ax)
%   sMid   = sg * edgeMid(:,ax)
%   sCutMain = prctile(sNodes, 80)
%   sCutTip  = prctile(sNodes, 92)
%   motorMask = (sMid < sCutMain)
%   tipMask   = (sMid >= sCutTip)
%
% Outputs:
%   (1) mean radius vs time (all/main/post)
%   (2) active tension ramp
%   (3) motor region sanity (subset of edges with motors ON)
%   (4) initial vs final node cloud
%   (5) rupture strain diagnostic (tip edges): max + 95th percentile vs eps_break
%   (6) rupture progression: fraction tip edges over eps_break and/or D > D_break
%       + tip opening proxy (RMS spread of tip nodes), right axis
%   (7) WHERE rupture happened: final-time shell with broken edges highlighted by break-time bin

if nargin < 1 || isempty(matFile)
    matFile = 'follicleShellSim_test2.mat';
end

close all;

%% ---------------- USER SETTINGS ----------------
eps_break = 0.25;
D_break   = 0.99;     % current sim saves damage D, where D=1 means failed

% Tip definition (MUST match sim)
mainPct = 80;         % posterior-most 20% has motors OFF
tipPct  = 92;         % posterior-most 8% defines tip (EDGE mask)

% Rupture map binning (static visualization)
nbins = 12;           % number of time bins to color broken edges (early->late)

%% ---------------- LOAD ----------------
S = load(matFile);

nodeCoords = S.nodeCoords;
edgeList   = S.edgeList;
tSol       = S.tSol(:);
Y          = S.Y;
param      = S.param;

N  = size(nodeCoords,1);
Nt = numel(tSol);
E  = edgeList;
Ne = size(E,1);

hasD = isfield(S,'D_time');
if hasD
    D_time = S.D_time;      % [Nt x Ne], 0=intact and 1=failed
else
    D_time = [];
    warning('D_time not found in %s; rupture plots will fall back to strain threshold.', matFile);
end

%% ---------------- POSITIONS OVER TIME ----------------
q   = reshape(Y, Nt, N, 3);                         % displacement
pos = q + reshape(nodeCoords, 1, N, 3);             % absolute positions
center = squeeze(median(pos,2));                    % robust center (Nt x 3)

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
sCutTip  = prctile(sNodes, tipPct);

motorEdge = (sMid <  sCutMain);     % motors ON in main-body (non-posterior)
tipEdge   = (sMid >= sCutTip);      % tip edges (posterior cap)

% Tip nodes for opening proxy + highlighting: take nodes in posterior-most tipPct%
tipNodeCut = prctile(sNodes, tipPct);
idxTipNodes = find(sNodes >= tipNodeCut);

fprintf('[MASK] motorEdge fraction = %.3f | tipEdge fraction = %.3f\n', ...
    mean(motorEdge), mean(tipEdge));

%% ---------------- MAIN vs POSTERIOR NODE GROUPS (for radius curves) ----------------
% Use same axis-consistent node mask (posterior-most 20% vs rest)
idxPost = find(sNodes >= sCutMain);     % posterior cap nodes
idxMain = find(sNodes <  sCutMain);     % main body nodes

%% ---------------- MEAN RADIUS CURVES ----------------
rMean = zeros(Nt,1);
rPost = zeros(Nt,1);
rMain = zeros(Nt,1);

for k = 1:Nt
    Pk = squeeze(pos(k,:,:));
    ck = center(k,:);
    r  = vecnorm(Pk - ck, 2, 2);
    rMean(k) = mean(r);
    rPost(k) = mean(r(idxPost));
    rMain(k) = mean(r(idxMain));
end

%% ---------------- EDGE STRAIN OVER TIME ----------------
i1 = E(:,1); i2 = E(:,2);
L_time = zeros(Nt, Ne);

for k = 1:Nt
    Pk = squeeze(pos(k,:,:));              % [N x 3]
    d  = Pk(i2,:) - Pk(i1,:);              % [Ne x 3]
    L_time(k,:) = sqrt(sum(d.^2,2)).';
end

if isfield(param,'L0_edge0')
    L0_init = param.L0_edge0(:).';
else
    L0_init = sqrt(sum((nodeCoords(i2,:) - nodeCoords(i1,:)).^2,2)).';
end
eps_time = (L_time - L0_init) ./ max(L0_init, 1e-18);

%% ---------------- RUPTURE METRICS (TIP ONLY) ----------------
eps_tip = eps_time(:, tipEdge);
maxEpsTip = max(eps_tip, [], 2);
p95EpsTip = prctile(eps_tip, 95, 2);
fracTipOverBreak = mean(eps_tip > eps_break, 2);

if hasD
    D_tip = D_time(:, tipEdge);
    fracTipBroken = mean(D_tip > D_break, 2);
else
    fracTipBroken = fracTipOverBreak;
end

% Tip opening proxy: RMS spread of tip nodes around their centroid
openR = zeros(Nt,1);
for k = 1:Nt
    Pk = squeeze(pos(k,:,:));
    Tn = Pk(idxTipNodes,:);
    ct = mean(Tn,1);
    openR(k) = sqrt(mean(sum((Tn-ct).^2,2)));
end

% Summary prints
if hasD
    breakMask_tip = (D_time(:,tipEdge) > D_break);
else
    breakMask_tip = (eps_time(:,tipEdge) > eps_break);
end
kFirst = find(any(breakMask_tip,2), 1, 'first');
if ~isempty(kFirst)
    fprintf('[RUPTURE] First TIP edge breaks at t = %.4f s\n', tSol(kFirst));
else
    fprintf('[RUPTURE] No TIP edges met rupture criterion.\n');
end
fprintf('[RUPTURE] End-state: fracTipOverBreak=%.3f | fracTipBroken=%.3f\n', ...
    fracTipOverBreak(end), fracTipBroken(end));

%% ---------------- BREAK TIME PER EDGE (WHERE) ----------------
tBreak = nan(1,Ne);
if hasD
    breakMask_t = (D_time > D_break);          % [Nt x Ne]
else
    breakMask_t = (eps_time > eps_break);      % [Nt x Ne]
end

for e = 1:Ne
    k0 = find(breakMask_t(:,e), 1, 'first');
    if ~isempty(k0)
        tBreak(e) = tSol(k0);
    end
end

%% =========================
% ===== FIGURE (3 x 2) =====
% =========================
figure('Color','w','Position',[100 100 1200 950]);

% (1) radius vs time
subplot(3,2,1); hold on; grid on; box on;
plot(tSol, rMean*1e6, 'LineWidth',2);
plot(tSol, rMain*1e6, 'LineWidth',2);
plot(tSol, rPost*1e6, 'LineWidth',2);
xlabel('Time (s)'); ylabel('Mean radius (\mum)');
title('Contraction check: mean radius vs time');
legend('All nodes','Main-body','Posterior cap','Location','best');
if isfield(param,'tActOn'), xline(param.tActOn,'--','tActOn'); end

% (2) active tension ramp
subplot(3,2,2); hold on; grid on; box on;
if isfield(param,'act_on') && param.act_on && isfield(param,'tActOn') && isfield(param,'tActTau')
    A = 1 ./ (1 + exp(-(tSol - param.tActOn)/max(param.tActTau,eps)));
    if isfield(param,'Tact_main')
        plot(tSol, A * param.Tact_main * 1e6, 'LineWidth',2);
    end
    if isfield(param,'Tact_tip') && param.Tact_tip ~= 0
        plot(tSol, A * param.Tact_tip * 1e6, 'LineWidth',2);
        legend('main active tension','tip active tension','Location','best');
    else
        legend('main active tension','Location','best');
    end
    if isfield(param,'tActOn'), xline(param.tActOn,'--','tActOn'); end
    xlabel('Time (s)'); ylabel('Active tension (\muN)');
    title('Actomyosin activation ramp');
else
    text(0.1,0.5,'Active tension parameters not saved.', ...
         'Units','normalized','FontSize',11);
    axis off;
end

% (3) initial geometry + motor edges (visual sanity)
subplot(3,2,3); hold on; grid on; box on;
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

% (4) initial vs final node cloud
subplot(3,2,4); hold on; grid on; box on;
Pk0 = squeeze(pos(1,:,:));
Pk1 = squeeze(pos(end,:,:));
scatter3(Pk0(:,1)*1e6,Pk0(:,2)*1e6,Pk0(:,3)*1e6,8,'filled');
scatter3(Pk1(:,1)*1e6,Pk1(:,2)*1e6,Pk1(:,3)*1e6,8,'filled');
axis equal;
xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
title('Initial vs final node cloud');
legend('t=0','t=end','Location','best');

% (5) rupture strain diagnostic (tip edges)
subplot(3,2,5); hold on; grid on; box on;
plot(tSol, maxEpsTip, 'LineWidth', 2);
plot(tSol, p95EpsTip, 'LineWidth', 2);
yline(eps_break, '--', 'eps\_break');
xlabel('Time (s)'); ylabel('Tip tensile strain');
title('Rupture diagnostic: tip-edge strain');
legend('max tip strain','95th % tip strain','eps\_break','Location','best');

% (6) rupture progression + opening proxy
subplot(3,2,6); hold on; grid on; box on;
yyaxis left
plot(tSol, fracTipOverBreak, 'LineWidth', 2);
plot(tSol, fracTipBroken, 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Fraction of tip edges');
title('Rupture progression + tip opening proxy');
if hasD
    legend('strain > eps\_break','D > D\_break','Location','best');
else
    legend('strain > eps\_break','(fallback) same curve','Location','best');
end
yyaxis right
plot(tSol, openR*1e6, 'LineWidth', 2);
ylabel('Tip opening proxy (\mum)');

%% =========================
% (7) WHERE rupture happened (STATIC MAP at final time)
% =========================
Pk = squeeze(pos(end,:,:));
ck = center(end,:);
Pk = Pk - ck;

% edge subset for clarity
rng(7);
showEdges = min(9000, Ne);
pickE = randperm(Ne, showEdges);

isBroken = ~isnan(tBreak);
brokenIdx = pickE(isBroken(pickE));
intactIdx = pickE(~isBroken(pickE));

figure('Color','w','Position',[150 120 1100 850]);
axh = axes; hold(axh,'on'); grid(axh,'on'); box(axh,'on');
axis(axh,'equal'); view(axh,3); rotate3d(axh,'on');
xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
title('Rupture location map (final time): broken edges grouped by break time');

% nodes + tip nodes
scatter3(Pk(:,1)*1e6,Pk(:,2)*1e6,Pk(:,3)*1e6,8,'filled');
Tn = Pk(idxTipNodes,:);
scatter3(Tn(:,1)*1e6,Tn(:,2)*1e6,Tn(:,3)*1e6,25,'filled');

% intact edges (batched)
[Xint,Yint,Zint] = batchEdges(Pk, E, intactIdx);
plot3(Xint*1e6, Yint*1e6, Zint*1e6, '-', 'LineWidth',0.5);

% broken edges binned by break time
if ~isempty(brokenIdx)
    tb = tBreak(brokenIdx);
    tmin = min(tb); tmax = max(tb);
    if tmax == tmin, tmax = tmin + eps; end
    bins = linspace(tmin, tmax, nbins+1);

    for b = 1:nbins
        inBin = brokenIdx(tb >= bins(b) & tb < bins(b+1));
        if isempty(inBin), continue; end
        [Xb,Yb,Zb] = batchEdges(Pk, E, inBin);
        plot3(Xb*1e6, Yb*1e6, Zb*1e6, '-', 'LineWidth',2.0);
    end
end

legend('All nodes','Tip nodes','Intact edges','Broken edges (time-binned)','Location','best');

% quick print: where did the earliest breaks occur?
if any(isBroken)
    tEarliest = min(tBreak(isBroken));
    earlyEdges = find(tBreak == tEarliest);
    fprintf('[WHERE] earliest break time = %.4f s | num edges = %d\n', tEarliest, numel(earlyEdges));
    fprintf('[WHERE] earliest breaks: tip fraction = %.3f\n', mean(tipEdge(earlyEdges)));
end

%% =========================
% (9) 3D node trajectories (spaghetti plot)
% =========================
% Controls
trajNumNodes = 300;        % how many nodes to plot (keep modest)
trajStrideT  = 3;          % plot every k-th timepoint (speed)
centerTraj   = true;       % subtract center(t) so you see shape motion
trajMicrons  = true;       % plot in microns

rng(10);
N = size(pos,2);
Nt = size(pos,1);

idxPlot = randperm(N, min(trajNumNodes, N));
tt = 1:trajStrideT:Nt;

figure('Color','w','Position',[120 120 1050 850]);
hAx = axes; hold(hAx,'on'); grid(hAx,'on'); box(hAx,'on');
axis(hAx,'equal'); view(hAx,3); rotate3d(hAx,'on');
xlabel(hAx,'x'); ylabel(hAx,'y'); zlabel(hAx,'z');
title(hAx, sprintf('3D node trajectories (N=%d nodes, stride=%d)', numel(idxPlot), trajStrideT));

for ii = 1:numel(idxPlot)
    nID = idxPlot(ii);

    P = squeeze(pos(tt, nID, :));              % [T x 3]
    if centerTraj
        C = center(tt,:);                      % [T x 3]
        P = P - C;
    end

    if trajMicrons
        P = P * 1e6;
    end

    plot3(hAx, P(:,1), P(:,2), P(:,3), '-', 'LineWidth', 0.7);
end


%% =========================
% (10) Trajectories with start/end markers (cleaner)
% =========================
trajNumNodes = 250;
trajStrideT  = 4;
centerTraj   = true;

rng(11);
idxPlot = randperm(N, min(trajNumNodes, N));
tt = 1:trajStrideT:Nt;

figure('Color','w','Position',[160 160 1100 850]);
hAx = axes; hold(hAx,'on'); grid(hAx,'on'); box(hAx,'on');
axis(hAx,'equal'); view(hAx,3); rotate3d(hAx,'on');
xlabel(hAx,'x (\mum)'); ylabel(hAx,'y (\mum)'); zlabel(hAx,'z (\mum)');
title(hAx,'3D trajectories (faint trails) + start/end markers');

for ii = 1:numel(idxPlot)
    nID = idxPlot(ii);

    P = squeeze(pos(tt, nID, :));
    if centerTraj
        P = P - center(tt,:);
    end
    P = P * 1e6;

    plot3(hAx, P(:,1), P(:,2), P(:,3), '-', 'LineWidth', 0.5);  % trail
end

% Start + end positions for ALL nodes (gives nice context)
P0 = squeeze(pos(1,:,:));
P1 = squeeze(pos(end,:,:));
if centerTraj
    P0 = P0 - center(1,:);
    P1 = P1 - center(end,:);
end
P0 = P0*1e6; P1 = P1*1e6;

scatter3(hAx, P0(:,1), P0(:,2), P0(:,3), 8, 'filled');   % start
scatter3(hAx, P1(:,1), P1(:,2), P1(:,3), 14, 'filled');  % end
legend(hAx, {'Trails (subset)','Start (all)','End (all)'}, 'Location','best');

end


%% -------- helper: batch edges into one NaN-separated line --------
function [X,Y,Z] = batchEdges(P, E, idx)
    % P: [N x 3], E: [Ne x 2], idx: edge indices
    if isempty(idx)
        X = nan; Y = nan; Z = nan; return;
    end
    ii = E(idx,1); jj = E(idx,2);
    X = [P(ii,1) P(jj,1) nan(numel(idx),1)].';
    Y = [P(ii,2) P(jj,2) nan(numel(idx),1)].';
    Z = [P(ii,3) P(jj,3) nan(numel(idx),1)].';
    X = X(:); Y = Y(:); Z = Z(:);
end
