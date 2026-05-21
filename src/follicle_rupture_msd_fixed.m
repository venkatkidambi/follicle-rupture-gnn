function follicle_rupture_msd_fixed
% Follicle rupture shell model
% State: X = [q; v; sigma; zeta; D]
%   q, v        node displacement and velocity
%   sigma,zeta  viscoelastic edge states
%   D           edge damage in [0,1], where 1 means failed

clc; close all;

rng(42);

%% USER-PARAMETERS
N = 200;
mVal = 1;

R = 125e-6;

% Baseline edge parameters
kOuter_range = [1, 1]; 
cOuter_range = [10, 10];

% Time
tEnd  = 60;                    % s
numSave = 300;                 % saved output frames; solver still adapts internally
runPostprocess = true;        % set true to generate diagnostic plots after saving
odeSolver = 'ode45';           % try 'ode15s' or 'ode23t' if this system behaves stiffly
progressEveryPhysSec = 2;      % console progress interval in physical seconds

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(repoRoot, 'results');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
simFile = fullfile(outputDir, 'follicleShellSim_test2.mat');

%% Nodes (ELLIPSOID SURFACE: "football")
a = 200e-6; b = 80e-6; c = 80e-6;
nodeCoords = genEllipsoidSurfaceNodes(N, a, b, c, [0 0 0]);

param.ell_a = a; param.ell_b = b; param.ell_c = c;

% Posterior axis convention
param.postAxis = 1;
param.postSign = -1;

% Surface penalty + normal damping
param.k_surf = 5e-6;     % [N] penalty strength
param.c_surf = 5e-2;     % [N*s/m] normal damping
param.c_drag = 1e-3;     % [N*s/m] whole-node damping
param.projectEdgeForcesToSurface = true;

%% Surface edges (shell-only)
DT = delaunayTriangulation(nodeCoords);
F = freeBoundary(DT);   % [nFaces x 3]
allEdges = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
allEdges = sort(allEdges,2);
edgeList = unique(allEdges,'rows');
Ne = size(edgeList,1);

if runPostprocess
    figure; scatter3(nodeCoords(:,1),nodeCoords(:,2),nodeCoords(:,3),10,'filled');
    axis equal; grid on; title('Ellipsoid surface nodes');
end

param.N = N;
param.nodeCoords = nodeCoords;
param.edgeList = edgeList;
param.m = mVal;
param.Ne = Ne;

% Fixed geometric rest lengths per edge
E = edgeList;
param.edge_i = E(:,1);
param.edge_j = E(:,2);
param.edgeAccumI = sparse(param.edge_i, (1:Ne)',  1, N, Ne);
param.edgeAccumJ = sparse(param.edge_j, (1:Ne)',  1, N, Ne);
param.L0_edge0 = sqrt(sum((nodeCoords(E(:,2),:) - nodeCoords(E(:,1),:)).^2, 2));

% Baseline per-edge Burgers parameters
param.k1_edge = kOuter_range(1) * ones(Ne,1);
param.k2_edge = kOuter_range(2) * ones(Ne,1);
param.c_edge  = cOuter_range(1) * ones(Ne,1);
param.c_series_edge = cOuter_range(2) * ones(Ne,1);

% Passive 1D connection model between cells.
% Options: 'spring', 'damper', 'maxwell', 'kelvin-voigt', 'zener', 'burgers'.
% You may also set param.connectionModel_edge to a per-edge cell/string array.
param.connectionModel = 'maxwell';
param.connectionModelKey = normalizeConnectionModel(param.connectionModel);
param.connectionModelId_edge = connectionModelIds(param, Ne);

%% Spatial localization masks (main-body vs posterior cap)
edgeMid = 0.5*(nodeCoords(E(:,1),:) + nodeCoords(E(:,2),:) );
ax = param.postAxis;
sg = param.postSign;

sNodes = sg * nodeCoords(:,ax);
sMid   = sg * edgeMid(:,ax);

sCutMain = prctile(sNodes, 80);    % posterior-most ~20% has motors OFF in main mask
sCutTip  = prctile(sNodes, 92);    % posterior-most ~8% defines tip

param.motorMask = double(sMid < sCutMain);    % motors ON in main-body
param.tipMask   = double(sMid >= sCutTip);    % posterior cap edges

param.damageMode = 'seam';    % 'seam' or 'tip'
param.seamWidthFrac = 0.04;   % fraction of posterior-score range

switch lower(param.damageMode)
    case 'seam'
        seamWidth = param.seamWidthFrac * range(sNodes);
        param.damageMask = double(abs(sMid - sCutTip) <= seamWidth);
    case 'tip'
        param.damageMask = param.tipMask;
    otherwise
        error('Unknown damageMode "%s". Use "seam" or "tip".', param.damageMode);
end

%% Active tension parameters
param.act_on      = true;
param.tActOn      = 0.5;     % [s]
param.tActTau     = 0.1;     % [s]

param.Tact_main   = 2e-6;   % [N] main-body active edge tension (TUNE)
param.Tact_tip    = 0;       % [N] usually 0

param.axisSelectiveActomyosin = true;
param.axisAlignmentMin = 0.65;  % 1=perfectly parallel, 0=no orientation filter

if param.axisSelectiveActomyosin
    edgeVec = nodeCoords(E(:,2),:) - nodeCoords(E(:,1),:);
    edgeLen = sqrt(sum(edgeVec.^2, 2));
    axisAlign = abs(edgeVec(:,param.postAxis)) ./ max(edgeLen, 1e-18);

    param.axisAlignment = axisAlign;
    param.motorMask = param.motorMask .* double(axisAlign >= param.axisAlignmentMin);
end

%% Damage / rupture parameters
param.tipDamage_on  = true;

param.sig_fail_tip  = 2e-8; % [N] tip failure threshold (TUNE)
param.tau_dmg       = 0.05;   % [s] damage timescale (TUNE)
param.dmg_pow       = 1;     % sharpness

param.eps_gate_on   = true;  % only damage if opening
param.eps_gate      = 0.01;

% Optional: gentle tip priming (adds small damage even before overload)
param.tipPriming_on   = true;
param.tipPrimingRate  = 0.02;  % [1/s]
param.tipPriming_on_t = 0.5;   % [s]
param.tipPriming_tau  = 0.3;   % [s]

% Minimum load-bearing fraction to avoid numerical singularities
param.g_floor_stiff = 0;        % allow fully failed
param.g_floor_damp  = 1e-6;     % keep a tiny damping floor

%% Nondimensionalization
L_scale  = R;
k1_scale = mean(param.k1_edge);
k2_scale = mean(param.k2_edge);
c0_scale = mean(param.c_edge);
k0_scale = k1_scale + k2_scale;
tau_scale = c0_scale / max(k2_scale,1e-18);
epsilon   = mVal * k0_scale / max(c0_scale^2,1e-18);

param.scale.L       = L_scale;
param.scale.k2      = k2_scale;
param.scale.tau     = tau_scale;
param.scale.epsilon = epsilon;
param.numSave = numSave;
param.runPostprocess = runPostprocess;
param.odeSolver = odeSolver;
param.progressEveryPhysSec = progressEveryPhysSec;

%% ================= State sizes =================
nPos = 3*N;
nVel = 3*N;

param.nPos = nPos;
param.nVel = nVel;

% X = [q; v; sigma; zeta; D]
nX = nPos + nVel + Ne + Ne + Ne;

% Damage block indices for extracting D_time after integration
iD1 = nPos + nVel + 2*Ne + 1;
iD2 = iD1 + Ne - 1;

%% Integrate
tFinal_dim = tEnd / tau_scale;
assert(isscalar(tFinal_dim) && isfinite(tFinal_dim));

X0_dim = zeros(nX,1);

opts = odeset( ...
    'RelTol', 1e-2, ...
    'AbsTol', 1e-4, ...
    'MaxStep', 5, ...
    'InitialStep', 1e-3, ...
    'OutputFcn', @(t_dim, y, flag) myProgressOutputFcn(t_dim, y, flag, tau_scale, tFinal_dim, param.progressEveryPhysSec), ...
    'Stats', 'on');

tDimQ_all = linspace(0, tFinal_dim, numSave);

odeRHS = @(t,X) nondimRHS_shell(t, X, param);
switch lower(param.odeSolver)
    case 'ode45'
        [tDimQ_all, XDimQ_all] = ode45(odeRHS, tDimQ_all, X0_dim, opts);
    case 'ode15s'
        [tDimQ_all, XDimQ_all] = ode15s(odeRHS, tDimQ_all, X0_dim, opts);
    case 'ode23t'
        [tDimQ_all, XDimQ_all] = ode23t(odeRHS, tDimQ_all, X0_dim, opts);
    otherwise
        error('Unknown odeSolver "%s". Use "ode45", "ode15s", or "ode23t".', param.odeSolver);
end

tSol  = tDimQ_all * tau_scale;
Y_dim = XDimQ_all(:, 1:nPos);
Y     = Y_dim * L_scale;

% extract damage over time
D_time = XDimQ_all(:, iD1:iD2);

save(simFile, ...
     'nodeCoords','edgeList','tSol','Y','D_time','param','-v7.3');

if param.runPostprocess
    try
        confirm_actomyosin_working(simFile);
    catch ME
        warning('follicle_rupture:postprocessFailed', ...
            'Postprocessing failed (%s): %s', ME.identifier, ME.message);
    end
end

end % main


%% ===================== Nondimensional wrapper RHS =====================
function dXdim = nondimRHS_shell(t_dim, X_dim, param)
nPos = param.nPos;
nVel = param.nVel;
Ne   = param.Ne;
scale = param.scale;

% indices in ND vector
iQ2  = nPos;
iV2  = iQ2 + nVel;
iS2  = iV2 + Ne;
iZ2  = iS2 + Ne;
iD2  = iZ2 + Ne;

q_dim   = X_dim(1:iQ2);
v_dim   = X_dim(iQ2+1:iV2);
sig_dim = X_dim(iV2+1:iS2);
zet_dim = X_dim(iS2+1:iZ2);
D_dim   = X_dim(iZ2+1:iD2);

% to physical
q_phys   = q_dim * scale.L;
v_phys   = v_dim * (scale.L/scale.tau);
sig_phys = sig_dim * (scale.k2*scale.L);
zet_phys = zet_dim * scale.L;
D_phys   = D_dim;
t_phys   = t_dim * scale.tau;

dX_phys = ShellRHS_physical(t_phys, [q_phys; v_phys; sig_phys; zet_phys; D_phys], param);

% back to ND
dqdt_dim   = (dX_phys(1:nPos)/scale.L) * scale.tau;
dvdt_dim   = scale.epsilon * ((dX_phys(nPos+1:nPos+nVel) / (scale.L/scale.tau)) * scale.tau);

dsig_dim   = (dX_phys(nPos+nVel+1 : nPos+nVel+Ne) / (scale.k2*scale.L)) * scale.tau;
dzet_dim   = (dX_phys(nPos+nVel+Ne+1 : nPos+nVel+2*Ne) / scale.L) * scale.tau;

off = nPos + nVel + 2*Ne;
dD_dim     = dX_phys(off+(1:Ne)) * scale.tau;

dXdim = [dqdt_dim; dvdt_dim; dsig_dim; dzet_dim; dD_dim];
end


%% ========================= PHYSICAL RHS =========================
function dxdt = ShellRHS_physical(t, x, param)
% X = [ q(3N); v(3N); sigma(Ne); zeta(Ne); D(Ne) ]

N    = param.N;
nPos = param.nPos;
nVel = param.nVel;
Ne   = param.Ne;

% Indices
iQ1 = 1;                 iQ2 = nPos;
iV1 = iQ2+1;             iV2 = iV1 + nVel - 1;
iS1 = iV2+1;             iS2 = iS1 + Ne - 1;
iZ1 = iS2+1;             iZ2 = iZ1 + Ne - 1;
iD1 = iZ2+1;             iD2 = iD1 + Ne - 1;

q     = x(iQ1:iQ2);
v     = x(iV1:iV2);
sigma = x(iS1:iS2);
zeta  = x(iZ1:iZ2);
D     = x(iD1:iD2);

% clamp D safely
D = min(max(D,0),1);

qMat = reshape(q,3,[])';
vMat = reshape(v,3,[])';
pos  = param.nodeCoords + qMat;
iEdge = param.edge_i;
jEdge = param.edge_j;

% Ellipsoid signed error and outward normals at current node positions
a = param.ell_a; b = param.ell_b; c = param.ell_c;
xx = pos(:,1); yy = pos(:,2); zz = pos(:,3);
fSurf = (xx./a).^2 + (yy./b).^2 + (zz./c).^2 - 1;
nxv = 2*xx/(a^2);
nyv = 2*yy/(b^2);
nzv = 2*zz/(c^2);
nrm = sqrt(nxv.^2 + nyv.^2 + nzv.^2) + 1e-18;
nSurf = [nxv./nrm, nyv./nrm, nzv./nrm];

% precompute active tension (additive)
TactVec = activeTension(t, param);

% Edge kinematics wrt fixed L0
rij = pos(jEdge,:) - pos(iEdge,:);
L = sqrt(sum(rij.^2, 2));
Lsafe = max(L, 1e-18);
n_hat = rij ./ Lsafe;
n_hat(L < 1e-18,:) = 0;

L0 = max(param.L0_edge0(:), 1e-18);
stretch = L - L0;
stretch_dot = sum((vMat(jEdge,:) - vMat(iEdge,:)) .* n_hat, 2);
eps_e = stretch ./ L0;

% Passive connection element forces
g  = max(1 - D, param.g_floor_stiff);
gd = max(1 - D, param.g_floor_damp);

k1 = g  .* param.k1_edge(:);
k2 = g  .* param.k2_edge(:);
cS = max(gd .* param.c_series_edge(:), 1e-12);
cM = max(gd .* param.c_edge(:), 1e-12);

[sig_passive, dsigma, dzeta] = connectionElement1DVectorized( ...
    param.connectionModelId_edge, stretch, stretch_dot, sigma, zeta, k1, k2, cM, cS);

% total tension = passive + active
sig_total = sig_passive + TactVec;

Fedge = sig_total .* n_hat;
Fi = Fedge;
Fj = -Fedge;

if isfield(param,'projectEdgeForcesToSurface') && param.projectEdgeForcesToSurface
    nI = nSurf(iEdge,:);
    nJ = nSurf(jEdge,:);
    Fi = Fi - sum(Fi .* nI, 2) .* nI;
    Fj = Fj - sum(Fj .* nJ, 2) .* nJ;
end

Fnet = zeros(N,3);
for dim = 1:3
    Fnet(:,dim) = param.edgeAccumI * Fi(:,dim) + param.edgeAccumJ * Fj(:,dim);
end

dD = zeros(Ne,1);

% ---- Localized damage from total tension (tension-only) ----
if isfield(param,'tipDamage_on') && param.tipDamage_on && isfield(param,'damageMask')
    if isfield(param,'eps_gate_on') && param.eps_gate_on
        gateOpen = (eps_e > param.eps_gate);
    else
        gateOpen = true(Ne,1);
    end

    eligible = (param.damageMask(:) > 0) & gateOpen;
    if any(eligible)
        sig_fail = param.sig_fail_tip;
        sig_tens = max(sig_total, 0);
        over = max(sig_tens - sig_fail, 0) ./ max(sig_fail, 1e-18);
        activeDamage = eligible & (over > 0);

        if any(activeDamage)
            pwr = max(param.dmg_pow, 1);
            tau = max(param.tau_dmg, 1e-6);
            dD(activeDamage) = dD(activeDamage) + ...
                (over(activeDamage).^pwr) .* (1 - D(activeDamage)) ./ tau;
        end
    end
end

% ---- Optional tip priming (gentle background weakening) ----
if isfield(param,'tipPriming_on') && param.tipPriming_on && isfield(param,'damageMask') && any(param.damageMask>0)
    A = 1 ./ (1 + exp(-(t - param.tipPriming_on_t)/max(param.tipPriming_tau,eps)));
    damageEdges = (param.damageMask > 0);
    dD(damageEdges) = dD(damageEdges) + (A * param.tipPrimingRate) .* (1 - D(damageEdges));
end

%% ===== Surface constraint: ellipsoid manifold via penalty + damping =====
F_surf  = -param.k_surf * fSurf .* nSurf;
vN      = sum(vMat .* nSurf, 2);
F_dampN = -param.c_surf * (vN .* nSurf);
F_drag  = -param.c_drag * vMat;

Fnet = Fnet + F_surf + F_dampN + F_drag;

% dq/dt
dqdt = reshape(vMat.', [], 1);

% dv/dt
dvdt = reshape((Fnet/param.m).',[],1);

% assemble
dxdt = [dqdt; dvdt; dsigma; dzeta; dD];
end


%% ---------------- Active tension helper ----------------
function T = activeTension(t, param)
Ne = numel(param.L0_edge0);
T = zeros(Ne,1);

if ~isfield(param,'act_on') || ~param.act_on
    return;
end

A = 1 ./ (1 + exp(-(t - param.tActOn)/max(param.tActTau,eps)));

if isfield(param,'motorMask')
    T = T + (A * param.Tact_main) .* param.motorMask(:);
else
    T = T + (A * param.Tact_main) * ones(Ne,1);
end

if isfield(param,'tipMask')
    T = T + (A * param.Tact_tip) .* param.tipMask(:);
end
end


%% ---------------- Passive connection elements ----------------
function [sig_passive, sig_dot, z_dot] = connectionElement1DVectorized(modelId, stretch, stretch_dot, sig, z, k1, k2, cM, cS)
sig_passive = zeros(size(stretch));
sig_dot = zeros(size(stretch));
z_dot = zeros(size(stretch));

isSpring = (modelId == 1);
if any(isSpring)
    sig_passive(isSpring) = k1(isSpring) .* stretch(isSpring);
end

isDamper = (modelId == 2);
if any(isDamper)
    sig_passive(isDamper) = cM(isDamper) .* stretch_dot(isDamper);
end

isMaxwell = (modelId == 3);
if any(isMaxwell)
    sig_passive(isMaxwell) = sig(isMaxwell);
    sig_dot(isMaxwell) = k1(isMaxwell) .* stretch_dot(isMaxwell) - ...
        (k1(isMaxwell) ./ cM(isMaxwell)) .* sig(isMaxwell);
end

isKelvinVoigt = (modelId == 4);
if any(isKelvinVoigt)
    sig_passive(isKelvinVoigt) = k1(isKelvinVoigt) .* stretch(isKelvinVoigt) + ...
        cM(isKelvinVoigt) .* stretch_dot(isKelvinVoigt);
end

isZener = (modelId == 5);
if any(isZener)
    sig_passive(isZener) = k1(isZener) .* stretch(isZener) + sig(isZener);
    sig_dot(isZener) = k2(isZener) .* stretch_dot(isZener) - ...
        (k2(isZener) ./ cM(isZener)) .* sig(isZener);
end

isBurgers = (modelId == 6);
if any(isBurgers)
    % Burgers = Kelvin-Voigt branch in series with Maxwell branch.
    z_dot(isBurgers) = (sig(isBurgers) - k1(isBurgers) .* z(isBurgers)) ./ cS(isBurgers);
    sig_dot(isBurgers) = k2(isBurgers) .* ...
        (stretch_dot(isBurgers) - z_dot(isBurgers)) - ...
        (k2(isBurgers) ./ cM(isBurgers)) .* sig(isBurgers);
    sig_passive(isBurgers) = sig(isBurgers);
end
end


function modelIds = connectionModelIds(param, Ne)
if isfield(param, 'connectionModel_edge') && ~isempty(param.connectionModel_edge)
    edgeModels = param.connectionModel_edge;
    modelIds = zeros(Ne,1);

    if iscell(edgeModels)
        for e = 1:Ne
            modelIds(e) = connectionModelId(normalizeConnectionModel(edgeModels{e}));
        end
    elseif isstring(edgeModels)
        for e = 1:Ne
            modelIds(e) = connectionModelId(normalizeConnectionModel(edgeModels(e)));
        end
    elseif iscategorical(edgeModels)
        for e = 1:Ne
            modelIds(e) = connectionModelId(normalizeConnectionModel(char(edgeModels(e))));
        end
    else
        error('param.connectionModel_edge must be a cell, string, or categorical array.');
    end
    return;
end

modelIds = connectionModelId(param.connectionModelKey) * ones(Ne,1);
end


function modelId = connectionModelId(modelKey)
switch modelKey
    case 'spring'
        modelId = 1;
    case 'damper'
        modelId = 2;
    case 'maxwell'
        modelId = 3;
    case 'kelvin-voigt'
        modelId = 4;
    case 'zener'
        modelId = 5;
    case 'burgers'
        modelId = 6;
    otherwise
        error('Unknown normalized connection model "%s".', modelKey);
end
end


function modelKey = normalizeConnectionModel(modelName)
if ~(ischar(modelName) || (isstring(modelName) && isscalar(modelName)))
    error('Connection model must be a character vector or scalar string.');
end

raw = lower(strtrim(char(modelName)));
raw = regexprep(raw, '^1d\s*', '');
compact = regexprep(raw, '[\s_\-/\(\)]', '');

switch compact
    case 'spring'
        modelKey = 'spring';
    case {'damper', 'dashpot'}
        modelKey = 'damper';
    case 'maxwell'
        modelKey = 'maxwell';
    case {'kelvinvoigt', 'kelvinvoight', 'voigt', 'kv'}
        modelKey = 'kelvin-voigt';
    case {'zener', 'standardlinearsolid', 'sls'}
        modelKey = 'zener';
    case 'burgers'
        modelKey = 'burgers';
    otherwise
        error(['Unknown connection model "%s". Use spring, damper, maxwell, ', ...
               'kelvin-voigt, zener, or burgers.'], modelName);
end
end


%% football nodes
function X = genEllipsoidSurfaceNodes(N, a, b, c, center)
if nargin < 5, center = [0 0 0]; end
i = (0:N-1)';
phi = (1 + sqrt(5)) / 2;
theta = 2*pi*i/phi;
z = 1 - 2*(i+0.5)/N;
r = sqrt(max(0, 1 - z.^2));
u = [r.*cos(theta), r.*sin(theta), z];
den = (u(:,1).^2)/(a^2) + (u(:,2).^2)/(b^2) + (u(:,3).^2)/(c^2);
t = 1 ./ sqrt(den);
X = u .* t + center;
end


%% --------- progress output function ----------
function status = myProgressOutputFcn(t_dim, ~, flag, tau_scale, tFinal_dim, progressEveryPhysSec)
persistent startWall nextPhysTick
status = 0;
progressEveryPhysSec = max(progressEveryPhysSec, eps);
switch flag
    case 'init'
        startWall    = tic;
        nextPhysTick = 0;
        fprintf('==> ODE start. Target t_phys = %.3f s (t_dim = %.3g)\n', ...
                tFinal_dim * tau_scale, tFinal_dim);
    case ''
        if isempty(startWall), startWall = tic; end
        if isempty(nextPhysTick), nextPhysTick = 0; end
        if ~isempty(t_dim)
            tdim  = t_dim(end);
            tphys = tdim * tau_scale;
            if tphys >= nextPhysTick
                frac    = max(0, min(tdim / tFinal_dim, 1));
                elapsed = toc(startWall);
                if frac > 0
                    remSec = elapsed * (1/frac - 1);
                    fprintf('t_phys = %6.1f s  (%3.0f%%)   ETA ~ %6.1fs\n', tphys, 100*frac, remSec);
                else
                    fprintf('t_phys = %6.1f s  (%3.0f%%)\n', tphys, 100*frac);
                end
                nextPhysTick = floor(tphys / progressEveryPhysSec) * progressEveryPhysSec + progressEveryPhysSec;
            end
        end
    case 'done'
        totalWall = toc(startWall);
        fprintf('==> ODE done. Elapsed wall time: %.2f s\n', totalWall);
        clear startWall nextPhysTick
end
end
