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

figure; scatter3(nodeCoords(:,1),nodeCoords(:,2),nodeCoords(:,3),10,'filled'); 
axis equal; grid on; title('Ellipsoid surface nodes');

param.N = N;
param.nodeCoords = nodeCoords;
param.edgeList = edgeList;
param.m = mVal;
param.Ne = Ne;

% Fixed geometric rest lengths per edge
E = edgeList;
param.L0_edge0 = sqrt(sum((nodeCoords(E(:,2),:) - nodeCoords(E(:,1),:)).^2, 2));

% Baseline per-edge Burgers parameters
param.k1_edge = kOuter_range(1) * ones(Ne,1);
param.k2_edge = kOuter_range(2) * ones(Ne,1);
param.c_edge  = cOuter_range(1) * ones(Ne,1);
param.c_series_edge = cOuter_range(2) * ones(Ne,1);

% Passive 1D connection model between cells.
% Options: 'spring', 'damper', 'maxwell', 'kelvin-voigt', 'zener', 'burgers'.
% You may also set param.connectionModel_edge to a per-edge cell/string array.
param.connectionModel = 'spring';
param.connectionModelKey = normalizeConnectionModel(param.connectionModel);

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
    'OutputFcn', @(t_dim, y, flag) myProgressOutputFcn(t_dim, y, flag, tau_scale, tFinal_dim), ...
    'Stats', 'on');

numSave = 1800;
tDimQ_all = linspace(0, tFinal_dim, numSave);

[tDimQ_all, XDimQ_all] = ode45(@(t,X) nondimRHS_shell(t, X, param), tDimQ_all, X0_dim, opts);

tSol  = tDimQ_all * tau_scale;
Y_dim = XDimQ_all(:, 1:nPos);
Y     = Y_dim * L_scale;

% extract damage over time
D_time = XDimQ_all(:, iD1:iD2);

save(simFile, ...
     'nodeCoords','edgeList','tSol','Y','D_time','param','-v7.3');

try
    confirm_actomyosin_working(simFile);
catch ME
    warning('follicle_rupture:postprocessFailed', ...
        'Postprocessing failed (%s): %s', ME.identifier, ME.message);
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
E    = param.edgeList;

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

% Edge strains wrt fixed L0
eps_e = zeros(Ne,1);
for e = 1:Ne
    i = E(e,1); j = E(e,2);
    rij  = pos(j,:) - pos(i,:);
    Lcur = norm(rij);
    L0e = max(param.L0_edge0(e), 1e-18);
    eps_e(e) = (Lcur - L0e)/L0e;
end

% Passive connection element forces
Fnet  = zeros(N,3);
dzeta = zeros(Ne,1);
dsigma= zeros(Ne,1);
dD    = zeros(Ne,1);

for e = 1:Ne
    i = E(e,1); j = E(e,2);

    rij = pos(j,:) - pos(i,:);
    L   = norm(rij);

    L0e = max(param.L0_edge0(e), 1e-18);
    stretch = L - L0e;

    if L < 1e-18
        n_hat = [0 0 0];
        stretch_dot = 0;
    else
        n_hat = rij / L;
        stretch_dot = dot(vMat(j,:) - vMat(i,:), n_hat);
    end

    % Damage -> load-bearing fraction
    g  = max(1 - D(e), param.g_floor_stiff);
    gd = max(1 - D(e), param.g_floor_damp);

    k1 = g  * param.k1_edge(e);
    k2 = g  * param.k2_edge(e);
    cS = gd * param.c_series_edge(e);
    cM = gd * param.c_edge(e);

    % avoid zeros
    cS = max(cS, 1e-12);
    cM = max(cM, 1e-12);

    % states
    z   = zeta(e);
    sig = sigma(e);

    modelKey = connectionModelForEdge(param, e);
    [sig_passive, sig_dot, z_dot] = connectionElement1D( ...
        modelKey, stretch, stretch_dot, sig, z, k1, k2, cM, cS);

    dzeta(e)  = z_dot;
    dsigma(e) = sig_dot;

    % total tension = passive + active
    sig_total = sig_passive + TactVec(e);

    Fi = sig_total * n_hat;
    Fj = -sig_total * n_hat;

    if isfield(param,'projectEdgeForcesToSurface') && param.projectEdgeForcesToSurface
        Fi = Fi - dot(Fi, nSurf(i,:)) * nSurf(i,:);
        Fj = Fj - dot(Fj, nSurf(j,:)) * nSurf(j,:);
    end

    Fnet(i,:) = Fnet(i,:) + Fi;
    Fnet(j,:) = Fnet(j,:) + Fj;

    % ---- Localized damage from total tension (tension-only) ----
    if isfield(param,'tipDamage_on') && param.tipDamage_on && isfield(param,'damageMask') && param.damageMask(e) > 0

        if ~isfield(param,'eps_gate_on') || ~param.eps_gate_on || (eps_e(e) > param.eps_gate)

            sig_fail = param.sig_fail_tip;
            sig_tens = max(sig_total, 0);

            over = max(sig_tens - sig_fail, 0) / max(sig_fail, 1e-18);
            if over > 0
                pwr = max(param.dmg_pow, 1);
                tau = max(param.tau_dmg, 1e-6);
                dD(e) = dD(e) + (over^pwr) * (1 - D(e)) / tau;
            end
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
function [sig_passive, sig_dot, z_dot] = connectionElement1D(modelKey, stretch, stretch_dot, sig, z, k1, k2, cM, cS)
sig_dot = 0;
z_dot = 0;

switch modelKey
    case 'spring'
        sig_passive = k1 * stretch;

    case 'damper'
        sig_passive = cM * stretch_dot;

    case 'maxwell'
        sig_passive = sig;
        sig_dot = k1 * stretch_dot - (k1 / cM) * sig;

    case 'kelvin-voigt'
        sig_passive = k1 * stretch + cM * stretch_dot;

    case 'zener'
        sig_passive = k1 * stretch + sig;
        sig_dot = k2 * stretch_dot - (k2 / cM) * sig;

    case 'burgers'
        % Burgers = Kelvin-Voigt branch in series with Maxwell branch.
        z_dot = (sig - k1 * z) / cS;
        sig_dot = k2 * (stretch_dot - z_dot) - (k2 / cM) * sig;
        sig_passive = sig;

    otherwise
        error('Unknown normalized connection model "%s".', modelKey);
end
end


function modelKey = connectionModelForEdge(param, e)
if isfield(param, 'connectionModel_edge') && ~isempty(param.connectionModel_edge)
    edgeModels = param.connectionModel_edge;

    if iscell(edgeModels)
        modelKey = normalizeConnectionModel(edgeModels{e});
    elseif isstring(edgeModels)
        modelKey = normalizeConnectionModel(edgeModels(e));
    elseif iscategorical(edgeModels)
        modelKey = normalizeConnectionModel(char(edgeModels(e)));
    else
        error('param.connectionModel_edge must be a cell, string, or categorical array.');
    end
    return;
end

if isfield(param, 'connectionModelKey')
    modelKey = param.connectionModelKey;
else
    modelKey = normalizeConnectionModel(param.connectionModel);
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
function status = myProgressOutputFcn(t_dim, ~, flag, tau_scale, tFinal_dim)
persistent startWall nextPhysTick
status = 0;
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
                nextPhysTick = floor(tphys) + 1;
            end
        end
    case 'done'
        totalWall = toc(startWall);
        fprintf('==> ODE done. Elapsed wall time: %.2f s\n', totalWall);
        clear startWall nextPhysTick
end
end
