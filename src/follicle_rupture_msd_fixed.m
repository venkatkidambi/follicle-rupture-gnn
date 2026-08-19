function follicle_rupture_msd_fixed
% Follicle rupture shell model
% State: X = [q; v; sigma; zeta; D]
%   q, v        node displacement and velocity
%   sigma,zeta  viscoelastic edge states
%   D           edge damage in [0,1], where 1 means failed

clc; close all;

rng(42);

%% USER-PARAMETERS
% Cell-scale coordinate and display configuration
cellDiameter_um = 10;          % positive cell diameter in microns
coordinateInputMode = 'microns'; % trackedNodes coordinates are treated as microns
displayMode = 'microns';       % 'cells' or 'microns'

% Generated ellipsoid semi-axes in coordinateInputMode units. Defaults are
% 50 x 15 x 15 cell diameters = 500 x 150 x 150 microns for 10 um cells.
generatedEllipsoidSemiAxes = [50 15 15];

% All coordinate inputs are converted to meters once at geometry ingestion.
% Mechanics and saved physical data remain SI: force N, stiffness N/m,
% damping N*s/m, mass kg, time s; strain and damage are dimensionless.
units = buildUnitSystem(cellDiameter_um, coordinateInputMode, displayMode);

% Simulation resolution / scale
N = 300;                    % overridden by tracked node count below
mVal = 1;                 % [kg] smaller mass reduces inertial coasting
useExperimentalNodeData = true;
experimentalNodeDataFileName = fullfile('tracked_nodes', 'trackedNodes_nuc_48_2_OA7_1.mat');
experimentalNodeDataVariable = 'trackedNodes';
experimentalNodeDataFrame = 1; % use initial tracked coordinates only
experimentalCenterMode = 'median';   % 'median', 'mean', or 'none'
experimentalOrientMode = 'principal-axes'; % 'principal-axes' aligns long axis to model x-axis
experimentalLongAxisSign = -1;        % set 1 to flip experimental anterior/posterior ends
experimentalSurfaceConstraintMode = 'reference-ellipsoid'; % preserve each measured node's initial shell level
experimentalAnchorOn = false;          % measured node clouds should not have a frozen anterior shell
experimentalEdgeMode = 'knn';          % 'knn' connects all measured nodes; 'hull' uses convex hull only
experimentalNeighborCount = 6;         % local neighbors per measured node for 'knn' edge mode

% Passive edge material parameters
kOuter_range = [1, 1];
cOuter_range = [1, 1];

% Time integration and saved output
tEnd  = 1800;                    % s
numSave = 300;                 % saved output frames; solver still adapts internally
ruptureSaveEnd = 120;          % [s], densely sample the rupture transient
ruptureSaveInterval = 0.25;    % [s], output spacing during the rupture window
odeSolver = 'ode45';           % instant damage is stiff; ode15s handles it better
progressEveryPhysSec = 5;      % console progress interval in physical seconds

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(repoRoot, 'results');
inputNodeDataDir = fullfile(repoRoot, 'input_node_data');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if ~exist(inputNodeDataDir, 'dir')
    mkdir(inputNodeDataDir);
end
simFile = fullfile(outputDir, 'follicleShellSim_test2.mat');
experimentalNodeDataFile = fullfile(inputNodeDataDir, experimentalNodeDataFileName);

%% Nodes
if useExperimentalNodeData
    [nodeCoords, experimentalCellIds] = loadExperimentalNodeCoords( ...
        experimentalNodeDataFile, units.coordinateToMeters, experimentalCenterMode, ...
        experimentalNodeDataVariable, experimentalNodeDataFrame);
    nodeCoords = orientExperimentalNodeCoords( ...
        nodeCoords, experimentalOrientMode, experimentalLongAxisSign);
    N = size(nodeCoords, 1);
    [a, b, c] = fitAxisAlignedEllipsoidSemiAxes(nodeCoords);
    L_scale = median([a b c]);  % representative experimental radius scale
else
    ellipsoidSemiAxes_m = generatedEllipsoidSemiAxes * units.coordinateToMeters;
    a = ellipsoidSemiAxes_m(1); b = ellipsoidSemiAxes_m(2); c = ellipsoidSemiAxes_m(3);
    L_scale = b;              % transverse follicle radius scale
    nodeCoords = genEllipsoidSurfaceNodes(N, a, b, c, [0 0 0]);
    experimentalCellIds = [];
end

param.ell_a = a; param.ell_b = b; param.ell_c = c;
param.surfaceLevel0 = (nodeCoords(:,1)./a).^2 + ...
    (nodeCoords(:,2)./b).^2 + (nodeCoords(:,3)./c).^2;
param.units = units;
param.generatedEllipsoidSemiAxes = generatedEllipsoidSemiAxes;
param.useExperimentalNodeData = useExperimentalNodeData;
param.experimentalNodeDataFile = experimentalNodeDataFile;
param.experimentalNodeDataVariable = experimentalNodeDataVariable;
param.experimentalNodeDataFrame = experimentalNodeDataFrame;
param.experimentalCoordinateScale = units.coordinateToMeters;
param.experimentalCenterMode = experimentalCenterMode;
param.experimentalOrientMode = experimentalOrientMode;
param.experimentalLongAxisSign = experimentalLongAxisSign;
param.experimentalSurfaceConstraintMode = experimentalSurfaceConstraintMode;
param.experimentalAnchorOn = experimentalAnchorOn;
param.experimentalEdgeMode = experimentalEdgeMode;
param.experimentalNeighborCount = experimentalNeighborCount;
param.experimentalCellIds = experimentalCellIds;

% Posterior axis convention
param.postAxis = 1;
param.postSign = 1; % tracked-node rupture/0-10% bands use the opposite physical tip

% Surface penalty + normal damping
param.k_surf = 5e-6;     % [N] penalty strength
param.c_surf = 1;        % [N*s/m] suppress motion away from the local shell tangent
param.c_drag = 0.1;        % [N*s/m] overdamped tissue drag prevents post-cut coasting
if useExperimentalNodeData
    param.surfaceConstraintMode = experimentalSurfaceConstraintMode;
    param.projectEdgeForcesToSurface = true;
else
    param.surfaceConstraintMode = 'ellipsoid';
    param.projectEdgeForcesToSurface = true;
end

%% Surface/local edges (shell-only)
if useExperimentalNodeData
    edgeList = buildExperimentalEdgeList(nodeCoords, experimentalEdgeMode, experimentalNeighborCount);
else
    DT = delaunayTriangulation(nodeCoords);
    F = freeBoundary(DT);   % [nFaces x 3]
    allEdges = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    allEdges = sort(allEdges,2);
    edgeList = unique(allEdges,'rows');
end
Ne = size(edgeList,1);

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
edgeNodeIncidence = sparse([param.edge_i; param.edge_j], [(1:Ne)'; (1:Ne)'], 1, N, Ne);
nodeDegree = full(sum(edgeNodeIncidence, 2));
fprintf('[MESH] nodes = %d | edges = %d | edge length median/p95/max = %.3g/%.3g/%.3g %s | isolated nodes = %d\n', ...
    N, Ne, median(param.L0_edge0) * units.displayFromMeters, ...
    prctile(param.L0_edge0, 95) * units.displayFromMeters, ...
    max(param.L0_edge0) * units.displayFromMeters, units.displayLengthUnit, ...
    nnz(nodeDegree == 0));
if any(nodeDegree == 0)
    warning('follicle_rupture:isolatedNodes', ...
        'Mesh contains %d isolated nodes that cannot move mechanically.', nnz(nodeDegree == 0));
end
% Baseline per-edge passive model parameters
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

%% Spatial localization masks (main-body vs tiny posterior tip)
edgeMid = 0.5*(nodeCoords(E(:,1),:) + nodeCoords(E(:,2),:) );
ax = param.postAxis;
sg = param.postSign;
param.edgeMid0 = edgeMid;

sNodes = sg * nodeCoords(:,ax);
sMid   = sg * edgeMid(:,ax);

param.mainPct = 98.5;      % posterior-most ~1.5% excluded from the main body mask
param.tipPct  = 98.5;    % posterior-most ~1.5% defines tiny autonomous rupture tip

sCutMain = prctile(sNodes, param.mainPct);
sCutTip  = prctile(sNodes, param.tipPct);

param.bodyMask = double(sMid < sCutMain);     % main-body edges
param.tipMask  = double(sMid >= sCutTip);     % tiny posterior cap edges

%% Damage region selection
param.damageMode = 'principal-0-10-tip-patch';
param.tipPatchEdgeCount = 8; % compact weak patch; increase for a larger opening
param.seamWidthFrac = 0.04;   % fraction of posterior-score range

% Compute the same initial long-axis coordinate used by the post-processing
% plot. Orient its sign by the existing posterior convention so percentile
% membership is deterministic rather than depending on EIG's arbitrary sign.
centeredNodes = nodeCoords - mean(nodeCoords, 1);
principalCovariance = (centeredNodes.' * centeredNodes) / max(N - 1, 1);
[principalVectors, principalValues] = eig(principalCovariance, 'vector');
[~, principalOrder] = sort(principalValues, 'descend');
principalLongAxis = principalVectors(:, principalOrder(1));
posteriorDirection = zeros(3,1);
posteriorDirection(param.postAxis) = param.postSign;
if dot(principalLongAxis, posteriorDirection) < 0
    principalLongAxis = -principalLongAxis;
end
principalNodeScore = centeredNodes * principalLongAxis;
principalEdgeScore = mean(principalNodeScore(E), 2);
principalTipCutoff = prctile(principalNodeScore, 10);

param.principalLongAxis0 = principalLongAxis;
param.principalNodeScore0 = principalNodeScore;
param.principalTipPercentile = 10;
param.principalTipCutoff = principalTipCutoff;

switch lower(param.damageMode)
    case 'seam'
        seamWidth = param.seamWidthFrac * range(sNodes);
        param.damageMask = double(abs(sMid - sCutTip) <= seamWidth);
    case 'tip'
        param.damageMask = param.tipMask;
    case 'tipmost-edge'
        param.damageMask = zeros(Ne,1);
        [~, tipmostEdgeIdx] = max(sMid);
        param.damageMask(tipmostEdgeIdx) = 1;
        param.tipmostDamageEdge = tipmostEdgeIdx;
    case 'principal-0-10-tip-edge'
        % Restrict candidates to edges fully contained in the plotted 0-10%
        % node band, then choose the one nearest its extreme tip. This keeps
        % rupture selection and visualization in the same coordinate system.
        principalTipCandidates = find( ...
            principalNodeScore(E(:,1)) <= principalTipCutoff & ...
            principalNodeScore(E(:,2)) <= principalTipCutoff);
        if isempty(principalTipCandidates)
            error('No mesh edge has both endpoints in the principal-axis 0-10%% tip band.');
        end
        [~, localTipIdx] = min(principalEdgeScore(principalTipCandidates));
        principalTipEdgeIdx = principalTipCandidates(localTipIdx);
        param.damageMask = zeros(Ne,1);
        param.damageMask(principalTipEdgeIdx) = 1;
        param.principalTipDamageEdge = principalTipEdgeIdx;
        param.principalTipCandidateEdges = principalTipCandidates;
    case 'principal-0-10-tip-patch'
        % Select a compact group of edges nearest the extreme of the same
        % 0-10% principal-axis band shown in post-processing. Each edge must
        % still reach its own mechanical failure threshold; this defines a
        % weak tissue patch, not a prescribed cut time.
        principalTipCandidates = find( ...
            principalNodeScore(E(:,1)) <= principalTipCutoff & ...
            principalNodeScore(E(:,2)) <= principalTipCutoff);
        if isempty(principalTipCandidates)
            error('No mesh edge has both endpoints in the principal-axis 0-10%% tip band.');
        end
        [~, patchOrder] = sort(principalEdgeScore(principalTipCandidates), 'ascend');
        patchCount = min(param.tipPatchEdgeCount, numel(principalTipCandidates));
        principalTipPatchEdges = principalTipCandidates(patchOrder(1:patchCount));
        param.damageMask = zeros(Ne,1);
        param.damageMask(principalTipPatchEdges) = 1;
        param.principalTipPatchEdges = principalTipPatchEdges;
        param.principalTipCandidateEdges = principalTipCandidates;
    otherwise
        error(['Unknown damageMode "%s". Use "seam", "tip", ', ...
               '"tipmost-edge", "principal-0-10-tip-edge", or ', ...
               '"principal-0-10-tip-patch".'], param.damageMode);
end

fprintf('[MASK] body edges = %d/%d | tip edges = %d | damage edges = %d\n', ...
    nnz(param.bodyMask), Ne, nnz(param.tipMask), nnz(param.damageMask));
if nnz(param.damageMask) == 0
    error('Rupture-region selection produced zero damage edges.');
end
selectedDamageEdge = find(param.damageMask > 0);
rupturePatchNodes = unique(E(selectedDamageEdge,:));
param.rupturePatchNodes = rupturePatchNodes;
if isscalar(selectedDamageEdge)
    selectedEndpointPercentiles = 100 * [ ...
        mean(principalNodeScore <= principalNodeScore(E(selectedDamageEdge,1))), ...
        mean(principalNodeScore <= principalNodeScore(E(selectedDamageEdge,2)))];
    fprintf('[RUPTURE LOCATION] edge = %d | principal-axis endpoint percentiles = %.3g, %.3g\n', ...
        selectedDamageEdge, selectedEndpointPercentiles(1), selectedEndpointPercentiles(2));
else
    patchPercentiles = 100 * arrayfun(@(nodeIdx) ...
        mean(principalNodeScore <= principalNodeScore(nodeIdx)), rupturePatchNodes);
    fprintf('[RUPTURE LOCATION] patch edges = %s | nodes = %d | percentile span = %.3g-%.3g\n', ...
        mat2str(selectedDamageEdge.'), numel(rupturePatchNodes), ...
        min(patchPercentiles), max(patchPercentiles));
end

%% Localized rim retraction
% Contract the intact shoulder edges that join the rupture-patch nodes to
% the surrounding shell. Contracting the rupture edges themselves collapses
% their endpoints into a cluster before failure. The shoulder arrangement
% instead loads the weak patch in tension and, after failure, keeps drawing
% its rim nodes back along the shell surface.
isPatchNodeI = ismember(E(:,1), rupturePatchNodes);
isPatchNodeJ = ismember(E(:,2), rupturePatchNodes);
% The patch is at the minimum principal score. Keep only links whose
% non-patch endpoint has the larger score, i.e. links directed back toward
% the follicle body instead of farther toward the extreme tip.
bodyFacingFromI = isPatchNodeI & ~isPatchNodeJ & ...
    principalNodeScore(E(:,2)) > principalNodeScore(E(:,1));
bodyFacingFromJ = isPatchNodeJ & ~isPatchNodeI & ...
    principalNodeScore(E(:,1)) > principalNodeScore(E(:,2));
param.retractionMask = double((bodyFacingFromI | bodyFacingFromJ) & param.damageMask == 0);
if ~any(param.retractionMask)
    error('Rupture patch has no intact shoulder edges available for retraction.');
end

param.contractionStartTime = 3;       % [s]
param.contractionRateFrac_edge = zeros(Ne,1);
param.contractionRateFrac_edge(param.retractionMask > 0) = 0.20; % [L0/s]
param.minRestLengthFrac_edge = ones(Ne,1);
param.minRestLengthFrac_edge(param.retractionMask > 0) = 0.30;

% The patch edges are weak and fail when stretched by the retracting rim;
% the shoulder edges remain intact so they continue pulling the rim back.
param.extensionCritFrac_edge = inf(Ne,1);
param.extensionCritFrac_edge(param.damageMask > 0) = 0.15;

if any(param.contractionRateFrac_edge < 0) || ...
        any(param.minRestLengthFrac_edge <= 0 | param.minRestLengthFrac_edge > 1)
    error('Contraction rates must be nonnegative and rest-length floors must lie in (0,1].');
end

% Optional weak anterior anchoring. This is a soft positional boundary
% condition to stabilize directional epithelial retraction, not a rigid pin.
if useExperimentalNodeData
    param.anchor_on = experimentalAnchorOn;
else
    param.anchor_on = false; % temporarily disabled to allow post-rupture recoil
end
param.anchorPct = 60;          % anterior-most percentile of nodes
param.k_anchor = 10;        % weak positional tether stiffness
param.c_anchor = 1;        % anchor damping

sCutAnchor = prctile(sNodes, param.anchorPct);
param.anchorNodeMask = double(sNodes <= sCutAnchor);
fprintf('[ANCHOR] enabled = %d | eligible anchor nodes = %d/%d\n', ...
    param.anchor_on, nnz(param.anchorNodeMask), N);

retractingEdges = find(param.retractionMask > 0);
fprintf('[RETRACTION] shoulder edges = %s | start = %.3g s | rate = %.3g L0/s | floor = %.3g L0\n', ...
    mat2str(retractingEdges.'), param.contractionStartTime, ...
    param.contractionRateFrac_edge(retractingEdges(1)), ...
    param.minRestLengthFrac_edge(retractingEdges(1)));
fprintf('[RUPTURE] weak patch edges = %s | break = %.3g L0 extension\n', ...
    mat2str(selectedDamageEdge.'), param.extensionCritFrac_edge(selectedDamageEdge(1)));

%% Damage / rupture parameters
param.ruptureDamageTau = 0.005; % [s], rapid continuous approximation to a break

% Minimum load-bearing fraction to avoid numerical singularities
param.g_floor_stiff = 0;        % allow fully failed
param.g_floor_damp  = 1e-6;     % keep a tiny damping floor

%% Nondimensionalization
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
param.ruptureSaveEnd = ruptureSaveEnd;
param.ruptureSaveInterval = ruptureSaveInterval;
param.odeSolver = odeSolver;
param.progressEveryPhysSec = progressEveryPhysSec;

%% ================= State sizes =================
nPos = 3*N;
nVel = 3*N;

param.nPos = nPos;
param.nVel = nVel;

% X = [q; v; sigma; zeta; D]
nX = nPos + nVel + Ne + Ne + Ne;

% Edge state block indices for extracting saved time courses
iD1 = nPos + nVel + 2*Ne + 1;
iD2 = iD1 + Ne - 1;

%% Integrate
tFinal_dim = tEnd / tau_scale;
assert(isscalar(tFinal_dim) && isfinite(tFinal_dim));

X0_dim = zeros(nX,1);

opts = odeset( ...
    'RelTol', 1e-5, ...
    'AbsTol', 1e-8, ...
    'MaxStep', 5, ...
    'InitialStep', 1e-3, ...
    'OutputFcn', @(t_dim, y, flag) myProgressOutputFcn(t_dim, y, flag, tau_scale, tFinal_dim, param.progressEveryPhysSec), ...
    'Stats', 'on');

ruptureSaveEnd = min(max(param.ruptureSaveEnd, 0), tEnd);
if ~(isscalar(param.ruptureSaveInterval) && isfinite(param.ruptureSaveInterval) && ...
        param.ruptureSaveInterval > 0)
    error('ruptureSaveInterval must be a positive finite scalar.');
end
tPhysQ_early = (0:param.ruptureSaveInterval:ruptureSaveEnd).';
tPhysQ_late = linspace(ruptureSaveEnd, tEnd, numSave).';
tPhysQ_all = unique([tPhysQ_early; tPhysQ_late]);
tDimQ_all = tPhysQ_all / tau_scale;

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
V_dim = XDimQ_all(:, nPos+1:nPos+nVel);
V     = V_dim * (L_scale/tau_scale);
L0_edge0 = param.L0_edge0;

% extract edge state time courses
D_time = min(max(XDimQ_all(:, iD1:iD2), 0), 1);

% Save the imposed natural-length history alongside displacement/damage so
% rupture can be checked directly in post-processing.
Lnatural_time = zeros(numel(tSol), Ne);
for it = 1:numel(tSol)
    Lnatural_time(it,:) = contractingNaturalLengths(tSol(it), param).';
end

save(simFile, ...
     'nodeCoords','edgeList','L0_edge0','Lnatural_time','tSol','Y','V','D_time','param','-v7.3');

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

% clamp bounded edge states safely
D = min(max(D,0),1);

qMat = reshape(q,3,[])';
vMat = reshape(v,3,[])';
pos  = param.nodeCoords + qMat;
iEdge = param.edge_i;
jEdge = param.edge_j;

surfaceConstraintMode = 'ellipsoid';
if isfield(param,'surfaceConstraintMode') && ~isempty(param.surfaceConstraintMode)
    surfaceConstraintMode = lower(strtrim(char(param.surfaceConstraintMode)));
end
if ~ismember(surfaceConstraintMode, {'ellipsoid', 'reference-ellipsoid', 'project-ellipsoid', 'none', 'off'})
    error(['Unknown surfaceConstraintMode "%s". Use "ellipsoid", ', ...
        '"reference-ellipsoid", "project-ellipsoid", or "none".'], surfaceConstraintMode);
end
useEllipsoidNormals = any(strcmp(surfaceConstraintMode, ...
    {'ellipsoid', 'reference-ellipsoid', 'project-ellipsoid'}));
useEllipsoidPenalty = any(strcmp(surfaceConstraintMode, {'ellipsoid', 'reference-ellipsoid'}));
useEllipsoidNormalDamping = useEllipsoidNormals;

% Ellipsoid signed error and outward normals at current node positions.
% Experimental node clouds can keep their measured shape by disabling this.
if useEllipsoidNormals
    a = param.ell_a; b = param.ell_b; c = param.ell_c;
    xx = pos(:,1); yy = pos(:,2); zz = pos(:,3);
    fSurf = (xx./a).^2 + (yy./b).^2 + (zz./c).^2 - 1;
    if strcmp(surfaceConstraintMode, 'reference-ellipsoid')
        fSurf = fSurf - (param.surfaceLevel0 - 1);
    end
    nxv = 2*xx/(a^2);
    nyv = 2*yy/(b^2);
    nzv = 2*zz/(c^2);
    nrm = sqrt(nxv.^2 + nyv.^2 + nzv.^2) + 1e-18;
    nSurf = [nxv./nrm, nyv./nrm, nzv./nrm];
else
    fSurf = zeros(N,1);
    nSurf = zeros(N,3);
end

% Edge kinematics relative to the imposed time-dependent natural lengths.
rij = pos(jEdge,:) - pos(iEdge,:);
L = sqrt(sum(rij.^2, 2));
Lsafe = max(L, 1e-18);
n_hat = rij ./ Lsafe;
n_hat(L < 1e-18,:) = 0;

[Lnatural, LnaturalDot] = contractingNaturalLengths(t, param);
L0_ref = max(param.L0_edge0(:), 1e-18);
geometricStretchDot = sum((vMat(jEdge,:) - vMat(iEdge,:)) .* n_hat, 2);
stretch = L - Lnatural;
stretch_dot = geometricStretchDot - LnaturalDot;

% Passive connection element forces
g  = max(1 - D, param.g_floor_stiff);
gd = max(1 - D, param.g_floor_damp);

k1 = g  .* param.k1_edge(:);
k2 = g  .* param.k2_edge(:);
cS = max(gd .* param.c_series_edge(:), 1e-12);
cM = max(gd .* param.c_edge(:), 1e-12);

[sig_passive, dsigma, dzeta] = connectionElement1DVectorized( ...
    param.connectionModelId_edge, stretch, stretch_dot, sigma, zeta, k1, k2, cM, cS);

% Ensure a fully damaged history-dependent edge carries no residual stored
% stress (Maxwell/Burgers sigma otherwise remains a state after failure).
sig_total = (1 - D) .* sig_passive;

Fedge = sig_total .* n_hat;
Fi = Fedge;
Fj = -Fedge;

if useEllipsoidNormals && isfield(param,'projectEdgeForcesToSurface') && param.projectEdgeForcesToSurface
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

% ---- Weakest-edge rupture from critical extension ----
% This is the direct 3D counterpart of ext >= delt_crit in the chain model.
criticalExtension = param.extensionCritFrac_edge(:) .* L0_ref;
breaking = (param.damageMask(:) > 0) & (stretch >= criticalExtension);
if any(breaking)
    tauBreak = max(param.ruptureDamageTau, 1e-8);
    dD(breaking) = (1 - D(breaking)) ./ tauBreak;
end

%% ===== Optional surface constraint: ellipsoid manifold via penalty + damping =====
F_surf = zeros(N,3);
F_dampN = zeros(N,3);
if useEllipsoidPenalty
    F_surf  = -param.k_surf * fSurf .* nSurf;
end
if useEllipsoidNormalDamping
    vN      = sum(vMat .* nSurf, 2);
    F_dampN = -param.c_surf * (vN .* nSurf);
end
F_drag  = -param.c_drag * vMat;

Fnet = Fnet + F_surf + F_dampN + F_drag;

% Soft anterior positional boundary condition. This weak tether biases
% directional retraction without rigidly pinning nodes or adding contact.
if isfield(param,'anchor_on') && param.anchor_on
    F_anchor = -param.k_anchor * qMat - param.c_anchor * vMat;
    Fnet = Fnet + param.anchorNodeMask .* F_anchor;
end

% dq/dt
dqdt = reshape(vMat.', [], 1);

% dv/dt
dvdt = reshape((Fnet/param.m).',[],1);

% assemble
dxdt = [dqdt; dvdt; dsigma; dzeta; dD];
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


function [Lnatural, LnaturalDot] = contractingNaturalLengths(t, param)
% Prescribed natural-length contraction, equivalent to
% Ln=max(L0-rate*(t-tStart), Lmin) in the 1D reference model.
L0 = param.L0_edge0(:);
rateFrac = param.contractionRateFrac_edge(:);
floorFrac = param.minRestLengthFrac_edge(:);
elapsed = max(t - param.contractionStartTime, 0);
unclampedFrac = 1 - rateFrac .* elapsed;
lengthFrac = max(unclampedFrac, floorFrac);
Lnatural = L0 .* lengthFrac;

isShortening = (t >= param.contractionStartTime) & ...
    (rateFrac > 0) & (unclampedFrac > floorFrac);
LnaturalDot = -L0 .* rateFrac .* isShortening;
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
oversampleFactor = 20;
M = max(oversampleFactor * N, 2000);
i = (0:M-1)';
phi = (1 + sqrt(5)) / 2;
theta = 2*pi*i/phi;
z = 1 - 2*(i+0.5)/M;
r = sqrt(max(0, 1 - z.^2));
u = [r.*cos(theta), r.*sin(theta), z];

% Area-weight candidates for the linear sphere-to-ellipsoid map. Uniform
% directions alone overpopulate visually compressed regions on elongated
% ellipsoids; equal-area quantiles give a more even surface point cloud.
areaWeight = a*b*c * sqrt((u(:,1)./a).^2 + (u(:,2)./b).^2 + (u(:,3)./c).^2);
cdf = cumsum(areaWeight);
targets = ((0:N-1)' + 0.5) * (cdf(end) / N);
idx = interp1(cdf, (1:M)', targets, 'nearest', 'extrap');
idx = min(max(idx, 1), M);

X = [a*u(idx,1), b*u(idx,2), c*u(idx,3)] + center;
end


function edgeList = buildExperimentalEdgeList(nodeCoords, edgeMode, neighborCount)
if nargin < 2 || isempty(edgeMode)
    edgeMode = 'knn';
end
if nargin < 3 || isempty(neighborCount)
    neighborCount = 6;
end

if isstring(edgeMode)
    edgeModeKey = lower(char(edgeMode));
else
    edgeModeKey = lower(edgeMode);
end

switch edgeModeKey
    case {'hull', 'convex-hull', 'convex_hull'}
        DT = delaunayTriangulation(nodeCoords);
        F = freeBoundary(DT);
        allEdges = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
        edgeList = unique(sort(allEdges, 2), 'rows');
    case {'knn', 'nearest', 'nearest-neighbor', 'nearest_neighbor'}
        N = size(nodeCoords, 1);
        k = max(1, min(round(neighborCount), N - 1));
        dist2 = squaredDistanceMatrix(nodeCoords);
        dist2(1:N+1:end) = inf;
        edges = zeros(N * k, 2);
        row0 = 0;
        for nodeIdx = 1:N
            [~, order] = sort(dist2(nodeIdx,:), 'ascend');
            nbr = order(1:k).';
            edges(row0 + (1:k), :) = [repmat(nodeIdx, k, 1), nbr];
            row0 = row0 + k;
        end
        edgeList = unique(sort(edges, 2), 'rows');
    otherwise
        error('Unknown experimentalEdgeMode "%s". Use "knn" or "hull".', edgeMode);
end

fprintf('[MESH] experimental edge mode = %s\n', edgeModeKey);
end


function dist2 = squaredDistanceMatrix(X)
rowNorm2 = sum(X.^2, 2);
dist2 = rowNorm2 + rowNorm2.' - 2 * (X * X.');
dist2 = max(dist2, 0);
end


function [nodeCoords, cellIds] = loadExperimentalNodeCoords(dataFile, coordinateScale, centerMode, variableName, frameIndex)
if nargin < 2 || isempty(coordinateScale)
    coordinateScale = 1;
end
if nargin < 3 || isempty(centerMode)
    centerMode = 'median';
end
if nargin < 4 || isempty(variableName)
    variableName = 'trackedNodes';
end
if nargin < 5 || isempty(frameIndex)
    frameIndex = 1;
end
if ~isfile(dataFile)
    error('Experimental node data file not found: %s', dataFile);
end

[~, ~, extension] = fileparts(dataFile);
switch lower(extension)
    case '.mat'
        fileData = load(dataFile, variableName);
        if ~isfield(fileData, variableName)
            error('MAT-file %s does not contain variable "%s".', dataFile, variableName);
        end
        tracked = fileData.(variableName);
        if ~isnumeric(tracked) || ndims(tracked) > 3 || ...
                size(tracked,2) ~= 3 || size(tracked,3) < frameIndex || frameIndex < 1
            error(['Variable "%s" must be numeric [N x 3 x T], and frameIndex ', ...
                   'must identify an available frame.'], variableName);
        end
        nodeCoords = double(tracked(:,:,frameIndex)) * coordinateScale;
        cellIds = (1:size(nodeCoords,1)).';
        placeholderRows = all(nodeCoords == 0, 2);
        if any(placeholderRows)
            fprintf('[GEOMETRY] omitted %d zero-placeholder tracks absent from frame %d\n', ...
                nnz(placeholderRows), frameIndex);
        end
        fprintf('[GEOMETRY] using frame %d only from %s (%s)\n', ...
            frameIndex, dataFile, variableName);
    otherwise
        raw = readmatrix(dataFile);
        if size(raw, 2) < 4
            error('Experimental node data must contain Cell ID, x, y, z in the first four columns.');
        end
        cellIds = raw(:,1);
        nodeCoords = raw(:,2:4) * coordinateScale;
        placeholderRows = false(size(nodeCoords,1),1);
end

validRows = all(isfinite(nodeCoords), 2) & ~placeholderRows;
cellIds = cellIds(validRows);
nodeCoords = nodeCoords(validRows, :);
[nodeCoords, uniqueIdx] = unique(nodeCoords, 'rows', 'stable');
cellIds = cellIds(uniqueIdx);

if size(nodeCoords, 1) < 4
    error('Experimental node data must contain at least four valid xyz rows.');
end

if isstring(centerMode)
    centerModeKey = lower(char(centerMode));
else
    centerModeKey = lower(centerMode);
end
switch centerModeKey
    case 'median'
        nodeCoords = nodeCoords - median(nodeCoords, 1, 'omitnan');
    case 'mean'
        nodeCoords = nodeCoords - mean(nodeCoords, 1, 'omitnan');
    case 'none'
        % Keep spreadsheet coordinates as-is.
    otherwise
        error('Unknown experimentalCenterMode "%s". Use "median", "mean", or "none".', centerMode);
end

fprintf('[GEOMETRY] loaded %d experimental nodes from %s\n', size(nodeCoords,1), dataFile);
end


function nodeCoords = orientExperimentalNodeCoords(nodeCoords, orientMode, longAxisSign)
if nargin < 2 || isempty(orientMode)
    orientMode = 'none';
end
if nargin < 3 || isempty(longAxisSign)
    longAxisSign = 1;
end

if isstring(orientMode)
    orientModeKey = lower(char(orientMode));
else
    orientModeKey = lower(orientMode);
end

switch orientModeKey
    case {'none', 'off'}
        return;
    case {'principal-axes', 'principal_axes', 'pca'}
        if ~isscalar(longAxisSign) || ~ismember(longAxisSign, [-1 1])
            error('experimentalLongAxisSign must be 1 or -1.');
        end

        [~, ~, principalAxes] = svd(nodeCoords, 'econ');

        % SVD axis signs are arbitrary. Make them deterministic before the
        % optional anterior/posterior flip so repeated runs orient the same.
        for axisIdx = 1:size(principalAxes, 2)
            [~, dominantCoord] = max(abs(principalAxes(:,axisIdx)));
            if principalAxes(dominantCoord, axisIdx) < 0
                principalAxes(:,axisIdx) = -principalAxes(:,axisIdx);
            end
        end
        if det(principalAxes) < 0
            principalAxes(:,3) = -principalAxes(:,3);
        end

        nodeCoords = nodeCoords * principalAxes;
        nodeCoords(:,1) = longAxisSign * nodeCoords(:,1);

        fprintf('[GEOMETRY] oriented experimental nodes with principal axes (long axis -> x)\n');
    otherwise
        error('Unknown experimentalOrientMode "%s". Use "principal-axes" or "none".', orientMode);
end
end


function [a, b, c] = fitAxisAlignedEllipsoidSemiAxes(nodeCoords)
axisWeights = [nodeCoords(:,1).^2, nodeCoords(:,2).^2, nodeCoords(:,3).^2];
validRows = all(isfinite(axisWeights), 2);
axisWeights = axisWeights(validRows, :);

if size(axisWeights, 1) < 4
    error('Could not fit ellipsoid: at least four finite xyz rows are required.');
end

invAxesSquared = axisWeights \ ones(size(axisWeights, 1), 1);
if any(~isfinite(invAxesSquared) | invAxesSquared <= 0)
    warning('follicle_rupture:ellipsoidFitFallback', ...
        'Least-squares ellipsoid fit was not positive; using coordinate-span semi-axes instead.');
    [a, b, c] = estimateAxisAlignedSemiAxes(nodeCoords);
    return;
end

semiAxes = sqrt(1 ./ invAxesSquared(:)).';

if any(~isfinite(semiAxes) | semiAxes <= 0)
    error('Could not estimate positive semi-axes from fitted experimental node coordinates.');
end

a = semiAxes(1);
b = semiAxes(2);
c = semiAxes(3);
fprintf('[GEOMETRY] fitted ellipsoid semi-axes: a=%.3g, b=%.3g, c=%.3g m\n', a, b, c);
end


function [a, b, c] = estimateAxisAlignedSemiAxes(nodeCoords)
semiAxes = max(abs(nodeCoords), [], 1, 'omitnan');
spanAxes = 0.5 * range(nodeCoords, 1);
semiAxes = max(semiAxes, spanAxes);
semiAxes(~isfinite(semiAxes) | semiAxes <= 0) = median(semiAxes(isfinite(semiAxes) & semiAxes > 0));

if any(~isfinite(semiAxes) | semiAxes <= 0)
    error('Could not estimate positive semi-axes from experimental node coordinates.');
end

a = semiAxes(1);
b = semiAxes(2);
c = semiAxes(3);
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


function units = buildUnitSystem(cellDiameter_um, coordinateInputMode, displayMode)
if ~(isnumeric(cellDiameter_um) && isreal(cellDiameter_um) && ...
        isscalar(cellDiameter_um) && isfinite(cellDiameter_um) && cellDiameter_um > 0)
    error('follicle_rupture:invalidCellDiameter', ...
        'cellDiameter_um must be a positive, finite numeric scalar.');
end

coordinateInputMode = normalizeUnitMode(coordinateInputMode, ...
    {'cells', 'microns', 'meters'}, 'coordinate input');
displayMode = normalizeUnitMode(displayMode, ...
    {'cells', 'microns'}, 'display');

cellDiameter_m = cellDiameter_um * 1e-6;
switch coordinateInputMode
    case 'cells'
        coordinateToMeters = cellDiameter_m;
    case 'microns'
        coordinateToMeters = 1e-6;
    case 'meters'
        coordinateToMeters = 1;
end

switch displayMode
    case 'cells'
        displayFromMeters = 1 / cellDiameter_m;
        displayLengthUnit = 'cell diameters';
        displayVelocityUnit = 'cells/s';
    case 'microns'
        displayFromMeters = 1e6;
        displayLengthUnit = 'µm';
        displayVelocityUnit = 'µm/s';
end

units.schemaVersion = 1;
units.internalLengthUnit = 'm';
units.internalTimeUnit = 's';
units.internalForceUnit = 'N';
units.coordinateInputMode = coordinateInputMode;
units.displayMode = displayMode;
units.cellDiameter_um = cellDiameter_um;
units.cellDiameter_m = cellDiameter_m;
units.coordinateToMeters = coordinateToMeters;
units.displayFromMeters = displayFromMeters;
units.displayLengthUnit = displayLengthUnit;
units.displayVelocityUnit = displayVelocityUnit;
units.stiffnessUnit = 'N/m';
units.dampingUnit = 'N*s/m';
units.massUnit = 'kg';
units.strainUnit = 'dimensionless';
units.damageUnit = 'dimensionless';
end


function mode = normalizeUnitMode(mode, allowedModes, modeDescription)
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('follicle_rupture:invalidUnitMode', ...
        '%s mode must be a character vector or scalar string.', modeDescription);
end

mode = lower(strtrim(char(mode)));
if ~ismember(mode, allowedModes)
    error('follicle_rupture:invalidUnitMode', ...
        'Unknown %s mode "%s". Valid modes are: %s.', ...
        modeDescription, mode, strjoin(allowedModes, ', '));
end
end
