function follicle_rupture_msd_fixed
% Follicle rupture shell model
% State: X = [q; v; sigma; zeta; D; A_edge]
%   q, v        node displacement and velocity
%   sigma,zeta  viscoelastic edge states
%   D           edge damage in [0,1], where 1 means failed
%   A_edge      active rest-length recruitment in [0,1]

clc; close all;

rng(42);

%% USER-PARAMETERS
% Simulation resolution / scale
N = 300;
mVal = 1;                 % smaller mass reduces inertial coasting
useExperimentalNodeData = false;
experimentalNodeDataFileName = 'nodeDataSample1.xlsx';
experimentalCoordinateScale = 1e-6;  % spreadsheet coordinates are in microns
experimentalCenterMode = 'median';   % 'median', 'mean', or 'none'
experimentalOrientMode = 'principal-axes'; % 'principal-axes' aligns long axis to model x-axis
experimentalLongAxisSign = -1;        % set 1 to flip experimental anterior/posterior ends
experimentalSurfaceConstraintMode = 'project-ellipsoid'; % project forces + damp normal velocity, but do not snap nodes
experimentalAnchorOn = false;          % measured node clouds should not have a frozen anterior shell
experimentalEdgeMode = 'knn';          % 'knn' connects all measured nodes; 'hull' uses convex hull only
experimentalNeighborCount = 6;         % local neighbors per measured node for 'knn' edge mode

% Passive edge material parameters
kOuter_range = [1, 1];
cOuter_range = [1, 1];

% Time integration and saved output
tEnd  = 1800;                    % s
numSave = 300;                 % saved output frames; solver still adapts internally
runPostprocess = true;         % set true to generate diagnostic plots after saving
makeRuptureMotionVideo = true;  % set true to export an MP4 after saving
odeSolver = 'ode45';           % instant damage is stiff; ode15s handles it better
progressEveryPhysSec = 5;      % console progress interval in physical seconds

ruptureVideoName = 'follicleRupture_motion.mp4';

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
ruptureVideoFile = fullfile(outputDir, ruptureVideoName);
experimentalNodeDataFile = fullfile(inputNodeDataDir, experimentalNodeDataFileName);

%% Nodes
if useExperimentalNodeData
    [nodeCoords, experimentalCellIds] = loadExperimentalNodeCoords( ...
        experimentalNodeDataFile, experimentalCoordinateScale, experimentalCenterMode);
    nodeCoords = orientExperimentalNodeCoords( ...
        nodeCoords, experimentalOrientMode, experimentalLongAxisSign);
    N = size(nodeCoords, 1);
    [a, b, c] = fitAxisAlignedEllipsoidSemiAxes(nodeCoords);
    L_scale = median([a b c]);  % representative experimental radius scale
    geometryTitle = 'Experimental surface nodes';
else
    a = 500e-6; b = 150e-6; c = 150e-6;
    L_scale = b;              % transverse follicle radius scale
    nodeCoords = genEllipsoidSurfaceNodes(N, a, b, c, [0 0 0]);
    experimentalCellIds = [];
    geometryTitle = 'Ellipsoid surface nodes';
end

param.ell_a = a; param.ell_b = b; param.ell_c = c;
param.useExperimentalNodeData = useExperimentalNodeData;
param.experimentalNodeDataFile = experimentalNodeDataFile;
param.experimentalCoordinateScale = experimentalCoordinateScale;
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
param.postSign = -1;

% Surface penalty + normal damping
param.k_surf = 5e-6;     % [N] penalty strength
param.c_surf = 5e-2;     % [N*s/m] normal damping
param.c_drag = 1e-4;     % [N*s/m] whole-node damping; higher = less post-cut coasting
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

if runPostprocess
    figure; scatter3(nodeCoords(:,1),nodeCoords(:,2),nodeCoords(:,3),10,'filled');
    axis equal; grid on; title(geometryTitle);
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
edgeNodeIncidence = sparse([param.edge_i; param.edge_j], [(1:Ne)'; (1:Ne)'], 1, N, Ne);
nodeDegree = full(sum(edgeNodeIncidence, 2));
fprintf('[MESH] nodes = %d | edges = %d | edge length median/p95/max = %.3g/%.3g/%.3g um | isolated nodes = %d\n', ...
    N, Ne, median(param.L0_edge0) * 1e6, prctile(param.L0_edge0, 95) * 1e6, ...
    max(param.L0_edge0) * 1e6, nnz(nodeDegree == 0));
if any(nodeDegree == 0)
    warning('follicle_rupture:isolatedNodes', ...
        'Mesh contains %d isolated nodes that cannot move mechanically.', nnz(nodeDegree == 0));
end
param.edgeNeighbor = spones(edgeNodeIncidence' * edgeNodeIncidence);
param.edgeNeighbor = param.edgeNeighbor - speye(Ne);
param.edgeNeighbor = spones(param.edgeNeighbor);
param.edgeNeighborCount = full(sum(param.edgeNeighbor, 2));

% Baseline per-edge passive model parameters
param.k1_edge = kOuter_range(1) * ones(Ne,1);
param.k2_edge = kOuter_range(2) * ones(Ne,1);
param.c_edge  = cOuter_range(1) * ones(Ne,1);
param.c_series_edge = cOuter_range(2) * ones(Ne,1);

% Passive 1D connection model between cells.
% Options: 'spring', 'damper', 'maxwell', 'kelvin-voigt', 'zener', 'burgers'.
% You may also set param.connectionModel_edge to a per-edge cell/string array.
param.connectionModel = 'burgers';
param.connectionModelKey = normalizeConnectionModel(param.connectionModel);
param.connectionModelId_edge = connectionModelIds(param, Ne);

%% Spatial localization masks (main-body vs tiny posterior tip)
edgeMid = 0.5*(nodeCoords(E(:,1),:) + nodeCoords(E(:,2),:) );
ax = param.postAxis;
sg = param.postSign;
param.edgeMid0 = edgeMid;

sNodes = sg * nodeCoords(:,ax);
sMid   = sg * edgeMid(:,ax);

param.mainPct = 98.5;      % posterior-most ~1.5% has motors OFF
param.tipPct  = 98.5;    % posterior-most ~1.5% defines tiny autonomous rupture tip

sCutMain = prctile(sNodes, param.mainPct);
sCutTip  = prctile(sNodes, param.tipPct);

param.motorMask = double(sMid < sCutMain);    % motors ON in main-body
param.tipMask   = double(sMid >= sCutTip);    % tiny posterior cap edges

%% Damage region selection
param.damageMode = 'tipmost-edge';     % 'seam', 'tip', or 'tipmost-edge'
param.seamWidthFrac = 0.04;   % fraction of posterior-score range

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
    otherwise
        error('Unknown damageMode "%s". Use "seam", "tip", or "tipmost-edge".', param.damageMode);
end

fprintf('[MASK] motor edges = %d/%d | tip edges = %d | damage edges = %d\n', ...
    nnz(param.motorMask), Ne, nnz(param.tipMask), nnz(param.damageMask));
if nnz(param.damageMask) == 0
    error('Tiny tip selection produced zero damage edges. Lower param.tipPct slightly.');
end

%% Active mechanics parameters
% Shared active timing
param.act_on      = false;    % force-based actomyosin off; rest-length shortening drives contraction
param.tActOn      = 0.5;     % [s]
param.tActTau     = 0.1;     % [s]

% Optional legacy force-based actomyosin mechanics.
% These remain available for comparison but are gated by param.act_on.
param.Tact_main   = 1e-7;   % [N] main-body active edge tension (TUNE)
param.Tact_tip    = 1e-8;    % [N] light tip-edge cohesion; keeps tip cap from overtaking
param.activeForceMode = 'posterior-traction'; % 'edge-tension' or 'posterior-traction'
param.balanceActiveTraction = true;           % remove net translation from polar traction
param.activeStrokeLimit = 100e-6;              % [m] finite motor stroke before active pull stalls
param.activeStrokeTau = 2e-6;                 % [m] smoothness of stroke cutoff

% Optional smooth phase-delayed wave for force-based actomyosin mechanics.
% Posterior motor edges turn on first; more anterior edges turn on later,
% giving a continuous traveling pull rather than discrete contractile rings.
param.activeWave_on = false;
param.activeWaveMode = 'continuous';  % 'continuous' or 'rings'
param.activeWaveStartPct = 98.5;        % posterior edge of motor domain
param.activeWaveDelayMax = 240;       % [s] delay from wave start to anterior end
param.activeWavePulse_on = true;      % each location turns off after its pulse
param.activeWavePulseDuration = 10;  % [s] active duration at each location
param.activeWaveOffTau = 3;           % [s] smooth shutoff timescale

% Optional orientation filter for active edge masks.
param.axisSelectiveActomyosin = false;
param.axisAlignmentMin = 0.65;  % 1=perfectly parallel, 0=no orientation filter

if param.axisSelectiveActomyosin
    edgeVec = nodeCoords(E(:,2),:) - nodeCoords(E(:,1),:);
    edgeLen = sqrt(sum(edgeVec.^2, 2));
    axisAlign = abs(edgeVec(:,param.postAxis)) ./ max(edgeLen, 1e-18);

    param.axisAlignment = axisAlign;
    param.motorMask = param.motorMask .* double(axisAlign >= param.axisAlignmentMin);
end

% Optional actomyosin-driven preferred junction shortening.
% This contracts selected epithelial edges by shortening their active rest length,
% while the passive connection model still supplies the stress response.
param.activeRestLength_on = true;
param.activeShorteningMax = 0.32;       % maximum fractional rest-length reduction
param.activeShorteningTau = 1;        % [s] contraction timescale
param.activeShorteningActivation = 'exponential'; % 'exponential', 'sigmoid', or 'accelerating'
param.activeShorteningAccelPower = 2; % >1 starts slow then speeds up; 2 = quadratic ramp
param.activeShorteningSigmoidDelay = 0; % [s] added to each edge's activation time for sigmoid midpoint
param.activeShorteningSigmoidTau = param.activeShorteningTau; % [s] sigmoid switch sharpness
param.activeShorteningMask = param.motorMask; % per-edge weights; can be any value in [0,1]

% Optional mechanics-only local recruitment for active rest-length shortening.
% When enabled, shortening spreads through edge-neighbor coupling instead of
% the prescribed activeRestLengthShortening(t,param) wave below.
param.activeRecruitment_on = true;
param.k_seed = 1/60;                         % [1/s] autonomous seed activation
param.k_neighbor = 1/90;                    % [1/s] neighbor recruitment
param.k_strainRateRecruit = 1/120;            % [1/s] combined neighbor/strain-rate recruitment
param.recruitmentThreshold = 0.03;             % suppress tiny neighbor activation leaks
param.recruitmentStrainRateThreshold = 1e-4; % [1/s] ignore tiny strain-rate noise
param.recruitmentStrainRateScale = 1e-3;     % [1/s] strain-rate signal saturation scale
param.recruitmentStrainRateGain = 2;         % strain-rate amplification near recruitment front
param.seedMask = param.tipMask;           % default posterior/tip seed
param.recruitmentMask = param.motorMask;     % recruitable active motor domain
param.recruitmentSignalMode = 'neighborA_plus_strainRate';   % 'neighborA' or 'neighborA_plus_strainRate'

% Optional wave for active rest-length shortening. This is independent of
% force-based activeWave_on, so preferred-length contraction can travel even
% when external/nodal active actomyosin forces are disabled.
param.activeRestLengthWave_on = false;
param.activeRestLengthWaveMode = 'continuous'; % 'continuous' or 'rings'
param.activeRestLengthWaveStartPct = 98.5;
param.activeRestLengthWaveDelayMax = 900;      % [s] delay from wave start to anterior end
param.activeRestLengthWavePulse_on = false;   % false = edges stay shortened after their wave reaches them
param.activeRestLengthWavePulseDuration = 1;  % [s] shortening duration at each location
param.activeRestLengthWaveOffTau = 1;          % [s] smooth recovery timescale

% Optional weak anterior anchoring. This is a soft positional boundary
% condition to stabilize directional epithelial retraction, not a rigid pin.
if useExperimentalNodeData
    param.anchor_on = experimentalAnchorOn;
else
    param.anchor_on = true;
end
param.anchorPct = 60;          % anterior-most percentile of nodes
param.k_anchor = 10;        % weak positional tether stiffness
param.c_anchor = 1;        % anchor damping

sCutAnchor = prctile(sNodes, param.anchorPct);
param.anchorNodeMask = double(sNodes <= sCutAnchor);
fprintf('[ANCHOR] anchor nodes = %d/%d\n', nnz(param.anchorNodeMask), N);

param = configureActiveForceWave(param, sNodes, sMid, Ne);
param = configureActiveRestLengthWave(param, sNodes, sMid, Ne);

fprintf('[ACTIVE REST LENGTH] enabled = %d | max shortening = %.3g | tau = %.3g s\n', ...
    isfield(param,'activeRestLength_on') && param.activeRestLength_on, ...
    param.activeShorteningMax, param.activeShorteningTau);

%% Damage / rupture parameters
param.tipDamage_on  = true;  % off because the tiny tip now ruptures autonomously

param.sig_fail_tip  = 1e-9; % [N] tip failure threshold (TUNE)
param.tau_dmg       = 0.005;   % [s] damage timescale (TUNE)
param.dmg_pow       = 2;     % sharpness

% Force the tiny tip to fail on a clock, independent of contractile pulling.
param.instantDamage_on   = true;
param.instantDamageTime  = 0;
param.instantDamageTau   = 1e-4;  % [s] smaller = closer to immediate cut

param.eps_gate_on   = true;  % only damage if opening
param.eps_gate      = 0.001;

% Optional: gentle tip priming (adds small damage even before overload)
param.tipPriming_on   = false; % off because instantDamage_on handles autonomous rupture
param.tipPrimingRate  = 0.05;  % [1/s]
param.tipPriming_on_t = 0.5;   % [s]
param.tipPriming_tau  = 0.3;   % [s]

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
param.runPostprocess = runPostprocess;
param.odeSolver = odeSolver;
param.progressEveryPhysSec = progressEveryPhysSec;

%% ================= State sizes =================
nPos = 3*N;
nVel = 3*N;

param.nPos = nPos;
param.nVel = nVel;

% X = [q; v; sigma; zeta; D; A_edge]
nX = nPos + nVel + Ne + Ne + Ne + Ne;

% Edge state block indices for extracting saved time courses
iD1 = nPos + nVel + 2*Ne + 1;
iD2 = iD1 + Ne - 1;
iA1 = iD2 + 1;
iA2 = iA1 + Ne - 1;

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

% extract edge state time courses
D_time = XDimQ_all(:, iD1:iD2);
A_time = min(max(XDimQ_all(:, iA1:iA2), 0), 1);

save(simFile, ...
     'nodeCoords','edgeList','tSol','Y','D_time','A_time','param','-v7.3');

if makeRuptureMotionVideo
    try
        export_follicle_rupture_motion_video(simFile, ruptureVideoFile);
    catch ME
        warning('follicle_rupture:videoExportFailed', ...
            'Video export failed (%s): %s', ME.identifier, ME.message);
    end
end

if param.runPostprocess
    try
        confirm_actomyosin_working(simFile);
    catch ME
        warning('follicle_rupture:postprocessFailed', ...
            'Postprocessing failed (%s): %s', ME.identifier, ME.message);
    end
end

end % main


%% ===================== Active wave setup helpers =====================
function param = configureActiveForceWave(param, sNodes, sMid, Ne)
if isfield(param,'activeWave_on') && param.activeWave_on
    sWaveStart = prctile(sNodes, param.activeWaveStartPct);
    motorEdgeForDelay = param.motorMask(:) > 0;
    if any(motorEdgeForDelay)
        sAnterior = min(sMid(motorEdgeForDelay));
    else
        sAnterior = min(sMid);
    end

    phaseContinuous = max(0, sWaveStart - sMid) ./ max(sWaveStart - sAnterior, eps);
    phaseContinuous = min(phaseContinuous, 1);

    if isfield(param,'activeWaveMode') && strcmpi(param.activeWaveMode, 'rings')
        if ~isfield(param,'activeWaveNumRings') || isempty(param.activeWaveNumRings)
            param.activeWaveNumRings = 2;
        end
        nRings = max(1, round(param.activeWaveNumRings));
        if nRings == 1
            ringId = ones(Ne,1);
        else
            ringId = floor(phaseContinuous * nRings) + 1;
            ringId = min(max(ringId, 1), nRings);
        end

        if isfield(param,'activeRingStartTimes') && numel(param.activeRingStartTimes) == nRings
            ringStartTimes = param.activeRingStartTimes(:);
        elseif nRings == 1
            ringStartTimes = param.tActOn;
        else
            ringStartTimes = linspace(param.tActOn, param.tActOn + param.activeWaveDelayMax, nRings).';
        end

        if isfield(param,'activeRingDuration') && numel(param.activeRingDuration) == nRings
            ringDurations = param.activeRingDuration(:);
        else
            ringDurations = param.activeWavePulseDuration * ones(nRings,1);
        end

        param.activeWaveNumRings = nRings;
        param.activeRingId_edge = ringId;
        param.activeRingStartTimes = ringStartTimes;
        param.activeRingStopTimes = ringStartTimes + ringDurations;
        param.tActOn_edge = ringStartTimes(ringId);
        param.tActOff_edge = param.activeRingStopTimes(ringId);

        motorCountsByRing = zeros(nRings,1);
        if any(motorEdgeForDelay)
            motorCountsByRing = accumarray(ringId(motorEdgeForDelay), 1, [nRings 1], @sum, 0);
        end
        fprintf('[ACTIVE WAVE] mode = rings | rings = %d | motor edges per ring: %s | starts: %s\n', ...
            nRings, sprintf('%d ', motorCountsByRing), sprintf('%.1f ', ringStartTimes));
    else
        param.activeWaveMode = 'continuous';
        param.activeRingId_edge = [];
        param.tActOn_edge = param.tActOn + param.activeWaveDelayMax * phaseContinuous;
        param.tActOff_edge = param.tActOn_edge + param.activeWavePulseDuration;

        if any(motorEdgeForDelay)
            fprintf('[ACTIVE WAVE] mode = continuous | delay range = %.1f-%.1f s | pulse duration = %.1f s\n', ...
                min(param.tActOn_edge(motorEdgeForDelay)), ...
                max(param.tActOn_edge(motorEdgeForDelay)), ...
                param.activeWavePulseDuration);
        else
            fprintf('[ACTIVE WAVE] mode = continuous | no motor edges selected\n');
        end
    end

    param.activeWavePhase_edge = phaseContinuous;
    param.activeWaveDelay_edge = param.tActOn_edge - param.tActOn;
else
    param.activeWaveDelay_edge = zeros(Ne,1);
    param.activeWavePhase_edge = zeros(Ne,1);
    param.activeRingId_edge = ones(Ne,1);
    param.tActOn_edge = param.tActOn * ones(Ne,1);
    param.tActOff_edge = inf(Ne,1);
end
end


function param = configureActiveRestLengthWave(param, sNodes, sMid, Ne)
if isfield(param,'activeRestLengthWave_on') && param.activeRestLengthWave_on
    sRestWaveStart = prctile(sNodes, param.activeRestLengthWaveStartPct);
    restMaskForDelay = param.activeShorteningMask(:);
    if isscalar(restMaskForDelay)
        restEdgeForDelay = repmat(restMaskForDelay > 0, Ne, 1);
    elseif numel(restMaskForDelay) == Ne
        restEdgeForDelay = restMaskForDelay > 0;
    else
        error('param.activeShorteningMask must be scalar or have one value per edge.');
    end

    if any(restEdgeForDelay)
        sAnteriorRest = min(sMid(restEdgeForDelay));
    else
        sAnteriorRest = min(sMid);
    end

    restPhaseContinuous = max(0, sRestWaveStart - sMid) ./ max(sRestWaveStart - sAnteriorRest, eps);
    restPhaseContinuous = min(restPhaseContinuous, 1);

    if isfield(param,'activeRestLengthWaveMode') && strcmpi(param.activeRestLengthWaveMode, 'rings')
        if ~isfield(param,'activeRestLengthWaveNumRings') || isempty(param.activeRestLengthWaveNumRings)
            param.activeRestLengthWaveNumRings = 2;
        end
        nRestRings = max(1, round(param.activeRestLengthWaveNumRings));
        if nRestRings == 1
            restRingId = ones(Ne,1);
        else
            restRingId = floor(restPhaseContinuous * nRestRings) + 1;
            restRingId = min(max(restRingId, 1), nRestRings);
        end

        if isfield(param,'activeRestLengthRingStartTimes') && numel(param.activeRestLengthRingStartTimes) == nRestRings
            restRingStartTimes = param.activeRestLengthRingStartTimes(:);
        elseif nRestRings == 1
            restRingStartTimes = param.tActOn;
        else
            restRingStartTimes = linspace(param.tActOn, param.tActOn + param.activeRestLengthWaveDelayMax, nRestRings).';
        end

        if isfield(param,'activeRestLengthRingDuration') && numel(param.activeRestLengthRingDuration) == nRestRings
            restRingDurations = param.activeRestLengthRingDuration(:);
        else
            restRingDurations = param.activeRestLengthWavePulseDuration * ones(nRestRings,1);
        end

        param.activeRestLengthWaveNumRings = nRestRings;
        param.activeRestLengthRingId_edge = restRingId;
        param.activeRestLengthRingStartTimes = restRingStartTimes;
        param.activeRestLengthRingStopTimes = restRingStartTimes + restRingDurations;
        param.tActiveRestLengthOn_edge = restRingStartTimes(restRingId);
        param.tActiveRestLengthOff_edge = param.activeRestLengthRingStopTimes(restRingId);

        restCountsByRing = zeros(nRestRings,1);
        if any(restEdgeForDelay)
            restCountsByRing = accumarray(restRingId(restEdgeForDelay), 1, [nRestRings 1], @sum, 0);
        end
        fprintf('[ACTIVE REST LENGTH WAVE] mode = rings | rings = %d | shortening edges per ring: %s | starts: %s\n', ...
            nRestRings, sprintf('%d ', restCountsByRing), sprintf('%.1f ', restRingStartTimes));
    else
        param.activeRestLengthWaveMode = 'continuous';
        param.activeRestLengthRingId_edge = [];
        param.tActiveRestLengthOn_edge = param.tActOn + param.activeRestLengthWaveDelayMax * restPhaseContinuous;
        param.tActiveRestLengthOff_edge = param.tActiveRestLengthOn_edge + param.activeRestLengthWavePulseDuration;

        if any(restEdgeForDelay)
            fprintf('[ACTIVE REST LENGTH WAVE] mode = continuous | delay range = %.1f-%.1f s | pulse duration = %.1f s\n', ...
                min(param.tActiveRestLengthOn_edge(restEdgeForDelay)), ...
                max(param.tActiveRestLengthOn_edge(restEdgeForDelay)), ...
                param.activeRestLengthWavePulseDuration);
        else
            fprintf('[ACTIVE REST LENGTH WAVE] mode = continuous | no shortening edges selected\n');
        end
    end

    if ~isfield(param,'activeRestLengthWavePulse_on') || ~param.activeRestLengthWavePulse_on
        param.tActiveRestLengthOff_edge = inf(Ne,1);
    end
    param.activeRestLengthWavePhase_edge = restPhaseContinuous;
    param.activeRestLengthWaveDelay_edge = param.tActiveRestLengthOn_edge - param.tActOn;
else
    param.activeRestLengthWaveDelay_edge = zeros(Ne,1);
    param.activeRestLengthWavePhase_edge = zeros(Ne,1);
    param.activeRestLengthRingId_edge = ones(Ne,1);
    param.tActiveRestLengthOn_edge = param.tActOn * ones(Ne,1);
    param.tActiveRestLengthOff_edge = inf(Ne,1);
end
end


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
iA2  = iD2 + Ne;

q_dim   = X_dim(1:iQ2);
v_dim   = X_dim(iQ2+1:iV2);
sig_dim = X_dim(iV2+1:iS2);
zet_dim = X_dim(iS2+1:iZ2);
D_dim   = X_dim(iZ2+1:iD2);
A_dim   = X_dim(iD2+1:iA2);

% to physical
q_phys   = q_dim * scale.L;
v_phys   = v_dim * (scale.L/scale.tau);
sig_phys = sig_dim * (scale.k2*scale.L);
zet_phys = zet_dim * scale.L;
D_phys   = D_dim;
A_phys   = A_dim;
t_phys   = t_dim * scale.tau;

dX_phys = ShellRHS_physical(t_phys, [q_phys; v_phys; sig_phys; zet_phys; D_phys; A_phys], param);

% back to ND
dqdt_dim   = (dX_phys(1:nPos)/scale.L) * scale.tau;
dvdt_dim   = scale.epsilon * ((dX_phys(nPos+1:nPos+nVel) / (scale.L/scale.tau)) * scale.tau);

dsig_dim   = (dX_phys(nPos+nVel+1 : nPos+nVel+Ne) / (scale.k2*scale.L)) * scale.tau;
dzet_dim   = (dX_phys(nPos+nVel+Ne+1 : nPos+nVel+2*Ne) / scale.L) * scale.tau;

off = nPos + nVel + 2*Ne;
dD_dim     = dX_phys(off+(1:Ne)) * scale.tau;
dA_dim     = dX_phys(off+Ne+(1:Ne)) * scale.tau;

dXdim = [dqdt_dim; dvdt_dim; dsig_dim; dzet_dim; dD_dim; dA_dim];
end


%% ========================= PHYSICAL RHS =========================
function dxdt = ShellRHS_physical(t, x, param)
% X = [ q(3N); v(3N); sigma(Ne); zeta(Ne); D(Ne); A_edge(Ne) ]

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
iA1 = iD2+1;             iA2 = iA1 + Ne - 1;

q     = x(iQ1:iQ2);
v     = x(iV1:iV2);
sigma = x(iS1:iS2);
zeta  = x(iZ1:iZ2);
D     = x(iD1:iD2);
A_edge = x(iA1:iA2);

% clamp bounded edge states safely
D = min(max(D,0),1);
A_edge = min(max(A_edge,0),1);

qMat = reshape(q,3,[])';
vMat = reshape(v,3,[])';
pos  = param.nodeCoords + qMat;
iEdge = param.edge_i;
jEdge = param.edge_j;

surfaceConstraintMode = 'ellipsoid';
if isfield(param,'surfaceConstraintMode') && ~isempty(param.surfaceConstraintMode)
    surfaceConstraintMode = lower(strtrim(char(param.surfaceConstraintMode)));
end
if ~ismember(surfaceConstraintMode, {'ellipsoid', 'project-ellipsoid', 'none', 'off'})
    error('Unknown surfaceConstraintMode "%s". Use "ellipsoid", "project-ellipsoid", or "none".', surfaceConstraintMode);
end
useEllipsoidNormals = any(strcmp(surfaceConstraintMode, {'ellipsoid', 'project-ellipsoid'}));
useEllipsoidPenalty = strcmp(surfaceConstraintMode, 'ellipsoid');
useEllipsoidNormalDamping = useEllipsoidNormals;

% Ellipsoid signed error and outward normals at current node positions.
% Experimental node clouds can keep their measured shape by disabling this.
if useEllipsoidNormals
    a = param.ell_a; b = param.ell_b; c = param.ell_c;
    xx = pos(:,1); yy = pos(:,2); zz = pos(:,3);
    fSurf = (xx./a).^2 + (yy./b).^2 + (zz./c).^2 - 1;
    nxv = 2*xx/(a^2);
    nyv = 2*yy/(b^2);
    nzv = 2*zz/(c^2);
    nrm = sqrt(nxv.^2 + nyv.^2 + nzv.^2) + 1e-18;
    nSurf = [nxv./nrm, nyv./nrm, nzv./nrm];
else
    fSurf = zeros(N,1);
    nSurf = zeros(N,3);
end

% precompute active edge tension (additive)
TactVec = activeTension(t, param);

% Edge kinematics with optional active rest-length shortening.
rij = pos(jEdge,:) - pos(iEdge,:);
L = sqrt(sum(rij.^2, 2));
Lsafe = max(L, 1e-18);
n_hat = rij ./ Lsafe;
n_hat(L < 1e-18,:) = 0;

L0_ref = max(param.L0_edge0(:), 1e-18);
geometricStretchDot = sum((vMat(jEdge,:) - vMat(iEdge,:)) .* n_hat, 2);
strainRate = geometricStretchDot ./ L0_ref;
dA = activeRecruitmentRHS(A_edge, param, strainRate);
[activeFrac, activeFracDot] = activeRestLengthShorteningForState(t, param, A_edge, dA);
L0_active = L0_ref .* (1 - activeFrac);
L0 = max(L0_active, 1e-18);
stretch = L - L0;
L0_dot = -L0_ref .* activeFracDot;
stretch_dot = geometricStretchDot - L0_dot;
eps_e = stretch ./ L0_ref;

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

FactiveTraction = activeTractionForces(t, param, N, pos);
if useEllipsoidNormals && isfield(param,'projectEdgeForcesToSurface') && param.projectEdgeForcesToSurface
    FactiveTraction = FactiveTraction - sum(FactiveTraction .* nSurf, 2) .* nSurf;
end
Fnet = Fnet + FactiveTraction;

dD = zeros(Ne,1);

% ---- Optional instantaneous cut of the damage ring ----
if isfield(param,'instantDamage_on') && param.instantDamage_on && ...
        isfield(param,'damageMask') && any(param.damageMask > 0)
    if t >= param.instantDamageTime
        damageEdges = (param.damageMask(:) > 0);
        tauInstant = max(param.instantDamageTau, 1e-8);
        dD(damageEdges) = dD(damageEdges) + (1 - D(damageEdges)) ./ tauInstant;
    end
end

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
dxdt = [dqdt; dvdt; dsigma; dzeta; dD; dA];
end


%% ---------------- Active tension helper ----------------
function T = activeTension(t, param)
Ne = numel(param.L0_edge0);
T = zeros(Ne,1);

if ~isfield(param,'act_on') || ~param.act_on
    return;
end

A = activeEdgeActivation(t, param);

if ~isfield(param,'activeForceMode') || ~strcmpi(param.activeForceMode, 'posterior-traction')
    if isfield(param,'motorMask')
        T = T + (param.Tact_main .* A) .* param.motorMask(:);
    else
        T = T + param.Tact_main .* A;
    end
end

if isfield(param,'tipMask')
    T = T + (param.Tact_tip .* A) .* param.tipMask(:);
end
end


function Fnode = activeTractionForces(t, param, N, pos)
Fnode = zeros(N,3);

if ~isfield(param,'act_on') || ~param.act_on
    return;
end
if ~isfield(param,'activeForceMode') || ~strcmpi(param.activeForceMode, 'posterior-traction')
    return;
end

A = activeEdgeActivation(t, param);
if isfield(param,'motorMask')
    edgeAmp = param.Tact_main .* A .* param.motorMask(:);
else
    edgeAmp = param.Tact_main .* A;
end

if isfield(param,'activeStrokeLimit') && isfinite(param.activeStrokeLimit) && param.activeStrokeLimit > 0 && ...
        isfield(param,'edgeMid0') && nargin >= 4
    currentEdgeMid = 0.5 * (pos(param.edge_i,:) + pos(param.edge_j,:));
    posteriorDisp = param.postSign * (currentEdgeMid(:,param.postAxis) - param.edgeMid0(:,param.postAxis));

    strokeTau = 2e-6;
    if isfield(param,'activeStrokeTau') && isfinite(param.activeStrokeTau) && param.activeStrokeTau > 0
        strokeTau = param.activeStrokeTau;
    end
    strokeFactor = 1 ./ (1 + exp((posteriorDisp - param.activeStrokeLimit) ./ strokeTau));
    edgeAmp = edgeAmp .* strokeFactor;
end

posteriorDir = zeros(1,3);
posteriorDir(param.postAxis) = param.postSign;

nodeAmp = 0.5 * (param.edgeAccumI * edgeAmp + param.edgeAccumJ * edgeAmp);
Fnode = nodeAmp .* posteriorDir;

if isfield(param,'balanceActiveTraction') && param.balanceActiveTraction
    Fnode = Fnode - mean(Fnode, 1);
end
end


function dA = activeRecruitmentRHS(A_edge, param, strainRate)
Ne = numel(param.L0_edge0);
dA = zeros(Ne,1);

if ~isfield(param,'activeRecruitment_on') || ~param.activeRecruitment_on
    return;
end
if nargin < 3 || isempty(strainRate)
    strainRate = zeros(Ne,1);
else
    strainRate = strainRate(:);
end

signalMode = 'neighborA';
if isfield(param,'recruitmentSignalMode') && ~isempty(param.recruitmentSignalMode)
    signalMode = lower(strtrim(char(param.recruitmentSignalMode)));
end

switch signalMode
    case 'neighbora'
        neighborA = (param.edgeNeighbor * A_edge) ./ max(param.edgeNeighborCount, 1);
        theta = param.recruitmentThreshold;
        neighborSignal = max(neighborA - theta, 0) ./ max(1 - theta, eps);
        recruitmentRate = param.k_neighbor;
    case 'neighbora_plus_strainrate'
        neighborA = (param.edgeNeighbor * A_edge) ./ max(param.edgeNeighborCount, 1);
        theta = param.recruitmentThreshold;
        neighborA = max(neighborA - theta, 0) ./ max(1 - theta, eps);

        threshold = param.recruitmentStrainRateThreshold;
        scale = param.recruitmentStrainRateScale;
        strainRateSignal = max(abs(strainRate) - threshold, 0);
        strainRateSignal = min(strainRateSignal ./ max(scale, eps), 1);

        alpha = param.recruitmentStrainRateGain;
        neighborSignal = neighborA .* (1 + alpha * strainRateSignal);
        neighborSignal = min(neighborSignal, 1);
        recruitmentRate = param.k_strainRateRecruit;
    otherwise
        error('Unknown recruitmentSignalMode "%s". Use "neighborA" or "neighborA_plus_strainRate".', signalMode);
end

seedMask = edgeMaskVector(param, 'seedMask', Ne, zeros(Ne,1));
recruitmentMask = edgeMaskVector(param, 'recruitmentMask', Ne, ones(Ne,1));

dA = param.k_seed .* seedMask .* (1 - A_edge) + ...
     recruitmentRate .* neighborSignal .* recruitmentMask .* (1 - A_edge);
end


function [activeFrac, activeFracDot] = activeRestLengthShorteningForState(t, param, A_edge, dA)
Ne = numel(param.L0_edge0);
activeFrac = zeros(Ne,1);
activeFracDot = zeros(Ne,1);

if ~isfield(param,'activeRestLength_on') || ~param.activeRestLength_on
    return;
end

if ~isfield(param,'activeRecruitment_on') || ~param.activeRecruitment_on
    [activeFrac, activeFracDot] = activeRestLengthShortening(t, param);
    return;
end

mask = edgeMaskVector(param, 'activeShorteningMask', Ne, ones(Ne,1));
rawFrac = param.activeShorteningMax .* A_edge .* mask;
rawFracDot = param.activeShorteningMax .* dA .* mask;

activeFrac = min(max(rawFrac, 0), 0.95);
activeFracDot = rawFracDot;
activeFracDot(rawFrac <= 0 | rawFrac >= 0.95) = 0;
end


function mask = edgeMaskVector(param, fieldName, Ne, defaultMask)
if isfield(param, fieldName) && ~isempty(param.(fieldName))
    mask = param.(fieldName)(:);
else
    mask = defaultMask(:);
end

if isscalar(mask)
    mask = mask * ones(Ne,1);
elseif numel(mask) ~= Ne
    error('param.%s must be scalar or have one value per edge.', fieldName);
end
end


function [activeFrac, activeFracDot] = activeRestLengthShortening(t, param)
Ne = numel(param.L0_edge0);
activeFrac = zeros(Ne,1);
activeFracDot = zeros(Ne,1);

if ~isfield(param,'activeRestLength_on') || ~param.activeRestLength_on
    return;
end

if isfield(param,'tActiveRestLengthOn_edge') && numel(param.tActiveRestLengthOn_edge) == Ne
    tOn = param.tActiveRestLengthOn_edge(:);
else
    tOn = param.tActOn * ones(Ne,1);
end

activationMode = 'exponential';
if isfield(param,'activeShorteningActivation') && ~isempty(param.activeShorteningActivation)
    activationMode = lower(strtrim(char(param.activeShorteningActivation)));
end

switch activationMode
    case 'exponential'
        tau = max(param.activeShorteningTau, eps);
        dtOn = t - tOn;
        isOn = dtOn >= 0;
        A = zeros(Ne,1);
        Adot = zeros(Ne,1);
        A(isOn) = 1 - exp(-dtOn(isOn) ./ tau);
        Adot(isOn) = exp(-dtOn(isOn) ./ tau) ./ tau;
    case {'accelerating', 'power', 'reverse-exponential'}
        tau = max(param.activeShorteningTau, eps);
        p = 2;
        if isfield(param,'activeShorteningAccelPower') && isfinite(param.activeShorteningAccelPower)
            p = max(param.activeShorteningAccelPower, 1);
        end
        dtOn = t - tOn;
        u = max(dtOn ./ tau, 0);
        isRamping = u < 1;
        A = min(u .^ p, 1);
        Adot = zeros(Ne,1);
        Adot(isRamping & u > 0) = (p ./ tau) .* u(isRamping & u > 0).^(p - 1);
    case {'sigmoid', 'delayed-sigmoid'}
        delay = 0;
        if isfield(param,'activeShorteningSigmoidDelay') && isfinite(param.activeShorteningSigmoidDelay)
            delay = param.activeShorteningSigmoidDelay;
        end
        sigmoidTau = param.activeShorteningTau;
        if isfield(param,'activeShorteningSigmoidTau') && isfinite(param.activeShorteningSigmoidTau) && param.activeShorteningSigmoidTau > 0
            sigmoidTau = param.activeShorteningSigmoidTau;
        end
        sigmoidTau = max(sigmoidTau, eps);
        tMid = tOn + delay;
        A = 1 ./ (1 + exp(-(t - tMid) ./ sigmoidTau));
        Adot = (A .* (1 - A)) ./ sigmoidTau;
        tau = sigmoidTau;
    otherwise
        error('Unknown activeShorteningActivation "%s". Use "exponential", "sigmoid", or "accelerating".', activationMode);
end

if isfield(param,'activeRestLengthWavePulse_on') && param.activeRestLengthWavePulse_on && ...
        isfield(param,'tActiveRestLengthOff_edge') && numel(param.tActiveRestLengthOff_edge) == Ne
    tOff = param.tActiveRestLengthOff_edge(:);
    offTau = tau;
    if isfield(param,'activeRestLengthWaveOffTau') && isfinite(param.activeRestLengthWaveOffTau) && param.activeRestLengthWaveOffTau > 0
        offTau = param.activeRestLengthWaveOffTau;
    end
    Aoff = 1 ./ (1 + exp((t - tOff) ./ max(offTau, eps)));
    AoffDot = -(Aoff .* (1 - Aoff)) ./ max(offTau, eps);
    Adot = Adot .* Aoff + A .* AoffDot;
    A = A .* Aoff;
end

if isfield(param,'activeShorteningMask') && ~isempty(param.activeShorteningMask)
    mask = param.activeShorteningMask(:);
else
    mask = ones(Ne,1);
end

if isscalar(mask)
    rawFrac = param.activeShorteningMax * A .* mask;
    rawFracDot = param.activeShorteningMax * Adot .* mask;
elseif numel(mask) == Ne
    rawFrac = param.activeShorteningMax * A .* mask;
    rawFracDot = param.activeShorteningMax * Adot .* mask;
else
    error('param.activeShorteningMask must be scalar or have one value per edge.');
end

activeFrac = min(max(rawFrac, 0), 0.95);
activeFracDot = rawFracDot;
activeFracDot(rawFrac <= 0 | rawFrac >= 0.95) = 0;
end


function A = activeEdgeActivation(t, param)
Ne = numel(param.L0_edge0);

if isfield(param,'tActOn_edge') && numel(param.tActOn_edge) == Ne
    tActOn = param.tActOn_edge(:);
else
    tActOn = param.tActOn * ones(Ne,1);
end

Aon = 1 ./ (1 + exp(-(t - tActOn)/max(param.tActTau,eps)));

if isfield(param,'activeWavePulse_on') && param.activeWavePulse_on && ...
        isfield(param,'tActOff_edge') && numel(param.tActOff_edge) == Ne
    tActOff = param.tActOff_edge(:);
    offTau = param.tActTau;
    if isfield(param,'activeWaveOffTau') && isfinite(param.activeWaveOffTau) && param.activeWaveOffTau > 0
        offTau = param.activeWaveOffTau;
    elseif isfield(param,'activeRingOffTau') && isfinite(param.activeRingOffTau) && param.activeRingOffTau > 0
        offTau = param.activeRingOffTau;
    end
    Aoff = 1 ./ (1 + exp((t - tActOff)/max(offTau,eps)));
    A = Aon .* Aoff;
else
    A = Aon;
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


function [nodeCoords, cellIds] = loadExperimentalNodeCoords(xlsxFile, coordinateScale, centerMode)
if nargin < 2 || isempty(coordinateScale)
    coordinateScale = 1;
end
if nargin < 3 || isempty(centerMode)
    centerMode = 'median';
end
if ~isfile(xlsxFile)
    error('Experimental node data file not found: %s', xlsxFile);
end

raw = readmatrix(xlsxFile);
if size(raw, 2) < 4
    error('Experimental node data must contain Cell ID, x, y, z in the first four columns.');
end

cellIds = raw(:,1);
nodeCoords = raw(:,2:4) * coordinateScale;
validRows = all(isfinite(nodeCoords), 2);
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

fprintf('[GEOMETRY] loaded %d experimental nodes from %s\n', size(nodeCoords,1), xlsxFile);
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

        fprintf(['[GEOMETRY] oriented experimental nodes with principal axes ', ...
                 '(long axis -> x, posterior model end is -x)\n']);
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
