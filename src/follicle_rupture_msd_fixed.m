function simFile = follicle_rupture_msd_fixed(mechanicsProfile, simFile, coordinateCalibration)
% Follicle rupture shell model
% Default geometry uses the source image calibration. Use 'cell-scale' for
% estimated cell mechanics in the overdamped limit; its output is separate.
% An optional second argument selects the output MAT-file.
% Third argument: 'source' (default), or 'input-units' to reproduce the old
% coordinate assumption. Mechanics and coordinate calibration are independent.
% State: X = [q; v; sigma; zeta; D; A]
%   q, v        node displacement and velocity
%   sigma,zeta  viscoelastic edge states
%   D           edge damage in [0,1], where 1 means failed
%   A           node adhesion release in [0,1], where 1 means detached
% Cell-scale integration replaces q with absolute position and eliminates v;
% the saved displacement/velocity schema is unchanged.

clc; close all;

if nargin < 1 || isempty(mechanicsProfile)
    mechanicsProfile = 'legacy';
end
mechanicsProfile = normalizeUnitMode(mechanicsProfile, ...
    {'legacy', 'cell-scale'}, 'mechanics profile');
if nargin < 3 || isempty(coordinateCalibration)
    coordinateCalibration = 'source';
end
coordinateCalibration = normalizeUnitMode(coordinateCalibration, ...
    {'source', 'input-units'}, 'coordinate calibration');

%% USER PARAMETERS
% Edit this section to configure a run. Derived arrays and timing values are
% constructed later from these scalar controls.

% Reproducibility, units, and generated geometry
randomSeed = 42;
cellDiameter_um = 10;            % positive cell diameter in microns
coordinateInputMode = 'microns'; % 'cells', 'microns', or 'meters' for generated/other inputs
% The known tracked dataset uses Step14 calibration when coordinateCalibration='source'.
displayMode = 'microns';         % 'cells' or 'microns'
generatedNodeCount = 300;         % used only for generated geometry
generatedEllipsoidSemiAxes = [50 15 15]; % coordinateInputMode units

% Experimental geometry and connectivity
useExperimentalNodeData = true;
experimentalNodeDataFileName = fullfile('tracked_nodes', 'trackedNodes_nuc_48_2_OA7_1.mat');
experimentalNodeDataVariable = 'trackedNodes';
experimentalNodeDataFrame = 1; % use initial tracked coordinates only
experimentalCenterMode = 'median';   % 'median', 'mean', or 'none'
experimentalOrientMode = 'principal-axes'; % 'principal-axes' aligns long axis to model x-axis
experimentalLongAxisSign = -1;        % set 1 to flip experimental anterior/posterior ends
experimentalSurfaceConstraintMode = 'reference-ellipsoid'; % preserve each measured node's initial shell level
experimentalAnchorOn = true;           % softly hold the far body end during directional retraction
experimentalEdgeMode = 'knn';          % 'knn' connects all measured nodes; 'hull' uses convex hull only
experimentalNeighborCount = 6;         % local neighbors per measured node for 'knn' edge mode

% Passive mechanics and surface coupling
mVal = 1;                          % legacy nominal mass; see README inertia audit
% Cell-scale estimates, independent of the image calibration/display units.
cellVolume_um3 = 500;               % assumed volume per cell, not measured
cellDensity_kg_m3 = 1000;           % water-like wet density
cellForceScale = 1e-3;              % scales ALL force/stiffness/damping terms
kEquilibrium = 2;                  % [N/m], equilibrium spring
kTransient = 2;                    % [N/m], transient spring
cMaxwell = 100;                    % [N*s/m], Maxwell branch
cKelvin = 100;                     % [N*s/m], Kelvin branch (Burgers only)
param.connectionModel = 'zener';    % persistent spring plus relaxing Maxwell branch
param.postAxis = 1;
param.postSign = 1;                  % tracked-node rupture tip convention
param.k_surf = 5e-6;                % [N] surface penalty strength
param.c_surf = 1;                   % [N*s/m] normal surface damping
param.c_drag = 0.1;                 % [N*s/m] background tissue drag
param.projectEdgeForcesToSurface = true;

% Rupture-region selection
param.damageMode = 'principal-0-10-tip-patch';
param.tipPct = 98.5;                % used by the alternate 'tip'/'seam' modes
param.seamWidthFrac = 0.04;         % used by the alternate 'seam' mode
param.principalTipPercentile = 10;
param.tipPatchEdgeCount = 22;

% Distributed contraction and directional traction
param.anchorPercentile = 90;
param.contractionMinAxialAlignment = 0.35;
param.prestressStartTime = 3;       % [s]
param.activationRiseTime = 900;     % [s], prolonged global slow-pull onset
param.contractionRateFrac = 4e-4;   % [L0/s], plateau rate; smooth ramps extend the duration
param.minRestLengthFrac = 0.50;     % maximum 50% natural-length shortening
param.activeTraction_on = false;    % post-rupture motion comes from stored energy
param.activeTractionForce = 1e-7;   % [N/node], away from the rupture tip

% Rupture-gated elastic adhesion clutch
param.adhesionClutch_on = true;
param.adhesionStiffness = 1;               % [N/m], reference-position tether
param.adhesionLoadThresholdTip = 5e-7;     % [N], tip releases most readily
param.adhesionLoadThresholdBody = 3e-6;    % [N] at the P90 anchor boundary
param.adhesionTransitionWidthFrac = 0.25;
param.adhesionReleaseTau = 5;              % [s], irreversible detachment time
param.adhesionRuptureGateDamage = 0.5;     % release starts only after rupture
param.adhesionDragStuck = 0.05;               % [N*s/m]
param.adhesionDragSliding = 0.01;           % [N*s/m]
param.adhesionLoadRegularization = 1e-12;  % [N]

% Far-body anchor and rupture law
param.k_anchor = 10;                % [N/m]
param.c_anchor = 1;                 % [N*s/m]
param.extensionCritFrac = 0.05;     % [L0] weak-edge extension to rupture
param.rupturePatchCascade_on = true; % first weak-edge break opens the local patch
param.ruptureDamageTau = 0.5;       % [s], fast relative to tissue relaxation
param.g_floor_stiff = 0;            % residual stiffness after failure
param.g_floor_damp = 1e-6;          % residual damping after failure

% Time integration and saved output
tEnd = 1800;                       % [s]
numSave = 300;                     % frames after the dense rupture window
ruptureSaveEnd = 900;              % [s], includes slow preload and recoil onset
ruptureSaveInterval = 1;           % [s]
odeSolver = 'ode45';               % vectorized explicit solve; rupture remains resolved
relativeTolerance = 1e-5;
absoluteTolerance = 1e-8;
maximumSolverStep = 5;
initialSolverStep = 1e-3;
progressEveryPhysSec = 5;          % [s]

%% Derived setup
% Inputs are converted to meters once at geometry ingestion. Mechanics and
% saved physical data remain SI; strain and damage are dimensionless.
rng(randomSeed);
units = buildUnitSystem(cellDiameter_um, coordinateInputMode, displayMode);

repoRoot = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(repoRoot, 'results');
inputNodeDataDir = fullfile(repoRoot, 'input_node_data');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if ~exist(inputNodeDataDir, 'dir')
    mkdir(inputNodeDataDir);
end
if nargin < 2 || isempty(simFile)
    if strcmp(mechanicsProfile, 'legacy')
        simFile = fullfile(outputDir, 'follicleShellSim_test2.mat');
    else
        simFile = fullfile(outputDir, 'follicleShellSim_cell_scale.mat');
    end
    if strcmp(coordinateCalibration, 'input-units')
        [folder, name, extension] = fileparts(simFile);
        simFile = fullfile(folder, [name '_input_units' extension]);
    end
end
if ~(ischar(simFile) || (isstring(simFile) && isscalar(simFile)))
    error('follicle_rupture:invalidOutputFile', 'simFile must be a path string.');
end
saveDir = fileparts(simFile);
if ~isempty(saveDir) && ~isfolder(saveDir)
    mkdir(saveDir);
end
experimentalNodeDataFile = fullfile(inputNodeDataDir, experimentalNodeDataFileName);

%% Nodes
if useExperimentalNodeData
    units = experimental_coordinate_units(experimentalNodeDataFile, units, coordinateCalibration);
    fprintf('[UNITS] input = %s | meters per raw xyz unit = %s | source = %s\n', ...
        units.coordinateInputMode, mat2str(units.coordinateToMeters, 12), units.calibrationSource);
    [nodeCoords, experimentalCellIds] = loadExperimentalNodeCoords( ...
        experimentalNodeDataFile, units.coordinateToMeters, experimentalCenterMode, ...
        experimentalNodeDataVariable, experimentalNodeDataFrame);
    nodeCoords = orientExperimentalNodeCoords( ...
        nodeCoords, experimentalOrientMode, experimentalLongAxisSign);
    N = size(nodeCoords, 1);
    [a, b, c] = fitAxisAlignedEllipsoidSemiAxes(nodeCoords);
    L_scale = median([a b c]);  % representative experimental radius scale
else
    N = generatedNodeCount;
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
param.coordinateCalibration = coordinateCalibration;
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

if useExperimentalNodeData
    param.surfaceConstraintMode = experimentalSurfaceConstraintMode;
else
    param.surfaceConstraintMode = 'ellipsoid';
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
param.k1_edge = kEquilibrium * ones(Ne,1);
param.k2_edge = kTransient * ones(Ne,1);
param.c_edge  = cMaxwell * ones(Ne,1);
param.c_series_edge = cKelvin * ones(Ne,1);
param.mechanicsProfile = mechanicsProfile;
if strcmp(mechanicsProfile, 'cell-scale')
    param = cell_scale_mechanics(param, cellVolume_um3, ...
        cellDensity_kg_m3, cellForceScale);
    odeSolver = 'ode15s'; % integrate after eliminating overdamped velocities
end

% You may set param.connectionModel_edge to a per-edge cell/string array.
param.connectionModelId_edge = connectionModelIds(param, Ne);

%% Spatial localization for alternate tip/seam rupture modes
edgeMid = 0.5*(nodeCoords(E(:,1),:) + nodeCoords(E(:,2),:) );
ax = param.postAxis;
sg = param.postSign;
param.edgeMid0 = edgeMid;

sNodes = sg * nodeCoords(:,ax);
sMid   = sg * edgeMid(:,ax);

sCutTip  = prctile(sNodes, param.tipPct);

param.tipMask  = double(sMid >= sCutTip);     % tiny posterior cap edges

%% Damage region selection
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
principalTipCutoff = prctile(principalNodeScore, ...
    param.principalTipPercentile);

param.principalLongAxis0 = principalLongAxis;
param.principalNodeScore0 = principalNodeScore;
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
        % 0-10% principal-axis band shown in post-processing. With the cascade
        % enabled, the first threshold crossing opens this entire patch.
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

fprintf('[MASK] tip edges = %d/%d | damage edges = %d\n', ...
    nnz(param.tipMask), Ne, nnz(param.damageMask));
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

%% Distributed bodyward-oriented prestress
% Form a chain of actively shortening, predominantly axial edges from the
% rupture tip toward a softly anchored far-body cap. With the far end held,
% shortening along this corridor draws each intervening percentile ring
% away from the tip. Rupture edges themselves never contract.
param.anchorCutoff = prctile(principalNodeScore, param.anchorPercentile);
param.anchorNodeMask = double(principalNodeScore >= param.anchorCutoff);

edgeLongAxisSpan = abs(principalNodeScore(E(:,2)) - ...
    principalNodeScore(E(:,1)));
edgeAxialAlignment = edgeLongAxisSpan ./ max(param.L0_edge0, 1e-18);
reachesTowardAnchor = principalEdgeScore <= param.anchorCutoff & ...
    ~(param.anchorNodeMask(E(:,1)) > 0 & param.anchorNodeMask(E(:,2)) > 0);
param.retractionMask = double(reachesTowardAnchor & ...
    edgeAxialAlignment >= param.contractionMinAxialAlignment & ...
    param.damageMask == 0);
if ~any(param.retractionMask)
    error('No intact axial edges are available for distributed bodyward prestress.');
end

% Start shortening everywhere in the axial corridor at the same time. A
% smooth, load-sensitive adhesion clutch below determines when each node can
% slide, so the apparent activation sequence is mechanical rather than a
% prescribed tip-to-body clock. The far-body anchor supplies the directional
% reference; optional external traction is disabled by default.
param.contractionRampTime = param.activationRiseTime;
param.contractionRateFrac_edge = zeros(Ne,1);
param.contractionRateFrac_edge(param.retractionMask > 0) = ...
    param.contractionRateFrac;
param.minRestLengthFrac_edge = ones(Ne,1);
param.minRestLengthFrac_edge(param.retractionMask > 0) = ...
    param.minRestLengthFrac;
param.contractionStartTime_edge = inf(Ne,1);
param.contractionStartTime_edge(param.retractionMask > 0) = ...
    param.prestressStartTime;
param.contractionDuration_edge = zeros(Ne,1);
param.contractionDuration_edge(param.retractionMask > 0) = (1 - ...
    param.minRestLengthFrac_edge(param.retractionMask > 0)) ./ ...
    param.contractionRateFrac_edge(param.retractionMask > 0) + ...
    param.contractionRampTime;
param.prestressDuration = max(param.contractionStartTime_edge(param.retractionMask > 0) + ...
    param.contractionDuration_edge(param.retractionMask > 0)) - ...
    param.prestressStartTime;
nodeScoreMin = min(principalNodeScore);
nodeScoreSpan = max(param.anchorCutoff - nodeScoreMin, 1e-18);
nodeBodywardFraction = min(max((principalNodeScore - nodeScoreMin) ./ ...
    nodeScoreSpan, 0), 1);
param.activeTractionStartTime_node = ...
    repmat(param.prestressStartTime, N, 1);
param.activeTractionStartTime_node(param.anchorNodeMask > 0) = inf;
param.activeTractionDuration = param.prestressDuration; % [s], through contraction
param.activeTractionRampTime = param.activationRiseTime;
% Adhesion release thresholds rise toward the far body. Every node begins with
% an elastic reference-position tether, and no tether may release until the
% weak patch ruptures. Subsequent detachment is driven only by local load.
param.adhesionLoadThreshold_node = param.adhesionLoadThresholdTip + ...
    (param.adhesionLoadThresholdBody - param.adhesionLoadThresholdTip) .* ...
    nodeBodywardFraction;
param.adhesionLoadThreshold_node(param.anchorNodeMask > 0) = ...
    param.adhesionLoadThresholdBody;

% The patch is eligible as soon as global contraction begins; its extension
% threshold, rather than a later clock, decides the actual rupture time.
param.ruptureEnableTime = param.prestressStartTime;

% The patch edges are weak and fail when stretched by the retracting tissue;
% the axial corridor remains intact so it continues drawing the rings back.
param.extensionCritFrac_edge = inf(Ne,1);
param.extensionCritFrac_edge(param.damageMask > 0) = ...
    param.extensionCritFrac;

if any(param.contractionRateFrac_edge < 0) || ...
        any(param.minRestLengthFrac_edge <= 0 | param.minRestLengthFrac_edge > 1)
    error('Contraction rates must be nonnegative and rest-length floors must lie in (0,1].');
end
if param.adhesionStiffness < 0 || param.adhesionReleaseTau <= 0 || ...
        param.adhesionLoadThresholdTip <= 0 || ...
        param.adhesionLoadThresholdBody <= 0
    error('Adhesion stiffness must be nonnegative; release time and load thresholds must be positive.');
end

% Soft far-body anchoring provides the directional reference for the
% distributed contraction. It is a positional tether, not a rigid pin.
if useExperimentalNodeData
    param.anchor_on = experimentalAnchorOn;
else
    param.anchor_on = true;
end
fprintf('[ANCHOR] enabled = %d | far-body nodes >= P%d = %d/%d\n', ...
    param.anchor_on, param.anchorPercentile, nnz(param.anchorNodeMask), N);

retractingEdges = find(param.retractionMask > 0);
fprintf('[PRESTRESS] axial corridor edges = %d | alignment >= %.2f | global start = %.3g s | rate = %.3g L0/s | floor = %.3g L0\n', ...
    numel(retractingEdges), param.contractionMinAxialAlignment, ...
    param.prestressStartTime, ...
    param.contractionRateFrac_edge(retractingEdges(1)), ...
    param.minRestLengthFrac_edge(retractingEdges(1)));
fprintf('[PRELOAD] global smoothstep rise = %.3g s; no spatial activation schedule\n', ...
    param.activationRiseTime);
fprintf('[TRACTION] enabled = %d | force = %.3g N/node | duration = %.3g s | ramp = %.3g s\n', ...
    param.activeTraction_on, param.activeTractionForce, param.activeTractionDuration, ...
    param.activeTractionRampTime);
fprintf(['[ADHESION CLUTCH] enabled = %d | stiffness = %.3g N/m | ', ...
    'threshold tip/body = %.3g/%.3g N | transition = %.3g threshold | ', ...
    'release tau = %.3g s | ', ...
    'drag stuck/sliding = %.3g/%.3g N*s/m\n'], ...
    param.adhesionClutch_on, param.adhesionStiffness, param.adhesionLoadThresholdTip, ...
    param.adhesionLoadThresholdBody, param.adhesionTransitionWidthFrac, ...
    param.adhesionReleaseTau, param.adhesionDragStuck, param.adhesionDragSliding);
fprintf(['[PRESTRESS/RUPTURE] corridor shortening duration = %.3g s | rupture eligibility = %.3g s | ', ...
    'weak patch edges = %s | break = %.3g L0 extension\n'], ...
    param.prestressDuration, param.ruptureEnableTime, ...
    mat2str(selectedDamageEdge.'), ...
    param.extensionCritFrac_edge(selectedDamageEdge(1)));

%% Nondimensionalization
k1_scale = mean(param.k1_edge);
k2_scale = mean(param.k2_edge);
c0_scale = mean(param.c_edge);
k0_scale = k1_scale + k2_scale;
tau_scale = c0_scale / max(k2_scale,1e-18);
epsilon   = param.m * k0_scale / max(c0_scale^2,1e-18);

param.scale.L       = L_scale;
param.scale.k2      = k2_scale;
param.scale.tau     = tau_scale;
param.scale.epsilon = epsilon;
% Physical inertia is negligible on the saved time scale. Record its size
% explicitly; cell-scale uses its overdamped limit, 0 = F_net, rather than
% increasing mass to make an explicit integration possible.
param.scale.inertia = param.m / (k2_scale * tau_scale^2);
if strcmp(mechanicsProfile, 'legacy')
    param.effectiveInertialMass = param.m / epsilon;
else
    param.effectiveInertialMass = 0; % overdamped approximation
    param.inertialRelaxationTime = param.m / param.c_drag;
end
fprintf('[MECHANICS] %s | node mass = %.3g kg | effective inertia = %.3g kg | k2 = %.3g N/m | relaxation = %.3g s\n', ...
    mechanicsProfile, param.m, param.effectiveInertialMass, k2_scale, tau_scale);
param.numSave = numSave;
param.ruptureSaveEnd = ruptureSaveEnd;
param.ruptureSaveInterval = ruptureSaveInterval;
param.odeSolver = odeSolver;
param.progressEveryPhysSec = progressEveryPhysSec;

%% ================= State sizes =================
nPos = 3*N;
nVel = 3*N;
nAdhesion = N;

param.nPos = nPos;
param.nVel = nVel;
param.nAdhesion = nAdhesion;

% X = [q; v; sigma; zeta; D; A]
nX = nPos + nVel + Ne + Ne + Ne + nAdhesion;

% Damage and adhesion-state block indices for saved time courses
iD1 = nPos + nVel + 2*Ne + 1;
iD2 = iD1 + Ne - 1;
iA1 = iD2 + 1;
iA2 = iA1 + nAdhesion - 1;

%% Integrate
tFinal_dim = tEnd / tau_scale;
assert(isscalar(tFinal_dim) && isfinite(tFinal_dim));

X0_dim = zeros(nX,1);
if strcmp(mechanicsProfile, 'cell-scale')
    % Absolute positions keep ode15s finite-difference perturbations above
    % coordinate roundoff during the very small initial displacements.
    X0_dim(1:nPos) = reshape(nodeCoords.', [], 1) / L_scale;
end

opts = odeset( ...
    'RelTol', relativeTolerance, ...
    'AbsTol', absoluteTolerance, ...
    'MaxStep', maximumSolverStep / tau_scale, ...
    'InitialStep', initialSolverStep / tau_scale, ...
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

odeRHS = @(t,X) follicle_shell_rhs(t, X, param);
if strcmp(mechanicsProfile, 'cell-scale')
    [XDimQ_all, param] = integrate_cell_scale(param, tDimQ_all, X0_dim, opts);
else
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
end

tSol  = tDimQ_all * tau_scale;
if strcmp(mechanicsProfile, 'cell-scale') && ...
        (abs(tSol(end) - tEnd) > 1e-8 * max(tEnd,1) || ...
         any(~isfinite(XDimQ_all), 'all'))
    error('follicle_rupture:incompleteIntegration', ...
        'Cell-scale integration did not complete; the output file was not replaced.');
end
Y_dim = XDimQ_all(:, 1:nPos);
if strcmp(mechanicsProfile, 'cell-scale')
    Y_dim = Y_dim - reshape(nodeCoords.', 1, []) / L_scale;
end
Y     = Y_dim * L_scale;
V_dim = XDimQ_all(:, nPos+1:nPos+nVel);
V     = V_dim * (L_scale/tau_scale);
L0_edge0 = param.L0_edge0;

% extract edge state time courses
D_time = min(max(XDimQ_all(:, iD1:iD2), 0), 1);
A_time = min(max(XDimQ_all(:, iA1:iA2), 0), 1);

% Save the imposed natural-length history alongside displacement/damage so
% rupture can be checked directly in post-processing.
Lnatural_time = zeros(numel(tSol), Ne);
% Evaluate forces through the same RHS used by the solver, including the
% viscoelastic history states, damage, and surface projection. Units are N.
edgeForce_time = zeros(numel(tSol), Ne);
nodeInternalForce_time = zeros(numel(tSol), nPos);
nodeNetForce_time = zeros(numel(tSol), nPos);
for it = 1:numel(tSol)
    Lnatural_time(it,:) = contracting_natural_lengths(tSol(it), param).';
    [~, forces] = follicle_shell_rhs(tDimQ_all(it), XDimQ_all(it,:).', param);
    edgeForce_time(it,:) = forces.edgeForce.';
    nodeInternalForce_time(it,:) = reshape(forces.nodeInternalForce.', 1, []);
    nodeNetForce_time(it,:) = reshape(forces.nodeNetForce.', 1, []);
end

save(simFile, ...
     'nodeCoords','edgeList','L0_edge0','Lnatural_time','tSol','Y','V', ...
     'D_time','A_time','edgeForce_time','nodeInternalForce_time', ...
     'nodeNetForce_time','param','-v7.3');

end % main


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

modelKey = normalizeConnectionModel(param.connectionModel);
modelIds = connectionModelId(modelKey) * ones(Ne,1);
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


%% Generated ellipsoid nodes
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
        nodeCoords = double(tracked(:,:,frameIndex)) .* coordinateScale;
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
        nodeCoords = raw(:,2:4) .* coordinateScale;
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
