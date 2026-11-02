%% ===================== Nondimensional wrapper RHS =====================
function [dXdim, forces] = follicle_shell_rhs(t_dim, X_dim, param)
%FOLLICLE_SHELL_RHS Dimensionless shell equations and physical force histories.
% 'cell-scale' uses absolute position states and velocity rows F/(k2*L).
% The overdamped helper eliminates velocities by setting those rows to zero.
% Dividing by mu=m/(k2*tau^2) recovers Newton's law for verification.
% Legacy files retain their original acceleration scaling.
nPos = param.nPos;
nVel = param.nVel;
Ne   = param.Ne;
N    = param.N;
scale = param.scale;

% indices in ND vector
iQ2  = nPos;
iV2  = iQ2 + nVel;
iS2  = iV2 + Ne;
iZ2  = iS2 + Ne;
iD2  = iZ2 + Ne;
iA2  = iD2 + N;

q_dim   = X_dim(1:iQ2);
if isfield(param, 'mechanicsProfile') && strcmp(param.mechanicsProfile, 'cell-scale')
    q_dim = q_dim - reshape(param.nodeCoords.', [], 1) / scale.L;
end
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

x_phys = [q_phys; v_phys; sig_phys; zet_phys; D_phys; A_phys];
if nargout > 1
    [dX_phys, forces] = ShellRHS_physical(t_phys, x_phys, param);
else
    dX_phys = ShellRHS_physical(t_phys, x_phys, param);
end

% back to ND
dqdt_dim   = (dX_phys(1:nPos)/scale.L) * scale.tau;
if isfield(param, 'mechanicsProfile') && strcmp(param.mechanicsProfile, 'cell-scale')
    dvdt_dim = scale.inertia * ((dX_phys(nPos+1:nPos+nVel) / (scale.L/scale.tau)) * scale.tau);
else
    % Compatibility only: this historically cancels m, giving m_eff=m/epsilon.
    dvdt_dim = scale.epsilon * ((dX_phys(nPos+1:nPos+nVel) / (scale.L/scale.tau)) * scale.tau);
end

dsig_dim   = (dX_phys(nPos+nVel+1 : nPos+nVel+Ne) / (scale.k2*scale.L)) * scale.tau;
dzet_dim   = (dX_phys(nPos+nVel+Ne+1 : nPos+nVel+2*Ne) / scale.L) * scale.tau;

off = nPos + nVel + 2*Ne;
dD_dim     = dX_phys(off+(1:Ne)) * scale.tau;
dA_dim     = dX_phys(off+Ne+(1:N)) * scale.tau;

dXdim = [dqdt_dim; dvdt_dim; dsig_dim; dzet_dim; dD_dim; dA_dim];
end


%% ========================= PHYSICAL RHS =========================
function [dxdt, forces] = ShellRHS_physical(t, x, param)
% X = [ q(3N); v(3N); sigma(Ne); zeta(Ne); D(Ne); A(N) ]

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
iA1 = iD2+1;             iA2 = iA1 + N - 1;

q     = x(iQ1:iQ2);
v     = x(iV1:iV2);
sigma = x(iS1:iS2);
zeta  = x(iZ1:iZ2);
D     = x(iD1:iD2);
A     = x(iA1:iA2);

% Clamp bounded damage/release states safely.
D = min(max(D,0),1);
A = min(max(A,0),1);

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
    [level, nSurf, surfaceNormalSpeed] = shell_surface_geometry(t,pos,param);
    fSurf = level - 1;
    if strcmp(surfaceConstraintMode, 'reference-ellipsoid')
        fSurf = fSurf - (param.surfaceLevel0 - 1);
    end
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

[Lnatural, LnaturalDot] = contracting_natural_lengths(t, param);
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

if nargout > 1
    forces.edgeForce = sig_total; % signed axial link force; tension is positive
    forces.nodeInternalForce = Fnet; % projected link forces before external loads
end

dD = zeros(Ne,1);

% ---- Weakest-edge rupture from critical extension ----
% This is the direct 3D counterpart of ext >= delt_crit in the chain model.
criticalExtension = param.extensionCritFrac_edge(:) .* L0_ref;
ruptureEnabled = ~isfield(param, 'ruptureEnableTime') || ...
    t >= param.ruptureEnableTime;
thresholdBreaking = ruptureEnabled & (param.damageMask(:) > 0) & ...
    (stretch >= criticalExtension);
breaking = thresholdBreaking;
if isfield(param, 'ruptureTime')
    % Cell-scale integration locates the first threshold crossing as an
    % event. Finite-difference Jacobian probes must not trigger the latch.
    breaking = (t >= param.ruptureTime) & (param.damageMask(:) > 0);
elseif isfield(param, 'rupturePatchCascade_on') && param.rupturePatchCascade_on
    patchAlreadyBreaking = any(D(param.damageMask(:) > 0) > 0);
    if patchAlreadyBreaking || any(thresholdBreaking)
        breaking = ruptureEnabled & (param.damageMask(:) > 0);
    end
end
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
    vN      = sum(vMat .* nSurf, 2) - surfaceNormalSpeed;
    F_dampN = -param.c_surf * (vN .* nSurf);
end
F_drag  = -param.c_drag * vMat;

% Directed active traction along the fitted surface. This represents the
% external reaction supplied by basal adhesion or surrounding tissue, which
% pairwise edge contraction cannot provide by itself.
F_active = zeros(N,3);
if isfield(param, 'activeTraction_on') && param.activeTraction_on
    tractionStart = param.activeTractionStartTime_node(:);
    rampTime = max(param.activeTractionRampTime, 1e-8);
    rampUp = smoothstep01((t - tractionStart) ./ rampTime);
    rampDown = smoothstep01((tractionStart + ...
        param.activeTractionDuration - t) ./ rampTime);
    tractionActivity = rampUp .* rampDown;

    bodywardDirection = repmat(param.principalLongAxis0(:).', N, 1);
    if useEllipsoidNormals
        bodywardDirection = bodywardDirection - ...
            sum(bodywardDirection .* nSurf, 2) .* nSurf;
    end
    directionNorm = vecnorm(bodywardDirection, 2, 2);
    validDirection = directionNorm > 1e-12;
    bodywardDirection(validDirection,:) = ...
        bodywardDirection(validDirection,:) ./ directionNorm(validDirection);
    bodywardDirection(~validDirection,:) = 0;
    F_active = param.activeTractionForce .* tractionActivity .* bodywardDirection;
end

% Elastic adhesion stores the global preload instead of merely slowing it.
% The clutch cannot detach before a selected weak edge actually ruptures.
% After rupture, its local stored load drives irreversible, threshold-dependent
% release; the tip-to-body threshold gradient can then produce sequential ring
% motion without assigning any node an activation time.
F_adhesion = zeros(N,3);
dA = zeros(N,1);
if isfield(param, 'adhesionClutch_on') && param.adhesionClutch_on
    tangentDisplacement = qMat;
    tangentVelocity = vMat;
    if useEllipsoidNormals
        tangentDisplacement = tangentDisplacement - ...
            sum(tangentDisplacement .* nSurf, 2) .* nSurf;
        tangentVelocity = tangentVelocity - ...
            sum(tangentVelocity .* nSurf, 2) .* nSurf;
    end

    elasticClutchForce = -param.adhesionStiffness .* tangentDisplacement;
    trialClutchForce = elasticClutchForce - ...
        param.adhesionDragStuck .* tangentVelocity;
    loadEpsilon = max(param.adhesionLoadRegularization, 1e-18);
    clutchLoad = sqrt(sum(trialClutchForce.^2, 2) + loadEpsilon^2) - ...
        loadEpsilon;
    loadThreshold = param.adhesionLoadThreshold_node(:);
    transitionWidth = max(param.adhesionTransitionWidthFrac .* ...
        loadThreshold, loadEpsilon);
    overloadActivity = 0.5 .* (1 + tanh((clutchLoad - ...
        loadThreshold) ./ transitionWidth));

    selectedDamage = D(param.damageMask(:) > 0);
    if isempty(selectedDamage)
        ruptureDamage = 0;
    else
        ruptureDamage = max(selectedDamage);
    end
    gateStart = min(max(param.adhesionRuptureGateDamage, 0), 1 - 1e-8);
    ruptureGate = smoothstep01((ruptureDamage - gateStart) ./ ...
        max(1 - gateStart, 1e-8));
    dA = ruptureGate .* overloadActivity .* (1 - A) ./ ...
        max(param.adhesionReleaseTau, 1e-8);
    dA(param.anchorNodeMask > 0) = 0;

    clutchEngagement = 1 - A;
    adhesionDrag = clutchEngagement .* param.adhesionDragStuck + ...
        A .* param.adhesionDragSliding;
    clutchEngagement(param.anchorNodeMask > 0) = 0;
    adhesionDrag(param.anchorNodeMask > 0) = 0;
    F_adhesion = clutchEngagement .* elasticClutchForce - ...
        adhesionDrag .* tangentVelocity;
end

Fnet = Fnet + F_surf + F_dampN + F_drag + F_active + F_adhesion;

% Soft anterior positional boundary condition. This weak tether biases
% directional retraction without rigidly pinning nodes or adding contact.
if isfield(param,'anchor_on') && param.anchor_on
    F_anchor = -param.k_anchor * qMat - param.c_anchor * vMat;
    Fnet = Fnet + param.anchorNodeMask .* F_anchor;
end

if nargout > 1
    forces.nodeNetForce = Fnet; % includes drag, adhesion, surface and anchor loads
end

% dq/dt
dqdt = reshape(vMat.', [], 1);

% dv/dt
dvdt = reshape((Fnet/param.m).',[],1);

% assemble
dxdt = [dqdt; dvdt; dsigma; dzeta; dD; dA];
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


function value = smoothstep01(coordinate)
coordinate = min(max(coordinate, 0), 1);
value = coordinate.^2 .* (3 - 2 * coordinate);
end


