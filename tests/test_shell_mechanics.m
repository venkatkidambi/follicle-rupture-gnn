function tests = test_shell_mechanics
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testNewtonUnitsAndMass(testCase)
[p, x] = twoCellFixture;
x(4) = x(4) + 0.1; % extend a 10 um link by 1 um
[dx, f] = follicle_shell_rhs(0, x, p);
physicalAcceleration = dx(7:12) / p.scale.inertia * p.scale.L / p.scale.tau^2;
verifyEqual(testCase, physicalAcceleration, ...
    reshape(f.nodeNetForce.', [], 1) / p.m, 'RelTol', 1e-12);
verifyEqual(testCase, f.edgeForce, 2e-9, 'RelTol', 1e-12);
verifyEqual(testCase, sum(f.nodeInternalForce,1), zeros(1,3), 'AbsTol', 1e-22);
p.m = 2*p.m;
p.scale.inertia = 2*p.scale.inertia;
dxHeavy = follicle_shell_rhs(0, x, p);
heavyAcceleration = dxHeavy(7:12) / p.scale.inertia * p.scale.L / p.scale.tau^2;
verifyEqual(testCase, heavyAcceleration, physicalAcceleration/2, 'RelTol', 1e-12);
end

function testFailedLinkCannotRetainHistoryForce(testCase)
[p, x] = twoCellFixture;
x(4) = x(4) + 0.1;
x(13) = 0.2; % nonzero Maxwell history
x(15) = 1;   % fully failed link
[~, f] = follicle_shell_rhs(0, x, p);
verifyEqual(testCase, f.edgeForce, 0);
verifyEqual(testCase, f.nodeInternalForce, zeros(2,3));
end

function testZenerRelaxation(testCase)
[p, x] = twoCellFixture;
x(13) = 0.2;
dx = follicle_shell_rhs(0, x, p);
% At fixed length, sigma_dot = -sigma/(c/k2), and tau=c/k2.
verifyEqual(testCase, dx(13), -0.2, 'AbsTol', 1e-14);
end

function testCommonForceScaling(testCase)
[p, x] = twoCellFixture;
x(4) = x(4) + 0.1;
x(7) = 0.01;
x(13) = 0.2;
[~, f] = follicle_shell_rhs(0, x, p);
scaled = cell_scale_mechanics(p, 500, 1000, 0.001);
scaled.scale.k2 = scaled.k2_edge;
scaled.scale.inertia = scaled.m / (scaled.scale.k2 * scaled.scale.tau^2);
[~, g] = follicle_shell_rhs(0, x, scaled);
verifyEqual(testCase, scaled.m, 5e-13, 'RelTol', 1e-14);
verifyEqual(testCase, scaled.c_edge/scaled.k2_edge, p.c_edge/p.k2_edge, 'RelTol',1e-14);
verifyEqual(testCase, g.nodeNetForce, 0.001*f.nodeNetForce, 'AbsTol', 1e-23);
verifyEqual(testCase, scaled.adhesionLoadThresholdTip, ...
    0.001*p.adhesionLoadThresholdTip);
verifyEqual(testCase, scaled.adhesionLoadRegularization, ...
    0.001*p.adhesionLoadRegularization);
end

function testContractionLengthAndDerivative(testCase)
[p, ~] = twoCellFixture;
p.contractionRateFrac_edge = 4e-4;
p.minRestLengthFrac_edge = 0.5;
p.contractionStartTime_edge = 3;
p.contractionRampTime = 900;
[initial, initialRate] = contracting_natural_lengths(3, p);
[final, finalRate] = contracting_natural_lengths(2153, p);
verifyEqual(testCase, initial, p.L0_edge0);
verifyEqual(testCase, final, 0.5*p.L0_edge0);
verifyEqual(testCase, [initialRate finalRate], [0 0]);
for t = [100 900 1100 1500 2100]
    [~, rate] = contracting_natural_lengths(t, p);
    h = 0.01;
    numericalRate = (contracting_natural_lengths(t+h,p) - ...
        contracting_natural_lengths(t-h,p)) / (2*h);
    verifyEqual(testCase, rate, numericalRate, 'AbsTol', 1e-17);
end
end

function testOverdampedAnalyticalRetraction(testCase)
[p, x] = twoCellFixture;
p.connectionModelId_edge = 1; % simple spring: exact relative motion is exponential
extension = 1e-6;
x(1) = x(1) - extension/(2*p.scale.L);
x(4) = x(4) + extension/(2*p.scale.L);
% Consistent algebraic initial velocities: drag balances link tension.
speed = p.k1_edge*extension/p.c_drag;
x(7) = speed*p.scale.tau/p.scale.L;
x(10) = -x(7);
massDiagonal = ones(size(x));
massDiagonal(7:12) = 0;
opts = odeset('RelTol',1e-8, 'AbsTol',1e-11, ...
    'Mass',spdiags(massDiagonal,0,numel(x),numel(x)), ...
    'MStateDependence','none', 'MassSingular','yes', ...
    'JPattern',shell_jacobian_pattern(p));
time = linspace(0,0.1,21).';
[~, states] = ode15s(@(t,y) follicle_shell_rhs(t,y,p), time/p.scale.tau, x, opts);
measured = (states(:,4) - states(:,1))*p.scale.L - p.L0_edge0;
expected = extension*exp(-2*p.k1_edge/p.c_drag*time);
verifyEqual(testCase, measured, expected, 'AbsTol',1e-12);
p.rupturePatchCascade_on = false;
reducedStates = integrate_cell_scale(p, time/p.scale.tau, x, opts);
reducedExtension = (reducedStates(:,4)-reducedStates(:,1))*p.scale.L-p.L0_edge0;
verifyEqual(testCase, reducedExtension, expected, 'AbsTol',1e-12);
for idx = [1 10 21]
    [~, forces] = follicle_shell_rhs(time(idx)/p.scale.tau,states(idx,:).',p);
    verifyLessThan(testCase, max(abs(forces.nodeNetForce),[],'all'), 1e-14);
end
end

function testOverdampedForceBalanceAllMaterials(testCase)
[p, x] = twoCellFixture;
p.nodeCoords = [0.8 0.3 0.2; 1.5 0.2 0.4]*1e-5;
p.L0_edge0 = norm(diff(p.nodeCoords));
p.ell_a = 2e-5; p.ell_b = 1e-5; p.ell_c = 1e-5;
p.surfaceConstraintMode = 'reference-ellipsoid';
p.surfaceLevel0 = sum((p.nodeCoords./[p.ell_a p.ell_b p.ell_c]).^2,2);
p.projectEdgeForcesToSurface = true;
p.anchor_on = true; p.anchorNodeMask = [0;1];
p.adhesionClutch_on = true;
p.adhesionTransitionWidthFrac = 0.25;
p.adhesionLoadThreshold_node = [p.adhesionLoadThresholdTip;p.adhesionLoadThresholdBody];
p.adhesionReleaseTau = 5; p.adhesionRuptureGateDamage = 0.5;
x(1:6) = reshape(p.nodeCoords.',[],1)/p.scale.L;
x(1) = x(1)-0.03;
x(13) = 0.2; x(14) = 0.01; x(15) = 0.1; x(16:17) = [0.2;0.4];
reduced = x([1:6 13:17]);
for model = 1:6
    p.connectionModelId_edge = model;
    [dx, fullState] = overdamped_shell_rhs(0,reduced,p);
    [~, forces] = follicle_shell_rhs(0,fullState,p);
    verifyTrue(testCase, all(isfinite(dx)));
    verifyLessThan(testCase, max(abs(forces.nodeNetForce),[],'all'), 1e-20);
end
% Repeat with a translating/deforming surface: algebraic drag elimination
% must include the prescribed normal surface speed for every material.
p.movingSurface.time=[0;10];
p.movingSurface.center=[0 0 0;1 2 3]*1e-6;
p.movingSurface.quadric(:,:,1)=diag(1./[p.ell_a p.ell_b p.ell_c].^2);
p.movingSurface.quadric(:,:,2)=p.movingSurface.quadric(:,:,1)/1.1^2;
for model = 1:6
    p.connectionModelId_edge = model;
    [dx,fullState]=overdamped_shell_rhs(5/p.scale.tau,reduced,p);
    [~,forces]=follicle_shell_rhs(5/p.scale.tau,fullState,p);
    verifyTrue(testCase,all(isfinite(dx)));
    verifyLessThan(testCase,max(abs(forces.nodeNetForce),[],'all'),1e-20);
end
% Frozen rotated/translated geometry must give the same dynamics as an
% independently constant surface history, for every supported material.
p.movingSurface.center(1,:)=[2 -1 3]*1e-6;
angle=.4; R=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
p.movingSurface.quadric(:,:,1)=R*p.movingSurface.quadric(:,:,1)*R.';
fixed=set_surface_mode(p,'fixed');
constant=p;
constant.movingSurface.center(2,:)=constant.movingSurface.center(1,:);
constant.movingSurface.quadric(:,:,2)=constant.movingSurface.quadric(:,:,1);
for model = 1:6
    fixed.connectionModelId_edge=model; constant.connectionModelId_edge=model;
    [dx,fullState]=overdamped_shell_rhs(5/p.scale.tau,reduced,fixed);
    expected=overdamped_shell_rhs(5/p.scale.tau,reduced,constant);
    [~,forces]=follicle_shell_rhs(5/p.scale.tau,fullState,fixed);
    verifyEqual(testCase,dx,expected,'AbsTol',1e-12);
    verifyLessThan(testCase,max(abs(forces.nodeNetForce),[],'all'),1e-20);
end
end

function testJacobianIncludesNonlocalRuptureGate(testCase)
[p, ~] = twoCellFixture;
p.damageMask = 1;
pattern = shell_jacobian_pattern(p);
verifyTrue(testCase, issparse(pattern));
verifyEqual(testCase, size(pattern), [17 17]);
verifyEqual(testCase, full(pattern(16:17,15)), ones(2,1));
verifyEqual(testCase, full(pattern(15,[1:6 15])), ones(1,7));
end

function testRuptureEventAndIrreversibleCascade(testCase)
[p, x] = twoCellFixture;
p.c_drag = 1; % almost stationary endpoints while the test link shortens
p.damageMask = 1;
p.extensionCritFrac_edge = 0.05;
p.contractionRateFrac_edge = 0.01;
p.minRestLengthFrac_edge = 0.5;
p.contractionStartTime_edge = 0;
p.contractionRampTime = 1;
massDiagonal = ones(size(x));
massDiagonal(7:12) = 0;
opts = odeset('RelTol',1e-7, 'AbsTol',1e-10, 'MaxStep',0.1/p.scale.tau, ...
    'Mass',spdiags(massDiagonal,0,numel(x),numel(x)), ...
    'MStateDependence','none', 'MassSingular','yes', ...
    'JPattern',shell_jacobian_pattern(p));
time = (0:0.1:10).';
[states, p] = integrate_cell_scale(p, time/p.scale.tau, x, opts);
verifyGreaterThan(testCase, p.ruptureTime, 5);
verifyLessThan(testCase, p.ruptureTime, 6);
expectedDamage = 1-exp(-max(time-p.ruptureTime,0)/p.ruptureDamageTau);
verifyEqual(testCase, states(:,15), expectedDamage, 'AbsTol',1e-6);
verifyEqual(testCase, states(time < p.ruptureTime,15), ...
    zeros(nnz(time < p.ruptureTime),1), 'AbsTol',1e-12);
end

function [p, x] = twoCellFixture
% Two cells, one initially unstressed 10 um link, no surface or anchor.
p.N = 2; p.Ne = 1; p.nPos = 6; p.nVel = 6;
p.nodeCoords = [0 0 0; 1e-5 0 0];
p.edge_i = 1; p.edge_j = 2;
p.edgeAccumI = sparse([1;0]); p.edgeAccumJ = sparse([0;1]);
p.L0_edge0 = 1e-5;
p.m = 5e-13;
p.k1_edge = 0.002; p.k2_edge = 0.002;
p.c_edge = 0.1; p.c_series_edge = 0.1;
p.connectionModelId_edge = 5;
p.surfaceConstraintMode = 'none';
p.c_drag = 1e-4; p.k_surf = 5e-9; p.c_surf = 0.001;
p.g_floor_stiff = 0; p.g_floor_damp = 1e-6;
p.damageMask = 0; p.extensionCritFrac_edge = inf;
p.ruptureEnableTime = 0; p.rupturePatchCascade_on = true;
p.ruptureDamageTau = 0.5;
p.contractionRateFrac_edge = 0;
p.minRestLengthFrac_edge = 1;
p.contractionStartTime_edge = inf;
p.contractionRampTime = 900;
p.activeTraction_on = false; p.activeTractionForce = 1e-10;
p.adhesionClutch_on = false;
p.adhesionStiffness = 0.001;
p.adhesionLoadThresholdTip = 5e-10; p.adhesionLoadThresholdBody = 3e-9;
p.adhesionDragStuck = 5e-5; p.adhesionDragSliding = 1e-5;
p.adhesionLoadRegularization = 1e-15;
p.anchor_on = false; p.k_anchor = 0.01; p.c_anchor = 0.001;
p.mechanicsProfile = 'cell-scale';
p.scale.L = 1e-5; p.scale.k2 = p.k2_edge; p.scale.tau = 50;
p.scale.inertia = p.m/(p.scale.k2*p.scale.tau^2);
x = zeros(17,1);
x(1:6) = reshape(p.nodeCoords.',[],1)/p.scale.L;
end
