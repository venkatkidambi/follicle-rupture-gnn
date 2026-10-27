function [states, param] = integrate_cell_scale(param, times, initialState, opts)
%INTEGRATE_CELL_SCALE Resolve the irreversible patch latch as a terminal event.
% Before rupture, damage remains zero. At the first eligible extension
% threshold crossing, restart with patch damage evolving as (1-D)/tauBreak.
% This implements the existing cascade law without letting numerical
% Jacobian perturbations of an intact D=0 state activate it prematurely.
keep = [1:param.nPos, param.nPos+param.nVel+1:numel(initialState)];
velocity = param.nPos+(1:param.nVel);
initialState = initialState(keep);
pattern = shell_jacobian_pattern(param);
% Eliminating local velocities introduces two-hop position/state couplings.
if any(ismember(param.connectionModelId_edge,[2 4]))
    reducedPattern = []; % coupled drag inversion can make the Jacobian dense
else
    localInverse = kron(speye(param.N),ones(3));
    reducedPattern = spones(pattern(keep,keep) + ...
        pattern(keep,velocity)*localInverse*pattern(velocity,keep));
end
opts = odeset(opts,'Mass',[],'MassSingular','no','JPattern',reducedPattern);
if ~param.rupturePatchCascade_on
    opts = odeset(opts, 'Events', []);
elseif times(1)*param.scale.tau >= param.ruptureEnableTime && ...
        ruptureThreshold(times(1), initialState, param) <= 0
    param.ruptureTime = times(1)*param.scale.tau;
    opts = odeset(opts, 'Events', []);
else
    param.ruptureTime = inf;
    opts = odeset(opts, 'Events', @(t,x) ruptureThreshold(t,x,param));
end
first = ode15s(@(t,x) overdamped_shell_rhs(t,x,param), ...
    [times(1) times(end)], initialState, opts);
if isfield(first, 'xe') && ~isempty(first.xe)
    ruptureTimeDim = first.xe(end);
    param.ruptureTime = ruptureTimeDim * param.scale.tau;
    fprintf('[RUPTURE EVENT] first extension threshold at %.9g s\n', param.ruptureTime);
    before = times <= ruptureTimeDim;
    states = zeros(numel(times), numel(initialState));
    states(before,:) = deval(first, times(before)).';
    if any(~before)
        opts = odeset(opts, 'Events', []);
        second = ode15s(@(t,x) overdamped_shell_rhs(t,x,param), ...
            [ruptureTimeDim times(end)], first.ye(:,end), opts);
        requireComplete(second, times(end));
        states(~before,:) = deval(second, times(~before)).';
    end
else
    requireComplete(first, times(end));
    states = deval(first, times).';
end
fullStates = zeros(numel(times),numel(keep)+param.nVel);
for idx = 1:numel(times)
    [~, fullState] = overdamped_shell_rhs(times(idx),states(idx,:).',param);
    fullStates(idx,:) = fullState.';
end
states = fullStates;
end

function [value, isterminal, direction] = ruptureThreshold(t, state, param)
pos = reshape(state(1:param.nPos), 3, []).' * param.scale.L;
edges = param.damageMask(:) > 0;
lengths = vecnorm(pos(param.edge_j(edges),:) - pos(param.edge_i(edges),:), 2, 2);
natural = contracting_natural_lengths(t*param.scale.tau, param);
margin = param.extensionCritFrac_edge(edges) - ...
    (lengths - natural(edges)) ./ param.L0_edge0(edges);
value = min(margin);
if isempty(value)
    value = 1;
elseif t*param.scale.tau < param.ruptureEnableTime
    value = max(value, eps);
end
isterminal = 1;
direction = -1;
end

function requireComplete(solution, finalTime)
if abs(solution.x(end) - finalTime) > 1e-10 * max(finalTime,1)
    error('follicle_rupture:incompleteIntegration', ...
        'Cell-scale integration stopped early; no output was saved.');
end
end
