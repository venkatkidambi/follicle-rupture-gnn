function [derivative, fullState] = overdamped_shell_rhs(t, state, param)
%OVERDAMPED_SHELL_RHS Eliminate velocity by solving the linear drag balance.
% The reduced state omits velocity; all other states retain their ordering.
% Forces and viscoelastic rates still come from the shared mechanical RHS.
nPos = param.nPos;
fullState = [state(1:nPos); zeros(param.nVel,1); state(nPos+1:end)];
[~, forces] = follicle_shell_rhs(t, fullState, param);
pos = reshape(state(1:nPos),3,[]).' * param.scale.L;
mode = lower(strtrim(char(param.surfaceConstraintMode)));
useNormals = ismember(mode, {'ellipsoid','reference-ellipsoid','project-ellipsoid'});
normal = zeros(param.N,3);
if useNormals
    [~, normal] = shell_surface_geometry(t*param.scale.tau,pos,param);
end
baseDrag = param.c_drag * ones(param.N,1);
if param.anchor_on
    baseDrag = baseDrag + param.c_anchor * param.anchorNodeMask;
end
adhesionDrag = zeros(param.N,1);
if param.adhesionClutch_on
    release = min(max(state(nPos+3*param.Ne+(1:param.N)),0),1);
    adhesionDrag = (1-release)*param.adhesionDragStuck + release*param.adhesionDragSliding;
    adhesionDrag(param.anchorNodeMask > 0) = 0;
end
tangentDrag = baseDrag + adhesionDrag;
normalDrag = baseDrag + double(useNormals)*param.c_surf;
if any(tangentDrag <= 0 | normalDrag <= 0)
    error('follicle_rupture:invalidDrag', 'Overdamped motion requires positive node drag.');
end
F = forces.nodeNetForce;
directDamping = ismember(param.connectionModelId_edge, [2 4]);
if ~any(directDamping)
    normalForce = sum(F.*normal,2).*normal;
    velocity = (F-normalForce)./tangentDrag + normalForce./normalDrag;
else
    % Pure dampers and Kelvin-Voigt links couple neighboring velocities.
    C = spdiags(repelem(tangentDrag,3),0,nPos,nPos);
    normalCoefficient = normalDrag-tangentDrag;
    for row = 1:3
        for col = 1:3
            C = C + sparse(row:3:nPos,col:3:nPos, ...
                normalCoefficient.*normal(:,row).*normal(:,col),nPos,nPos);
        end
    end
    edges = find(directDamping);
    i = param.edge_i(edges); j = param.edge_j(edges);
    separation = pos(j,:)-pos(i,:);
    length = vecnorm(separation,2,2);
    direction = separation./max(length,1e-18);
    direction(length < 1e-18,:) = 0;
    projectedI = direction; projectedJ = direction;
    if useNormals && param.projectEdgeForcesToSurface
        projectedI = direction-sum(direction.*normal(i,:),2).*normal(i,:);
        projectedJ = direction-sum(direction.*normal(j,:),2).*normal(j,:);
    end
    damage = min(max(state(nPos+2*param.Ne+edges),0),1);
    damping = (1-damage).*max(max(1-damage,param.g_floor_damp).*param.c_edge(edges),1e-12);
    row = repelem((1:numel(edges)).',3);
    colsI = reshape((3*i-2+(0:2)).',[],1);
    colsJ = reshape((3*j-2+(0:2)).',[],1);
    directions = reshape(direction.',[],1);
    B = sparse([row;row],[colsI;colsJ],[directions;-directions],numel(edges),nPos);
    G = sparse([colsI;colsJ],[row;row], ...
        [reshape(projectedI.',[],1);-reshape(projectedJ.',[],1)],nPos,numel(edges));
    C = C + G*spdiags(damping,0,numel(edges),numel(edges))*B;
    velocity = reshape(C\reshape(F.',[],1),3,[]).';
end
fullState(nPos+(1:param.nVel)) = reshape(velocity.',[],1)*param.scale.tau/param.scale.L;
fullDerivative = follicle_shell_rhs(t, fullState, param);
derivative = [fullDerivative(1:nPos); fullDerivative(nPos+param.nVel+1:end)];
end
