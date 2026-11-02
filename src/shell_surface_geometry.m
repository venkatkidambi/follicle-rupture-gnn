function [level, normal, normalSpeed] = shell_surface_geometry(t, pos, param)
%SHELL_SURFACE_GEOMETRY Shared fixed/moving surface and its normal velocity.
% Linear interpolation of positive-definite quadratic forms stays positive
% definite and avoids eigenvector flips in nearly axisymmetric ellipsoids.
if isfield(param,'movingSurface')
    s = param.movingSurface;
    interval = find(s.time<=t,1,'last');
    if isempty(interval), interval=1; end
    interval = min(interval,numel(s.time)-1);
    dt = s.time(interval+1)-s.time(interval);
    weight = min(max((t-s.time(interval))/dt,0),1);
    Qdot = (s.quadric(:,:,interval+1)-s.quadric(:,:,interval))/dt;
    centerDot = (s.center(interval+1,:)-s.center(interval,:))/dt;
    Q = s.quadric(:,:,interval)+weight*dt*Qdot;
    center = s.center(interval,:)+weight*dt*centerDot;
    if t<s.time(1) || t>s.time(end), Qdot(:)=0; centerDot(:)=0; end
    r = pos-center;
    gradient = 2*r*Q;
    level = sum((r*Q).*r,2);
    partialTime = sum((r*Qdot).*r,2)-sum(gradient.*centerDot,2);
else
    Q = diag(1./[param.ell_a param.ell_b param.ell_c].^2);
    gradient = 2*pos*Q;
    level = sum((pos*Q).*pos,2);
    partialTime = zeros(size(level));
end
gradientNorm = vecnorm(gradient,2,2)+1e-18;
normal = gradient./gradientNorm;
normalSpeed = -partialTime./gradientNorm;
end
