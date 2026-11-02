function [surface, diagnostics] = fit_experimental_surface(data)
%FIT_EXPERIMENTAL_SURFACE Full ellipsoid at every registered measured frame.
% Only training cells enter the fits. All frame times enter the prescribed
% geometry, so subsequent temporal validation is conditional, not forecasting.
Nt = numel(data.time);
surface.time = data.time;
surface.center = zeros(Nt,3);
surface.quadric = zeros(3,3,Nt);
diagnostics.semiAxes_um = zeros(Nt,3);
diagnostics.rmsDistance_um = zeros(Nt,1);
diagnostics.cellCount = zeros(Nt,1);
unit = 1e-4;
for frame = 1:Nt
    mask = data.trainNodes & data.valid(:,frame);
    points = data.position(mask,:,frame);
    assert(size(points,1)>=10,'At least ten valid training cells are required.');
    origin = mean(points,1);
    x = (points-origin)/unit;
    design = [x.^2,2*x(:,1).*x(:,2),2*x(:,1).*x(:,3), ...
        2*x(:,2).*x(:,3),x];
    assert(rank(design)==9,'Experimental ellipsoid fit is rank deficient.');
    coefficient = design\ones(size(x,1),1);
    Q = [coefficient(1) coefficient(4) coefficient(5); ...
         coefficient(4) coefficient(2) coefficient(6); ...
         coefficient(5) coefficient(6) coefficient(3)];
    center = -0.5*(Q\coefficient(7:9));
    Q = Q/(1+center.'*Q*center);
    [~,flag] = chol(Q);
    assert(flag==0,'Frame %d did not yield a positive-definite ellipsoid.',frame);
    surface.center(frame,:) = origin+unit*center.';
    surface.quadric(:,:,frame) = Q/unit^2;
    diagnostics.semiAxes_um(frame,:) = sort(1./sqrt(eig(Q)),'descend').'*unit*1e6;
    r = x-center.';
    residual = sum((r*Q).*r,2)-1;
    distance = residual./vecnorm(2*r*Q,2,2)*unit;
    diagnostics.rmsDistance_um(frame) = sqrt(mean(distance.^2))*1e6;
    diagnostics.cellCount(frame) = nnz(mask);
end
diagnostics.volume_um3 = 4*pi/3*prod(diagnostics.semiAxes_um,2);
surface.description = ['Full positive-definite ellipsoid fitted to training cells at every ' ...
    'pose-corrected frame; centers and quadratic forms interpolated linearly'];
end
