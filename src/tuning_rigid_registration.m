function [aligned, rotation, translation] = tuning_rigid_registration(points, reference, mask)
%TUNING_RIGID_REGISTRATION Robust proper rotation/translation, no scale fit.
mask = mask & all(isfinite(points),2) & all(isfinite(reference),2);
if nnz(mask)<4
    error('follicle_rupture:registrationData','At least four valid reference cells are needed.');
end
X = points(mask,:); Y = reference(mask,:);
weights = ones(size(X,1),1);
for iteration = 1:5
    centerX = sum(weights.*X,1)/sum(weights);
    centerY = sum(weights.*Y,1)/sum(weights);
    [U,~,V] = svd((weights.*(X-centerX)).'*(Y-centerY));
    rotation = U*diag([1 1 det(U*V.')])*V.';
    translation = centerY-centerX*rotation;
    residual = vecnorm(X*rotation+translation-Y,2,2);
    weights = min(1,3e-6./max(residual,eps));
end
aligned = points*rotation+translation;
end
