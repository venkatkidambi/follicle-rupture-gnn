function fit = tune_moving_network_simulation(maxEvaluations, resumeSearch)
%TUNE_MOVING_NETWORK_SIMULATION Reopen all previously fitted degrees of freedom.
% Optimize graph density/anisotropy, local link swaps, weak-patch location and
% size, seven mechanical coefficients, and three resistance controls. The
% measured ellipsoid history remains fixed; node motion follows force balance.
if nargin<1, maxEvaluations=160; end
if nargin<2, resumeSearch=false; end
validateattributes(maxEvaluations,{'numeric'},{'scalar','integer','>=',80});
fit=tune_network_simulation(maxEvaluations,resumeSearch,'moving');
end
