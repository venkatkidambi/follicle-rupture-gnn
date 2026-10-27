function param = tuning_parameters(base, values)
%TUNING_PARAMETERS Update a few global mechanical/activation coefficients.
% No observed trajectories, fitted node positions, or displacement feedback.
% values: contraction rate, global activation onset, Maxwell relaxation time,
% adhesion load multiplier, adhesion detachment time, rupture strain, and
% optionally the global contraction rise time (six-value presets retain it).
validateattributes(values,{'numeric'},{'real','finite','vector'});
assert(ismember(numel(values),[6 7]),'Expected six or seven fitted coefficients.');
param = base;
if isfield(param,'ruptureTime'), param = rmfield(param,'ruptureTime'); end
if numel(values)==7
    param.activationRiseTime=values(7);
    param.contractionRampTime=values(7);
    param.activeTractionRampTime=values(7);
end
param.contractionRateFrac = values(1);
param.prestressStartTime = values(2);
param.c_edge(:) = values(3)*param.k2_edge;
param.adhesionLoadThresholdTip = base.adhesionLoadThresholdTip*values(4);
param.adhesionLoadThresholdBody = base.adhesionLoadThresholdBody*values(4);
param.adhesionLoadThreshold_node = base.adhesionLoadThreshold_node*values(4);
param.adhesionReleaseTau = values(5);
param.extensionCritFrac = values(6);
active = param.retractionMask>0;
param.contractionRateFrac_edge(active) = values(1);
param.contractionStartTime_edge(active) = values(2);
param.contractionDuration_edge(active) = (1-param.minRestLengthFrac_edge(active))/values(1) + ...
    param.contractionRampTime;
param.prestressDuration = max(param.contractionDuration_edge(active));
param.activeTractionStartTime_node(:) = values(2);
param.activeTractionStartTime_node(param.anchorNodeMask>0) = inf;
param.activeTractionDuration = param.prestressDuration;
param.ruptureEnableTime = values(2);
param.extensionCritFrac_edge(param.damageMask>0) = values(6);
param.scale.tau = mean(param.c_edge)/mean(param.k2_edge);
param.scale.inertia = param.m/(param.scale.k2*param.scale.tau^2);
param.scale.epsilon = param.m*(mean(param.k1_edge)+mean(param.k2_edge))/mean(param.c_edge)^2;
% Reject overlapping contraction ramps rather than exploit that unsupported case.
if any((1-param.minRestLengthFrac_edge(active))/values(1)<param.contractionRampTime)
    error('follicle_rupture:tuningRamp','Contraction plateau duration must exceed ramp duration.');
end
end
