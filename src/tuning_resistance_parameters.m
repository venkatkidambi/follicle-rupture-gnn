function param = tuning_resistance_parameters(base, values)
%TUNING_RESISTANCE_PARAMETERS Change exactly three resistance controls.
% values = [body/tip release-load ratio, tether stiffness multiplier,
%           shared background/sliding drag multiplier], relative to base.
validateattributes(values,{'numeric'},{'real','finite','positive','vector','numel',3});
param = base;
if isfield(param,'ruptureTime'), param=rmfield(param,'ruptureTime'); end
tip = base.adhesionLoadThresholdTip;
oldRatio = base.adhesionLoadThresholdBody/tip;
if values(1) ~= oldRatio
    fraction = (base.adhesionLoadThreshold_node-tip) / ...
        (base.adhesionLoadThresholdBody-tip);
    assert(all(fraction >= -1e-12 & fraction <= 1+1e-12), ...
        'Expected the existing linear tip-to-body threshold profile.');
    param.adhesionLoadThresholdBody = tip*values(1);
    param.adhesionLoadThreshold_node = tip + ...
        (param.adhesionLoadThresholdBody-tip)*fraction;
end
param.adhesionStiffness = base.adhesionStiffness*values(2);
param.c_drag = base.c_drag*values(3);
param.adhesionDragSliding = base.adhesionDragSliding*values(3);
end
