function [Lnatural, LnaturalDot] = contracting_natural_lengths(t, param)
% Prescribed natural-length contraction, equivalent to
% Ln=max(L0-rate*(t-tStart), Lmin) in the 1D reference model.
L0 = param.L0_edge0(:);
rateFrac = param.contractionRateFrac_edge(:);
floorFrac = param.minRestLengthFrac_edge(:);
startTime = param.contractionStartTime_edge(:);
riseTime = max(param.contractionRampTime, 1e-8);
isActiveEdge = rateFrac > 0;
rawDuration = zeros(size(rateFrac));
rawDuration(isActiveEdge) = (1 - floorFrac(isActiveEdge)) ./ ...
    rateFrac(isActiveEdge);
elapsed = max(t - startTime, 0);
effectiveElapsed = zeros(size(rateFrac));
rateActivation = zeros(size(rateFrac));

isRampUp = isActiveEdge & elapsed < riseTime;
if any(isRampUp)
    rampCoordinate = elapsed(isRampUp) ./ riseTime;
    rampIntegral = rampCoordinate.^3 - 0.5 * rampCoordinate.^4;
    effectiveElapsed(isRampUp) = riseTime .* rampIntegral;
    rateActivation(isRampUp) = smoothstep01(rampCoordinate);
end

isPlateau = isActiveEdge & elapsed >= riseTime & elapsed <= rawDuration;
effectiveElapsed(isPlateau) = elapsed(isPlateau) - 0.5 * riseTime;
rateActivation(isPlateau) = 1;

isRampDown = isActiveEdge & elapsed > rawDuration & ...
    elapsed < rawDuration + riseTime;
if any(isRampDown)
    rampCoordinate = (elapsed(isRampDown) - rawDuration(isRampDown)) ./ riseTime;
    rampIntegral = rampCoordinate.^3 - 0.5 * rampCoordinate.^4;
    effectiveElapsed(isRampDown) = rawDuration(isRampDown) - ...
        0.5 * riseTime + riseTime .* (rampCoordinate - rampIntegral);
    rateActivation(isRampDown) = 1 - smoothstep01(rampCoordinate);
end

isComplete = isActiveEdge & elapsed >= rawDuration + riseTime;
effectiveElapsed(isComplete) = rawDuration(isComplete);

lengthFrac = max(1 - rateFrac .* effectiveElapsed, floorFrac);
Lnatural = L0 .* lengthFrac;
LnaturalDot = -L0 .* rateFrac .* rateActivation;
end


function value = smoothstep01(coordinate)
coordinate = min(max(coordinate, 0), 1);
value = coordinate.^2 .* (3 - 2 * coordinate);
end


