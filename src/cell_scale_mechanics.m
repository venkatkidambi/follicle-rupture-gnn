function param = cell_scale_mechanics(param, volume_um3, density_kg_m3, forceScale)
%CELL_SCALE_MECHANICS Apply a common force scale and a wet cell mass estimate.
% These are order-of-magnitude estimates, not a fit to measured Drosophila
% rupture forces. Scaling every force coefficient together preserves elastic
% force ratios, strain thresholds and c/k relaxation times. Correcting the
% artificial legacy inertia with overdamped motion changes transient motion.
validateattributes(volume_um3, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(density_kg_m3, {'numeric'}, {'real','scalar','finite','positive'});
validateattributes(forceScale, {'numeric'}, {'real','scalar','finite','positive'});

forceFields = {'k1_edge', 'k2_edge', 'c_edge', 'c_series_edge', ...
    'k_surf', 'c_surf', 'c_drag', 'activeTractionForce', ...
    'adhesionStiffness', 'adhesionLoadThresholdTip', ...
    'adhesionLoadThresholdBody', 'adhesionDragStuck', ...
    'adhesionDragSliding', 'adhesionLoadRegularization', ...
    'k_anchor', 'c_anchor'};
for idx = 1:numel(forceFields)
    name = forceFields{idx};
    param.(name) = forceScale * param.(name);
end
param.cellVolume_um3 = volume_um3;
param.cellDensity_kg_m3 = density_kg_m3;
param.forceScaleFromLegacy = forceScale;
param.m = density_kg_m3 * volume_um3 * 1e-18;
end
