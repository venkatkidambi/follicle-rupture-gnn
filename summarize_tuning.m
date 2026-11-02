function summarize_tuning(fitFile)
%SUMMARIZE_TUNING Export a comparison figure and auditable numerical report.
repo=fileparts(mfilename('fullpath'));
if nargin<1
    fitFile=fullfile(repo,'results','network_tuning','best_fit.mat');
    if ~isfile(fitFile), fitFile=fullfile(repo,'results','parameter_tuning','best_fit.mat'); end
end
loaded=load(fitFile,'fit'); fit=loaded.fit; data=fit.data;
isMoving=isfield(fit.param,'movingSurface');
cellTable=table(fit.param.experimentalCellIds(:),data.trainNodes, ...
    sqrt(mean(fit.metrics.error_um.^2,2)),fit.metrics.error_um(:,end), ...
    fit.metrics.observedAxial(:,end)*1e6,fit.metrics.predictedAxial(:,end)*1e6, ...
    'VariableNames',{'CellID','UsedForTraining','TrajectoryRMSE_um','FinalError_um', ...
    'ObservedFinalAxial_um','PredictedFinalAxial_um'});
writetable(cellTable,fullfile(fileparts(fitFile),'cell_errors.csv'));
fig=plot_tuning_comparison(fitFile); close(fig);
if isfield(fit,'structure'), fig=plot_network_fit(fitFile); close(fig); end
if isMoving, fig=plot_moving_surface_fit(fitFile); close(fig); end
reportFile=fullfile(fileparts(fitFile),'tuning_report.txt');
fid=fopen(reportFile,'w');
assert(fid>=0,'Cannot open tuning report.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'FOLLICLE TRAJECTORY PARAMETER FIT\n\n');
fprintf(fid,'%d cells, %d frames; 0 to %.0f s, assuming uniform %.9f s intervals.\n', ...
    fit.param.N,numel(data.time),data.time(end),data.frameInterval);
fprintf(fid,'Source: %s\n',fit.param.experimentalNodeDataFile);
fprintf(fid,'Raw xyz spacing: %s micrometers.\n',mat2str(fit.param.units.coordinateToMeters*1e6,12));
fprintf(fid,'Use only cells present in frame 1, with original row identities.\n');
fprintf(fid,'Pose correction: %s.\n',data.registrationDescription);
fprintf(fid,'No registration is applied to the simulated node positions.\n\n');
fprintf(fid,'METHOD\n');
fprintf(fid,'Fit %d global parameters with bounded, cached coordinate search.\n',numel(fit.values));
if isMoving
    fprintf(fid,'Full ellipsoid fitted to training cells at every registered experimental frame.\n');
    fprintf(fid,'Interpolate centers and positive-definite quadratic forms linearly; preserve initial node surface levels.\n');
    fprintf(fid,'Normal damping uses velocity relative to the moving surface.\n');
    if isfield(fit,'structure')
        fprintf(fid,'Network connections, weak-patch location/size, and all ten mechanical/resistance controls are refitted.\n');
        referenceFile=fit.multiplierReferenceFile;
    else
        fprintf(fid,'Network and weak patch are frozen; ten mechanical/resistance controls are refitted.\n');
        referenceFile=fit.sourceFitFile;
    end
    fprintf(fid,'Multipliers reference the previous resistance fit: %s\n',referenceFile);
    fprintf(fid,'Geometry at ALL frames is prescribed external forcing; validation is conditional on that geometry.\n');
end
if isfield(fit,'fitType') && strcmp(fit.fitType,'resistance')
    fprintf(fid,'Three-control resistance extension of the saved network fit.\n');
    fprintf(fid,'Network, weak patch, geometry, seven fitted coefficients, and data split are frozen.\n');
    fprintf(fid,'Tip release threshold is fixed; body/tip ratio changes the existing linear threshold gradient.\n');
    fprintf(fid,'One multiplier scales background and sliding drag; attached drag is fixed.\n');
    fprintf(fid,'Bounds were fixed before search: ratio 1-12, stiffness and drag multipliers 0.2-5.\n');
    fprintf(fid,'Multipliers reference the frozen network fit, not a previous trial.\n');
    fprintf(fid,'Source fit: %s\n',fit.sourceFitFile);
end
if isfield(fit,'structure')
    fprintf(fid,'Also fit a static local graph and connected weak patch using training error only.\n');
    fprintf(fid,'Candidate links: symmetric %d-neighbor pool plus original links; connected graphs only.\n',fit.geometry.neighborPoolSize);
    fprintf(fid,'Search: anisotropic local neighbors, surface patch centers/sizes, and local edge swaps.\n');
    fprintf(fid,'Connections: %d; weak links: %d; rupture-patch cell IDs: %s.\n', ...
        fit.param.Ne,nnz(fit.param.damageMask),mat2str(fit.param.experimentalCellIds(fit.param.rupturePatchNodes).'));
    fprintf(fid,'Explicit connections and weak-link labels are saved in connections.csv.\n');
    if isfield(fit,'searchCoverage')
        fprintf(fid,'Every global control received an evaluated proposal: %d.\n',all(fit.searchCoverage.evaluatedControls));
        fprintf(fid,'Selected connections changed: %d; selected weak-link pairs changed: %d.\n', ...
            fit.searchCoverage.connectionsChanged,fit.searchCoverage.weakPairsChanged);
        proposals=fit.searchCoverage.proposals;
        fprintf(fid,'Evaluated graph / patch-center / patch-size / scalar / link-swap proposals: %d / %d / %d / %d / %d.\n', ...
            nnz(startsWith(proposals,'local graph')),nnz(startsWith(proposals,'patch center')), ...
            nnz(startsWith(proposals,'patch size')),nnz(startsWith(proposals,'mechanics')), ...
            nnz(strcmp(proposals,'local edge swap')));
    end
end
fprintf(fid,'Training objective: mean squared 3D position error (um^2).\n');
fprintf(fid,'%d training cells x %d training frames; %d held-out cells x %d held-out frames.\n', ...
    nnz(data.trainNodes),nnz(data.trainFrames),nnz(data.validationNodes),nnz(data.validationFrames));
fprintf(fid,'Every fifth spatially ordered cell and alternating frames are excluded from the mechanical fitting objective.\n');
fprintf(fid,'Validation cells do not participate in rigid registration. Held-out frames use\n');
fprintf(fid,'training cap cells at that frame for pose correction; this is interpolation\n');
fprintf(fid,'validation within one specimen, not independent biological validation or forecasting.\n');
if isMoving
    fprintf(fid,'Training trajectories also determine global surface geometry at all frames. No prescribed node displacements,\n');
else
    fprintf(fid,'Measured trajectories enter only the objective. No prescribed displacements,\n');
end
fprintf(fid,'per-node target-position forces, node-specific fitted schedules, or imposed rupture time.\n');
fprintf(fid,'The existing single global natural-length contraction law supplies active loading.\n');
fprintf(fid,'Velocities follow overdamped mass-spring-damper force balance.\n\n');
fprintf(fid,'RESULTS (micrometers, 3D Euclidean errors)\n');
fprintf(fid,'Model                    Training RMSE   Joint holdout RMSE\n');
if isMoving
    fprintf(fid,'Holdout scores below are conditional on measured surface geometry, not temporal forecasts.\n');
    fprintf(fid,'Prior fixed surface      %13.6f   %18.6f\n',fit.fixedSurfaceMetrics.trainRMSE_um,fit.fixedSurfaceMetrics.validationRMSE_um);
    fprintf(fid,'Moving before refit      %13.6f   %18.6f\n',fit.movingSurfaceInitialMetrics.trainRMSE_um,fit.movingSurfaceInitialMetrics.validationRMSE_um);
end
fprintf(fid,'Original cell-scale      %13.6f   %18.6f\n',fit.baselineMetrics.trainRMSE_um,fit.baselineMetrics.validationRMSE_um);
fprintf(fid,'No-motion benchmark      %13.6f   %18.6f\n',fit.stationaryMetrics.trainRMSE_um,fit.stationaryMetrics.validationRMSE_um);
if isfield(fit,'parameterOnlyMetrics')
    fprintf(fid,'Previous parameter fit   %13.6f   %18.6f\n',fit.parameterOnlyMetrics.trainRMSE_um,fit.parameterOnlyMetrics.validationRMSE_um);
end
if isfield(fit,'networkMetrics')
    fprintf(fid,'Frozen network benchmark %13.6f   %18.6f\n',fit.networkMetrics.trainRMSE_um,fit.networkMetrics.validationRMSE_um);
end
fprintf(fid,'Fitted mechanics         %13.6f   %18.6f\n',fit.metrics.trainRMSE_um,fit.metrics.validationRMSE_um);
fprintf(fid,'Training mean squared error: %.8f um^2\n',fit.metrics.objective);
fprintf(fid,'Held-out median error: %.6f um\n',fit.metrics.validationMedian_um);
late=data.validationNodes & data.validationFrames.' & (data.time.'>=1200) & data.valid;
fprintf(fid,'Last 10 min held-out RMSE: fitted %.6f um; no motion %.6f um\n', ...
    sqrt(mean(fit.metrics.error_um(late).^2)),sqrt(mean(fit.stationaryMetrics.error_um(late).^2)));
fprintf(fid,'Raw unregistered all-cell/all-frame RMSE (diagnostic only): %.6f um\n',fit.metrics.rawRMSE_um);
fprintf(fid,'Predicted extension-triggered rupture: %.6f s (not an observed rupture label).\n',fit.prediction.param.ruptureTime);
fprintf(fid,'Final band mean axial displacements, observed / predicted (um):\n');
for band=1:3
    mask=data.bandMasks(:,band);
    fprintf(fid,'  %d-%d%%: %.6f / %.6f\n',data.percentileRanges(band,:), ...
        mean(fit.metrics.observedAxial(mask,end))*1e6,mean(fit.metrics.predictedAxial(mask,end))*1e6);
end
fprintf(fid,'\nPARAMETERS (SI units; rate in 1/s, time constants in s, multipliers/strain dimensionless)\n');
for idx=1:numel(fit.values)
    fprintf(fid,'%s = %.12g; bounds [%.12g, %.12g]\n', ...
        fit.parameterNames{idx},fit.values(idx),fit.lowerBounds(idx),fit.upperBounds(idx));
end
fprintf(fid,'Fixed k1 = k2 = %.9g N/m; mass %.9g kg/cell.\n',mean(fit.param.k1_edge),fit.param.m);
fprintf(fid,'Fitted edge dashpot %.9g N*s/m; adhesion tip/body loads %.9g / %.9g N.\n', ...
    mean(fit.param.c_edge),fit.param.adhesionLoadThresholdTip,fit.param.adhesionLoadThresholdBody);
if isfield(fit,'fitType') && strcmp(fit.fitType,'resistance')
    fprintf(fid,'Adhesion stiffness %.9g N/m; background/sliding drag %.9g / %.9g N*s/m.\n', ...
        fit.param.adhesionStiffness,fit.param.c_drag,fit.param.adhesionDragSliding);
    fprintf(fid,'Fixed attached drag %.9g N*s/m.\n',fit.param.adhesionDragStuck);
    fprintf(fid,'Seven coefficients held fixed (adhesion multiplier is the original tip-scale reference):\n');
    for idx=1:numel(fit.fixedParameterValues)
        fprintf(fid,'  %s = %.12g\n',fit.fixedParameterNames{idx},fit.fixedParameterValues(idx));
    end
    fprintf(fid,'\nLOCAL SENSITIVITY AT THE NETWORK BENCHMARK\n');
    fprintf(fid,'Central log perturbations x1.2 and /1.2; training positions only.\n');
    for idx=1:numel(fit.values)
        fprintf(fid,'%s: position derivative RMS %.6g um/log-unit; train RMSE (-,+) %.6f / %.6f um.\n', ...
            fit.parameterNames{idx},fit.sensitivity.derivativeRMS_um(idx),fit.sensitivity.trainRMSE_um(idx,:));
    end
    fprintf(fid,'Pairwise derivative cosine matrix (near +/-1 indicates local similarity):\n');
    fprintf(fid,'  % .6f % .6f % .6f\n',fit.sensitivity.derivativeCosine.');
    fprintf(fid,'These local diagnostics do not establish parameter identifiability.\n');
    atBound=abs(log(fit.values./fit.lowerBounds))<1e-6 | abs(log(fit.values./fit.upperBounds))<1e-6;
    fprintf(fid,'Selected parameters at a search bound: %s\n',strjoin(fit.parameterNames(atBound),', '));
end
fprintf(fid,'\nCOMPUTATION AND LIMITS\n');
if isMoving
    d=fit.surfaceDiagnostics;
    fprintf(fid,'Initial/final fitted ellipsoid volumes %.9g / %.9g um^3 (%.4f%% change).\n', ...
        d.volume_um3(1),d.volume_um3(end),100*(d.volume_um3(end)/d.volume_um3(1)-1));
    fprintf(fid,'Maximum approximate surface-fit RMS residual %.6f um.\n',max(d.rmsDistance_um));
    fprintf(fid,'Ellipsoid volume is a shape-fit proxy; an open shell does not establish true enclosed volume.\n');
end
fprintf(fid,'%d coarse evaluations (%d completed; %d censored/failed), %d final fine solves.\n', ...
    numel(fit.history),nnz(strcmp({fit.history.status},'ok')),nnz(~strcmp({fit.history.status},'ok')),numel(fit.fineCandidates));
unresolved=0;
for idx=find(~strcmp({fit.history.status},'ok'))
    if isfield(fit,'structure')
        failed=fit.history(idx).candidate;
        resolved=arrayfun(@(entry) strcmp(entry.status,'ok') && ...
            isequal(entry.candidate.values,failed.values) && ...
            isequal(entry.candidate.edges,failed.edges) && ...
            isequal(entry.candidate.weakEdges,failed.weakEdges),fit.history);
    else
        resolved=arrayfun(@(entry) strcmp(entry.status,'ok') && ...
            isequal(entry.values,fit.history(idx).values),fit.history);
    end
    unresolved=unresolved+~any(resolved);
end
fprintf(fid,'%d recorded censored/failed attempts have no later successful solve at the same parameters.\n',unresolved);
fprintf(fid,'Final search pass wall time %.2f s.\n',fit.totalWallSeconds);
fprintf(fid,'Recorded coarse-evaluation wall time %.2f s; final fine solves %.2f s.\n', ...
    sum([fit.history.wallSeconds]),sum([fit.fineCandidates.wallSeconds]));
fprintf(fid,'Search selects by training loss only, including final fine-solver ranking.\n');
fprintf(fid,'Uniform timing and rigid motion of the reference cap are assumptions.\n');
fprintf(fid,'Absolute stiffness/force scale is not identifiable from these trajectories:\n');
fprintf(fid,'common scaling of forces and damping leaves overdamped motion unchanged.\n');
fprintf(fid,'The mass and spring scale remain assumptions, not fitted measurements.\n');
fprintf(fid,'This is a local parameter fit with no confidence intervals or uniqueness claim.\n');
if isfield(fit,'structure')
    fprintf(fid,'The inferred graph/weak patch are fitted hypotheses, not identified biological contacts.\n');
    fprintf(fid,'Search is bounded and heuristic, not an exhaustive optimum over all graphs.\n');
    fprintf(fid,'Rupture time remains an emergent extension-threshold event.\n');
end
fprintf(fid,'Residual trajectory mismatch can reflect missing mechanics, not just parameter error.\n');
fprintf('Wrote %s\n',reportFile);
end
