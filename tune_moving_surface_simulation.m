function fit = tune_moving_surface_simulation(maxEvaluations, resumeSearch)
%TUNE_MOVING_SURFACE_SIMULATION Refit ten global controls with measured geometry.
% Keep the previously fitted network/weak patch fixed. No node trajectory,
% target-position force, or observed rupture time enters the forward solve.
if nargin<1, maxEvaluations=120; end
if nargin<2, resumeSearch=false; end
validateattributes(maxEvaluations,{'numeric'},{'scalar','integer','>=',21});
repo=fileparts(mfilename('fullpath')); addpath(fullfile(repo,'src'));
sourceFile=fullfile(repo,'results','resistance_tuning','best_fit.mat');
loaded=load(sourceFile,'fit'); previous=loaded.fit;
base=previous.param; data=previous.data;
[surface,diagnostics]=fit_experimental_surface(data);
base.movingSurface=surface;
base.surfaceConstraintMode='reference-ellipsoid';
base.surfaceLevel0=shell_surface_geometry(data.time(1),base.nodeCoords,base);
names={'contractionRateFrac','prestressStartTime','maxwellRelaxationTime', ...
    'adhesionThresholdMultiplier','adhesionReleaseTau','extensionCritFrac', ...
    'activationRiseTime','bodyTipAdhesionRatio','adhesionStiffnessMultiplier', ...
    'backgroundSlidingDragMultiplier'};
initial=[base.contractionRateFrac,base.prestressStartTime, ...
    mean(base.c_edge./base.k2_edge),1,base.adhesionReleaseTau, ...
    base.extensionCritFrac,base.contractionRampTime, ...
    base.adhesionLoadThresholdBody/base.adhesionLoadThresholdTip,1,1];
lower=[2e-5 0 5 .15 2 .01 60 1 .2 .2];
upper=[4e-4 1200 300 8 180 .15 900 12 5 5];
logarithmic=true(1,10); logarithmic(2)=false;
low=lower; high=upper;
low(logarithmic)=log(low(logarithmic)); high(logarithmic)=log(high(logarithmic));
out=fullfile(repo,'results','moving_surface_tuning');
if ~isfolder(out), mkdir(out); end
configuration=struct('version',1,'base',base,'data',data,'lower',lower,'upper',upper);
checkpoint=fullfile(out,'search_history.mat');
history=struct('values',{},'objective',{},'trainRMSE_um',{},'wallSeconds',{},'status',{});
if resumeSearch && isfile(checkpoint)
    saved=load(checkpoint);
    assert(isequaln(saved.configuration,configuration),'Checkpoint inputs changed.');
    history=saved.history;
end
cache=containers.Map('KeyType','char','ValueType','double');
bestValues=initial; bestObjective=inf; started=tic;
for idx=1:numel(history)
    cache(key(history(idx).values))=idx;
    if history(idx).objective<bestObjective
        bestObjective=history(idx).objective; bestValues=history(idx).values;
    end
end
incumbent=tuning_simulate(parameters(initial),data.time,'fine',300,false);
incumbentMetrics=tuning_score(incumbent,data);
fprintf('Moving-surface initial train %.5f um; fixed-surface prior %.5f um.\n', ...
    incumbentMetrics.trainRMSE_um,previous.metrics.trainRMSE_um);
evaluate(initial);
for step=[.16 .08 .04 .02]
    for sweep=1:3
        before=bestObjective; count=numel(history);
        for dim=[1 2 7 4 8 9 10 5 6 3]
            center=encode(bestValues);
            for direction=[-1 1]
                z=center; z(dim)=min(max(z(dim)+direction*step,0),1);
                evaluate(decode(z));
            end
        end
        if numel(history)>=maxEvaluations || numel(history)==count || bestObjective>=before-1e-8, break; end
    end
    if numel(history)>=maxEvaluations, break; end
end
bestResult=incumbent; bestMetrics=incumbentMetrics; selected=initial;
fine=struct('values',initial,'objective',incumbentMetrics.objective, ...
    'metrics',incumbentMetrics,'wallSeconds',incumbent.wallSeconds);
[~,rank]=sort([history.objective]); rank=rank(isfinite([history(rank).objective]));
for idx=1:min(5,numel(rank))
    values=history(rank(idx)).values;
    if strcmp(key(values),key(initial)), continue; end
    result=tuning_simulate(parameters(values),data.time,'fine',300,false);
    metrics=tuning_score(result,data);
    fine(end+1)=struct('values',values,'objective',metrics.objective, ...
        'metrics',metrics,'wallSeconds',result.wallSeconds); %#ok<AGROW>
    if metrics.objective<bestMetrics.objective
        selected=values; bestResult=result; bestMetrics=metrics;
    end
end
fit=struct('fitType','moving-surface','parameterNames',{names},'values',selected, ...
    'lowerBounds',lower,'upperBounds',upper,'param',parameters(selected), ...
    'data',data,'prediction',bestResult,'metrics',bestMetrics,'history',history, ...
    'fineCandidates',fine,'baselineMetrics',previous.baselineMetrics, ...
    'stationaryMetrics',previous.stationaryMetrics,'fixedSurfaceMetrics',previous.metrics, ...
    'movingSurfaceInitialMetrics',incumbentMetrics,'surfaceDiagnostics',diagnostics, ...
    'sourceFitFile',sourceFile,'frozenBase',base,'totalWallSeconds',toc(started), ...
    'forwardModel','Measured moving ellipsoid; fixed network; ten global force-law controls', ...
    'validationDescription','Held-out cells at non-training frames, conditional on measured geometry at every frame');
save(fullfile(out,'best_fit.mat'),'fit','-v7.3');
writetable(table(names.',selected.',lower.',upper.', ...
    'VariableNames',{'Parameter','FittedValue','LowerBound','UpperBound'}),fullfile(out,'parameters.csv'));
writetable(table(data.time,surface.center*1e6,diagnostics.semiAxes_um, ...
    diagnostics.volume_um3,diagnostics.rmsDistance_um,diagnostics.cellCount, ...
    'VariableNames',{'Time_s','Center_um','SemiAxes_um','EllipsoidVolume_um3', ...
    'ApproxSurfaceResidualRMS_um','TrainingCellCount'}),fullfile(out,'ellipsoid_fits.csv'));
fid=fopen(fullfile(out,'tuning_report.txt'),'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'EXPERIMENTALLY CONDITIONED MOVING-ELLIPSOID FIT\n\n');
fprintf(fid,'45 measured frames across 1800 s; training cells only fit full ellipsoids at every frame.\n');
fprintf(fid,'Rigid far-cap pose correction; no scale correction. Positive-definite quadratic forms and centers interpolate linearly.\n');
fprintf(fid,'Initial node-specific surface levels are preserved; normal damping is relative to moving surface speed.\n');
fprintf(fid,'No node target positions, trajectory forcing, or observed rupture time. Surface geometry IS prescribed external forcing.\n');
fprintf(fid,'Fixed prior network/weak patch; ten global controls reoptimized using training position MSE only.\n');
fprintf(fid,'Multipliers are relative to the previous resistance fit. Local bounded search; no uniqueness claim.\n');
fprintf(fid,'Validation is conditional on geometry at all frames, not independent forecasting or biological validation.\n\n');
fprintf(fid,'Model                       Training RMSE_um  Conditional holdout RMSE_um\n');
fprintf(fid,'Prior fixed surface         %.6f  %.6f\n',previous.metrics.trainRMSE_um,previous.metrics.validationRMSE_um);
fprintf(fid,'Moving surface before refit %.6f  %.6f\n',incumbentMetrics.trainRMSE_um,incumbentMetrics.validationRMSE_um);
fprintf(fid,'Moving surface after refit  %.6f  %.6f\n',bestMetrics.trainRMSE_um,bestMetrics.validationRMSE_um);
fprintf(fid,'Predicted rupture: %.9g s (extension threshold, not measured label).\n',bestResult.param.ruptureTime);
fprintf(fid,'Search trials: %d; fine candidates: %d; wall time: %.1f s.\n',numel(history),numel(fine),toc(started));
fprintf(fid,'Initial/final ellipsoid volumes: %.9g / %.9g um^3 (%.3f%% change).\n', ...
    diagnostics.volume_um3(1),diagnostics.volume_um3(end),100*(diagnostics.volume_um3(end)/diagnostics.volume_um3(1)-1));
fprintf(fid,'Ellipsoid volume is a shape-fit proxy; an open retracting shell is not a closed volume measurement.\n');
fprintf('MOVING SURFACE COMPLETE: train %.5f -> %.5f um; conditional holdout %.5f -> %.5f um.\n', ...
    previous.metrics.trainRMSE_um,bestMetrics.trainRMSE_um,previous.metrics.validationRMSE_um,bestMetrics.validationRMSE_um);

    function p=parameters(values)
        p=tuning_parameters(base,values(1:7));
        p=tuning_resistance_parameters(p,values(8:10));
    end
    function evaluate(values)
        if numel(history)>=maxEvaluations || isKey(cache,key(values)), return; end
        timer=tic;
        try
            result=tuning_simulate(parameters(values),data.time,'coarse',90,false);
            metrics=tuning_score(result,data);
            objective=metrics.objective; rmse=metrics.trainRMSE_um; status='ok';
        catch exception
            if ~ismember(exception.identifier,{'follicle_rupture:tuningTimeout', ...
                    'follicle_rupture:incompleteIntegration','MATLAB:deval:SolOutsideInterval', ...
                    'follicle_rupture:tuningRamp'}), rethrow(exception); end
            objective=inf; rmse=inf; status=exception.identifier;
        end
        history(end+1)=struct('values',values,'objective',objective,'trainRMSE_um',rmse, ...
            'wallSeconds',toc(timer),'status',status);
        cache(key(values))=numel(history);
        if objective<bestObjective, bestObjective=objective; bestValues=values; end
        fprintf('Moving surface %d/%d: train %.5f um; best %.5f; %.1f s; %s\n', ...
            numel(history),maxEvaluations,rmse,sqrt(bestObjective),toc(timer),status);
        save(checkpoint,'configuration','history');
    end
    function z=encode(values)
        z=values; z(logarithmic)=log(z(logarithmic)); z=(z-low)./(high-low);
    end
    function values=decode(z)
        values=low+z.*(high-low); values(logarithmic)=exp(values(logarithmic));
    end
end
function value=key(values)
value=sprintf('%.12g,',values);
end
