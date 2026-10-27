function fit = tune_simulation(maxEvaluations, resumeSearch)
%TUNE_SIMULATION Bounded, deterministic fit of seven global model parameters.
% All measured displacements are used exclusively in tuning_score. Forward
% predictions use the original force law and fixed initial geometry/mesh.
% No Optimization Toolbox, parallel pool, or videos are needed during search.
if nargin<1, maxEvaluations=48; end
if nargin<2, resumeSearch=false; end
validateattributes(maxEvaluations,{'numeric'},{'scalar','integer','positive'});
repo=fileparts(mfilename('fullpath'));
addpath(fullfile(repo,'src'));
out=fullfile(repo,'results','parameter_tuning');
if ~isfolder(out), mkdir(out); end
baselineFile=fullfile(repo,'results','follicleShellSim_cell_scale.mat');
if ~isfile(baselineFile)
    follicle_rupture_msd_fixed('cell-scale');
end
baseline=load(baselineFile,'param','Y','tSol');
base=baseline.param;
assert(strcmp(base.mechanicsProfile,'cell-scale') && strcmp(base.units.coordinateInputMode,'voxels'), ...
    'Regenerate the calibrated cell-scale baseline before tuning.');
data=tuning_prepare_data(base,1800);
baselinePrediction.Y=interp1(baseline.tSol,baseline.Y,data.time,'pchip');
baselineMetrics=tuning_score(baselinePrediction,data);
stationary.Y=zeros(size(baselinePrediction.Y));
stationaryMetrics=tuning_score(stationary,data);
names={'contractionRateFrac','prestressStartTime','maxwellRelaxationTime', ...
    'adhesionThresholdMultiplier','adhesionReleaseTau','extensionCritFrac','activationRiseTime'};
lower=[2e-5,0,5,0.15,2,0.01,60];
upper=[4e-4,1200,300,8,180,0.15,900];
logarithmic=[true false true true true true true];
low=lower; high=upper;
low(logarithmic)=log(lower(logarithmic));
high(logarithmic)=log(upper(logarithmic));
configuration=struct('base',base,'data',data,'lowerBounds',lower,'upperBounds',upper, ...
    'objective','mean squared 3D position error in um^2');
history=struct('values',{},'objective',{},'trainRMSE_um',{},'wallSeconds',{},'status',{});
cache=containers.Map('KeyType','char','ValueType','double');
bestObjective=inf; bestValues=[];
started=tic;
% Reuse completed solves after interruption or when extending the budget.
% Timeouts are computationally censored, not evidence of poor fit; retry
% them with the normal budget before drawing conclusions about parameters.
checkpointFile=fullfile(out,'search_history.mat');
if resumeSearch && isfile(checkpointFile)
    saved=load(checkpointFile);
    % A six-parameter checkpoint fixed the rise at the baseline's 900 s.
    % Embed those same predictions in the expanded search without new solves.
    if isfield(saved,'configuration') && numel(saved.configuration.lowerBounds)==6
        previous=saved.configuration;
        assert(isequaln(previous.base,base) && isequaln(previous.data,data) && ...
            isequal(previous.lowerBounds,lower(1:6)) && isequal(previous.upperBounds,upper(1:6)) && ...
            strcmp(previous.objective,configuration.objective), ...
            'Cannot extend a checkpoint with different geometry, targets, or bounds.');
        for row=1:numel(saved.history)
            saved.history(row).values(7)=base.contractionRampTime;
        end
        saved.configuration=configuration;
    end
    assert(isfield(saved,'configuration') && isequaln(saved.configuration,configuration), ...
        'Checkpoint geometry, observations, objective, or bounds changed; start a fresh search.');
    history=saved.history;
    % RMSE retains all information needed to rescore older screening runs
    % by squared error without repeating their forward integrations.
    for row=1:numel(history)
        history(row).objective=history(row).trainRMSE_um^2;
    end
    for row=1:numel(history)
        if strcmp(history(row).status,'ok')
            cache(sprintf('%.12g,',history(row).values))=row;
            if history(row).objective<bestObjective
                bestObjective=history(row).objective;
                bestValues=history(row).values;
            end
        end
    end
    retry=history(~strcmp({history.status},'ok'));
    for row=1:numel(retry)
        if numel(history)>=maxEvaluations, break; end
        evaluate(retry(row).values);
    end
end
fprintf('Training: %d cells x %d frames. Joint holdout: %d cells x %d frames.\n', ...
    nnz(data.trainNodes),nnz(data.trainFrames),nnz(data.validationNodes),nnz(data.validationFrames));
fprintf('Unfitted model training RMSE %.4f um; stationary %.4f um.\n', ...
    baselineMetrics.trainRMSE_um,stationaryMetrics.trainRMSE_um);
% Sparse screening of amplitude/onset before spending evaluations on all
% dimensions. These are parameter guesses, never trajectory forcing.
seeds=[ ...
    6e-5 3 50 1 5 .05; ...
    1e-4 3 50 1 5 .05; ...
    1.5e-4 450 50 1 5 .05; ...
    2e-4 750 50 1 5 .05; ...
    3e-4 1050 50 1 5 .05; ...
    1e-4 300 100 .5 30 .03; ...
    1.5e-4 600 100 .5 30 .03; ...
    2e-4 900 100 .5 30 .03; ...
    1e-4 300 20 2 60 .08; ...
    2e-4 600 20 2 60 .08; ...
    2e-4 0 50 1 5 .15; ...
    1.5e-4 0 50 .5 30 .12; ...
    2e-4 300 50 .5 30 .12; ...
    3e-4 600 50 .5 30 .12];
seeds(:,7)=base.contractionRampTime;
for row=1:size(seeds,1)
    if numel(history)>=maxEvaluations, break; end
    evaluate(seeds(row,:));
end
if isempty(bestValues), error('No screening candidate completed. Inspect search_history.mat.'); end
% Resolve late loading with the existing global ramp, without fitting a
% spatial activation schedule or assigning any cell's displacement.
timingSeeds=repmat(bestValues,8,1);
timingSeeds(:,[1 2 7])=[3e-4 1000 300;4e-4 1100 300;3e-4 1100 100; ...
    2e-4 900 300;4e-4 1200 100;3e-4 900 600;2e-4 1000 100;4e-4 1000 600];
for row=1:size(timingSeeds,1)
    if numel(history)>=maxEvaluations, break; end
    evaluate(timingSeeds(row,:));
end
% Bounded coordinate pattern search; cache avoids repeated solves at bounds.
steps=[.18 .09 .045 .09];
dimensions={[1 2 7],[1 2 7],[1 2 7],[6 4 3 5]};
for stage=1:numel(steps)
    step=steps(stage);
    for dimension=dimensions{stage}
        center=encode(bestValues);
        for direction=[-1 1]
            if numel(history)>=maxEvaluations, break; end
            trial=center;
            trial(dimension)=min(max(trial(dimension)+direction*step,0),1);
            evaluate(decode(trial));
        end
    end
end

% Rank by training loss only. Validation observations never select a model.
[~,rank]=sort([history.objective]);
rank=rank(isfinite([history(rank).objective]));
fine=struct('values',{},'objective',{},'metrics',{},'wallSeconds',{});
bestFineObjective=inf;
for idx=1:min(3,numel(rank))
    values=history(rank(idx)).values;
    param=tuning_parameters(base,values);
    result=tuning_simulate(param,data.time,'fine',180,false);
    metrics=tuning_score(result,data);
    fine(idx)=struct('values',values,'objective',metrics.objective,'metrics',metrics, ...
        'wallSeconds',result.wallSeconds);
    fprintf('Fine candidate %d: training RMSE %.4f um, objective %.6f\n',idx,metrics.trainRMSE_um,metrics.objective);
    if metrics.objective<bestFineObjective
        bestFineObjective=metrics.objective;
        bestValues=values; bestResult=result; bestMetrics=metrics;
    end
end
fit.parameterNames=names;
fit.values=bestValues;
fit.lowerBounds=lower; fit.upperBounds=upper;
fit.param=tuning_parameters(base,bestValues);
fit.metrics=bestMetrics;
fit.baselineMetrics=baselineMetrics;
fit.stationaryMetrics=stationaryMetrics;
fit.data=data;
fit.prediction=bestResult;
fit.history=history; fit.fineCandidates=fine;
fit.totalWallSeconds=toc(started);
fit.forwardSolveSeconds=sum([history.wallSeconds])+sum([fine.wallSeconds]);
fit.forwardModel='Unchanged overdamped mass-spring-damper force law; no observation feedback';
fit.fixedForceScale=mean(base.k1_edge);
save(fullfile(out,'best_fit.mat'),'fit','-v7.3');
parameterTable=table(names.',bestValues.',lower.',upper.', ...
    'VariableNames',{'Parameter','FittedValue','LowerBound','UpperBound'});
writetable(parameterTable,fullfile(out,'parameters.csv'));
fprintf('FIT COMPLETE: training RMSE %.4f um; held-out RMSE %.4f um; stationary holdout %.4f um.\n', ...
    bestMetrics.trainRMSE_um,bestMetrics.validationRMSE_um,stationaryMetrics.validationRMSE_um);
disp(parameterTable);
fprintf('%d coarse solves, %d fine solves, %.1f s total.\n',numel(history),numel(fine),fit.totalWallSeconds);

    function evaluate(values)
        key=sprintf('%.12g,',values);
        if isKey(cache,key), return; end
        trialStarted=tic;
        try
            param=tuning_parameters(base,values);
            result=tuning_simulate(param,data.time,'coarse',60,false);
            metrics=tuning_score(result,data);
            objective=metrics.objective; rmse=metrics.trainRMSE_um; status='ok';
        catch exception
            if ~ismember(exception.identifier,{'follicle_rupture:tuningTimeout', ...
                    'follicle_rupture:incompleteIntegration','MATLAB:deval:SolOutsideInterval'})
                rethrow(exception);
            end
            objective=inf; rmse=inf; status=exception.identifier;
        end
        entry=struct('values',values,'objective',objective,'trainRMSE_um',rmse, ...
            'wallSeconds',toc(trialStarted),'status',status);
        history(end+1)=entry;
        cache(key)=numel(history);
        if objective<bestObjective
            bestObjective=objective; bestValues=values;
        end
        fprintf('Trial %02d: loss %.5f, train RMSE %.3f um, %.2f s, best %.5f | %s\n', ...
            numel(history),objective,rmse,entry.wallSeconds,bestObjective,mat2str(values,4));
        save(checkpointFile,'history','bestValues','bestObjective','configuration');
    end

    function z=encode(values)
        z=values;
        z(logarithmic)=log(z(logarithmic));
        z=(z-low)./(high-low);
    end

    function values=decode(z)
        values=low+z.*(high-low);
        values(logarithmic)=exp(values(logarithmic));
    end
end
