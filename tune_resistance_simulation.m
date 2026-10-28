function fit = tune_resistance_simulation(maxEvaluations, resumeSearch)
%TUNE_RESISTANCE_SIMULATION Three-control extension of the frozen network fit.
% Fixed bounds/budget and training-only selection; holdout is reported last.
if nargin<1, maxEvaluations=64; end
if nargin<2, resumeSearch=false; end
validateattributes(maxEvaluations,{'numeric'},{'scalar','integer','>=',8});
repo=fileparts(mfilename('fullpath')); addpath(fullfile(repo,'src'));
sourceFile=fullfile(repo,'results','network_tuning','best_fit.mat');
loaded=load(sourceFile,'fit'); previous=loaded.fit;
base=previous.param; data=previous.data;
out=fullfile(repo,'results','resistance_tuning');
if ~isfolder(out), mkdir(out); end
names={'bodyTipAdhesionRatio','adhesionStiffnessMultiplier','backgroundSlidingDragMultiplier'};
initial=[base.adhesionLoadThresholdBody/base.adhesionLoadThresholdTip,1,1];
lower=[1 .2 .2]; upper=[12 5 5];
configuration=struct('version',1,'base',base,'data',data,'lowerBounds',lower, ...
    'upperBounds',upper,'initial',initial,'selection','training 3D position MSE only');
checkpoint=fullfile(out,'search_history.mat');
history=struct('values',{},'objective',{},'trainRMSE_um',{},'wallSeconds',{},'status',{});
if resumeSearch && isfile(checkpoint)
    saved=load(checkpoint);
    assert(isequaln(saved.configuration,configuration),'Checkpoint inputs changed.');
    history=saved.history;
end
started=tic;
% Fine central differences in log parameters, measured on training positions.
% This diagnostic cannot change bounds or use held-out errors for selection.
sensitivityFile=fullfile(out,'sensitivity.mat');
if resumeSearch && isfile(sensitivityFile)
    saved=load(sensitivityFile);
    assert(isequaln(saved.configuration,configuration),'Sensitivity inputs changed.');
    sensitivity=saved.sensitivity; incumbent=saved.incumbent;
else
    incumbent=tuning_simulate(tuning_resistance_parameters(base,initial),data.time,'fine',300,false);
    incumbentMetrics=tuning_score(incumbent,data);
    mask=data.trainNodes & data.trainFrames.' & data.valid;
    h=log(1.2); J=zeros(3*nnz(mask),3);
    sensitivity.parameterNames=names;
    sensitivity.logStep=h;
    sensitivity.trainRMSE_um=zeros(3,2);
    sensitivity.positionChangeRMS_um=zeros(3,2);
    for dim=1:3
        vectors=cell(1,2);
        for side=1:2
            values=initial; values(dim)=values(dim)*exp((2*side-3)*h);
            result=tuning_simulate(tuning_resistance_parameters(base,values),data.time,'fine',300,false);
            metrics=tuning_score(result,data);
            vectors{side}=trainingVector(result.Y,mask);
            sensitivity.trainRMSE_um(dim,side)=metrics.trainRMSE_um;
            delta=reshape(vectors{side}-trainingVector(incumbent.Y,mask),3,[]);
            sensitivity.positionChangeRMS_um(dim,side)=sqrt(mean(sum(delta.^2,1)));
        end
        J(:,dim)=(vectors{2}-vectors{1})/(2*h);
        fprintf('Sensitivity %s: train RMSE [-,+] = %.5f, %.5f um; motion change %.5f, %.5f um\n', ...
            names{dim},sensitivity.trainRMSE_um(dim,:),sensitivity.positionChangeRMS_um(dim,:));
    end
    sensitivity.derivativeRMS_um=sqrt(sum(J.^2,1)/nnz(mask));
    normalized=J./max(vecnorm(J),eps);
    sensitivity.derivativeCosine=normalized.'*normalized;
    sensitivity.singularValues=svd(J,0);
    sensitivity.incumbentTrainRMSE_um=incumbentMetrics.trainRMSE_um;
    save(sensitivityFile,'sensitivity','incumbent','configuration');
    writetable(table(names.',sensitivity.derivativeRMS_um.', ...
        sensitivity.trainRMSE_um(:,1),sensitivity.trainRMSE_um(:,2), ...
        'VariableNames',{'Parameter','PositionSensitivity_um_per_log_unit','TrainRMSE_divide1p2','TrainRMSE_multiply1p2'}), ...
        fullfile(out,'sensitivity.csv'));
end
incumbentMetrics=tuning_score(incumbent,data);
cache=containers.Map('KeyType','char','ValueType','double');
bestValues=initial; bestObjective=inf;
for idx=1:numel(history)
    cache(key(history(idx).values))=idx;
    if history(idx).objective<bestObjective
        bestObjective=history(idx).objective; bestValues=history(idx).values;
    end
end
evaluate(initial);
% Broad deterministic log-space probes, followed by bounded pattern search.
for dim=1:3
    for factor=[.5 2]
        values=initial; values(dim)=min(max(initial(dim)*factor,lower(dim)),upper(dim));
        evaluate(values);
    end
end
for idx=1:12
    z=mod(idx*[sqrt(2) sqrt(3) sqrt(5)],1);
    evaluate(decode(z));
end
for step=[.2 .1 .05 .025]
    for sweep=1:4
        before=bestObjective; count=numel(history);
        for dim=1:3
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
% Keep the unchanged incumbent eligible; select fine solutions by training only.
bestResult=incumbent; bestMetrics=incumbentMetrics; selected=initial;
fine=struct('values',initial,'objective',incumbentMetrics.objective, ...
    'metrics',incumbentMetrics,'wallSeconds',incumbent.wallSeconds);
[~,rank]=sort([history.objective]); rank=rank(isfinite([history(rank).objective]));
for idx=1:min(5,numel(rank))
    values=history(rank(idx)).values;
    if strcmp(key(values),key(initial)), continue; end
    result=tuning_simulate(tuning_resistance_parameters(base,values),data.time,'fine',300,false);
    metrics=tuning_score(result,data);
    fine(end+1)=struct('values',values,'objective',metrics.objective,'metrics',metrics,'wallSeconds',result.wallSeconds); %#ok<AGROW>
    if metrics.objective<bestMetrics.objective
        selected=values; bestResult=result; bestMetrics=metrics;
    end
    fprintf('Fine candidate %d: train RMSE %.6f um\n',idx,metrics.trainRMSE_um);
end
fit=struct('fitType','resistance','parameterNames',{names},'values',selected, ...
    'lowerBounds',lower,'upperBounds',upper,'param',tuning_resistance_parameters(base,selected), ...
    'data',data,'prediction',bestResult,'metrics',bestMetrics,'history',history, ...
    'fineCandidates',fine,'baselineMetrics',previous.baselineMetrics, ...
    'stationaryMetrics',previous.stationaryMetrics,'networkMetrics',incumbentMetrics, ...
    'sensitivity',sensitivity,'sourceFitFile',sourceFile,'frozenBase',base, ...
    'fixedParameterNames',{previous.parameterNames},'fixedParameterValues',previous.values, ...
    'totalWallSeconds',toc(started),'forwardModel','Frozen network and seven coefficients; three resistance controls only');
allowed={'adhesionLoadThresholdBody','adhesionLoadThreshold_node','adhesionStiffness', ...
    'c_drag','adhesionDragSliding','ruptureTime'};
for field=setdiff(fieldnames(base),allowed).'
    assert(isequaln(fit.param.(field{1}),base.(field{1})),'Changed frozen field: %s',field{1});
end
save(fullfile(out,'best_fit.mat'),'fit','-v7.3');
writetable(table(names.',selected.',lower.',upper.', ...
    'VariableNames',{'Parameter','FittedValue','LowerBound','UpperBound'}),fullfile(out,'parameters.csv'));
fprintf('RESISTANCE FIT COMPLETE: training %.6f -> %.6f um; joint holdout %.6f -> %.6f um.\n', ...
    incumbentMetrics.trainRMSE_um,bestMetrics.trainRMSE_um,incumbentMetrics.validationRMSE_um,bestMetrics.validationRMSE_um);
disp(selected);

    function evaluate(values)
        if numel(history)>=maxEvaluations || isKey(cache,key(values)), return; end
        timer=tic;
        try
            result=tuning_simulate(tuning_resistance_parameters(base,values),data.time,'coarse',90,false);
            metrics=tuning_score(result,data);
            objective=metrics.objective; rmse=metrics.trainRMSE_um; status='ok';
        catch exception
            if ~ismember(exception.identifier,{'follicle_rupture:tuningTimeout', ...
                    'follicle_rupture:incompleteIntegration','MATLAB:deval:SolOutsideInterval'}), rethrow(exception); end
            objective=inf; rmse=inf; status=exception.identifier;
        end
        history(end+1)=struct('values',values,'objective',objective,'trainRMSE_um',rmse, ...
            'wallSeconds',toc(timer),'status',status);
        cache(key(values))=numel(history);
        if objective<bestObjective, bestObjective=objective; bestValues=values; end
        fprintf('Resistance %d/%d: train %.5f um; best %.5f; values %s; %s\n', ...
            numel(history),maxEvaluations,rmse,sqrt(bestObjective),mat2str(values,5),status);
        save(checkpoint,'history','configuration');
    end
    function z=encode(values)
        z=(log(values)-log(lower))./(log(upper)-log(lower));
    end
    function values=decode(z)
        values=exp(log(lower)+z.*(log(upper)-log(lower)));
    end
end
function value=key(values)
value=sprintf('%.12g,',values);
end
function vector=trainingVector(Y,mask)
positions=reshape(Y.',3,size(mask,1),[])*1e6;
positions=reshape(positions,3,[]);
vector=reshape(positions(:,mask(:)),[],1);
end
