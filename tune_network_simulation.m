function fit = tune_network_simulation(maxEvaluations, resumeSearch, surfaceMode)
%TUNE_NETWORK_SIMULATION Fit static connections, weak patch, and global mechanics.
% Motion is predicted by force balance. Only training position error selects
% candidates; validation data never propose edges, patches, or coefficients.
% The original parameter-only fit is retained in results/parameter_tuning.
if nargin<1, maxEvaluations=160; end
if nargin<2, resumeSearch=false; end
if nargin<3, surfaceMode='fixed'; end
assert(ismember(surfaceMode,{'fixed','moving'}),'Use fixed or moving surface mode.');
movingSearch=strcmp(surfaceMode,'moving');
validateattributes(maxEvaluations,{'numeric'},{'scalar','integer','>=',8});
repo=fileparts(mfilename('fullpath')); addpath(fullfile(repo,'src'));
if movingSearch
    outputName='moving_network_tuning';
else
    outputName='network_tuning';
end
out=fullfile(repo,'results',outputName);
if ~isfolder(out), mkdir(out); end
if movingSearch
    sourceFile=fullfile(repo,'results','moving_surface_tuning','best_fit.mat');
    loaded=load(sourceFile,'fit'); previous=loaded.fit;
    % Retain the original multiplier reference so the incumbent is reproduced
    % exactly and changing the graph does not compound parameter multipliers.
    base=previous.frozenBase; data=previous.data;
    names=previous.parameterNames;
    lower=previous.lowerBounds; upper=previous.upperBounds;
    initialValues=previous.values;
    baselineMetrics=previous.baselineMetrics;
    stationaryMetrics=previous.stationaryMetrics;
    initialEdges=previous.param.edgeList;
    initialWeakEdges=find(previous.param.damageMask>0);
    scalarOrder=[1 2 7 4 8 9 10 5 6 3];
    trialSeconds=180;
else
baselineFile=fullfile(repo,'results','follicleShellSim_cell_scale.mat');
if ~isfile(baselineFile), follicle_rupture_msd_fixed('cell-scale'); end
loaded=load(baselineFile,'param','Y','tSol'); base=loaded.param;
assert(strcmp(base.mechanicsProfile,'cell-scale'),'Network search requires the cell-scale model.');
data=tuning_prepare_data(base,1800);
original.Y=interp1(loaded.tSol,loaded.Y,data.time,'pchip');
baselineMetrics=tuning_score(original,data);
stationary.Y=zeros(size(original.Y)); stationaryMetrics=tuning_score(stationary,data);
names={'contractionRateFrac','prestressStartTime','maxwellRelaxationTime', ...
    'adhesionThresholdMultiplier','adhesionReleaseTau','extensionCritFrac','activationRiseTime'};
lower=[2e-5 0 5 .15 2 .01 60]; upper=[4e-4 1200 300 8 180 .15 900];
preset=readtable(fullfile(repo,'tuned_parameters.csv'));
initialValues=preset.FittedValue.';
seedFile=fullfile(repo,'results','parameter_tuning','best_fit.mat');
if isfile(seedFile)
    seed=load(seedFile,'fit');
    assert(isequal(seed.fit.param.nodeCoords,base.nodeCoords) && ...
        isequal(seed.fit.data,data),'Parameter-only seed uses different geometry or observations.');
    initialValues=seed.fit.values;
end
if numel(initialValues)==6, initialValues(7)=base.contractionRampTime; end
initialEdges=base.edgeList; initialWeakEdges=find(base.damageMask>0);
scalarOrder=[1 2 7 4 5 6 3]; trialSeconds=60;
end
logarithmic=true(size(lower)); logarithmic(2)=false;
low=lower; high=upper; low(logarithmic)=log(low(logarithmic)); high(logarithmic)=log(high(logarithmic));
geometry=tuning_network_geometry(base);
initial=struct('values',initialValues,'edges',initialEdges, ...
    'weakEdges',initialWeakEdges,'centerNode',0,'patchSize',numel(initialWeakEdges), ...
    'proposal','original fitted network');
configuration=struct('version',1,'base',base,'data',data,'lower',lower,'upper',upper, ...
    'geometry',geometry,'initial',initial,'objective','training mean squared 3D position error');
if movingSearch, configuration.version=2; configuration.surfaceMode=surfaceMode; end
history=struct('candidate',{},'objective',{},'trainRMSE_um',{},'wallSeconds',{},'status',{});
checkpoint=fullfile(out,'search_history.mat');
if resumeSearch && isfile(checkpoint)
    saved=load(checkpoint);
    assert(isequaln(saved.configuration,configuration),'Checkpoint inputs changed; start a fresh search.');
    history=saved.history;
end
cache=containers.Map('KeyType','char','ValueType','double');
best=initial; bestObjective=inf; started=tic;
for idx=1:numel(history)
    if strcmp(history(idx).status,'ok')
        cache(candidateKey(history(idx).candidate))=idx;
        if history(idx).objective<bestObjective
            best=history(idx).candidate; bestObjective=history(idx).objective;
        end
    end
end
fprintf('Static-network fit: %d training cells x %d frames; %d joint holdout cells x %d frames.\n', ...
    nnz(data.trainNodes),nnz(data.trainFrames),nnz(data.validationNodes),nnz(data.validationFrames));
evaluate(initial);
% Keep an independent fine solve of the previous fit as an incumbent. A
% coarse-ranking artifact cannot replace it with a worse fine training fit.
incumbent=tuning_simulate(tuning_rebuild_network(base,initial.edges,initial.weakEdges,initial.values), ...
    data.time,'fine',300,false);
incumbentMetrics=tuning_score(incumbent,data);
if movingSearch
    assert(max(abs(incumbent.Y-previous.prediction.Y),[],'all')<1e-10, ...
        'Rebuilding the incumbent changed its prediction; check multiplier references.');
end
stream=RandStream('mt19937ar','Seed',72419);
% Alternate topology, patch, scalar, and local edge searches. Budgets reserve
% work for every family, instead of letting a broad first stage consume it.
% Resuming starts a new refinement round rather than repeatedly resetting
% to broad surface screening and the largest scalar step.
roundNumber=ceil(max(0,numel(history)-1)/80);
while numel(history)<maxEvaluations
    roundNumber=roundNumber+1; before=numel(history);
    remaining=maxEvaluations-numel(history);
    roundBudget=min(remaining,80);
    stopGraph=numel(history)+max(1,floor(.18*roundBudget));
    for count=[4 6 8 10]
        for weight=[.5 1 2]
            if numel(history)>=stopGraph, break; end
            trial=best;
            trial.edges=tuning_local_network(base,geometry,count,weight);
            trial=refreshPatch(trial);
            trial.proposal=sprintf('local graph k=%d axial weight=%g',count,weight);
            evaluate(trial);
        end
    end
    stopPatch=min(maxEvaluations,numel(history)+max(1,floor(.25*roundBudget)));
    % First broad surface coverage, then neighbors of the best patch center.
    centers=geometry.patchCenters;
    if roundNumber>1
        centerNode=best.centerNode;
        if centerNode==0, centerNode=nearestPatchCenter(best); end
        [~,order]=sort(sum((base.nodeCoords-base.nodeCoords(centerNode,:)).^2,2));
        centers=order(base.anchorNodeMask(order)==0).';
    end
    for node=centers
        if numel(history)>=stopPatch, break; end
        trial=best; trial.centerNode=node;
        trial.weakEdges=tuning_rupture_patch(base,trial.edges,node,trial.patchSize);
        trial.proposal=sprintf('patch center cell %d',base.experimentalCellIds(node));
        evaluate(trial);
    end
    for count=[1 3 8 16 32]
        if numel(history)>=maxEvaluations, break; end
        trial=best; trial.patchSize=count;
        if trial.centerNode==0, trial.centerNode=nearestPatchCenter(trial); end
        trial=refreshPatch(trial); trial.proposal=sprintf('patch size %d',count);
        evaluate(trial);
    end
    stopScalars=min(maxEvaluations,numel(history)+max(2,floor(.30*roundBudget)));
    for dimension=scalarOrder
        center=encode(best.values);
        for direction=[-1 1]
            if numel(history)>=stopScalars, break; end
            z=center; z(dimension)=min(max(z(dimension)+direction*.12/sqrt(roundNumber),0),1);
            trial=best; trial.values=decode(z); trial.proposal=['mechanics ' names{dimension}];
            evaluate(trial);
        end
    end
    stopEdges=min(maxEvaluations,numel(history)+max(2,roundBudget-(numel(history)-before)));
    attempts=0;
    while numel(history)<stopEdges && attempts<100
        attempts=attempts+1;
        trial=best;
        % Geometry-local changes to actual connection placement. Protect the
        % spanning tree and minimum degree 3; never use motion residuals to
        % prescribe which cells should move or assign per-cell actuation.
        missing=setdiff(geometry.pool,trial.edges,'rows');
        protected=ismember(trial.edges,geometry.backbone,'rows');
        if trial.centerNode==0, protected(trial.weakEdges)=true; end
        degree=accumarray(trial.edges(:),1,[base.N 1]);
        removable=find(~protected & all(degree(trial.edges)>3,2));
        if isempty(missing) || isempty(removable), break; end
        addition=missing(randi(stream,size(missing,1)),:);
        local=removable(any(ismember(trial.edges(removable,:),addition),2));
        if isempty(local), continue; end
        remove=local(randi(stream,numel(local)));
        oldWeakPairs=trial.edges(trial.weakEdges,:);
        trial.edges(remove,:)=[];
        trial.edges=sortrows([trial.edges;addition]);
        if trial.centerNode==0
            [~,trial.weakEdges]=ismember(oldWeakPairs,trial.edges,'rows');
        else
            trial=refreshPatch(trial);
        end
        trial.proposal='local edge swap'; evaluate(trial);
    end
    if numel(history)==before, break; end
end
assert(isfinite(bestObjective),'No candidate completed; inspect the checkpoint.');
% Final selection uses fine training error, including the previous fit.
bestResult=incumbent; bestMetrics=incumbentMetrics; selected=initial;
fine=struct('values',initial.values,'objective',incumbentMetrics.objective, ...
    'metrics',incumbentMetrics,'wallSeconds',incumbent.wallSeconds,'candidate',initial);
[~,rank]=sort([history.objective]); rank=rank(isfinite([history(rank).objective]));
for idx=1:min(5,numel(rank))
    trial=history(rank(idx)).candidate;
    if strcmp(candidateKey(trial),candidateKey(initial)), continue; end
    try
        result=tuning_simulate(tuning_rebuild_network(base,trial.edges,trial.weakEdges,trial.values), ...
            data.time,'fine',300,false);
    catch exception
        if isCensored(exception), fprintf('Fine candidate skipped: %s\n',exception.identifier); continue; end
        rethrow(exception);
    end
    metrics=tuning_score(result,data);
    fine(end+1)=struct('values',trial.values,'objective',metrics.objective,'metrics',metrics, ...
        'wallSeconds',result.wallSeconds,'candidate',trial); %#ok<AGROW>
    if metrics.objective<bestMetrics.objective
        selected=trial; bestMetrics=metrics; bestResult=result;
    end
end
fit.parameterNames=names; fit.values=selected.values;
fit.lowerBounds=lower; fit.upperBounds=upper;
fit.param=tuning_rebuild_network(base,selected.edges,selected.weakEdges,selected.values);
fit.metrics=bestMetrics; fit.baselineMetrics=baselineMetrics; fit.stationaryMetrics=stationaryMetrics;
fit.data=data; fit.prediction=bestResult;
if movingSearch
    fit.fitType='moving-surface-network';
    fit.sourceFitFile=sourceFile;
    fit.multiplierReferenceFile=previous.sourceFitFile;
    fit.frozenBase=base; % coefficient reference, not a frozen graph
    fit.previousMovingSurfaceMetrics=incumbentMetrics;
    fit.fixedSurfaceMetrics=previous.fixedSurfaceMetrics;
    fit.movingSurfaceInitialMetrics=incumbentMetrics;
    fit.surfaceDiagnostics=previous.surfaceDiagnostics;
    fit.validationDescription=previous.validationDescription;
else
    fit.parameterOnlyMetrics=incumbentMetrics;
end
fit.history=history; fit.fineCandidates=fine; fit.structure=selected;
fit.geometry=geometry; fit.searchVersion=configuration.version;
fit.totalWallSeconds=toc(started);
fit.forwardSolveSeconds=sum([history.wallSeconds])+sum([fine.wallSeconds]);
fit.finalCoarseSolveBudgetSeconds=trialSeconds;
fit.forwardModel='Static fitted local network and weak patch; unchanged force law; no observation feedback';
if movingSearch
    fit.forwardModel='Measured moving ellipsoid; fitted local graph, rupture patch, and all ten global controls';
    % Audit search coverage, including trials that retained incumbent values.
    proposals={history.candidate};
    fit.searchCoverage.proposals=cellfun(@(x) x.proposal,proposals,'UniformOutput',false);
    fit.searchCoverage.globalControls=names;
    fit.searchCoverage.evaluatedControls=cellfun(@(name) ...
        any(strcmp(fit.searchCoverage.proposals,['mechanics ' name])),names);
    assert(all(fit.searchCoverage.evaluatedControls), ...
        'Search budget did not cover every global control.');
    fit.searchCoverage.connectionsChanged=~isequal(selected.edges,initial.edges);
    fit.searchCoverage.weakPairsChanged=~isequal(selected.edges(selected.weakEdges,:), ...
        initial.edges(initial.weakEdges,:));
end
fit.fixedForceScale=mean(base.k1_edge);
save(fullfile(out,'best_fit.mat'),'fit','-v7.3');
writetable(table(names.',selected.values.',lower.',upper.', ...
    'VariableNames',{'Parameter','FittedValue','LowerBound','UpperBound'}),fullfile(out,'parameters.csv'));
ids=base.experimentalCellIds;
writetable(table(ids(selected.edges(:,1)),ids(selected.edges(:,2)),fit.param.L0_edge0*1e6, ...
    logical(fit.param.damageMask),'VariableNames',{'CellID1','CellID2','InitialLength_um','WeakLink'}), ...
    fullfile(out,'connections.csv'));
fprintf('NETWORK FIT COMPLETE: training %.4f -> %.4f um; joint holdout %.4f -> %.4f um.\n', ...
    incumbentMetrics.trainRMSE_um,bestMetrics.trainRMSE_um,incumbentMetrics.validationRMSE_um,bestMetrics.validationRMSE_um);

    function evaluate(trial)
        if numel(history)>=maxEvaluations, return; end
        key=candidateKey(trial);
        if isKey(cache,key), return; end
        timer=tic;
        try
            param=tuning_rebuild_network(base,trial.edges,trial.weakEdges,trial.values);
            result=tuning_simulate(param,data.time,'coarse',trialSeconds,false);
            metrics=tuning_score(result,data);
            objective=metrics.objective; rmse=metrics.trainRMSE_um; status='ok';
        catch exception
            if ~isCensored(exception) && ~ismember(exception.identifier, ...
                    {'follicle_rupture:invalidNetwork','follicle_rupture:tuningRamp'})
                rethrow(exception);
            end
            objective=inf; rmse=inf; status=exception.identifier;
        end
        history(end+1)=struct('candidate',trial,'objective',objective,'trainRMSE_um',rmse, ...
            'wallSeconds',toc(timer),'status',status);
        if strcmp(status,'ok'), cache(key)=numel(history); end
        if objective<bestObjective, best=trial; bestObjective=objective; end
        fprintf('Network trial %d/%d: RMSE %.4f um, %.1f s, best %.4f | %s | %s\n', ...
            numel(history),maxEvaluations,rmse,history(end).wallSeconds,sqrt(bestObjective),trial.proposal,status);
        save(checkpoint,'configuration','history');
    end
    function trial=refreshPatch(trial)
        if trial.centerNode==0, trial.centerNode=nearestPatchCenter(trial); end
        trial.weakEdges=tuning_rupture_patch(base,trial.edges,trial.centerNode,trial.patchSize);
    end
    function node=nearestPatchCenter(~)
        center=mean(base.nodeCoords(base.rupturePatchNodes,:),1);
        distance=sum((base.nodeCoords-center).^2,2); distance(base.anchorNodeMask>0)=inf;
        [~,node]=min(distance);
    end
    function z=encode(values)
        z=values; z(logarithmic)=log(z(logarithmic)); z=(z-low)./(high-low);
    end
    function values=decode(z)
        values=low+z.*(high-low); values(logarithmic)=exp(values(logarithmic));
    end
end

function key=candidateKey(candidate)
key=[sprintf('%.12g,',candidate.values) '|' sprintf('%d,',candidate.edges.') ...
    '|' sprintf('%d,',sort(candidate.weakEdges))];
end

function answer=isCensored(exception)
answer=ismember(exception.identifier,{'follicle_rupture:tuningTimeout', ...
    'follicle_rupture:incompleteIntegration','MATLAB:deval:SolOutsideInterval'});
end
