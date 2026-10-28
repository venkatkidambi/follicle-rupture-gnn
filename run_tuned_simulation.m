function simFile = run_tuned_simulation(fitFile)
%RUN_TUNED_SIMULATION Run the fitted mechanics with no observed trajectories.
% Prefer the local network fit, then the parameter-only fit, then the CSV
% preset. Only parameters and initial geometry enter the forward solver.
repo=fileparts(mfilename('fullpath'));
addpath(fullfile(repo,'src'));
if nargin<1
    fitFile=fullfile(repo,'results','network_tuning','best_fit.mat');
    if ~isfile(fitFile), fitFile=fullfile(repo,'results','parameter_tuning','best_fit.mat'); end
    if ~isfile(fitFile), fitFile=fullfile(repo,'tuned_parameters.csv'); end
end
[~,~,extension]=fileparts(fitFile);
if strcmpi(extension,'.csv')
    preset=readtable(fitFile);
    expected={'contractionRateFrac';'prestressStartTime';'maxwellRelaxationTime'; ...
        'adhesionThresholdMultiplier';'adhesionReleaseTau';'extensionCritFrac';'activationRiseTime'};
    assert(ismember(height(preset),[6 7]) && ...
        isequal(cellstr(preset.Parameter),expected(1:height(preset))), ...
        'Unexpected parameter names/order in preset.');
    baselineFile=fullfile(repo,'results','follicleShellSim_cell_scale.mat');
    if ~isfile(baselineFile), follicle_rupture_msd_fixed('cell-scale'); end
    baseline=load(baselineFile,'param');
    assert(strcmp(baseline.param.mechanicsProfile,'cell-scale') && ...
        strcmp(baseline.param.units.coordinateInputMode,'voxels'), ...
        'The tuned preset requires the source-calibrated cell-scale geometry.');
    param=tuning_parameters(baseline.param,preset.FittedValue.');
else
    loaded=load(fitFile,'fit');
    param=loaded.fit.param;
    clear loaded;
end
time=unique([linspace(0,1800,181),linspace(0,1800,45)]).';
result=tuning_simulate(param,time,'fine',300,true);
simFile=fullfile(repo,'results','follicleShellSim_tuned.mat');
if isfield(param,'damageMode') && strcmp(param.damageMode,'fitted-static-connected-patch')
    simFile=fullfile(repo,'results','follicleShellSim_network_tuned.mat');
end
save(simFile,'-struct','result','-v7.3');
fprintf('Saved independently simulated fitted mechanics: %s\n',simFile);
end
