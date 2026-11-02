function simFile = run_tuned_simulation(fitFile)
%RUN_TUNED_SIMULATION Run fitted mechanics without per-node observed trajectories.
% Prefer the local network fit, then the parameter-only fit, then the CSV
% preset. Moving-surface fits also carry the prescribed global geometry
% history; no per-node measured trajectories enter the forward solver.
repo=fileparts(mfilename('fullpath'));
addpath(fullfile(repo,'src'));
if nargin<1
    fitFile=fullfile(repo,'results','network_tuning','best_fit.mat');
    if ~isfile(fitFile), fitFile=fullfile(repo,'results','parameter_tuning','best_fit.mat'); end
    if ~isfile(fitFile), fitFile=fullfile(repo,'tuned_parameters.csv'); end
end
[~,~,extension]=fileparts(fitFile);
isResistanceFit=false;
isMovingSurfaceFit=false;
isMovingNetworkFit=false;
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
    isResistanceFit=isfield(loaded.fit,'fitType') && strcmp(loaded.fit.fitType,'resistance');
    isMovingSurfaceFit=isfield(loaded.fit,'fitType') && strcmp(loaded.fit.fitType,'moving-surface');
    isMovingNetworkFit=isfield(loaded.fit,'fitType') && strcmp(loaded.fit.fitType,'moving-surface-network');
    clear loaded;
end
time=unique([linspace(0,1800,181),linspace(0,1800,45)]).';
result=tuning_simulate(param,time,'fine',300,true);
simFile=fullfile(repo,'results','follicleShellSim_tuned.mat');
if isfield(param,'damageMode') && strcmp(param.damageMode,'fitted-static-connected-patch')
    simFile=fullfile(repo,'results','follicleShellSim_network_tuned.mat');
end
if isResistanceFit
    simFile=fullfile(repo,'results','follicleShellSim_resistance_tuned.mat');
end
if isMovingSurfaceFit
    simFile=fullfile(repo,'results','follicleShellSim_moving_surface_tuned.mat');
end
if isMovingNetworkFit
    simFile=fullfile(repo,'results','follicleShellSim_moving_network_tuned.mat');
end
save(simFile,'-struct','result','-v7.3');
fprintf('Saved independently simulated fitted mechanics: %s\n',simFile);
end
