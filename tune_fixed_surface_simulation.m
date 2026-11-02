function fit = tune_fixed_surface_simulation(maxEvaluations, resumeSearch, sourceFitFile)
%TUNE_FIXED_SURFACE_SIMULATION Refit the latest network with expansion disabled.
% Freeze the exact first measured ellipsoid; optimize the graph, weak patch,
% and all ten mechanical/resistance controls. Save separately from prior fits.
if nargin<1, maxEvaluations=160; end
if nargin<2, resumeSearch=false; end
validateattributes(maxEvaluations,{'numeric'},{'scalar','integer','>=',80});
if nargin<3
    repo=fileparts(mfilename('fullpath'));
    sourceFitFile=fullfile(repo,'results','moving_network_tuning','best_fit.mat');
    if ~isfile(sourceFitFile)
        sourceFitFile=fullfile(repo,'results','moving_surface_tuning','best_fit.mat');
    end
end
assert(isfile(sourceFitFile),'follicle_rupture:missingSurfaceFit', ...
    'A saved moving-surface fit is required: %s',sourceFitFile);
fit=tune_network_simulation(maxEvaluations,resumeSearch,'fixed',sourceFitFile);
end
