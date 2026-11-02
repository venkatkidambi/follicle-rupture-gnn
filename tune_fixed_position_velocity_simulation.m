function fit = tune_fixed_position_velocity_simulation(maxEvaluations, resumeSearch, sourceFitFile)
%TUNE_FIXED_POSITION_VELOCITY_SIMULATION Full fixed-surface position/velocity fit.
% Reopen graph density/anisotropy, local edges, weak-patch location/size, and
% all ten global controls. Each loss uses training-only observed motion scales.
if nargin<1, maxEvaluations=240; end
if nargin<2, resumeSearch=false; end
if nargin<3
    repo=fileparts(mfilename('fullpath'));
    sourceFitFile=fullfile(repo,'results','fixed_surface_network_tuning','best_fit.mat');
    if ~isfile(sourceFitFile)
        sourceFitFile=fullfile(repo,'results','archive','02_fixed_surface_resistance','best_fit.mat');
    end
end
assert(isfile(sourceFitFile),'follicle_rupture:missingSurfaceFit', ...
    'Provide the completed fixed-surface position-only fit: %s',sourceFitFile);
source=load(sourceFitFile,'fit');
assert(~isfield(source.fit.param,'movingSurface'), ...
    'follicle_rupture:fixedSurfaceRequired','The source fit must already use a fixed surface.');
fit=tune_network_simulation(maxEvaluations,resumeSearch,'fixed',sourceFitFile,'position-velocity');
end
