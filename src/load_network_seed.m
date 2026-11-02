function previous = load_network_seed(sourceFile)
%LOAD_NETWORK_SEED Recover the original reference for all ten fitted controls.
loaded=load(sourceFile,'fit'); previous=loaded.fit;
if isfield(previous,'fitType') && strcmp(previous.fitType,'resistance')
    % The resistance stage freezes seven network-fit controls. Recover their
    % pre-fit reference rather than applying their multipliers twice.
    networkFile=previous.sourceFitFile;
    if ~isfile(networkFile)
        networkFile=fullfile(fileparts(mfilename('fullpath')),'..','results', ...
            'archive','01_fixed_surface_network','best_fit.mat');
    end
    checkpoint=fullfile(fileparts(networkFile),'search_history.mat');
    network=load(networkFile,'fit'); original=load(checkpoint,'configuration');
    assert(isequal(network.fit.values,previous.fixedParameterValues), ...
        'follicle_rupture:seedReference','Resistance fit and original network controls differ.');
    previous.frozenBase=original.configuration.base;
    previous.parameterNames=[reshape(previous.fixedParameterNames,1,[]), ...
        reshape(previous.parameterNames,1,[])];
    previous.values=[reshape(previous.fixedParameterValues,1,[]),reshape(previous.values,1,[])];
    previous.lowerBounds=[reshape(network.fit.lowerBounds,1,[]),reshape(previous.lowerBounds,1,[])];
    previous.upperBounds=[reshape(network.fit.upperBounds,1,[]),reshape(previous.upperBounds,1,[])];
    previous.multiplierReferenceFile=checkpoint;
end
assert(isfield(previous,'frozenBase') && numel(previous.values)==10, ...
    'follicle_rupture:invalidSurfaceSeed', ...
    'The source must supply all ten controls and their original coefficient reference.');
end
