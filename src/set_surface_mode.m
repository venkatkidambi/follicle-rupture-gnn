function param = set_surface_mode(param, mode)
%SET_SURFACE_MODE Select saved geometry, measured motion, or its initial shape.
% Fixed mode preserves the exact first-frame center, orientation, quadratic
% form, and node reference levels. Later measured frames are removed entirely.
mode = validatestring(mode, {'saved','fixed','moving'}, mfilename, 'mode');
if strcmp(mode,'saved'), return; end
if strcmp(mode,'moving')
    assert(isfield(param,'movingSurface'), 'follicle_rupture:missingMovingSurface', ...
        'Moving mode requires a fit containing measured surface history.');
    if isfield(param,'fixedSurface'), param=rmfield(param,'fixedSurface'); end
elseif isfield(param,'movingSurface')
    surface=param.movingSurface;
    param.fixedSurface=struct('center',surface.center(1,:), ...
        'quadric',surface.quadric(:,:,1));
    param=rmfield(param,'movingSurface');
end
end
