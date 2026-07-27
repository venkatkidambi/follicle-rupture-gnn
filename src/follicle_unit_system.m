function units = follicle_unit_system(cellDiameter_um, coordinateInputMode, displayMode)
%FOLLICLE_UNIT_SYSTEM Validate and centralize simulation/display length units.
%
% New simulations call:
%   units = follicle_unit_system(cellDiameter_um, coordinateInputMode, displayMode)
%
% Plotters may call:
%   units = follicle_unit_system(param)
%
% The struct form resolves saved metadata and supplies backward-compatible
% defaults for older MAT files: coordinates are interpreted as meters
% internally and plots are displayed in microns.

if nargin == 1 && isstruct(cellDiameter_um)
    param = cellDiameter_um;
    if isfield(param, 'units') && isstruct(param.units)
        savedUnits = param.units;
        cellDiameter_um = getUnitField(savedUnits, 'cellDiameter_um', 10);
        coordinateInputMode = getUnitField(savedUnits, 'coordinateInputMode', 'meters');
        displayMode = getUnitField(savedUnits, 'displayMode', 'microns');
    else
        cellDiameter_um = 10;
        coordinateInputMode = 'meters';
        displayMode = 'microns';
    end
elseif nargin ~= 3
    error('follicle_rupture:invalidUnitConfiguration', ...
        ['Call follicle_unit_system with (cellDiameter_um, coordinateInputMode, ', ...
         'displayMode) or with a saved param struct.']);
end

if ~(isnumeric(cellDiameter_um) && isreal(cellDiameter_um) && ...
        isscalar(cellDiameter_um) && isfinite(cellDiameter_um) && cellDiameter_um > 0)
    error('follicle_rupture:invalidCellDiameter', ...
        'cellDiameter_um must be a positive, finite numeric scalar.');
end

coordinateInputMode = normalizeMode(coordinateInputMode, ...
    {'cells', 'microns', 'meters'}, 'coordinate input');
displayMode = normalizeMode(displayMode, ...
    {'cells', 'microns'}, 'display');

cellDiameter_m = cellDiameter_um * 1e-6;
switch coordinateInputMode
    case 'cells'
        coordinateToMeters = cellDiameter_m;
    case 'microns'
        coordinateToMeters = 1e-6;
    case 'meters'
        coordinateToMeters = 1;
end

switch displayMode
    case 'cells'
        displayFromMeters = 1 / cellDiameter_m;
        displayLengthUnit = 'cell diameters';
        displayVelocityUnit = 'cells/s';
    case 'microns'
        displayFromMeters = 1e6;
        displayLengthUnit = 'µm';
        displayVelocityUnit = 'µm/s';
end

units.schemaVersion = 1;
units.internalLengthUnit = 'm';
units.internalTimeUnit = 's';
units.internalForceUnit = 'N';
units.coordinateInputMode = coordinateInputMode;
units.displayMode = displayMode;
units.cellDiameter_um = cellDiameter_um;
units.cellDiameter_m = cellDiameter_m;
units.coordinateToMeters = coordinateToMeters;
units.displayFromMeters = displayFromMeters;
units.displayLengthUnit = displayLengthUnit;
units.displayVelocityUnit = displayVelocityUnit;

% Mechanical parameters remain SI and are not rescaled by coordinate input.
units.stiffnessUnit = 'N/m';
units.dampingUnit = 'N*s/m';
units.massUnit = 'kg';
units.strainUnit = 'dimensionless';
units.damageUnit = 'dimensionless';
end


function value = getUnitField(units, fieldName, defaultValue)
if isfield(units, fieldName) && ~isempty(units.(fieldName))
    value = units.(fieldName);
else
    value = defaultValue;
end
end


function mode = normalizeMode(mode, allowedModes, modeDescription)
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('follicle_rupture:invalidUnitMode', ...
        '%s mode must be a character vector or scalar string.', modeDescription);
end
mode = lower(strtrim(char(mode)));
if ~ismember(mode, allowedModes)
    error('follicle_rupture:invalidUnitMode', ...
        'Unknown %s mode "%s". Valid modes are: %s.', ...
        modeDescription, mode, strjoin(allowedModes, ', '));
end
end
