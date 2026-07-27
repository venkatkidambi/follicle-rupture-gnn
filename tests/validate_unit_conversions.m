% Lightweight validation for the centralized follicle coordinate unit system.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'src'));

tol = 100 * eps;

units10 = follicle_unit_system(10, 'cells', 'microns');
metersFrom50Cells10 = 50 * units10.coordinateToMeters;
assert(abs(metersFrom50Cells10 - 500e-6) < tol, ...
    '50 cells at 10 microns/cell must equal 500 microns.');

units12 = follicle_unit_system(12, 'cells', 'microns');
metersFrom50Cells12 = 50 * units12.coordinateToMeters;
assert(abs(metersFrom50Cells12 - 600e-6) < tol, ...
    '50 cells at 12 microns/cell must equal 600 microns.');

micronUnits = follicle_unit_system(10, 'microns', 'microns');
meterUnits = follicle_unit_system(10, 'meters', 'microns');
coordsFromMicrons = [500 150 150] * micronUnits.coordinateToMeters;
coordsFromMeters = [500e-6 150e-6 150e-6] * meterUnits.coordinateToMeters;
assert(max(abs(coordsFromMicrons - coordsFromMeters)) < tol, ...
    'Micron and meter coordinate modes must produce equivalent SI coordinates.');

cellDisplayUnits = follicle_unit_system(10, 'meters', 'cells');
micronDisplayUnits = follicle_unit_system(10, 'meters', 'microns');
assert(abs(500e-6 * cellDisplayUnits.displayFromMeters - 50) < tol, ...
    'Cell display conversion must recover 50 cell diameters.');
assert(abs(500e-6 * micronDisplayUnits.displayFromMeters - 500) < tol, ...
    'Micron display conversion must recover 500 microns.');

legacyUnits = follicle_unit_system(struct());
assert(strcmp(legacyUnits.coordinateInputMode, 'meters'), ...
    'Legacy files must default to meter coordinate interpretation.');
assert(strcmp(legacyUnits.displayMode, 'microns'), ...
    'Legacy files must default to micron plotting.');

assertUnitError(@() follicle_unit_system(0, 'cells', 'microns'));
assertUnitError(@() follicle_unit_system(NaN, 'cells', 'microns'));
assertUnitError(@() follicle_unit_system(Inf, 'cells', 'microns'));
assertUnitError(@() follicle_unit_system(10, 'millimeters', 'microns'));
assertUnitError(@() follicle_unit_system(10, 'cells', 'meters'));

fprintf('[UNIT TEST] All follicle unit conversions passed.\n');


function assertUnitError(testFcn)
didError = false;
try
    testFcn();
catch
    didError = true;
end
assert(didError, 'Expected invalid unit configuration to raise an error.');
end
