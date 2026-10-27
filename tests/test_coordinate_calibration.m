function tests = test_coordinate_calibration
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'src'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testSourceSpacingAndAnisotropy(testCase)
units = experimental_coordinate_units( ...
    fullfile('tracked_nodes','trackedNodes_nuc_48_2_OA7_1.mat'), baseUnits, 'source');
% A displacement of two raw units has different physical length along z.
steps_um = (2*eye(3) .* units.coordinateToMeters)*1e6;
verifyEqual(testCase, diag(steps_um), [0.929832008623325;0.929832008623325;2], ...
    'AbsTol',1e-14);
verifyEqual(testCase, units.coordinateInputMode, 'voxels');
verifyEqual(testCase, units.schemaVersion, 2);
verifyEqual(testCase, units.displayFromMeters, 1e6);
verifyTrue(testCase, contains(units.calibrationSource, 'Step14_findingAxis_v2.m'));
end

function testCalibrationDoesNotLeakToOtherFiles(testCase)
for name = {'different_tracks.mat','nodeDataSample1.xlsx'}
    units = experimental_coordinate_units(name{1}, baseUnits, 'source');
    verifyEqual(testCase, units.coordinateToMeters, 1e-6);
    verifyEqual(testCase, units.coordinateInputMode, 'microns');
end
end

function testLegacyCoordinateOptOut(testCase)
units = experimental_coordinate_units('trackedNodes_nuc_48_2_OA7_1.mat', ...
    baseUnits, 'input-units');
verifyEqual(testCase, units.coordinateToMeters, 1e-6);
verifyEqual(testCase, units.coordinateInputMode, 'microns');
end

function units = baseUnits
units.schemaVersion = 1;
units.coordinateInputMode = 'microns';
units.coordinateToMeters = 1e-6;
units.displayFromMeters = 1e6;
end
