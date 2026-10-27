function units = experimental_coordinate_units(dataFile, units, calibrationMode)
%EXPERIMENTAL_COORDINATE_UNITS Resolve dataset-specific raw xyz spacing.
% Step14_findingAxis_v2.m loads this exact tracked-node filename, sets the
% meter-per-unit spacings on lines 33-37, and applies them on lines 99-101.
% Its z expression ScalingX*zscale is ScalingZ, not an additional z scaling.
% Other files retain coordinateInputMode; never apply this acquisition's
% calibration to a different spreadsheet or tracked dataset automatically.
units.calibrationSource = 'coordinateInputMode configuration';
if strcmp(calibrationMode, 'input-units')
    return;
end
if ~strcmp(calibrationMode, 'source')
    error('follicle_rupture:invalidCalibration', ...
        'Calibration must be source or input-units.');
end
[~, name, extension] = fileparts(dataFile);
if strcmp([name extension], 'trackedNodes_nuc_48_2_OA7_1.mat')
    units.schemaVersion = 2;
    units.coordinateInputMode = 'voxels';
    units.coordinateToMeters = [2.3245800215583136e-7 / 0.5, ...
        2.3245800215583136e-7 / 0.5, 9.9999999999999995e-7];
    units.calibrationSource = 'Step14_findingAxis_v2.m:33-37,99-101';
end
end
