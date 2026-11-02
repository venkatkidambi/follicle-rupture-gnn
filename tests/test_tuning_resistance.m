function tests = test_tuning_resistance
tests = functiontests(localfunctions);
end
function setupOnce(testCase)
testCase.TestData.oldPath=path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
end
function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end
function testIncumbentIsExactlyPreserved(testCase)
base=fixture;
actual=tuning_resistance_parameters(base,[6 1 1]);
verifyEqual(testCase,actual,rmfield(base,'ruptureTime'));
end
function testOnlyAuthorizedResistanceFieldsChange(testCase)
base=fixture;
actual=tuning_resistance_parameters(base,[2 .5 3]);
verifyEqual(testCase,actual.adhesionLoadThresholdTip,base.adhesionLoadThresholdTip);
verifyEqual(testCase,actual.adhesionLoadThreshold_node,[1;1.5;2]*1e-10,'RelTol',1e-14);
verifyEqual(testCase,actual.adhesionLoadThresholdBody,2e-10);
verifyEqual(testCase,actual.adhesionStiffness,.5e-3);
verifyEqual(testCase,actual.c_drag,3e-4,'RelTol',1e-14);
verifyEqual(testCase,actual.adhesionDragSliding,3e-5,'RelTol',1e-14);
changed={'adhesionLoadThresholdBody','adhesionLoadThreshold_node', ...
    'adhesionStiffness','c_drag','adhesionDragSliding','ruptureTime'};
for field=setdiff(fieldnames(base),changed).'
    verifyEqual(testCase,actual.(field{1}),base.(field{1}));
end
end
function testMultipliersAlwaysReferenceFrozenBase(testCase)
base=fixture;
a=tuning_resistance_parameters(base,[3 2 .4]);
b=tuning_resistance_parameters(base,[3 2 .4]);
verifyEqual(testCase,a,b);
verifyEqual(testCase,a.c_drag/a.adhesionDragSliding,base.c_drag/base.adhesionDragSliding);
end
function testCanRefitAfterUniformThresholdSelected(testCase)
base=fixture;
base.adhesionLoadThresholdBody=base.adhesionLoadThresholdTip;
base.adhesionLoadThreshold_node(:)=base.adhesionLoadThresholdTip;
base.principalNodeScore0=[-2;0;2]*1e-5;
base.anchorCutoff=2e-5; base.anchorNodeMask=[0;0;1];
actual=tuning_resistance_parameters(base,[3 1 1]);
verifyEqual(testCase,actual.adhesionLoadThreshold_node,[1;2;3]*1e-10,'RelTol',1e-14);
verifyEqual(testCase,actual.adhesionLoadThresholdTip,base.adhesionLoadThresholdTip);
end
function base=fixture
base=struct('adhesionLoadThresholdTip',1e-10,'adhesionLoadThresholdBody',6e-10, ...
    'adhesionLoadThreshold_node',[1;3.5;6]*1e-10,'adhesionStiffness',1e-3, ...
    'c_drag',1e-4,'adhesionDragSliding',1e-5,'adhesionDragStuck',5e-5, ...
    'ruptureTime',853,'edgeList',[1 2;2 3],'damageMask',[1;0], ...
    'contractionRateFrac',1.7e-4,'k1_edge',[.002;.002],'k_anchor',.01);
end
