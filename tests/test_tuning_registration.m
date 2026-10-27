function tests=test_tuning_registration
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath=path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testOnlyRigidPoseIsRemoved(testCase)
reference=[0 0 0;1 0 0;0 1 0;0 0 1;1 1 1;2 1 0]*1e-5;
rotation=[cos(.3) -sin(.3) 0;sin(.3) cos(.3) 0;0 0 1];
observed=reference*rotation+[3 -2 4]*1e-6;
anchors=[true(4,1);false(2,1)];
observed(5,:)=observed(5,:)+[2 -1 3]*1e-6;
[aligned,R,~]=tuning_rigid_registration(observed,reference,anchors);
verifyEqual(testCase,aligned(anchors,:),reference(anchors,:),'AbsTol',1e-18);
verifyGreaterThan(testCase,norm(aligned(5,:)-reference(5,:)),3e-6);
verifyEqual(testCase,det(R),1,'AbsTol',1e-13);
verifyEqual(testCase,norm(aligned(5,:)-aligned(6,:)), ...
    norm(observed(5,:)-observed(6,:)),'AbsTol',1e-18);
end

function testHeldOutCellsCannotChangeRegistration(testCase)
reference=[0 0 0;1 0 0;0 1 0;0 0 1;1 1 1]*1e-5;
observed=reference+[2 3 -4]*1e-6;
mask=[true(4,1);false];
[~,R1,T1]=tuning_rigid_registration(observed,reference,mask);
observed(5,:)=[1 2 3];
[~,R2,T2]=tuning_rigid_registration(observed,reference,mask);
verifyEqual(testCase,R1,R2);
verifyEqual(testCase,T1,T2);
end

function testHeldOutTargetsCannotChangeTrainingScore(testCase)
data.position=zeros(6,3,5);
data.displacement=zeros(6,3,5);
data.rawDisplacement=data.displacement;
data.axis=[1;0;0];
data.trainNodes=[true(4,1);false(2,1)];
data.validationNodes=~data.trainNodes;
data.trainFrames=[false;true;false;true;false];
data.validationFrames=[false;false;true;false;true];
data.valid=true(6,5);
data.huberDelta_um=5;
result.Y=zeros(5,18);
before=tuning_score(result,data);
data.displacement(data.validationNodes,1,:)=10e-6;
data.displacement(:,1,data.validationFrames)=10e-6;
after=tuning_score(result,data);
verifyEqual(testCase,after.objective,before.objective);
verifyEqual(testCase,after.trainRMSE_um,0);
verifyEqual(testCase,after.validationRMSE_um,10,'AbsTol',1e-12);
end
