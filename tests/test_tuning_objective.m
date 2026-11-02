function tests=test_tuning_objective
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
testCase.TestData.originalPath=path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
end
function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end
function testExactAndStationaryPredictions(testCase)
data=fixture;
data=configure_tuning_objective(data,'position-velocity');
result=prediction(data.displacement);
metrics=tuning_score(result,data);
verifyEqual(testCase,metrics.objective,0,'AbsTol',1e-25);
verifyEqual(testCase,metrics.trainVelocityRMSE_um_s,0,'AbsTol',1e-25);
stationary=prediction(zeros(size(data.displacement)));
metrics=tuning_score(stationary,data);
verifyEqual(testCase,metrics.trainPositionTerm,1,'AbsTol',1e-14);
verifyEqual(testCase,metrics.trainVelocityTerm,1,'AbsTol',1e-14);
verifyEqual(testCase,metrics.objective,2,'AbsTol',1e-14);
end
function testHoldoutCannotChangeTrainingObjectiveOrScales(testCase)
data=fixture; original=configure_tuning_objective(data,'position-velocity');
result=prediction(.8*data.displacement);
reference=tuning_score(result,original);
data.displacement(:,:,data.validationFrames)=data.displacement(:,:,data.validationFrames)+100;
data.displacement(data.validationNodes,:,:)=data.displacement(data.validationNodes,:,:)-200;
changed=configure_tuning_objective(data,'position-velocity');
verifyEqual(testCase,changed.objective,original.objective);
metrics=tuning_score(result,changed);
verifyEqual(testCase,metrics.objective,reference.objective);
verifyEqual(testCase,metrics.trainVelocityRMSE_um_s,reference.trainVelocityRMSE_um_s);
verifyGreaterThan(testCase,metrics.validationRMSE_um,reference.validationRMSE_um);
end
function testIrregularIntervalsAndConstantOffset(testCase)
data=configure_tuning_objective(fixture,'position-velocity');
offset=2e-6;
values=data.displacement; values(:,1,:)=values(:,1,:)+offset;
metrics=tuning_score(prediction(values),data);
verifyEqual(testCase,metrics.trainRMSE_um,2,'AbsTol',1e-12);
verifyEqual(testCase,metrics.trainVelocityRMSE_um_s,0,'AbsTol',1e-14);
verifyEqual(testCase,metrics.objective,4/data.objective.positionScale_um^2,'AbsTol',1e-12);
end
function testInvalidEndpointsExcludedAndLegacyUnchanged(testCase)
data=fixture; result=prediction(.8*data.displacement);
legacy=tuning_score(result,data);
position=configure_tuning_objective(data,'position');
verifyEqual(testCase,tuning_score(result,position),legacy);
data.valid(1,4)=false;
data.displacement(1,:,4)=NaN;
data=configure_tuning_objective(data,'position-velocity');
metrics=tuning_score(result,data);
verifyEqual(testCase,metrics.trainVelocitySamples,2);
verifyTrue(testCase,isfinite(metrics.objective));
end
function testVelocityDistinguishesEqualPositionErrors(testCase)
data=configure_tuning_objective(fixture,'position-velocity');
constant=data.displacement; changing=data.displacement;
constant(:,1,:)=constant(:,1,:)+1e-6;
changing(:,1,[2 6])=changing(:,1,[2 6])+1e-6;
changing(:,1,4)=changing(:,1,4)-1e-6;
a=tuning_score(prediction(constant),data);
b=tuning_score(prediction(changing),data);
verifyEqual(testCase,a.trainRMSE_um,b.trainRMSE_um,'AbsTol',1e-12);
verifyGreaterThan(testCase,b.trainVelocityRMSE_um_s,a.trainVelocityRMSE_um_s);
verifyGreaterThan(testCase,b.objective,a.objective);
end
function data=fixture
data.time=[0;1;2;4;7;11];
data.trainNodes=[true;true;false]; data.validationNodes=~data.trainNodes;
data.trainFrames=[false;true;false;true;false;true];
data.validationFrames=~data.trainFrames;
data.valid=true(3,6); data.huberDelta_um=5; data.axis=[1;0;0];
data.displacement=zeros(3,3,6);
for t=1:6
    data.displacement(:,:,t)=[1 2 3;2 -1 1;3 2 -1]*data.time(t)*1e-6;
end
data.position=data.displacement; data.rawDisplacement=data.displacement;
end
function result=prediction(values)
result.Y=reshape(permute(values,[2 1 3]),size(values,1)*3,[]).';
end
