function verification = verify_tuning(fitFile)
%VERIFY_TUNING Check convergence and force balance of the selected fit.
repo=fileparts(mfilename('fullpath'));
addpath(fullfile(repo,'src'));
if nargin<1, fitFile=fullfile(repo,'results','parameter_tuning','best_fit.mat'); end
loaded=load(fitFile,'fit'); fit=loaded.fit;
% Tenfold tighter error tolerances and half the maximum time step. The
% measured positions still do not enter this independent forward solve.
reference=tuning_simulate(fit.param,fit.data.time,'verification',300,true);
delta=reshape((reference.Y-fit.prediction.Y).',3,[]);
verification.maxPositionDifference_um=max(vecnorm(delta,2,1))*1e6;
if isinf(reference.param.ruptureTime) && isinf(fit.prediction.param.ruptureTime)
    verification.ruptureTimeDifference_s=0;
else
    verification.ruptureTimeDifference_s=abs(reference.param.ruptureTime-fit.prediction.param.ruptureTime);
end
verification.maxNetForceResidual_N=max(vecnorm(reshape(reference.nodeNetForce_time.',3,[]),2,1));
verification.maxAbsoluteEdgeForce_nN=max(abs(reference.edgeForce_time),[],'all')*1e9;
verification.wallSeconds=reference.wallSeconds;
assert(verification.maxPositionDifference_um<0.1,'Fitted trajectory is not numerically converged.');
assert(verification.ruptureTimeDifference_s<1,'Rupture time is not numerically converged.');
assert(verification.maxNetForceResidual_N<1e-18,'Overdamped force balance failed.');
verification.passed=true;
save(fullfile(fileparts(fitFile),'verification.mat'),'verification');
reportFile=fullfile(fileparts(fitFile),'tuning_report.txt');
fid=fopen(reportFile,'a');
assert(fid>=0,'Cannot append verification report.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'\nNUMERICAL VERIFICATION\n');
fprintf(fid,'10x tighter tolerances, maximum step reduced from 5 to 2.5 s.\n');
fprintf(fid,'Maximum node-position difference: %.9g um\n',verification.maxPositionDifference_um);
fprintf(fid,'Rupture-time difference: %.9g s\n',verification.ruptureTimeDifference_s);
fprintf(fid,'Maximum net-force residual: %.9g N\n',verification.maxNetForceResidual_N);
fprintf(fid,'Maximum absolute edge force at observed frame times: %.9g nN\n',verification.maxAbsoluteEdgeForce_nN);
fprintf(fid,'Passed; numerical consistency does not establish biological validity.\n');
disp(verification);
end
