function tests = test_moving_surface
tests = functiontests(localfunctions);
end
function setupOnce(testCase)
testCase.TestData.originalPath=path;
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'src'));
end
function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end
function testRecoverRotatedTranslatedEllipsoidAndExcludeHoldout(testCase)
n=150; z=1-2*((1:n).'-0.5)/n; phi=(1:n).'*pi*(3-sqrt(5));
sphere=[sqrt(1-z.^2).*cos(phi),sqrt(1-z.^2).*sin(phi),z];
angle=.43; R=[cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
axes=[250 90 70]*1e-6; center=[20 -15 8]*1e-6;
data.time=[0;10]; data.trainNodes=true(n,1); data.trainNodes(5:5:end)=false;
data.valid=true(n,2);
data.position(:,:,1)=(sphere.*axes)*R.'+center;
data.position(:,:,2)=(sphere.*(axes*1.1))*R.'+center+[3 2 1]*1e-6;
data.position(~data.trainNodes,:,:)=100; % must not influence geometry
[surface,d]=fit_experimental_surface(data);
verifyEqual(testCase,surface.center(1,:),center,'AbsTol',1e-16);
verifyEqual(testCase,d.semiAxes_um(1,:),axes*1e6,'AbsTol',1e-8);
verifyEqual(testCase,d.volume_um3(2)/d.volume_um3(1),1.1^3,'AbsTol',1e-10);
verifyLessThan(testCase,max(d.rmsDistance_um),1e-8);
end
function testSurfaceVelocityMatchesImplicitDerivative(testCase)
p.movingSurface.time=[0;10];
p.movingSurface.center=[0 0 0;3 1 -2]*1e-6;
p.movingSurface.quadric(:,:,1)=diag(1./([20 10 8]*1e-6).^2);
p.movingSurface.quadric(:,:,2)=diag(1./([22 11 9]*1e-6).^2);
pos=[10 5 4;15 -3 2]*1e-6; t=4; h=1e-4;
[~,normal,speed]=shell_surface_geometry(t,pos,p);
% A point moving at this normal speed must keep its implicit level.
plus=shell_surface_geometry(t+h,pos+h*speed.*normal,p);
minus=shell_surface_geometry(t-h,pos-h*speed.*normal,p);
verifyEqual(testCase,(plus-minus)/(2*h),zeros(2,1),'AbsTol',1e-10);
verifyEqual(testCase,vecnorm(normal,2,2),ones(2,1),'AbsTol',1e-14);
end
function testStaticLimitAndClampedEndpoints(testCase)
p.ell_a=20e-6; p.ell_b=10e-6; p.ell_c=8e-6;
pos=[10 4 -3;18 1 1]*1e-6;
[level,normal,speed]=shell_surface_geometry(5,pos,p);
p.movingSurface.time=[0;10]; p.movingSurface.center=zeros(2,3);
p.movingSurface.quadric=repmat(diag(1./[p.ell_a p.ell_b p.ell_c].^2),1,1,2);
[movingLevel,movingNormal,movingSpeed]=shell_surface_geometry(5,pos,p);
verifyEqual(testCase,movingLevel,level); verifyEqual(testCase,movingNormal,normal);
verifyEqual(testCase,movingSpeed,speed);
p.movingSurface.center(2,1)=1e-6;
[~,~,speed]=shell_surface_geometry(11,pos,p);
verifyEqual(testCase,speed,zeros(2,1));
end
