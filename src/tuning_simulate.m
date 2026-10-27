function result = tuning_simulate(param, time, fidelity, maxSeconds, saveForces)
%TUNING_SIMULATE Forward prediction using only parameters and initial geometry.
% No measurements, target displacements, or time-varying position constraints
% are accepted. The existing natural-length actuation drives the mechanics.
if nargin<5, saveForces=false; end
if nargin<4, maxSeconds=inf; end
if nargin<3, fidelity='fine'; end
switch fidelity
    case 'coarse'
        relativeTolerance=3e-4; absoluteTolerance=3e-7; maximumStep=20;
    case 'fine'
        relativeTolerance=1e-5; absoluteTolerance=1e-8; maximumStep=5;
    case 'verification'
        relativeTolerance=1e-6; absoluteTolerance=1e-9; maximumStep=2.5;
    otherwise
        error('follicle_rupture:tuningFidelity','Use coarse, fine, or verification.');
end
initial = zeros(param.nPos+param.nVel+3*param.Ne+param.N,1);
initial(1:param.nPos) = reshape(param.nodeCoords.',[],1)/param.scale.L;
started = tic;
opts=odeset('RelTol',relativeTolerance,'AbsTol',absoluteTolerance, ...
    'MaxStep',maximumStep/param.scale.tau,'InitialStep',1e-3/param.scale.tau, ...
    'Stats','off','OutputFcn',@checkBudget);
[states,param]=integrate_cell_scale(param,time/param.scale.tau,initial,opts);
result.wallSeconds=toc(started);
result.param=param;
result.tSol=time;
result.nodeCoords=param.nodeCoords;
result.edgeList=param.edgeList;
result.L0_edge0=param.L0_edge0;
result.Y=(states(:,1:param.nPos)-reshape(param.nodeCoords.',1,[])/param.scale.L)*param.scale.L;
result.V=states(:,param.nPos+(1:param.nVel))*param.scale.L/param.scale.tau;
off=param.nPos+param.nVel+2*param.Ne;
result.D_time=min(max(states(:,off+(1:param.Ne)),0),1);
result.A_time=min(max(states(:,off+param.Ne+(1:param.N)),0),1);
if saveForces
    Nt=numel(time);
    result.Lnatural_time=zeros(Nt,param.Ne);
    result.edgeForce_time=zeros(Nt,param.Ne);
    result.nodeInternalForce_time=zeros(Nt,param.nPos);
    result.nodeNetForce_time=zeros(Nt,param.nPos);
    for idx=1:Nt
        result.Lnatural_time(idx,:)=contracting_natural_lengths(time(idx),param).';
        [~,forces]=follicle_shell_rhs(time(idx)/param.scale.tau,states(idx,:).',param);
        result.edgeForce_time(idx,:)=forces.edgeForce.';
        result.nodeInternalForce_time(idx,:)=reshape(forces.nodeInternalForce.',1,[]);
        result.nodeNetForce_time(idx,:)=reshape(forces.nodeNetForce.',1,[]);
    end
end

    function status=checkBudget(~,~,flag)
        status=0;
        if isempty(flag) && toc(started)>maxSeconds
            error('follicle_rupture:tuningTimeout','Trial exceeded %.1f s computational budget.',maxSeconds);
        end
    end
end
