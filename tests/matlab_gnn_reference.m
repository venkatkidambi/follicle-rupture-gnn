function matlab_gnn_reference(outputPath)
% Deterministic fixtures invoking unchanged MATLAB mechanics, all five families.
addpath('src');
p.N=3; p.Ne=2; p.nPos=9; p.nVel=9;
p.nodeCoords=[0 0 0;1 .2 .1;2 -.1 .3];
p.edge_i=[1;2]; p.edge_j=[2;3];
p.edgeAccumI=sparse(p.edge_i,(1:2)',1,3,2);
p.edgeAccumJ=sparse(p.edge_j,(1:2)',1,3,2);
p.scale=struct('L',1,'tau',1,'k2',1,'inertia',1);
p.mechanicsProfile='cell-scale'; p.surfaceConstraintMode='none';
p.L0_edge0=vecnorm(p.nodeCoords(p.edge_j,:)-p.nodeCoords(p.edge_i,:),2,2);
p.contractionRateFrac_edge=[.2;.15]; p.minRestLengthFrac_edge=[.5;.5];
p.contractionStartTime_edge=[0;0]; p.contractionRampTime=.1;
p.k1_edge=[.7;.9]; p.k2_edge=[.8;1.1];
p.c_edge=[.3;.4]; p.c_series_edge=[.25;.35];
p.g_floor_stiff=0; p.g_floor_damp=1e-6;
p.damageMask=[1;0]; p.extensionCritFrac_edge=[.02;.02];
p.ruptureDamageTau=.08; p.rupturePatchCascade_on=false;
p.c_drag=1; p.c_surf=0; p.k_surf=0; p.m=1;
p.anchor_on=true; p.anchorNodeMask=[0;0;1]; p.k_anchor=1; p.c_anchor=.2;
p.activeTraction_on=false; p.adhesionClutch_on=true;
p.adhesionStiffness=.7; p.adhesionDragStuck=.05; p.adhesionDragSliding=.01;
p.adhesionLoadRegularization=1e-12; p.adhesionLoadThreshold_node=.008*ones(3,1);
p.adhesionTransitionWidthFrac=.25; p.adhesionRuptureGateDamage=.5; p.adhesionReleaseTau=.15;
names={'elastic','maxwell','kelvin-voigt','zener','burgers'}; modelIds=[1 3 4 5 6];
output=struct('generator','Unmodified follicle_shell_rhs and overdamped_shell_rhs', 'matlab_version',version);
pos=p.nodeCoords+[0 0 0;.08 -.01 .02;.05 .01 -.01];
state=[reshape(pos.',[],1);.12;-.03;.02;.01;.7;0;.1;.2;0];
for index=1:numel(names)
    p.connectionModelId_edge=repmat(modelIds(index),2,1);
    [derivative,fullState]=overdamped_shell_rhs(.4,state,p);
    [~,forces]=follicle_shell_rhs(.4,fullState,p);
    output.cases(index)=struct('family',names{index},'derivative',derivative, ...
        'net_force',forces.nodeNetForce,'edge_force',forces.edgeForce,'times',[],'history',[]);
    % Smooth interval with no threshold crossing, for full integrator comparison.
    p.extensionCritFrac_edge=[Inf;Inf];
    times=linspace(.4,.5,6)';
    options=odeset('RelTol',1e-11,'AbsTol',1e-13,'MaxStep',.0005);
    [~,history]=ode15s(@(t,x)overdamped_shell_rhs(t,x,p),times,state,options);
    output.cases(index).times=times;
    output.cases(index).history=history;
    p.extensionCritFrac_edge=[.02;.02];
end
sampleTimes=[-.1 0 .05 .1 1 2.5 2.55 2.6 3.5];
for index=1:numel(sampleTimes)
    [L,Ldot]=contracting_natural_lengths(sampleTimes(index),p);
    output.contraction(index)=struct('time',sampleTimes(index),'length',L,'rate',Ldot);
end
fid=fopen(outputPath,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(output));
end
