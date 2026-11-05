function matlab_zener_reference(outputPath)
% New fixtures from the unchanged original MATLAB RHS; no model files edited.
addpath('src');
p.N=2; p.Ne=1; p.nPos=6; p.nVel=6;
p.nodeCoords=[-.5 0 0; .5 0 0]; p.edge_i=1; p.edge_j=2;
p.edgeAccumI=sparse(1,1,1,2,1); p.edgeAccumJ=sparse(2,1,1,2,1);
p.scale=struct('L',1,'tau',1,'k2',1,'inertia',1);
p.mechanicsProfile='cell-scale'; p.surfaceConstraintMode='none';
p.L0_edge0=1; p.contractionRateFrac_edge=.2; p.minRestLengthFrac_edge=.5;
p.contractionStartTime_edge=0; p.contractionRampTime=.1;
p.k1_edge=.7; p.k2_edge=1; p.c_series_edge=.3;
p.g_floor_stiff=0; p.g_floor_damp=0;
p.damageMask=1; p.extensionCritFrac_edge=Inf; p.ruptureDamageTau=.1;
p.c_drag=1; p.c_surf=0; p.k_surf=0; p.m=1;
p.anchor_on=false; p.anchorNodeMask=zeros(2,1); p.activeTraction_on=false;
p.adhesionClutch_on=false;
times=linspace(0,.6,7).';
output=struct('generator','Unchanged src/overdamped_shell_rhs.m and follicle_shell_rhs.m', ...
    'matlab_version',version,'times',times);
names={'maxwell','zener','zener_damaged'};
for index=1:3
    q=1; p.connectionModelId_edge=5; p.c_edge=.3;
    if index==1, p.connectionModelId_edge=3; p.c_edge=.7*.3; end
    if index==3, q=.2; end
    initial=[-.55 0 0; .55 0 0];
    state=[reshape(initial.',[],1); .1; 0; 1-q; 0; 0];
    options=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.001);
    [~,states]=ode15s(@(t,x)overdamped_shell_rhs(t,x,p),times,state,options);
    output.(names{index})=struct('q',q,'positions',states(:,1:6), ...
        'transient_force',q*states(:,7));
end
% Analytical clamped step-strain response with sigma(0)=k2*epsilon.
output.step_strain=struct('strain',.1,'k_equilibrium',.7,'k_transient',1,'tau',.3, ...
    'maxwell_force',.1*exp(-times/.3),'zener_force',.7*.1+.1*exp(-times/.3));
fid=fopen(outputPath,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(output));
end
