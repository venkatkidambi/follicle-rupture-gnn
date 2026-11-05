function matlab_maxwell_reference(outputPath)
% Generate small reference cases using the ORIGINAL MATLAB mechanical RHS.
% Run from repository root, then compare with tests/python/test_viscoelastic.py.
addpath('src');
p.N=2; p.Ne=1; p.nPos=6; p.nVel=6;
p.nodeCoords=[-.5 0 0; .5 0 0]; p.edge_i=1; p.edge_j=2;
p.edgeAccumI=sparse(1,1,1,2,1); p.edgeAccumJ=sparse(2,1,1,2,1);
p.scale=struct('L',1,'tau',1,'k2',1,'inertia',1);
p.mechanicsProfile='cell-scale'; p.surfaceConstraintMode='none';
p.L0_edge0=1; p.contractionRateFrac_edge=.2; p.minRestLengthFrac_edge=.5;
p.contractionStartTime_edge=0; p.contractionRampTime=.1;
p.connectionModelId_edge=3; p.k1_edge=1; p.k2_edge=1;
p.c_edge=.3; p.c_series_edge=.3; p.g_floor_stiff=0; p.g_floor_damp=0;
p.damageMask=1; p.extensionCritFrac_edge=Inf; p.ruptureDamageTau=.1;
p.c_drag=1; p.c_surf=0; p.k_surf=0; p.m=1;
p.anchor_on=false; p.anchorNodeMask=zeros(2,1); p.activeTraction_on=false;
p.adhesionStiffness=2; p.adhesionDragStuck=.5; p.adhesionDragSliding=.1;
p.adhesionLoadRegularization=1e-12; p.adhesionLoadThreshold_node=.001*ones(2,1);
p.adhesionTransitionWidthFrac=.25; p.adhesionRuptureGateDamage=.5; p.adhesionReleaseTau=.1;
times=linspace(0,.6,7).';
output=struct('generator','Original src/contracting_natural_lengths.m and src/overdamped_shell_rhs.m', ...
    'matlab_version',version,'times',times);
for caseIndex=1:2
    q=1; p.adhesionClutch_on=false;
    if caseIndex==2, q=.2; p.adhesionClutch_on=true; end
    state=[reshape(p.nodeCoords.',[],1); .1; 0; 1-q; 0; 0];
    options=odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.001);
    [~,states]=ode15s(@(t,x)overdamped_shell_rhs(t,x,p),times,state,options);
    result=struct('q',q,'positions',states(:,1:6),'stress',q*states(:,7),'release',states(:,10:11));
    output.(sprintf('case%d',caseIndex))=result;
end
sampleTimes=[-.1 0 .05 .1 1 2.5 2.55 2.6 3];
lengths=zeros(size(sampleTimes)); rates=lengths;
for index=1:numel(sampleTimes)
    [lengths(index),rates(index)]=contracting_natural_lengths(sampleTimes(index),p);
end
output.contraction=struct('times',sampleTimes,'lengths',lengths,'rates',rates);
fid=fopen(outputPath,'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(output));
end
