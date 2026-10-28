function tests=test_tuning_network
tests=functiontests(localfunctions);
end
function setupOnce(testCase)
testCase.TestData.originalPath=path;
testCase.TestData.repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(testCase.TestData.repo,'src'));
end
function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end
function testLocalConnectedGraphsAndPatches(testCase)
base=fixture;
geometry=tuning_network_geometry(base);
for k=[4 6 10]
    for weight=[.5 1 2]
        edges=tuning_local_network(base,geometry,k,weight);
        verifyTrue(testCase,all(ismember(edges,geometry.pool,'rows')));
        verifyEqual(testCase,max(conncomp(graph(edges(:,1),edges(:,2),[],base.N))),1);
        verifyTrue(testCase,all(degree(graph(edges(:,1),edges(:,2),[],base.N))>=k));
        weak=tuning_rupture_patch(base,edges,geometry.patchCenters(1),8);
        nodes=unique(edges(weak,:));
        verifyEqual(testCase,numel(weak),8);
        verifyFalse(testCase,any(base.anchorNodeMask(nodes)>0));
        G=graph(edges(weak,1),edges(weak,2),[],base.N);
        verifyEqual(testCase,max(conncomp(subgraph(G,nodes))),1);
    end
end
end
function testRebuildChangesStateSizeAndPreservesInitialPositions(testCase)
base=fixture; geometry=tuning_network_geometry(base);
edges=tuning_local_network(base,geometry,4,.5);
weak=tuning_rupture_patch(base,edges,1,8);
param=tuning_rebuild_network(base,edges,weak,[1e-4 300 100 .5 30 .03 500]);
verifyEqual(testCase,param.nodeCoords,base.nodeCoords);
verifyEqual(testCase,param.Ne,size(edges,1));
verifyEqual(testCase,size(param.edgeAccumI),[base.N param.Ne]);
verifyEqual(testCase,full(sum(param.edgeAccumI,1)),ones(1,param.Ne));
verifyEqual(testCase,full(sum(param.edgeAccumJ,1)),ones(1,param.Ne));
verifyEqual(testCase,param.L0_edge0,sqrt(sum((base.nodeCoords(edges(:,1),:)-base.nodeCoords(edges(:,2),:)).^2,2)));
verifyEqual(testCase,param.contractionRateFrac_edge(weak),zeros(numel(weak),1));
verifyTrue(testCase,all(isinf(param.extensionCritFrac_edge(param.damageMask==0))));
verifyFalse(testCase,isfield(param,'ruptureTime'));
verifyEqual(testCase,param.c_edge,100*param.k2_edge);
verifyEqual(testCase,param.adhesionLoadThreshold_node,.5*base.adhesionLoadThreshold_node);
end
function testInvalidGraphsRejected(testCase)
base=fixture;
verifyError(testCase,@() tuning_rebuild_network(base,[1 2;2 3],1,[1e-4 300 100 .5 30 .03 500]), ...
    'follicle_rupture:invalidNetwork');
verifyError(testCase,@() tuning_rebuild_network(base,[1 2;1 2],1,[1e-4 300 100 .5 30 .03 500]), ...
    'follicle_rupture:invalidNetwork');
end
function testOriginalNetworkPreservesForces(testCase)
file=fullfile(testCase.TestData.repo,'results','follicleShellSim_cell_scale.mat');
assumeTrue(testCase,isfile(file),'Local calibrated baseline required for integration regression.');
loaded=load(file,'param'); base=loaded.param;
values=[1e-4 300 100 .5 30 .03 500];
a=tuning_parameters(base,values);
b=tuning_rebuild_network(base,base.edgeList,find(base.damageMask),values);
fields={'L0_edge0','k1_edge','k2_edge','c_edge','damageMask','retractionMask', ...
    'contractionStartTime_edge','contractionDuration_edge','extensionCritFrac_edge','edgeAccumI','edgeAccumJ'};
for idx=1:numel(fields), verifyEqual(testCase,a.(fields{idx}),b.(fields{idx})); end
x=zeros(a.nPos+a.nVel+3*a.Ne+a.N,1);
x(1:a.nPos)=reshape(a.nodeCoords.',[],1)/a.scale.L;
x(1:3:a.nPos)=x(1:3:a.nPos)+.001*sin((1:a.N)');
[dxA,fA]=follicle_shell_rhs(10,x,a);
[dxB,fB]=follicle_shell_rhs(10,x,b);
verifyEqual(testCase,dxA,dxB);
verifyEqual(testCase,fA,fB);
end
function base=fixture
[x,y,z]=ndgrid(0:4,0:2,0:1);
base.nodeCoords=[x(:)*2,y(:),z(:)]*1e-5; base.N=numel(x);
base.edgeList=[(1:base.N-1)',(2:base.N)'];
base.principalNodeScore0=base.nodeCoords(:,1); base.principalLongAxis0=[1;0;0];
base.anchorCutoff=7e-5; base.anchorNodeMask=double(base.principalNodeScore0>=base.anchorCutoff);
base.postSign=1; base.postAxis=1; base.tipPct=90;
base.k1_edge=.002; base.k2_edge=.002; base.c_edge=.1;
base.c_series_edge=.1; base.connectionModelId_edge=3;
base.contractionMinAxialAlignment=.3; base.minRestLengthFrac=.5;
base.adhesionLoadThresholdTip=1e-10; base.adhesionLoadThresholdBody=1e-9;
base.adhesionLoadThreshold_node=ones(base.N,1)*1e-9;
base.activeTractionStartTime_node=zeros(base.N,1);
base.scale.k2=.002; base.m=5e-13; base.ruptureTime=42;
end
