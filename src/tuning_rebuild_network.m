function param = tuning_rebuild_network(base, edges, weakEdges, values)
%TUNING_REBUILD_NETWORK Build a static mechanical network, then set coefficients.
% Inputs contain initial geometry, connections, weak links, and global scalars
% only. Observed trajectories never enter network construction or dynamics.
N=base.N;
validateattributes(edges,{'numeric'},{'integer','positive','ncols',2,'<=',N});
assert(all(edges(:,1)<edges(:,2)) && size(unique(edges,'rows'),1)==size(edges,1), ...
    'follicle_rupture:invalidNetwork','Edges must be unique ordered pairs.');
Ne=size(edges,1);
validateattributes(weakEdges,{'numeric'},{'vector','integer','positive','<=',Ne,'nonempty'});
assert(max(conncomp(graph(edges(:,1),edges(:,2),[],N)))==1, ...
    'follicle_rupture:invalidNetwork','The initial network must be connected.');
param=base;
stale={'principalTipDamageEdge','tipmostDamageEdge','principalTipCandidateEdges','principalTipPatchEdges','ruptureTime'};
for k=1:numel(stale)
    if isfield(param,stale{k}), param=rmfield(param,stale{k}); end
end
param.edgeList=edges; param.Ne=Ne;
param.experimentalEdgeMode='fitted-static-local-network';
param.experimentalNeighborCount=[]; % explicit links, not one fixed k
param.edge_i=edges(:,1); param.edge_j=edges(:,2);
param.edgeAccumI=sparse(edges(:,1),(1:Ne)',1,N,Ne);
param.edgeAccumJ=sparse(edges(:,2),(1:Ne)',1,N,Ne);
param.L0_edge0=sqrt(sum((base.nodeCoords(edges(:,2),:)-base.nodeCoords(edges(:,1),:)).^2,2));
assert(all(param.L0_edge0>0),'follicle_rupture:invalidNetwork','Zero-length connection.');
fields={'k1_edge','k2_edge','c_edge','c_series_edge','connectionModelId_edge'};
for k=1:numel(fields)
    original=base.(fields{k});
    assert(all(original==original(1)), 'Network fitting requires uniform passive link coefficients.');
    param.(fields{k})=repmat(original(1),Ne,1);
end
param.edgeMid0=(base.nodeCoords(edges(:,1),:)+base.nodeCoords(edges(:,2),:))/2;
param.tipMask=double(base.postSign*param.edgeMid0(:,base.postAxis)>= ...
    prctile(base.postSign*base.nodeCoords(:,base.postAxis),base.tipPct));
param.damageMode='fitted-static-connected-patch';
param.damageMask=zeros(Ne,1); param.damageMask(weakEdges)=1;
param.rupturePatchNodes=unique(edges(weakEdges,:));
param.tipPatchEdgeCount=numel(unique(weakEdges));
score=base.principalNodeScore0;
mid=mean(score(edges),2);
alignment=abs(score(edges(:,2))-score(edges(:,1)))./param.L0_edge0;
reaches=mid<=base.anchorCutoff & ...
    ~(base.anchorNodeMask(edges(:,1))>0 & base.anchorNodeMask(edges(:,2))>0);
param.retractionMask=double(reaches & alignment>=base.contractionMinAxialAlignment & param.damageMask==0);
assert(any(param.retractionMask),'follicle_rupture:invalidNetwork','No active intact connections.');
active=param.retractionMask>0;
param.contractionRateFrac_edge=zeros(Ne,1);
param.minRestLengthFrac_edge=ones(Ne,1);
param.minRestLengthFrac_edge(active)=base.minRestLengthFrac;
param.contractionStartTime_edge=inf(Ne,1);
param.contractionDuration_edge=zeros(Ne,1);
param.extensionCritFrac_edge=inf(Ne,1);
param=tuning_parameters(param,values);
end
