function geometry = tuning_network_geometry(base)
%TUNING_NETWORK_GEOMETRY Fixed local candidate links from initial positions.
% A symmetric 12-neighbor pool plus the original links bounds all proposals.
% A Euclidean minimum spanning tree keeps every proposed graph connected.
xyz=base.nodeCoords; N=base.N;
d2=zeros(N);
for dim=1:3, d2=d2+(xyz(:,dim)-xyz(:,dim).').^2; end
d2(1:N+1:end)=inf;
[~,order]=sort(d2,2);
k=min(12,N-1);
pool=unique(sort([repelem((1:N)',k),reshape(order(:,1:k).',[],1);base.edgeList],2),'rows');
lengths=vecnorm(xyz(pool(:,2),:)-xyz(pool(:,1),:),2,2);
G=graph(pool(:,1),pool(:,2),lengths,N);
assert(max(conncomp(G))==1,'Candidate local-link pool is disconnected.');
tree=minspantree(G);
geometry.pool=pool;
geometry.backbone=sortrows(sort(tree.Edges.EndNodes,2));
geometry.lengths=lengths;
geometry.axialSpan=abs(base.principalNodeScore0(pool(:,2))-base.principalNodeScore0(pool(:,1)));
geometry.neighborPoolSize=k;
% Farthest-point centers cover the nonanchored initial surface without
% inspecting future motion or treating node IDs as an ordered variable.
eligible=find(base.anchorNodeMask==0);
[~,tip]=min(base.principalNodeScore0(eligible));
centers=eligible(tip); distance=inf(N,1);
for idx=2:min(20,numel(eligible))
    distance=min(distance,sum((xyz-xyz(centers(end),:)).^2,2));
    [~,next]=max(distance(eligible));
    centers(end+1)=eligible(next); %#ok<AGROW>
end
geometry.patchCenters=centers;
end
