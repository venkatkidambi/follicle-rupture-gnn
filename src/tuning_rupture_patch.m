function weakEdges = tuning_rupture_patch(base, edges, centerNode, edgeCount)
%TUNING_RUPTURE_PATCH Grow a compact connected patch on the initial graph.
validateattributes(centerNode,{'numeric'},{'scalar','integer','positive','<=',base.N});
validateattributes(edgeCount,{'numeric'},{'scalar','integer','positive'});
assert(base.anchorNodeMask(centerNode)==0,'Patch center cannot be an anchored node.');
mid=(base.nodeCoords(edges(:,1),:)+base.nodeCoords(edges(:,2),:))/2;
distance=sum((mid-base.nodeCoords(centerNode,:)).^2,2);
eligible=~any(base.anchorNodeMask(edges)>0,2);
selected=false(size(edges,1),1); nodes=false(base.N,1); nodes(centerNode)=true;
for k=1:edgeCount
    frontier=find(eligible & ~selected & any(nodes(edges),2));
    if isempty(frontier), break; end
    [~,idx]=min(distance(frontier)); edge=frontier(idx);
    selected(edge)=true; nodes(edges(edge,:))=true;
end
weakEdges=find(selected);
assert(~isempty(weakEdges),'follicle_rupture:invalidNetwork','No eligible patch connections.');
end
