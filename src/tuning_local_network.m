function edges = tuning_local_network(base, geometry, neighborCount, axialWeight)
%TUNING_LOCAL_NETWORK Anisotropic local nearest-neighbor graph plus backbone.
validateattributes(neighborCount,{'numeric'},{'scalar','integer','positive','<=',geometry.neighborPoolSize});
validateattributes(axialWeight,{'numeric'},{'scalar','positive','finite'});
pool=geometry.pool;
cost=sqrt(max(geometry.lengths.^2+(axialWeight^2-1)*geometry.axialSpan.^2,0));
selected=false(size(cost));
for node=1:base.N
    incident=find(any(pool==node,2));
    [~,order]=sort(cost(incident));
    selected(incident(order(1:min(neighborCount,numel(order)))))=true;
end
edges=unique([pool(selected,:);geometry.backbone],'rows');
end
