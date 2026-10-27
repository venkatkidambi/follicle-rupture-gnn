function data = tuning_prepare_data(param, duration)
%TUNING_PREPARE_DATA Observations for scoring only; never used by an ODE RHS.
% Row identities and the initial geometry are fixed. Register each observed
% frame using only training cells in the initially anchored cap. Registration
% removes rigid pose, never scale or deformation, and never moves model nodes.
source = load(param.experimentalNodeDataFile, param.experimentalNodeDataVariable);
raw = double(source.(param.experimentalNodeDataVariable));
raw = raw(param.experimentalCellIds,:,:);
Nt = size(raw,3);
data.time = linspace(0,duration,Nt).';
affine = [raw(:,:,1),ones(param.N,1)] \ param.nodeCoords;
data.rawPosition = zeros(size(raw));
data.valid = reshape(all(isfinite(raw),2) & ~all(raw==0,2),param.N,Nt);
for frame = 1:Nt
    data.rawPosition(:,:,frame) = [raw(:,:,frame),ones(param.N,1)]*affine;
end
[~, order] = sort(param.principalNodeScore0);
data.trainNodes = true(param.N,1);
data.trainNodes(order(5:5:end)) = false;
data.validationNodes = ~data.trainNodes;
data.trainFrames = mod((1:Nt).',2)==0;
data.validationFrames = ~data.trainFrames;
data.validationFrames(1) = false; % initial conditions are not predictions
data.registrationNodes = data.trainNodes & param.anchorNodeMask>0;
data.position = zeros(size(raw));
data.rotation = zeros(3,3,Nt);
data.translation = zeros(Nt,3);
for frame = 1:Nt
    mask = data.registrationNodes & data.valid(:,frame);
    [data.position(:,:,frame),data.rotation(:,:,frame),data.translation(frame,:)] = ...
        tuning_rigid_registration(data.rawPosition(:,:,frame),param.nodeCoords,mask);
end
% Preserve exact zero initial displacement, avoiding an irrelevant fit offset.
data.position(:,:,1) = param.nodeCoords;
data.displacement = data.position-param.nodeCoords;
data.rawDisplacement = data.rawPosition-param.nodeCoords;
data.axis = param.principalLongAxis0;
data.percentileRanges = [0 10;25 35;45 55];
data.bandMasks = false(param.N,3);
for band = 1:3
    bounds = prctile(param.principalNodeScore0,data.percentileRanges(band,:));
    data.bandMasks(:,band) = param.principalNodeScore0>=bounds(1) & ...
        param.principalNodeScore0<=bounds(2);
end
data.frameInterval = duration/(Nt-1);
data.registrationDescription = 'Rigid pose fitted on training nodes in initial P90 cap; no scaling';
data.huberDelta_um = 5;
end
