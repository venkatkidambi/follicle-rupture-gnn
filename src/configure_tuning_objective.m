function data = configure_tuning_objective(data, mode)
%CONFIGURE_TUNING_OBJECTIVE Fix objective scales using training observations only.
% Velocities are interval averages between successive frames in the SAME
% split, using both valid endpoints. No smoothing crosses held-out frames.
mode=validatestring(mode,{'position','position-velocity'});
if isfield(data,'objective'), data=rmfield(data,'objective'); end
if strcmp(mode,'position'), return; end
trainFrames=find(data.trainFrames);
validationFrames=find(data.validationFrames);
assert(numel(trainFrames)>=2 && numel(validationFrames)>=2, ...
    'follicle_rupture:velocityFrames','Each split needs at least two frames.');
assert(all(diff(data.time)>0),'Time must be strictly increasing.');
objective.mode=mode;
objective.velocityWeight=1;
objective.trainPairs=[trainFrames(1:end-1),trainFrames(2:end)];
objective.validationPairs=[validationFrames(1:end-1),validationFrames(2:end)];
train=data.trainNodes & data.trainFrames.' & data.valid;
positionSquared=reshape(sum(data.displacement.^2,2),size(data.valid))*1e12;
objective.positionScale_um=sqrt(mean(positionSquared(train)));
pairs=objective.trainPairs;
dt=reshape(data.time(pairs(:,2))-data.time(pairs(:,1)),1,1,[]);
velocity=diffEndpoints(data.displacement,pairs)./dt;
valid=data.trainNodes & data.valid(:,pairs(:,1)) & data.valid(:,pairs(:,2));
velocitySquared=reshape(sum(velocity.^2,2),size(valid))*1e12;
objective.velocityScale_um_s=sqrt(mean(velocitySquared(valid)));
assert(any(valid,'all') && isfinite(objective.positionScale_um) && ...
    objective.positionScale_um>0 && isfinite(objective.velocityScale_um_s) && ...
    objective.velocityScale_um_s>0,'follicle_rupture:objectiveScale', ...
    'Joint fitting requires nonzero finite training displacement and velocity scales.');
objective.description=['MSE(position)/training observed displacement RMS^2 + ' ...
    'MSE(interval 3D velocity)/training observed interval velocity RMS^2; weight 1'];
data.objective=objective;
end

function delta=diffEndpoints(values,pairs)
delta=values(:,:,pairs(:,2))-values(:,:,pairs(:,1));
end
