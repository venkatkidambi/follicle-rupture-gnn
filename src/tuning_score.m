function metrics = tuning_score(result, data)
%TUNING_SCORE Compare predictions with measurements, outside the forward model.
prediction = permute(reshape(result.Y.',3,size(data.position,1),[]),[2 1 3]);
error_um = reshape(vecnorm(prediction-data.displacement,2,2),size(data.valid))*1e6;
train = data.trainNodes & data.trainFrames.' & data.valid;
validation = data.validationNodes & data.validationFrames.' & data.valid;
delta=data.huberDelta_um;
loss=delta^2*(sqrt(1+(error_um/delta).^2)-1);
metrics.objective=mean(loss(train));
metrics.trainRMSE_um=sqrt(mean(error_um(train).^2));
% Use squared error for selection: the genuine large late retraction must
% not be downweighted as an outlier. Keep robust loss only as a diagnostic.
metrics.robustObjective=metrics.objective;
metrics.objective=metrics.trainRMSE_um^2;
metrics.positionObjective=metrics.objective;
metrics.validationRMSE_um=sqrt(mean(error_um(validation).^2));
metrics.trainMedian_um=median(error_um(train));
metrics.validationMedian_um=median(error_um(validation));
metrics.error_um=error_um;
metrics.predictedAxial=reshape(sum(prediction.*reshape(data.axis,1,3,1),2),size(prediction,1),[]);
metrics.observedAxial=reshape(sum(data.displacement.*reshape(data.axis,1,3,1),2),size(prediction,1),[]);
metrics.rawRMSE_um=sqrt(mean(sum((prediction-data.rawDisplacement).^2,2),'all'))*1e6;
if isfield(data,'objective') && strcmp(data.objective.mode,'position-velocity')
    o=data.objective;
    [metrics.trainVelocityRMSE_um_s,metrics.trainVelocitySamples]= ...
        velocityError(prediction,data,o.trainPairs,data.trainNodes);
    [metrics.validationVelocityRMSE_um_s,metrics.validationVelocitySamples]= ...
        velocityError(prediction,data,o.validationPairs,data.validationNodes);
    metrics.trainPositionTerm=metrics.positionObjective/o.positionScale_um^2;
    metrics.trainVelocityTerm=metrics.trainVelocityRMSE_um_s^2/o.velocityScale_um_s^2;
    metrics.objective=metrics.trainPositionTerm+o.velocityWeight*metrics.trainVelocityTerm;
    metrics.validationJointObjective=metrics.validationRMSE_um^2/o.positionScale_um^2 + ...
        o.velocityWeight*metrics.validationVelocityRMSE_um_s^2/o.velocityScale_um_s^2;
end
end

function [rmse,count]=velocityError(prediction,data,pairs,nodes)
dt=reshape(data.time(pairs(:,2))-data.time(pairs(:,1)),1,1,[]);
predicted=(prediction(:,:,pairs(:,2))-prediction(:,:,pairs(:,1)))./dt;
observed=(data.displacement(:,:,pairs(:,2))-data.displacement(:,:,pairs(:,1)))./dt;
valid=nodes & data.valid(:,pairs(:,1)) & data.valid(:,pairs(:,2));
errorSquared=reshape(sum((predicted-observed).^2,2),size(valid))*1e12;
count=nnz(valid);
rmse=sqrt(mean(errorSquared(valid)));
end
