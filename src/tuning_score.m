function metrics = tuning_score(result, data)
%TUNING_SCORE Compare predictions with measurements, outside the forward model.
prediction = permute(reshape(result.Y.',3,size(data.position,1),[]),[2 1 3]);
error_um = squeeze(vecnorm(prediction-data.displacement,2,2))*1e6;
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
metrics.validationRMSE_um=sqrt(mean(error_um(validation).^2));
metrics.trainMedian_um=median(error_um(train));
metrics.validationMedian_um=median(error_um(validation));
metrics.error_um=error_um;
metrics.predictedAxial=reshape(sum(prediction.*reshape(data.axis,1,3,1),2),size(prediction,1),[]);
metrics.observedAxial=reshape(sum(data.displacement.*reshape(data.axis,1,3,1),2),size(prediction,1),[]);
metrics.rawRMSE_um=sqrt(mean(sum((prediction-data.rawDisplacement).^2,2),'all'))*1e6;
end
