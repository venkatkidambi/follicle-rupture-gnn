function figureHandle=plot_velocity_fit(fitFile)
%PLOT_VELOCITY_FIT Matched interval velocities, separated by training/holdout.
loaded=load(fitFile,'fit'); fit=loaded.fit; data=fit.data;
assert(isfield(data,'objective'),'Fit has no interval-velocity objective.');
figureHandle=figure('Color','w','Visible','off','Position',[100 100 1450 800]);
layout=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for split=1:2
    if split==1
        pairs=data.objective.trainPairs; nodes=data.trainNodes; label='Training';
    else
        pairs=data.objective.validationPairs; nodes=data.validationNodes; label='Held-out';
    end
    dt=(data.time(pairs(:,2))-data.time(pairs(:,1))).';
    time=mean([data.time(pairs(:,1)),data.time(pairs(:,2))],2)/60;
    observed=diffPairs(fit.metrics.observedAxial,pairs)./dt*1e6;
    predicted=diffPairs(fit.metrics.predictedAxial,pairs)./dt*1e6;
    source=diffPairs(fit.sourceMetrics.predictedAxial,pairs)./dt*1e6;
    for band=1:3
        valid=nodes & data.bandMasks(:,band) & ...
            data.valid(:,pairs(:,1)) & data.valid(:,pairs(:,2));
        a=observed; b=predicted; c=source;
        a(~valid)=NaN; b(~valid)=NaN; c(~valid)=NaN;
        nexttile; hold on;
        plot(time,mean(a,1,'omitnan'),'ko-','MarkerSize',3);
        plot(time,mean(c,1,'omitnan'),'--','Color',[.5 .5 .5],'LineWidth',1.4);
        plot(time,mean(b,1,'omitnan'),'Color',[.05 .4 .8],'LineWidth',1.8);
        yline(0,':','HandleVisibility','off'); grid on;
        title(sprintf('%s cells: initial %d-%d%% band',label,data.percentileRanges(band,:)));
        xlabel('Interval midpoint (min)'); ylabel('Mean axial velocity (\mum/s)');
        legend('Experiment','Position-only source','Position + velocity fit','Location','best');
    end
end
title(layout,'Matched interval velocities; objective uses all cell-wise 3D components');
exportgraphics(figureHandle,fullfile(fileparts(fitFile),'velocity_comparison.png'),'Resolution',180);
end

function values=diffPairs(values,pairs)
values=values(:,pairs(:,2))-values(:,pairs(:,1));
end
