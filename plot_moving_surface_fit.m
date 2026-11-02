function fig = plot_moving_surface_fit(fitFile)
%PLOT_MOVING_SURFACE_FIT Geometry diagnostics and conditional fit comparison.
repo=fileparts(mfilename('fullpath'));
if nargin<1, fitFile=fullfile(repo,'results','moving_surface_tuning','best_fit.mat'); end
loaded=load(fitFile,'fit'); fit=loaded.fit;
time=fit.data.time/60; d=fit.surfaceDiagnostics;
fig=figure('Visible','off','Color','w','Position',[100 100 1350 760]);
layout=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; plot(time,2*d.semiAxes_um,'LineWidth',1.5); grid on;
xlabel('Time (min)'); ylabel('Full axis length (\mum)');
legend('Longest','Middle','Shortest','Location','best'); title('Experimental ellipsoid axes');
nexttile; plot(time,d.volume_um3/d.volume_um3(1),'k-','LineWidth',1.5); grid on;
xlabel('Time (min)'); ylabel('Fitted volume / initial fitted volume');
title('Shape-fit proxy, not enclosed-volume measurement');
nexttile; plot(time,d.rmsDistance_um,'LineWidth',1.5); grid on;
xlabel('Time (min)'); ylabel('Approximate surface residual RMS (\mum)');
title('Ellipsoid fit residual on training cells');
for band=1:3
    nexttile; mask=fit.data.bandMasks(:,band);
    plot(time,mean(fit.metrics.observedAxial(mask,:),1)*1e6,'ko-','MarkerSize',3); hold on;
    plot(time,mean(fit.fixedSurfaceMetrics.predictedAxial(mask,:),1)*1e6,'--','LineWidth',1.5);
    if isfield(fit,'previousMovingSurfaceMetrics')
        plot(time,mean(fit.previousMovingSurfaceMetrics.predictedAxial(mask,:),1)*1e6, ...
            ':','Color',[.2 .5 .2],'LineWidth',1.8);
    end
    plot(time,mean(fit.metrics.predictedAxial(mask,:),1)*1e6,'LineWidth',1.8); grid on;
    xlabel('Time (min)'); ylabel('Mean axial displacement (\mum)');
    title(sprintf('Initial %d-%d%% band',fit.data.percentileRanges(band,:)));
    legend('Measured','Prior fixed surface','Refitted moving surface','Location','northwest');
    if isfield(fit,'previousMovingSurfaceMetrics')
        legend('Measured','Prior fixed surface','Prior moving surface','Full moving-surface refit','Location','northwest');
    end
end
title(layout,sprintf('Measured geometry at all frames | training RMSE %.2f µm | conditional holdout %.2f µm', ...
    fit.metrics.trainRMSE_um,fit.metrics.validationRMSE_um));
exportgraphics(fig,fullfile(fileparts(fitFile),'moving_surface_comparison.png'),'Resolution',180);
end
