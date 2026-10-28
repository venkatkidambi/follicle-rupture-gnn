function figureHandle=plot_network_fit(fitFile)
%PLOT_NETWORK_FIT Compare the original and inferred static mechanical graphs.
repo=fileparts(mfilename('fullpath'));
if nargin<1, fitFile=fullfile(repo,'results','network_tuning','best_fit.mat'); end
loaded=load(fitFile,'fit'); fit=loaded.fit;
source=load(fullfile(fileparts(fitFile),'search_history.mat'),'configuration');
original=source.configuration.base;
figureHandle=figure('Color','w','Visible','off','Position',[100 100 1400 620]);
layout=tiledlayout(1,2,'TileSpacing','compact');
for idx=1:2
    if idx==1, param=original; label='Original network'; else, param=fit.param; label='Fitted network'; end
    nexttile; hold on;
    xyz=param.nodeCoords*1e6; edges=param.edgeList;
    lines(edges,[.75 .78 .81],.4);
    lines(edges(param.damageMask>0,:),[.85 .1 .15],2.5);
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),8,[.2 .25 .3],'filled');
    nodes=param.rupturePatchNodes;
    scatter3(xyz(nodes,1),xyz(nodes,2),xyz(nodes,3),30,[.85 .1 .15],'filled');
    axis equal; view(25,20); grid on;
    xlabel('x (\mum)'); ylabel('y (\mum)'); zlabel('z (\mum)');
    title(sprintf('%s: %d links, %d weak links',label,param.Ne,nnz(param.damageMask)));
end
title(layout,'Static connections and weak rupture patch (red); positions are initial measurements');
exportgraphics(figureHandle,fullfile(fileparts(fitFile),'network_comparison.png'),'Resolution',180);
    function lines(edges,color,width)
        X=[xyz(edges(:,1),1),xyz(edges(:,2),1),nan(size(edges,1),1)].';
        Y=[xyz(edges(:,1),2),xyz(edges(:,2),2),nan(size(edges,1),1)].';
        Z=[xyz(edges(:,1),3),xyz(edges(:,2),3),nan(size(edges,1),1)].';
        plot3(X(:),Y(:),Z(:),'Color',color,'LineWidth',width);
    end
end
