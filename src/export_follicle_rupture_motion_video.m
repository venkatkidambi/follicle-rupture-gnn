function export_follicle_rupture_motion_video(matFile, outFile, opts)
% Standalone MP4 renderer for follicleShellSim_*.mat output.
% No arrows. Opaque markers to reduce video artifacts / scan lines.

if nargin < 1 || isempty(matFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    matFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end
if nargin < 2 || isempty(outFile)
    [matDir, matBase] = fileparts(matFile);
    outFile = fullfile(matDir, [matBase '_motion_noArrows.mp4']);
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

opts = video_defaults(opts);

S = load(matFile);
nodeCoords = S.nodeCoords;
edgeList = S.edgeList;
tSol = S.tSol(:);
Y = S.Y;
param = S.param;

Nt = numel(tSol);
E = edgeList;
Ne = size(E,1);

if isfield(S,'D_time')
    D_time = S.D_time;
else
    D_time = [];
end

[motorEdge, ~, damageEdge, postNodeMask, tipNodeMask] = follicle_video_masks(nodeCoords, E, param);

i1 = E(:,1);
i2 = E(:,2);

if isfield(param,'L0_edge0')
    L0 = param.L0_edge0(:);
else
    L0 = sqrt(sum((nodeCoords(i2,:) - nodeCoords(i1,:)).^2, 2));
end
L0 = max(L0, 1e-18);

[xLimFix, yLimFix, zLimFix] = fixed_video_limits(nodeCoords, Y, Nt);

fprintf('\n=== Exporting follicle rupture video: %s ===\n', outFile);

frameSizePx = [900 760];

fig = figure( ...
    'Color','w', ...
    'Units','pixels', ...
    'Position',[100 100 frameSizePx], ...
    'Renderer','painters', ...
    'Resize','off');

ax = axes('Parent',fig, ...
    'Units','pixels', ...
    'Position',[1 1 frameSizePx], ...
    'ActivePositionProperty','position');
hold(ax,'on');
axis(ax,'equal');
grid(ax,'off');
box(ax,'on');
configure_video_axes_geometry(fig, ax, frameSizePx);

vw = VideoWriter(outFile, 'MPEG-4');
vw.FrameRate = opts.fps;
vw.Quality = opts.quality;
open(vw);

frameSkip = max(1, floor(Nt / opts.targetFrames));

for k = 1:frameSkip:Nt
    cla(ax);

    U = reshape(Y(k,:),3,[]).';
    P = nodeCoords + U;
    C = median(P, 1, 'omitnan');
    Pc_um = (P - C) * 1e6;
    disp_um = sqrt(sum(U.^2, 2)) * 1e6;
    disp_um_finite = disp_um(isfinite(disp_um));

    d = P(i2,:) - P(i1,:);
    L = sqrt(sum(d.^2, 2));
    epsEdge = (L - L0) ./ L0;

    if ~isempty(D_time)
        brokenEdge = D_time(k,:).' > opts.damageBreak;
        damageBroken = D_time(k,damageEdge) > opts.damageBreak;
    else
        brokenEdge = epsEdge > opts.strainBreak;
        damageBroken = epsEdge(damageEdge) > opts.strainBreak;
    end

    nDamageEdges = nnz(damageEdge);
    nDamageBroken = nnz(damageBroken);
    damageFrac = nDamageBroken / max(nDamageEdges, 1);

    draw_video_edges(ax, Pc_um, E, true(Ne,1), opts.showBaseEdges, [0.82 0.82 0.82], 0.25);
    draw_video_edges(ax, Pc_um, E, motorEdge, opts.showMotorEdges, [0.18 0.45 0.85], 0.65);
    draw_video_edges(ax, Pc_um, E, damageEdge, opts.showDamageEdges, [0.95 0.55 0.10], 1.15);
    draw_video_edges(ax, Pc_um, E, brokenEdge, opts.showBrokenEdges, [0.85 0.05 0.05], 2.2);

    mainNodeMask = ~(postNodeMask | tipNodeMask);
    postOnlyNodeMask = postNodeMask & ~tipNodeMask;

    scatter3(ax, Pc_um(mainNodeMask,1), Pc_um(mainNodeMask,2), Pc_um(mainNodeMask,3), ...
        14, [0.35 0.35 0.35], 'filled', 'MarkerEdgeColor','none');

    scatter3(ax, Pc_um(postOnlyNodeMask,1), Pc_um(postOnlyNodeMask,2), Pc_um(postOnlyNodeMask,3), ...
        18, [0.05 0.55 0.15], 'filled', 'MarkerEdgeColor','none');

    scatter3(ax, Pc_um(tipNodeMask,1), Pc_um(tipNodeMask,2), Pc_um(tipNodeMask,3), ...
        24, [0.85 0.05 0.05], 'filled', 'MarkerEdgeColor','none');

    xlim(ax, xLimFix);
    ylim(ax, yLimFix);
    zlim(ax, zLimFix);

    view(ax, opts.viewAz, opts.viewEl);
    daspect(ax,[1 1 1]);
    axis(ax,'vis3d');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    xlabel(ax,''); ylabel(ax,''); zlabel(ax,'');

    [maxEps, p95Eps] = damage_strain_summary(epsEdge, damageEdge);
    if isempty(disp_um_finite)
        medianDisp = NaN;
        maxDisp = NaN;
    else
        medianDisp = median(disp_um_finite);
        maxDisp = max(disp_um_finite);
    end

    titleText = sprintf(['Follicle rupture shell | t = %.2f s | ', ...
        'damage broken = %.2f (%d/%d) | max/p95 strain = %.3g / %.3g | ', ...
        'median/max disp = %.2f / %.2f um'], ...
        tSol(k), damageFrac, nDamageBroken, nDamageEdges, maxEps, p95Eps, ...
        medianDisp, maxDisp);

    title(ax, titleText, ...
        'FontWeight','bold', 'FontSize', 10);
    text(ax, 0.5, 0.98, titleText, ...
        'Units','normalized', ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontWeight','bold', ...
        'FontSize',10, ...
        'Color','k', ...
        'BackgroundColor','w', ...
        'Margin',2, ...
        'Clipping','off', ...
        'Interpreter','none');

    configure_video_axes_geometry(fig, ax, frameSizePx);
    drawnow;
    frame = getframe(ax);
    writeVideo(vw, frame);
end

close(vw);
close(fig);

fprintf('Saved %s\n', outFile);
end


function configure_video_axes_geometry(fig, ax, frameSizePx)
set(fig, ...
    'Units','pixels', ...
    'Position',[100 100 frameSizePx], ...
    'Renderer','painters', ...
    'Resize','off');

set(ax, ...
    'Units','pixels', ...
    'Position',[1 1 frameSizePx], ...
    'ActivePositionProperty','position');
end


function opts = video_defaults(opts)
defaults = struct( ...
    'fps', 30, ...
    'quality', 95, ...
    'targetFrames', 200, ...
    'viewAz', -35, ...
    'viewEl', 18, ...
    'strainBreak', 0.25, ...
    'damageBreak', 0.99, ...
    'showBaseEdges', false, ...
    'showMotorEdges', false, ...
    'showDamageEdges', false, ...
    'showBrokenEdges', false);

names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
end


function [xLimFix, yLimFix, zLimFix] = fixed_video_limits(nodeCoords, Y, Nt)
PallMin = [inf inf inf];
PallMax = [-inf -inf -inf];

for k = 1:Nt
    U = reshape(Y(k,:),3,[]).';
    P = nodeCoords + U;
    C = median(P, 1, 'omitnan');
    Pc_um = (P - C) * 1e6;

    PallMin = min(PallMin, min(Pc_um, [], 1, 'omitnan'));
    PallMax = max(PallMax, max(Pc_um, [], 1, 'omitnan'));
end

span = PallMax - PallMin;
span(~isfinite(span) | span <= 0) = 1;

pad = 0.10 * span;

xLimFix = [PallMin(1)-pad(1), PallMax(1)+pad(1)];
yLimFix = [PallMin(2)-pad(2), PallMax(2)+pad(2)];
zLimFix = [PallMin(3)-pad(3), PallMax(3)+pad(3)];
end


function [motorEdge, tipEdge, damageEdge, postNodeMask, tipNodeMask] = follicle_video_masks(nodeCoords, E, param)
edgeMid0 = 0.5 * (nodeCoords(E(:,1),:) + nodeCoords(E(:,2),:));

if isfield(param,'postAxis')
    postAxis = param.postAxis;
else
    postAxis = 1;
end

if isfield(param,'postSign')
    postSign = param.postSign;
else
    postSign = -1;
end

sNodes = postSign * nodeCoords(:,postAxis);
sMid = postSign * edgeMid0(:,postAxis);

sCutMain = prctile(sNodes, getParamOrDefault(param, 'mainPct', 80));
sCutTip = prctile(sNodes, getParamOrDefault(param, 'tipPct', 92));

if isfield(param,'motorMask')
    motorEdge = logical(param.motorMask(:));
else
    motorEdge = sMid < sCutMain;
end

if isfield(param,'tipMask')
    tipEdge = logical(param.tipMask(:));
else
    tipEdge = sMid >= sCutTip;
end

if isfield(param,'damageMask')
    damageEdge = logical(param.damageMask(:));
else
    damageEdge = tipEdge;
end

postNodeMask = sNodes >= sCutMain;
tipNodeMask = sNodes >= sCutTip;
end


function val = getParamOrDefault(param, fieldName, defaultVal)
if isfield(param, fieldName) && ~isempty(param.(fieldName))
    val = param.(fieldName);
else
    val = defaultVal;
end
end


function draw_video_edges(ax, P, E, edgeMask, shouldDraw, colorRGB, lineWidth)
if ~shouldDraw
    return;
end

idx = find(edgeMask(:));
if isempty(idx)
    return;
end

[X,Y,Z] = video_batch_edges(P, E, idx);

plot3(ax, X, Y, Z, '-', ...
    'Color', colorRGB, ...
    'LineWidth', lineWidth, ...
    'Clipping','off', ...
    'HandleVisibility','off');
end


function [X,Y,Z] = video_batch_edges(P, E, idx)
ii = E(idx,1);
jj = E(idx,2);

X = [P(ii,1) P(jj,1) nan(numel(idx),1)].';
Y = [P(ii,2) P(jj,2) nan(numel(idx),1)].';
Z = [P(ii,3) P(jj,3) nan(numel(idx),1)].';

X = X(:);
Y = Y(:);
Z = Z(:);
end


function [maxEps, p95Eps] = damage_strain_summary(epsEdge, damageEdge)
epsDamage = epsEdge(damageEdge);
epsDamageFinite = epsDamage(isfinite(epsDamage));

if isempty(epsDamageFinite)
    maxEps = NaN;
    p95Eps = NaN;
else
    maxEps = max(epsDamageFinite);
    p95Eps = prctile(epsDamageFinite, 95);
end
end
