function export_follicle_rupture_motion_video(matFile, outFile, opts)
% export_follicle_rupture_motion_video
% Standalone MP4 renderer for follicleShellSim_*.mat output.

if nargin < 1 || isempty(matFile)
    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    matFile = fullfile(repoRoot, 'results', 'follicleShellSim_test2.mat');
end
if nargin < 2 || isempty(outFile)
    [matDir, matBase] = fileparts(matFile);
    outFile = fullfile(matDir, [matBase '_motion.mp4']);
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

[xLimFix, yLimFix, zLimFix, globalCap] = fixed_video_limits_and_arrow_cap( ...
    nodeCoords, Y, Nt, opts.arrowClipPct);

fprintf('\n=== Exporting follicle rupture video: %s ===\n', outFile);

fig = figure('Color','w','Units','inches','Position',[1 1 6.2 5.2]);
ax = axes('Parent',fig);
hold(ax,'on');
axis(ax,'equal');
grid(ax,'off');
box(ax,'on');

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

    d = P(i2,:) - P(i1,:);
    L = sqrt(sum(d.^2, 2));
    epsEdge = (L - L0) ./ L0;

    if ~isempty(D_time)
        brokenEdge = D_time(k,:).' > opts.damageBreak;
        damageBroken = D_time(k,damageEdge) > opts.damageBreak;
        nDamageEdges = nnz(damageEdge);
        nDamageBroken = nnz(damageBroken);
        damageFrac = nDamageBroken / max(nDamageEdges, 1);
    else
        brokenEdge = epsEdge > opts.strainBreak;
        damageBroken = epsEdge(damageEdge) > opts.strainBreak;
        nDamageEdges = nnz(damageEdge);
        nDamageBroken = nnz(damageBroken);
        damageFrac = nDamageBroken / max(nDamageEdges, 1);
    end

    draw_video_edges(ax, Pc_um, E, true(Ne,1), opts.showBaseEdges, [0.82 0.82 0.82], 0.25);
    draw_video_edges(ax, Pc_um, E, motorEdge, opts.showMotorEdges, [0.18 0.45 0.85], 0.65);
    draw_video_edges(ax, Pc_um, E, damageEdge, opts.showDamageEdges, [0.95 0.55 0.10], 1.15);
    draw_video_edges(ax, Pc_um, E, brokenEdge, opts.showBrokenEdges, [0.85 0.05 0.05], 2.2);

    mainNodeMask = ~(postNodeMask | tipNodeMask);
    scatter3(ax, Pc_um(mainNodeMask,1), Pc_um(mainNodeMask,2), Pc_um(mainNodeMask,3), ...
        14, [0.35 0.35 0.35], 'filled', 'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor','none');
    scatter3(ax, Pc_um(postNodeMask,1), Pc_um(postNodeMask,2), Pc_um(postNodeMask,3), ...
        18, [0.05 0.55 0.15], 'filled', 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor','none');
    scatter3(ax, Pc_um(tipNodeMask,1), Pc_um(tipNodeMask,2), Pc_um(tipNodeMask,3), ...
        24, [0.85 0.05 0.05], 'filled', 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','none');

    uPlot = clipped_relative_motion_um(U, globalCap);
    quiver3(ax, Pc_um(:,1), Pc_um(:,2), Pc_um(:,3), ...
        opts.arrowScale*uPlot(:,1), opts.arrowScale*uPlot(:,2), opts.arrowScale*uPlot(:,3), ...
        0, 'Color', [0 0.75 0.78], 'LineWidth', 1.0, 'MaxHeadSize', 0.45);

    xlim(ax, xLimFix);
    ylim(ax, yLimFix);
    zlim(ax, zLimFix);
    view(ax, opts.viewAz, opts.viewEl);
    daspect(ax,[1 1 1]);
    axis(ax,'vis3d');
    set(ax,'XTick',[],'YTick',[],'ZTick',[]);
    xlabel(ax,''); ylabel(ax,''); zlabel(ax,'');

    [maxEps, p95Eps] = damage_strain_summary(epsEdge, damageEdge);
    title(ax, sprintf(['Follicle rupture shell | t = %.2f s | ', ...
        'damage broken = %.2f (%d/%d) | max/p95 strain = %.3g / %.3g'], ...
        tSol(k), damageFrac, nDamageBroken, nDamageEdges, maxEps, p95Eps), ...
        'FontWeight','bold', 'FontSize', 10);

    drawnow;
    writeVideo(vw, getframe(fig));
end

close(vw);
close(fig);
fprintf('Saved %s\n', outFile);
end


function opts = video_defaults(opts)
defaults = struct( ...
    'fps', 30, ...
    'quality', 95, ...
    'targetFrames', 200, ...
    'arrowScale', 0.75, ...
    'arrowClipPct', 85, ...
    'viewAz', -35, ...
    'viewEl', 18, ...
    'scaleBarUm', 50, ...
    'scaleBarLabel', '50 \mum', ...
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


function [xLimFix, yLimFix, zLimFix, globalCap] = fixed_video_limits_and_arrow_cap(nodeCoords, Y, Nt, clipPct)
PallMin = [inf inf inf];
PallMax = [-inf -inf -inf];
magAll = [];

for k = 1:Nt
    U = reshape(Y(k,:),3,[]).';
    P = nodeCoords + U;
    C = median(P, 1, 'omitnan');
    Pc_um = (P - C) * 1e6;
    PallMin = min(PallMin, min(Pc_um, [], 1, 'omitnan'));
    PallMax = max(PallMax, max(Pc_um, [], 1, 'omitnan'));

    uRel_um = (U - median(U, 1, 'omitnan')) * 1e6;
    mag = vecnorm(uRel_um, 2, 2);
    magAll = [magAll; mag(isfinite(mag) & mag > 0)]; %#ok<AGROW>
end

span = PallMax - PallMin;
span(~isfinite(span) | span <= 0) = 1;
pad = 0.10 * span;
xLimFix = [PallMin(1)-pad(1), PallMax(1)+pad(1)];
yLimFix = [PallMin(2)-pad(2), PallMax(2)+pad(2)];
zLimFix = [PallMin(3)-pad(3), PallMax(3)+pad(3)];

if isempty(magAll)
    globalCap = [];
else
    globalCap = prctile(magAll, clipPct);
end
end


function uPlot = clipped_relative_motion_um(U, globalCap)
uRel_um = (U - median(U, 1, 'omitnan')) * 1e6;
mag = vecnorm(uRel_um, 2, 2);
good = isfinite(mag) & mag > 0;

if any(good) && ~isempty(globalCap) && isfinite(globalCap) && globalCap > 0
    scale = ones(size(mag));
    scale(good) = min(1, globalCap ./ mag(good));
    uPlot = uRel_um .* scale;
else
    uPlot = uRel_um;
end
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


function add_video_scale_bar_3D(ax, L, varargin)
p = inputParser;
p.addParameter('Label','', @(s)ischar(s)||isstring(s));
p.addParameter('Axis','x', @(s)ischar(s)||isstring(s));
p.addParameter('Frac',0.07, @(v)isnumeric(v)&&isscalar(v)&&v>0&&v<0.5);
p.addParameter('lw',2.0, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('Color',[0 0 0], @(v)isnumeric(v)&&numel(v)>=3);
p.parse(varargin{:});
S = p.Results;

xl = xlim(ax); yl = ylim(ax); zl = zlim(ax);
dx = diff(xl); dy = diff(yl); dz = diff(zl);
if dx <= 0, dx = 1; end
if dy <= 0, dy = 1; end
if dz <= 0, dz = 1; end

x0 = xl(1) + S.Frac*dx;
y0 = yl(1) + S.Frac*dy;
z0 = zl(1) + S.Frac*dz;
x1 = x0; y1 = y0; z1 = z0;

switch lower(string(S.Axis))
    case "x"
        x1 = x0 + L;
    case "y"
        y1 = y0 + L;
    case "z"
        z1 = z0 + L;
    otherwise
        error('Axis must be x, y, or z.');
end

plot3(ax, [x0 x1], [y0 y1], [z0 z1], '-', ...
    'Color', S.Color(1:3), 'LineWidth', S.lw, 'Clipping','off', 'HandleVisibility','off');

if strlength(string(S.Label)) > 0
    text(ax, (x0+x1)/2, y0 + 0.03*dy, z0, char(S.Label), ...
        'Color', S.Color(1:3), 'FontSize', 10, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'Clipping','off', 'Interpreter','tex', 'HandleVisibility','off');
end
end
