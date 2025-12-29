function plotMNIslice(lnm, brain, varargin)
% PLOTMNISLICE  Display one slice from a volumetric lesion network map (LNM)
% on top of an anatomical underlay. Supports sagittal, coronal, and transverse 
% planes, with masking and independent colormaps for underlay and overlay.
%
% Input:
%   lnm       : 3-D numeric array (overlay volume; e.g., LNM values).
%   brain     : 3-D numeric array (underlay volume; e.g., MNI152 T1).
%               Both volumes should be in the same space/grid/orientation.
%
% Optional:
%   'slice_index' : Slice index along the chosen plane (scalar). Default = mid-slice.
%   'brainMask'   : 3-D logical mask; outside voxels are hidden. Default = brain~=0.
%   'plane'       : 'sagittal' | 'coronal' | 'transverse'. Default 'transverse'.
%   'cmap_LNM'    : 256×3 colormap for overlay. Default = cyan→green→blue→red→yellow.
%   'cmap_brain'  : 256×3 colormap for underlay. Default = gray(256).
%   'title'       : Figure title (char). Default ''.
%   'vmin'        : Lower color limit for overlay. Default = min(lnm(:)).
%   'vmax'        : Upper color limit for overlay. Default = max(lnm(:)).
%   'ax'          : Target axes handle to draw into (for subplots). Default = [] → new figure.
%
% Output
%   Creates a figure with brain underlay and semi-transparent LNM overlay.
%
% Example
%   [cmap_lnm, cmap_brain, ~] = get_custom_colormaps();
%   plotMNIslice(LNM_vol, brain_vol, ...
%       'plane','sagittal', 'slice_index', 42, ...
%       'brainMask', brain_mask, ...
%       'cmap_LNM', cmap_lnm, 'cmap_brain', cmap_brain, ...
%       'vmin', -5, 'vmax', 5, 'title', 'LNM — sagittal slice');


%% parse inputs
p = inputParser;
p.addRequired('lnm',    @(x) isnumeric(x) && ndims(x)==3);
p.addRequired('brain',  @(x) isnumeric(x) && ndims(x)==3);
p.addParameter('slice_index', round(size(lnm,3)/2), @(x) isnumeric(x) && isscalar(x));
p.addParameter('brainMask',    [], @(x) islogical(x) || isempty(x));
p.addParameter('plane',        'transverse', @(x) ischar(x));
p.addParameter('cmap_LNM',     [], @(x) isnumeric(x) && size(x,2)==3);
p.addParameter('cmap_brain',   [], @(x) isnumeric(x) && size(x,2)==3);
p.addParameter('title',        '', @ischar);
p.addParameter('vmin',         [], @(x) isnumeric(x) && isscalar(x));
p.addParameter('vmax',         [], @(x) isnumeric(x) && isscalar(x));
p.addParameter('ax',           [], @(h) isempty(h) || (ishandle(h) && strcmp(get(h,'Type'),'axes')));

p.parse(lnm, brain, varargin{:});

slice_idx = p.Results.slice_index;
mask      = p.Results.brainMask;
plane     = lower(p.Results.plane);
cmap_LNM  = p.Results.cmap_LNM;
cmap_brain= p.Results.cmap_brain;
Title     = p.Results.title;
vmin      = p.Results.vmin;
vmax      = p.Results.vmax;
ax        = p.Results.ax;

% defaults
if isempty(mask),      mask = brain~=0;                   end
if isempty(cmap_LNM)   % cyan→green→blue→red→yellow
    base = [0 1 1; 0 1 0; 0 0 1; 1 0 0; 1 1 0];
    cmap_LNM = interp1(linspace(0,1,5)', base, linspace(0,1,256)');
end
if isempty(cmap_brain), cmap_brain = gray(256);            end
if isempty(vmin),      vmin = min(lnm(:));                 end
if isempty(vmax),      vmax = max(lnm(:));                 end

% extract slice according to plane
switch plane
  case 'sagittal'
    % slice at x = slice_idx
    bslice = squeeze(brain(slice_idx,:,:))';
    mslice = squeeze(mask(slice_idx,:,:))';
    lslice = squeeze(lnm(slice_idx,:,:))';
  case {'coronal','frontal'}
    % slice at y = slice_idx
    bslice = squeeze(brain(:,slice_idx,:))';
    mslice = squeeze(mask(:,slice_idx,:))';
    lslice = squeeze(lnm(:,slice_idx,:))';
  otherwise  % transverse
    bslice = squeeze(brain(:,:,slice_idx))';
    mslice = squeeze(mask(:,:,slice_idx))';
    lslice = squeeze(lnm(:,:,slice_idx))';
end


%% plot mni brain with LNM overlay
% mask outside
bslice(~mslice) = NaN;
lslice(~mslice) = NaN;

if isempty(ax)
    fig = figure('Color','w','Position',[100 100 400 400]);
    ax  = axes('Parent',fig,'Position',[0 0 1 1]);
    overlayInTile = false;
else
    fig = ancestor(ax,'figure');
    overlayInTile = isa(get(ax,'Parent'),'matlab.graphics.layout.TiledChartLayout');
end

% Underlay in provided/new axes
imagesc(ax, bslice, 'AlphaData', ~isnan(bslice));
colormap(ax, cmap_brain);
axis(ax,'off','equal');  ax.YDir = 'normal';
hold(ax,'on');

% Overlay axes: same tile if tiledlayout, else same figure/position
if overlayInTile
    ax2 = axes('Parent', get(ax,'Parent'));          % same tiledlayout
    ax2.Layout.Tile     = ax.Layout.Tile;            % SAME TILE
    ax2.Layout.TileSpan = ax.Layout.TileSpan;        % SAME SPAN
else
    ax2 = axes('Parent', fig, 'Units', get(ax,'Units'), 'Position', get(ax,'Position'));
end

imagesc(ax2, lslice, 'AlphaData', 0.7*~isnan(lslice));
colormap(ax2, cmap_LNM);
caxis(ax2, [vmin vmax]);
axis(ax2,'off','equal');  ax2.YDir = 'normal';
set([ax,ax2],'Color','none');

% title
if ~isempty(Title)
    text(ax2, 0.5, 0.97, Title, ...
         'Units','normalized','HorizontalAlignment','center', ...
         'FontSize',12,'Color','k','Interpreter','none');
end
end