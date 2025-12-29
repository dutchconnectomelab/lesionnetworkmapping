function corr_degree = plotCorrelation(x, y, varargin)
% PLOTCORRELATION Plots a scatter of paired observations (x vs y), computes the 
% Pearson correlation coefficient, and annotates the figure with r. Any pair where 
% x or y is NaN or exactly 0 is removed.
%
% Inputs:
%   x, y      : Numeric vectors or arrays of equal length .
%
% Optional:
%   'xlabel'  : X-axis label (char). Default ''.
%   'ylabel'  : Y-axis label (char). Default ''.
%   'title'   : Figure title (char). Default 'Correlation scatter plot'.
%   's'       : Size of points in scatter plot
%   'ax'      : Target axes handle to draw into (for subplots). Default = [] → new figure.
%
% Output
%   Creates a scatter plot of x versus y with the correlation value in the corner.

%% parse inputs
p = inputParser;
p.addRequired('x', @isnumeric);
p.addRequired('y', @isnumeric);
p.addParameter('s', 1, @isnumeric);
p.addParameter('xlabel', '', @ischar);
p.addParameter('ylabel', '', @ischar);
p.addParameter('title',  'Correlation scatter plot', @ischar);
p.addParameter('ax', [], @(h) isempty(h) || (ishandle(h) && strcmp(get(h,'Type'),'axes')));
p.parse(x, y, varargin{:});

s         = p.Results.s;
labelX    = p.Results.xlabel;
labelY    = p.Results.ylabel;
Title     = p.Results.title;
ax        = p.Results.ax;

assert(numel(x) == numel(y), 'plotCorrelation:InputSizeMismatch', ...
        'x and y must have the same length.');

%% compute & plot correlation
x = x(:);  y = y(:);
mask = isnan(x) | isnan(y) | x==0 | y==0;
x = x(~mask); y = y(~mask);

if isempty(x) || isempty(y)
    fprintf('[INFO] All entries were masked out → no valid data left. Please check the input data\n');
    corr_degree = NaN;
    return;
end

corr_degree = corr(x, y);

% choose/create target axes; keep old behavior when no 'ax' is provided
if isempty(ax)
    fig = figure('Color','w', 'Position',[100 100 750 750]);
    ax  = axes('Parent',fig);
end
axes(ax); %#ok<LAXES>

% scatter plot
h = scatter(x, y, s, 'k', 'filled'); hold(ax,'on');
h.SizeData = 20;

% customize axes appearance
set(ax,'Box','off','FontSize',12,'TickLength',[0.01 0.01]);
set(ax,'XAxisLocation','bottom','YAxisLocation','left');

% annotate correlation
text(0.05, 0.90, sprintf('r = %.2f', corr_degree), ...
    'Units','normalized', 'FontSize', 12, 'VerticalAlignment', 'top', 'Parent', ax);

% labels and title
if ~isempty(labelX), xlabel(ax,labelX, 'FontSize', 12); end
if ~isempty(labelY), ylabel(ax,labelY, 'FontSize', 12); end
if ~isempty(Title),  title(ax, Title, 'FontSize', 20, 'FontWeight', 'normal'); end


end
