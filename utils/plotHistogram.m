function plotHistogram(data, varargin)
% fancyHistogram - Creates an enhanced histogram with custom styling.
%
% Usage:
%   plotHistogram(data)
%   plotHistogram(data, 'NumBins', 20, 'Color', [0.2 0.6 0.8], 'ShowDensity', true)
%
% Inputs:
%   data         : numeric vector
%   'NumBins'    : number of bins (default 10)
%   'Color'      : RGB vector for bar color (default [0 0.447 0.741])
%   'ShowDensity': true/false, overlay kernel density (default true)

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'NumBins', 10, @isnumeric);
    addParameter(p, 'Color', [0 0.447 0.741], @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'ShowDensity', true, @islogical);
    addParameter(p, 'Title', 'Title', @ischar);
    addParameter(p, 'xLabel', 'Count', @ischar);
    addParameter(p, 'yLabel', 'Frequency', @ischar);
    
    parse(p, varargin{:});
    
    numBins = p.Results.NumBins;
    colorBar = p.Results.Color;
    showDensity = p.Results.ShowDensity;
    Title = p.Results.Title;
    xLabel = p.Results.xLabel;
    yLabel = p.Results.yLabel;
    
    % Create histogram
    h = histogram(data, numBins, 'FaceColor', colorBar, 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'LineWidth', 1.5);
    hold on

    % Overlay kernel density estimate if requested
    if showDensity
        [f, xi] = ksdensity(data);  % Kernel density
        f = f * numel(data) * (h.BinWidth); % Scale to match histogram counts
        plot(xi, f, 'r-', 'LineWidth', 2);
    end

    % Fancy grid and labels
    grid on
    set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12)
    xlabel(xLabel)
    ylabel(yLabel)
    title(Title)
    hold off
end