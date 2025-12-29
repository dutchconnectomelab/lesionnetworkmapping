function [custom_cmap, cmap_brain, cmap_reds, cmap_hot] = get_custom_colormaps()
% GET_CUSTOM_COLORMAPS  Returns three 256×3 colormaps tailored for this project:
%     • custom_cmap (LNM): multi-hue palette (cyan→green→blue→red→yellow) interpolated to 256 steps.
%     • cmap_brain (structural): grayscale for anatomical backgrounds.
%     • cmap_reds (prevalence/lesions): sequential white→red for proportions/rates.
%
% Input:
%       None
%
% Optional:
%       None
%
% Output:
%   custom_cmap : 256×3 double, multi-hue colormap for LNM surfaces/volumes.
%   cmap_brain  : 256×3 double, grayscale colormap for structural images.
%   cmap_reds   : 256×3 double, white-to-red sequential colormap for prevalence.
%
% Example
%   [cmap_lnm, cmap_brain, cmap_reds] = get_custom_colormaps();

% Define the custom colormap
custom_cmap = [
    0   1   1   % cyan
    0   1   0   % green
    0   0   1   % blue
    1   0   0   % red
    1   1   0   % yellow
];

nColors = 256;
custom_cmap = interp1(linspace(0,1,size(custom_cmap,1)), custom_cmap, linspace(0,1,nColors));

% define brain colormap 
cmap_brain = gray(256);

% simple red‐white diverging
% cmap_reds = [ones(nColors,1) linspace(1,0,nColors)' linspace(1,0,nColors)'];
ctrl = [ ...
    255 245 240
    254 224 210
    252 187 161
    252 146 114
    251 106  74
    239  59  44
    203  24  29
    165  15  21
    103   0  13] / 255;

t_ctrl = linspace(0,1,size(ctrl,1));
t      = linspace(0,1,nColors);

cmap_reds = [ ...
    interp1(t_ctrl, ctrl(:,1), t, 'linear')' ...
    interp1(t_ctrl, ctrl(:,2), t, 'linear')' ...
    interp1(t_ctrl, ctrl(:,3), t, 'linear')' ];

cmap_hot = hot(256);

end