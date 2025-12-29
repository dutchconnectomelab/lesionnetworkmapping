function plot_rh_surface(data, varargin)
    % PLOT_RH_SURFACE Displays parcel-wise data on the fsaverage right hemisphere 
    % by mapping Schaefer1000 parcel values to vertices and rendering them on
    % rh.inflated. Accepts vectors of length 500 (RH only), 1000 (LH+RH; RH part is used), 
    % or 1054 (Melbourne54 + LH+RH; RH part is used).
    %
    % Input:
    %   data        : Numeric or logical vector of parcel values. Expected sizes:
    %                 • 500  → RH Schaefer1000
    %                 • 1000 → LH+RH Schaefer1000 (uses entries 501:1000)
    %                 • 1054 → Melbourne54+Schaefer1000 (uses last 500)
    %
    % Optional:
    %   'title'     : Figure title (char). Default ''.
    %   'vmin'      : Color lower bound (scalar). Default = min(data).
    %   'vmax'      : Color upper bound (scalar). Default = max(data).
    %   'vcenter'   : Diverging midpoint (scalar). Default = 0 if range spans 0,
    %                 otherwise (vmin+vmax)/2.
    %   'cmap'      : Colormap (N×3). Default = simple red–white diverging (256×3).
    %   'ax'        : Target axes handle to draw into (for subplots). Default = [] → new figure.
    %
    % Output:
    %   Creates a figure with the right-hemisphere surface colored
    %   by the mapped overlay. Optionally displays a colorbar (see plot_surface).
    %
    % Dependencies
    %   • read_surf, read_annotation (FreeSurfer I/O)
    %   • schaefer2fsaverage (parcel→vertex mapping on fsaverage)
    %   • Required files under <projectDir>/data/freesurfer_files/:
    %       - fsaverage/surf/rh.inflated
    %       - annotation/rh.schaefer1000-yeo7.annot
    %
    % Example
    %   plot_rh_surface(LNM(501:end), 'title','LNM — RH', 'vmin',-4, 'vmax',4, ...
    %                   'vcenter',0, 'projectDir', projectDir)

    %% parse inputs
    p = inputParser;
    p.addRequired('data', @(x) (isnumeric(x) || islogical(x)) && isvector(x));
    p.addParameter('title',   '',        @ischar);
    p.addParameter('vmin',    [],        @(x) isempty(x)||isscalar(x));
    p.addParameter('vmax',    [],        @(x) isempty(x)||isscalar(x));
    p.addParameter('vcenter', [],        @(x) isempty(x)||isscalar(x));
    p.addParameter('cmap',    [],        @(x) isnumeric(x)&&size(x,2)==3);
    p.addParameter('cbar',    true,      @islogical);
    p.addParameter('projectDir','',      @ischar);
    p.addParameter('ax',        [],      @(h) isempty(h) || (ishandle(h) && strcmp(get(h,'Type'),'axes')));
    p.parse(data, varargin{:});
    
    Title      = p.Results.title;
    vmin       = p.Results.vmin;
    vmax       = p.Results.vmax;
    vcenter    = p.Results.vcenter;
    cmap       = p.Results.cmap;
    projectDir = p.Results.projectDir;
    ax         = p.Results.ax;
    cbar       = p.Results.cbar;
    % set defaults
    if isempty(vmin);     vmin = min(data(:)); end
    if isempty(vmax);     vmax = max(data(:)); end
    if isempty(vcenter)
        if vmin<0 && vmax>0
            vcenter = 0;
        else
            vcenter = (vmin+vmax)/2;
        end
    end
    
    if isempty(cmap)
        % simple red‐white cmap
        nColors = 256;
        % cmap = [ones(n,1) linspace(1,0,n)' linspace(1,0,n)'];

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
        
        cmap = [ ...
            interp1(t_ctrl, ctrl(:,1), t, 'linear')' ...
            interp1(t_ctrl, ctrl(:,2), t, 'linear')' ...
            interp1(t_ctrl, ctrl(:,3), t, 'linear')' ];
    end

    data = data(:);
    nParcels = size(data, 1);
    if nParcels==500
        fprintf('[INFO] %d parcels → assuming Schaefer-1000 rh for plotting rh surface.\n',nParcels);
    elseif nParcels==1000
        fprintf('[INFO] %d parcels → assuming lh+rh (Schaefer-1000) for plotting rh surface.\n',nParcels);
        data = data(501:end);
    elseif nParcels==1054
    %    fprintf('[INFO] %d parcels → assuming Melbourne54+lh+rh (Schaefer-1000) for plotting rh surface.\n',nParcels);
        data = data(end-499:end);
    else
        warning('[WARNING] Got %d parcels; expected 500, 1000, or 1054.',nParcels);
    end

    %% plot surface
    % load surfaces
    [rh_ver, rh_faces] = read_surf(fpath(projectDir,'rh.inflated', 'mustContain', 'fsaverage7'));
    
    % load annotations
    rh_annot_path = fpath(projectDir,'rh.schaefer1000-yeo7_fs7.annot');
    
    [~,~,rh_ct] = read_annotation(rh_annot_path, 0);
    
    % skip medial wall
    rh_rois = rh_ct.struct_names(2:end);
    
    % project to surface
    rh_data = schaefer2fsaverage(data, rh_annot_path, rh_rois, 'rh');
    
    % plot surface (pass ax through)
    plot_surface(rh_ver, rh_faces, rh_data, 'cmap', cmap, 'vmin', vmin, 'vmax', vmax, ...
                 'vcenter', vcenter, 'azimuth', 0, 'elevation', 90, 'ax', ax, 'cbar', cbar)
    
    if ~isempty(Title)
        if ~isempty(ax) && ishandle(ax)
            title(ax, Title, 'FontSize',12);
        else
            title(Title, 'FontSize',12);
        end
    end
end

function plot_surface(v, f, o, varargin)
    % medial wall view azimuth=0, elevation=90
    % lateral view azimuth=0, elevation=270
    
    %% parse inputs
    p = inputParser;
    p.addParameter('cmap',      'jet', @(x) (ischar(x) || (isnumeric(x)&&size(x,2)==3)));
    p.addParameter('cbar',      true,  @islogical);
    p.addParameter('hemi',      'lh',  @ischar);
    p.addParameter('vmin',      [],    @(x) isempty(x)||isscalar(x));
    p.addParameter('vmax',      [],    @(x) isempty(x)||isscalar(x));
    p.addParameter('vcenter',   [],    @(x) isempty(x)||isscalar(x));
    p.addParameter('reduce_res',false, @islogical);
    p.addParameter('azimuth',   0,     @isscalar);
    p.addParameter('elevation', 270,   @isscalar);
    p.addParameter('ax',        [],    @(x) isempty(x)||ishandle(x));
    p.addParameter('nan',false, @islogical);
    p.parse(varargin{:});
    args = p.Results;

    % set defaults (and note use of 'o' instead of nonexistent 'values')
    if isempty(args.vmin),       args.vmin    = min(o);       end
    if isempty(args.vmax),       args.vmax    = max(o);       end
    if isempty(args.vcenter)
        if args.vmin<0 && args.vmax>0
            args.vcenter = 0;
        else
            args.vcenter = (args.vmin + args.vmax)/2;
        end
    end

    % prepare axes
    if isempty(args.ax)
        fig = figure('Color','w');
        ax = axes('Parent', fig);
    else
        ax = args.ax;
        axes(ax);
    end
        
    % set range
    v = (v - min(v)) ./ (max(v) - min(v));
    if isempty(args.vmin), args.vmin = min(o); end
    if isempty(args.vmax), args.vmax = max(o); end
    if isempty(args.vcenter), args.vcenter = (args.vmin + args.vmax) / 2; end

    %%  plot surface with the modified overlay
    surf = trisurf(f+1, v(:, 1), v(:, 2), v(:, 3), 'EdgeColor', 'none');

    % Replace NaN values with a specific value beyond the range of the data
    if args.nan
        range_add = abs(0.1 * (args.vmax - args.vmin));
        nan_mask = isnan(o);
        o(nan_mask) = args.vmax + range_add; % Setting NaN values to a value beyond the vmax
        
        % Set vertex colors
        set(surf, 'FaceVertexCData', o, 'FaceColor', 'interp');
        colormap(ax, [args.cmap; 1 1 1]); % Append white color for NaN values
        caxis(ax, [args.vmin, args.vmax + range_add]);

    else
        set(surf, 'FaceVertexCData', o, 'FaceColor', 'interp');
        colormap(ax, args.cmap);
        caxis(ax, [args.vmin, args.vmax]); 
    end
    
    % Create a LightSource
    shading interp
    lightangle(args.elevation, args.azimuth)
    surf.FaceLighting = 'gouraud';
    surf.AmbientStrength = 0.4;
    surf.DiffuseStrength = 0.8;
    surf.SpecularStrength = 0.1;
    
    % Set the limits of the plot from the data
    if args.cbar
        colorbar(ax);
    end
    
    xlim([min(v(:, 1)), max(v(:, 1))]);
    ylim([min(v(:, 2)), max(v(:, 2))]);
    zlim([min(v(:, 3)), max(v(:, 3))]);
    axis off;
    view(ax, args.elevation, args.azimuth);
end
