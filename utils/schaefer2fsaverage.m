function x = schaefer2fsaverage(value, annot_path, ROIs, hemi, varargin)
    % SCHAEFER2FSAVERAGE Assigns parcel-wise values to the vertices of an fsaverage 
    % surface using a FreeSurfer *.annot file. Only parcels matching the requested 
    % hemisphere are mapped; all other vertices (e.g., medial wall / non-parcel) are 
    % set to NaN by default or to 0 if 'fill_zeros' is true.
    %
    % Inputs
    %   value       : Column vector of parcel values for the requested hemisphere
    %                 (length 500 for Schaefer-1000 per hemi).
    %   annot_path  : Path to the hemisphere-specific annotation file, e.g.
    %                 .../annotation/lh.schaefer1000-yeo7_fs5.annot  or
    %                 .../annotation/rh.schaefer1000-yeo7_fs5.annot
    %   ROIs        : Cell array of Schaefer-1000 parcel names (including hemi tags,
    %                 e.g., '..._LH_...' or '..._RH_...'), in the same order as VALUE.
    %   hemi        : 'lh' or 'rh' (hemisphere to map).
    %
    % Optional
    %   'fill_zeros': Logical. If true, initialize unmapped vertices to 0; otherwise
    %                 to NaN. Default false.
    %
    % Output
    %   x           : [Nvert × 1] vertex-wise array on the given hemisphere
    %                 (Nvert = 10242 for fsaverage5). Vertices of parcels present in
    %                 the annotation receive the corresponding VALUE; all other
    %                 vertices are NaN (or 0 if 'fill_zeros' is true).
    %
    % Example
    %   Right hemisphere: map 500 Schaefer values to fsaverage5 vertices
    %   rh_vals = some_vector_500x1;
    %   rh_annot = fullfile(projectDir,'data','freesurfer_files','annotation',...
    %                       'rh.schaefer1000-yeo7_fs5.annot');
    %   x_rh = schaefer2fsaverage(rh_vals, rh_annot, rh_roi_names, 'rh', 'fill_zeros', true);
    
    %% parse inputs
    p = inputParser;
    addParameter(p, 'fill_zeros', false, @(x) islogical(x) || isempty(x));
    parse(p, varargin{:});
    fill_zeros = p.Results.fill_zeros;
    
    %% convert schaefer parcel to fsaverage
    [vertex, label, colortable] = read_annotation(annot_path, 0);
    num_vertices = length(vertex);

    ROIs = ROIs(:); % Flatten the ROIs array
    if isvector(value)
        value = value(:); % Ensure value is a column vector
    end

    annot_roi_names = colortable.struct_names;

    if fill_zeros
        x = zeros(num_vertices, 1); % Initialize x with zeros values
    else
        x = nan(num_vertices, 1); % Initialize x with NaN values
    end
    
    % val_idx = 1;
    % for roi_indx = 1:length(ROIs)
    %     roi = ROIs{roi_indx};
    % 
    %     if contains(roi, ['_' upper(hemi) '_'])
    %         region_indx = find(strcmp(colortable.struct_names, {roi}));
    %         vertex_indices = find(label == colortable.table(region_indx, 5));
    %         x(vertex_indices) = value(val_idx);
    %         val_idx = val_idx + 1;
    %     end
    % end

    % ACCELERATED VERSION : Precompute ROI hemisphere match and struct name indices
    numROIs = length(ROIs);
    roi_mask = contains(ROIs, ['_' upper(hemi) '_']); % Logical array
    roi_names_filtered = ROIs(roi_mask);              % Only matching ROIs
    
    % Precompute mapping from ROI names to indices in colortable
    [~, region_indices] = ismember(roi_names_filtered, colortable.struct_names);
    
    % Precompute lookup values from colortable.table
    table_values = colortable.table(region_indices, 5);
    
    % Preallocate or ensure 'x' is the correct size
    if ~exist('x', 'var') || length(x) ~= length(label)
        x = zeros(size(label));
    end
    
    % Main loop: lighter now
    for i = 1:length(region_indices)
        vertex_indices = (label == table_values(i));
        x(vertex_indices) = value(i);
    end


end