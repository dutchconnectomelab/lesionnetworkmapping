function mat = map_voxel2schaeferMelbourne(data, varargin)
    % MAP_VOXEL2SCHAEFERMELBOURNE Computes parcel-wise mean values from a voxelwise image using 
    % a combined Melbourne-54 + Schaefer-1000 parcellation in MNI152 2 mm RAS space. 
    % The output vector follows the parcel order defined by the CTAB (filtered to IDs
    % present in the atlas; medial wall is excluded).
    %
    % Input:
    %   data        : 3-D numeric array (voxelwise image) already aligned to the
    %                 atlas grid (MNI152 2 mm RAS). Orientation must match the atlas
    %                 file used below; resample/reorient beforehand if needed.
    %
    % Optional:
    %   'projectDir' : Path to repository root containing the parcellation and ctab files. 
    %                  Default '' (relative paths).
    %
    % Output:
    %   mat         : 1 × N vector of parcel means (double), where N is the number
    %                 of nonzero parcel IDs present in the atlas. Parcels with no
    %                 in-atlas voxels receive NaN.
    %
    % Dependencies
    %   • load_nifti (to read the atlas NIfTI)
    %
    % Example
    %   vol = load_nifti(fullfile(projectDir,'data','degree','GSP1000_degree.nii.gz')).vol;
    %   x   = map_voxel2atlas(vol, 'projectDir', projectDir);   % 1×N parcel means
    
    %% parse inputs
    p = inputParser;
    addParameter(p, 'projectDir', '', @ischar);
    addParameter(p, 'ctab', '', @istable);
    addParameter(p, 'atlas_vol', '', @isnumeric);
    
    parse(p, varargin{:});
    projectDir = p.Results.projectDir;
    
    parse(p, varargin{:});
    ctab = p.Results.ctab;
    
    parse(p, varargin{:});
    atlas_vol = p.Results.atlas_vol;
    
    % Sanity check: input data must match atlas grid
    assert(ndims(data) == 3, 'map_voxel2atlas:InputNot3D', ...
        'Input "data" must be a 3-D volume. Got ndims(data) = %d.', ndims(data));
    
    assert(isequal(size(data), size(atlas_vol)), 'map_voxel2atlas:SizeMismatch', ...
        ['Input "data" size (%s) does not match atlas size (%s). ' ...
         'Resample/reorient data to the atlas grid (MNI152 2mm RAS) before calling map_voxel2atlas.'], ...
        mat2str(size(data)), mat2str(size(atlas_vol)));
    
    % Warn if no projectDir was provided
    if isempty(projectDir)
        warning(...
          'load_LNM:NoProjectDir', ...
          ['No projectDir given, using relative paths.\n' ...
           'If this fails, please call map_voxel2schaeferMelbourne(''projectDir'',<yourPath>)']);
    end
    
    %% map voxel-wise map to atlas regions
    % Identify parcels present in the atlas (exclude ID=0)
    parcelIDs = unique(atlas_vol(:));
    parcelIDs(parcelIDs==0) = [];
    
    % Filter CTAB rows to those IDs
    mask = ismember(ctab.ID,parcelIDs);
    ctab = ctab(mask,:);
    
    % Prepare output
    nRegions = height(ctab);
    mat      = nan(1,nRegions);
    
    % % Compute mean for each parcel, preserving CTAB order
    try % we notice that some older MATLAB version on WINDOWS do not have accumarray, use catch then
        % fast implementation
        ids  = ctab.ID(:);                   % region IDs (nRegions×1)
        idx  = atlas_vol(:);                 % flatten atlas
        vals = data(:);                      % flatten data
        
        % remove NaNs
        nanmask = ~isnan(vals);
        idx  = idx(nanmask);
        vals = vals(nanmask);
        
        % mask out background (IDs <= 0)
        mask = idx > 0;
        idx  = idx(mask);
        vals = vals(mask);
        
        % now idx contains only positive integers
        vals = double(vals);
        meansAll = accumarray(idx, vals, [], @mean, NaN);
        
        % extract the regions you want in ctab
        mat = meansAll(ctab.ID)';   % mat(i) = mean over region i
    catch % we notice that some older MATLAB version on WINDOWS do not have accumarray
        for i = 1:nRegions
            rid = ctab.ID(i);
            mask_vox = (atlas_vol == rid);
            if any(mask_vox(:))
                vals      = data(mask_vox);
                mat(i)    = mean(vals(~isnan(vals)));
            else
                mat(i) = NaN;
            end
        end
    end

end