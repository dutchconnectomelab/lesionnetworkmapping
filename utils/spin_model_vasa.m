function X_spin = spin_model_vasa(X, varargin)

    % SPIN_MODEL_VASA Generates spin-permuted Schaefer1000 maps using the
    % parcel centroid rotation method of František Váša
    % (https://github.com/frantisekvasa/rotate_parcellation).
    %
    % Each map (Schaefer1000 parcels: LH 1–500, RH 501–1000) is rotated on the fsaverage5
    % spherical surface via hemisphere-specific rotations derived from parcel centroids.
    % The same rotation is applied across subjects for each permutation (default behavior),
    % or independent rotations are used if 'method' is set to 'spin_independent'.
    %
    % Input:
    %   X        : [nVectors × 1000] stacked Schaefer1000 vectors (LH parcels 1:500, RH parcels 
    %                                501:1000).
    % Optional:
    %   'nPerm'      : Number of spin permutations (default 1).
    %   'random_number_generator' : Scalar seed. If NaN, uses rng('shuffle') (default NaN).
    %   'method'     : 'spin_same_angle' (default) or 'spin_independent'.
    %
    % Output:
    %   X_spin   : [nPerm × nVectors × 1000] array of spin-permuted maps.
    %              If nPerm == 1, the output is squeezed to [nVectors × 1000].
    %
    % Dependencies
    %   • rotate_parcellation (František Váša): github.com/frantisekvasa/rotate_parcellation
    %   • centroid_extraction_sphere (František Váša): github.com/frantisekvasa/rotate_parcellation
    %   • read_surf, read_annotation (FreeSurfer-compatible)
    %   • Schaefer1000 annotation files on fsaverage5 surface
    %
    % Reference
    %   Váša F. rotate_parcellation. GitHub repository.
    %   https://github.com/frantisekvasa/rotate_parcellation
    %
    % Example
    %   X_spin = spin_model_vasa(X, 'nPerm', 1000, 'random_number_generator', 42);
    
    %% parse inputs
    p = inputParser;
    addParameter(p, 'projectDir', nan);
    addParameter(p, 'random_number_generator', nan);
    addParameter(p, 'method', 'spin_same_angle');  
    addParameter(p, 'nPerm', 1,  @isnumeric);
    addParameter(p, 'verbose', false,  @islogical);
    addParameter(p, 'verbose_spin', false,  @islogical);

    % set defaults
    parse(p, varargin{:});
    
    if ~isnan(p.Results.random_number_generator)
        rng('default');
        random_number_generator = p.Results.random_number_generator;
        rng(random_number_generator);
    else
        rng('shuffle');
        random_number_generator = NaN;
    end

    projectDir=p.Results.projectDir;
    nPerm = p.Results.nPerm;
    verbose = p.Results.verbose;
    verbose_spin = p.Results.verbose_spin;
    method = validatestring(p.Results.method, {'spin_same_angle','spin_independent'});


    if size(X, 2) ~= 1000
        error('spin_model:InputSizeMismatch', ...
            ['Input X has %d columns, but this function expects 1000 ' ...
             '(Schaefer1000 parcels). Please make sure to only select ' ...
             'the Schaefer1000 parcels before calling this function.'], ...
             size(X, 2));
    end


    %% permutation test

    % get number of subjects and parcels
    [nVectors, nParcels] = size(X);

    lh_annot_path = fpath(projectDir, 'lh.schaefer1000-yeo7_fs5.annot');
    rh_annot_path = fpath(projectDir, 'rh.schaefer1000-yeo7_fs5.annot');
    
    lh_sphere_path = fpath(projectDir,'lh.sphere', 'mustContain', 'fsaverage5');
    rh_sphere_path = fpath(projectDir,'rh.sphere', 'mustContain', 'fsaverage5');
    
    coordsL = centroid_extraction_sphere(lh_sphere_path, lh_annot_path);  % 500×3
    coordsR = centroid_extraction_sphere(rh_sphere_path, rh_annot_path);  % 500×3
    
    % initialize permuted lesions matrix
    X_spin = zeros(nPerm, nVectors, nParcels);

    if verbose
        fprintf('\nStart spin permuting %d atlas-based vectors.', nVectors);
    end

    if strcmp(method, 'spin_independent')

        for i = 1:nVectors
            if verbose
                if mod(i, max(1,floor(nVectors/10)))==0
                    fprintf(' %d%%', round(100*i/nVectors));
                end
            end
                
            % parcel centroid rotation
            rot = rotate_parcellation(coordsL(2:end, :), coordsR(2:end, :), nPerm, 'verbose_spin', verbose_spin); 

            % get schaefer1000 vector
            x = X(i, :);
            
            % apply rotation
            X_spin(:, i, :) = x(rot)'; 

        end

    elseif strcmp(method, 'spin_same_angle')

        % parcel centroid rotation
        rot = rotate_parcellation(coordsL(2:end, :), coordsR(2:end, :), nPerm, 'verbose_spin', verbose_spin); 

        for i = 1:nVectors
            if verbose
                if mod(i, max(1,floor(nVectors/10)))==0
                    fprintf(' %d%%', round(100*i/nVectors));
                end
            end

            % get schaefer1000 vector
            x = X(i, :);
            
            % apply rotation
            X_spin(:, i, :) = x(rot)'; 

        end

    end

    X_spin = squeeze(X_spin);
end