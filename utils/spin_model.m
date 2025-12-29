function X_spin = spin_model(X, varargin)

    % SPIN_MODEL Generates spin-permuted lesion maps while preserving lesion 
    % prevalence. Each subject’s lesion vector (Schaefer1000 parcels: LH 1–500, RH 501–1000)
    % is mapped to the fsaverage5 surface, left and right hemispheres are rotated 
    % independently via the BrainSpace spherical spin procedure using the SAME rotation for
    % all subjects per permutation, and values are averaged back into Schaefer1000 parcels. 
    % Output maps are not binarized; parcel values reflect the mean of rotated vertex values 
    %
    % Input:
    %   X        : [nVectors × 1000] stacked Schaefer1000 vectors (LH parcels 1:500, RH parcels 
    %                                501:1000).
    %
    % Optional:
    %   'nPerm'      : Number of spin permutations (default 1).
    %   'random_number_generator' : Scalar seed for RNG. If NaN, RNG state is
    %                  left unchanged (default NaN).
    %
    % Output:
    %   M_spin   : [nPerm × nVectors × 1000] array of spin-permuted lesion maps.
    %              Note: if nPerm == 1, the output is squeezed to [nVectorsects × 1000].
    %              Values are parcel-wise means in [0,1]; no post-hoc binarization is applied.
    %
    % Dependencies
    %   • BrainSpace toolbox: spin_permutations, convert_surface
    %   • MATLAB statistics and machine learning toolbox: KDTreeSearcher
    %   • FreeSurfer-compatible readers: read_surf, read_annotation.
    %   • Helper: schaefer2fsaverage (maps Schaefer1000 parcels to fsaverage5 vertices).
    %   • Required files under <projectDir>/lnm_toolbox_data/freesurfer_files/:
    %       - fsaverage5/surf/lh.sphere, rh.sphere
    %       - annotation/lh.schaefer1000-yeo7_fs5.annot, rh.schaefer1000-yeo7_fs5.annot
    %
    % Reference
    %   Vos de Wael et al., BrainSpace: a toolbox for the analysis of macroscale gradients
    %   in neuroimaging and connectomics. Communications Biology (2020).
    %   https://doi.org/10.1038/s42003-020-0794-7
    %
    %
    % Example
    %   X_spin = spin_permute_lesions(X, 'nPerm', 1000, 'random_number_generator', 42);
    %
    
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
    method = validatestring(p.Results.method, {'spin_same_angle','spin_independent', 'fix_seed_count'});

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

    % load sphere & create struc for brainspace spin permute function
    [lh_ver, lh_faces] = read_surf(fpath(projectDir,'lh.sphere', 'mustContain', 'fsaverage5'));
    [rh_ver, rh_faces] = read_surf(fpath(projectDir,'rh.sphere', 'mustContain', 'fsaverage5'));

    lh_sphere.vertices = lh_ver;
    lh_sphere.faces = lh_faces;

    rh_sphere.vertices = rh_ver;
    rh_sphere.faces = rh_faces;

    % read annotation files
    lh_annot_path = fpath(projectDir, 'lh.schaefer1000-yeo7_fs5.annot');
    rh_annot_path = fpath(projectDir, 'rh.schaefer1000-yeo7_fs5.annot');

    [~, lh_annot_labels, lh_ct] = read_annotation(lh_annot_path, 0);
    [~, rh_annot_labels, rh_ct] = read_annotation(rh_annot_path, 0);

    % get freesurfer identifiers
    lh_annot_id = lh_ct.table(2:end, 5); % skip medial wall
    rh_annot_id = rh_ct.table(2:end, 5); % skip medial wall
    
    % get ROIs name + skip medial wall
    lh_rois = lh_ct.struct_names(2:end);
    rh_rois = rh_ct.struct_names(2:end);

    % precompute mapping vertex to schaefer
    [~,grpL] = ismember(lh_annot_labels, lh_annot_id);
    maskL = grpL>0; colsL = find(maskL); rowsL = grpL(maskL);
    cntL  = accumarray(rowsL,1,[numel(lh_annot_id),1]);
    Wlh   = sparse(rowsL, colsL, 1./cntL(rowsL), numel(lh_annot_id), numel(lh_annot_labels));
    
    [~,grpR] = ismember(rh_annot_labels, rh_annot_id);
    maskR = grpR>0; colsR = find(maskR); rowsR = grpR(maskR);
    cntR  = accumarray(rowsR,1,[numel(rh_annot_id),1]);
    Wrh   = sparse(rowsR, colsR, 1./cntR(rowsR), numel(rh_annot_id), numel(rh_annot_labels));
    
    % initialize permuted lesions matrix
    X_spin = zeros(nPerm, nVectors, nParcels);

    if verbose
        fprintf('\nStart spin permuting %d atlas-based vectors for %d spins.\n', nVectors, nPerm);
    end
    
    if strcmp(method, 'spin_independent')
       
        % generate random seed per vector
        vectorSeeds = randi(intmax('int32'), [nVectors 1]);

        for i = 1:nVectors

            if verbose
                if mod(i, max(1,floor(nVectors/10)))==0
                    fprintf(' %d%%', round(100*i/nVectors));
                end
            end
   
            % get schaefer1000 vector
            x = X(i, :);
            
            % project schaefer1000 values to fsaverage5 surface

            x_lh = schaefer2fsaverage(x(1:500), lh_annot_path, lh_rois, 'lh', 'fill_zeros', true);
            x_rh = schaefer2fsaverage(x(501:end), rh_annot_path, rh_rois, 'rh', 'fill_zeros', true);
    
            % apply rotation
            x_spin = spin_permutations({x_lh, x_rh}, {lh_sphere, rh_sphere}, nPerm, 'random_state', vectorSeeds(i));
            
            x_spin{1} = squeeze(x_spin{1});
            x_spin{2} = squeeze(x_spin{2});

            for iPerm=1:nPerm
                x_lh_rot = x_spin{1}(:, iPerm);
                x_rh_rot = x_spin{2}(:, iPerm);
                
                % map back to schaefer1000 parcels
                x_vals = [Wlh * x_lh_rot; Wrh * x_rh_rot];
                X_spin(iPerm, i, :) = x_vals; 
            end
        end

    elseif strcmp(method, 'spin_same_angle')

        % define fsaverage vertex indices
        idx_lh = (1:10242)'; % 10242 vertices on fsaverage 5 surface
        idx_rh = (1:10242)';

        % spin permute the vertex indices
        x_spin = spin_permutations({idx_lh, idx_rh}, {lh_sphere, rh_sphere}, nPerm, 'random_state', random_number_generator, 'verbose', verbose_spin);
        
        x_spin{1} = squeeze(x_spin{1});
        x_spin{2} = squeeze(x_spin{2});
        
        for i = 1:nVectors
   
            % get schaefer1000 vector
            x = X(i, :);
            
            % project schaefer1000 values to fsaverage5 surface
            x_lh = schaefer2fsaverage(x(1:500), lh_annot_path, lh_rois, 'lh', 'fill_zeros', true);
            x_rh = schaefer2fsaverage(x(501:end), rh_annot_path, rh_rois, 'rh', 'fill_zeros', true);
    
            for iPerm=1:nPerm
    
                % get rotated indices
                idx_lh = x_spin{1}(:, iPerm);
                idx_rh = x_spin{2}(:, iPerm);
           
                % apply rotation
                x_lh_rot = x_lh(idx_lh);
                x_rh_rot = x_rh(idx_rh);
                
                % map back to schaefer1000 parcels
                x_vals = [Wlh * x_lh_rot; Wrh * x_rh_rot];
        
                X_spin(iPerm, i, :) = x_vals; 
            end
        end
        
    elseif strcmp(method, 'fix_seed_count')
            
            % make sure X is binary
            n_multiplier = 5;
            max_iters = n_multiplier*nPerm;
            idx_lh = (1:10242)'; % 10242 vertices on fsaverage 5 surface
            idx_rh = (1:10242)';

            % generate random seed per vector
            vectorSeeds = randi(intmax('int32'), [nVectors 1]);
        
            for i = 1:nVectors

                if verbose
                    if mod(i, max(1,floor(nVectors/10)))==0
                        fprintf(' %d%%', round(100*i/nVectors));
                    end
                end

                x_spin = spin_permutations({idx_lh, idx_rh}, {lh_sphere, rh_sphere}, max_iters, 'random_state', vectorSeeds(i));
                x_spin{1} = squeeze(x_spin{1});
                x_spin{2} = squeeze(x_spin{2});
        
                % get lesion and lesion size
                x = X(i, :);
                num_seeds = sum(x);
                
                % convert schaefer to fsaverage5 surface
                x_lh = schaefer2fsaverage(x(1:500), lh_annot_path, lh_rois, 'lh', 'fill_zeros', true);
                x_rh = schaefer2fsaverage(x(501:end), rh_annot_path, rh_rois, 'rh', 'fill_zeros', true);
                
                jPerm = 1;
                for iPerm=1:nPerm

                    % spin permute lesion untill the spun lesion matches the seed‐count
                    for iter = 1:max_iters
                        
                        % get rotated indices
                        idx_lh = x_spin{1}(:, jPerm);
                        idx_rh = x_spin{2}(:, jPerm);
            
                        % apply rotation
                        x_lh_rot = x_lh(idx_lh);
                        x_rh_rot = x_rh(idx_rh);
                        
                        % map back to schaefer1000 parcels
                        x_vals = [Wlh * x_lh_rot; Wrh * x_rh_rot];
                        
                        % binarize & match seed‐count
                        sel = x_vals > 0;
                        cnt = sum(sel);
                        
                        jPerm = jPerm + 1;

                        if cnt > num_seeds
                          % too many → keep top num_seeds highest‐valued parcels
                          thr = sort(x_vals, 'descend');
                          thr = thr(num_seeds);
                          sel = (x_vals >= thr);
                        
                        elseif cnt < num_seeds 
                          % too few → grow one‐per‐contiguous‐run
                          % continue
                          need = num_seeds - cnt;
                          padded = [0; sel; 0];
                          starts = find(padded(1:end-1)==0 & padded(2:end)==1);
                          ends   = find(padded(1:end-1)==1 & padded(2:end)==0) - 1;
                        
                          candidates = [];
                          for j = 1:numel(starts)
                            s = starts(j); e = ends(j);
                            neigh = [];
                            if s-1 >= 1        && ~sel(s-1),    neigh(end+1)=s-1; end
                            if e+1 <= 1000     && ~sel(e+1),    neigh(end+1)=e+1; end
                            if ~isempty(neigh)
                              % pick the neighbor with highest m_vals
                              [~, best] = max(x_vals(neigh));
                              candidates(end+1) = neigh(best);
                            end
                          end
                        
                          if numel(candidates) < need
                            continue  % skip this spin, try next
                          end
                        
                          % add the top 'need' candidates
                          [~, ord] = sort(x_vals(candidates), 'descend');
                          to_add = candidates(ord(1:need));
                          sel(to_add) = true;
                        end
                        
                        % check for success
                        if sum(sel) == num_seeds
                          X_spin(iPerm, i, :) = double(sel);
                          break
                        end
                    end
                    if iter == max_iters
                        error( ...
                          'spin_model:NoConvergence', ...
                          ['Subject %d: failed to match %d lesion parcels ', ...
                           'in %d spins. ', ...
                           'Consider increasing ''max_iters'', or pick different method.'], ...
                          iSub, num_seeds, max_iters);
                    end
                end
            end

    else
         error('spin_model:UnknownMethod', ...
            ['Unknown method "%s". Valid options are: "spin_same_angle" (all subjects share ' ...
             'the same rotation per permutation) and "spin_independent" (each subject has ' ...
             'independent rotations per permutation).'], method);

    end

    X_spin = squeeze(X_spin);
end
