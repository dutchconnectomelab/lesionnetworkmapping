function M_random = random_permute_lesions(M, varargin)
    % RANDOM_PERMUTE_LESIONS  Generates a lesion null by, for each subject, 
    % selecting uniformly at random the same number of parcels as in the original 
    % lesion vector (seed count), without replacement, over all Schaefer-1000 parcels
    % (N=1000). Spatial structure is not preserved; only the per-subject count is matched.
    %
    % Input:
    %   M        : [num_subjects × 1000] lesion matrix over Schaefer1000 parcels
    %              (binary 0/1, LH parcels 1:500, RH parcels 501:1000).
    %
    % Optional:
    %   'nPerm'      : Number of spin permutations (default 1).
    %   'random_number_generator' : Scalar seed for RNG. If NaN, RNG state is
    %                  left unchanged (default NaN).
    % Output:
    %   M_random  : [num_subjects × 1000] binary matrix of randomly permuted
    %               lesions with exactly the original seed count per subject.
    %
    % Example
    %   M_null = random_permute_lesions(M, 'random_number_generator', 123);

    %% parse inputs
    p = inputParser;
    addParameter(p, 'nPerm', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(p, 'random_number_generator', nan, ...
        @(x) (isnumeric(x) && isscalar(x)) || (isnan(x)));
    addParameter(p, 'verbose', false, @islogical);
    parse(p, varargin{:});

    nPerm = double(p.Results.nPerm);
    seed  = p.Results.random_number_generator;
    verbose = p.Results.verbose;

    %% random permute lesion
    [num_subj, n_parcels] = size(M);

    % Precompute per-subject seed counts (treat any nonzero as seeded)
    % If M is non-binary, this still counts "active" parcels as >0.
    seed_counts = sum(M > 0, 2);

    if any(seed_counts > n_parcels)
        error('random_permute_lesions:InvalidSeedCount', ...
              'A subject requests %d seeds but only %d parcels exist.', ...
              max(seed_counts), n_parcels);
    end

    % preallocate
    M_random = false(nPerm, num_subj, n_parcels);

    % generate permutations
    for ip = 1:nPerm

        % one permutation
        for iSub = 1:num_subj
            k = seed_counts(iSub);
            if k == 0
                continue; 
            end

            % sample k unique parcel indices
            idx = randperm(n_parcels, k);
            M_random(ip, iSub, idx) = true;
        end

    end

    % squeeze if single permutation to match your spin_model behavior
    if nPerm == 1
        M_random = squeeze(M_random);
    end
end
