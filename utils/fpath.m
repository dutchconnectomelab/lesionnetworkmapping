function fpath = fpath(projectDir,fname, varargin)
    % FPATH Recursively locate a file inside a projectDir and return its full path.
    %
    %   fpath = fpath(projectDir, fname)
    %   fpath = fpath(projectDir, fname, 'mustContain', subString)
    %
    %   This helper searches recursively under projectDir for a file named fname.
    %   If multiple matches are found, the first match is returned (with a warning).
    %   Optionally, results can be filtered by requiring that the containing folder
    %   path includes a given substring (e.g., 'fsaverage5', 'fsaverage7', etc.).
    %
    % Inputs:
    %   projectDir   : Path to repository root (char).
    %   fname        : File name to search for (char), e.g., 'rh.inflated'.
    %
    % Optional:
    %   'mustContain': Substring that must appear in the matched file's folder path.
    %                 Default '' (no filtering).
    %
    % Output:
    %   fpath        : Full path (char) to the matched file.
    %
    % Behavior:
    %   • If no match is found → throws error('fpath:NotFound', ...).
    %   • If multiple matches are found → warns and returns the first match.
    %
    % Example
    %   surfPath = fpath(projectDir, 'rh.inflated', 'mustContain', 'fsaverage7');
    %   annotPath = fpath(projectDir, 'rh.schaefer1000-yeo7_fs7.annot');
    %
    %% Parse inputs
    p = inputParser;
    addParameter(p, 'mustContain', '', @ischar);
    parse(p, varargin{:});
    mustContain = p.Results.mustContain;

    % Find project root
    %[projectDir, ~] = setup_project(projectDir);

    % Search recursively
    f = dir(fullfile(projectDir, '**', fname));

    % Filter by mustContain if provided
    if ~isempty(mustContain)
        keep = contains({f.folder}, mustContain);
        f = f(keep);
    end

    % Handle no matches
    if isempty(f)
        error('fpath:NotFound', ...
              'File "%s" not found in repository.', fname);
    elseif numel(f) > 1
        warning('fpath:MultipleFound', ...
                'Multiple matches for "%s". Using first: %s', fname, f(1).folder);
    end

    % Build full path
    fpath = fullfile(f(1).folder, f(1).name);
end