function [lh_vertex_val, rh_vertex_val] = mri_vol2surf(input_nii, freesurferDir, varargin)
    % MRI_VOL2SURF  Project a volume in MNI152 space onto FreeSurfer fsaverage surface.
    %
    %   [lh_vertex_val, rh_vertex_val] = mri_vol2surf(input_nii, freesurferDir, ...)
    %
    % Inputs:
    %   input_nii      : Path to input .nii or .nii.gz file (in MNI152 space)
    %   freesurferDir  : Path to FreeSurfer installation (should contain 'bin' directory)
    %
    % Optional:
    %   'fsaverage'  : Which fsaverage template to use ('5' or '7', or the full 
    %                  names 'fsaverage5', 'fsaverage7', 'fsaverage'). 
    %                  Default '5'.
    %   'projfrac'   : Projection fraction (default = 0.5)
    %
    % Outputs:
    %   lh_vertex_val  : Vector of surface values (left hemisphere)
    %   rh_vertex_val  : Vector of surface values (right hemisphere)
    %
    
    %% Parse input
    p = inputParser;
    addRequired(p, 'input_nii', @(x) ischar(x) || isstring(x));
    addRequired(p, 'freesurferDir', @(x) ischar(x) || isstring(x));
    addParameter(p, 'fsaverage', '5', @ischar);
    addParameter(p, 'projfrac', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    parse(p, input_nii, freesurferDir, varargin{:});
    
    if strcmp(p.Results.fsaverage, '7') || strcmpi(p.Results.fsaverage, 'fsaverage7') || strcmpi(p.Results.fsaverage, 'fsaverage')
        fsaverage = 'fsaverage';
    elseif strcmp(p.Results.fsaverage, '5') || strcmpi(p.Results.fsaverage, 'fsaverage5')
        fsaverage = 'fsaverage5';
    else
        error('mri_vol2surf:InvalidFsaverage', ...
            'Only fsaverage5 and fsaverage7 are supported. Got: %s', p.Results.fsaverage);
    end
    
    projfrac = num2str(p.Results.projfrac);
    
    %%
    
    mri_vol2surf_bin = fullfile(freesurferDir, 'bin', 'mri_vol2surf');
    
    % make a tmp directory
    tmpDir = fullfile(tempdir, ['mri_vol2surf_' dec2hex(randi(2^20))]);
    if ~exist(tmpDir, 'dir'), mkdir(tmpDir); end
    out_lh = fullfile(tmpDir, 'out_lh.mgz');
    out_rh = fullfile(tmpDir, 'out_rh.mgz');
    
    % set env
    setenv('FREESURFER_HOME', freesurferDir);
    SUBJECTDIR = fullfile(freesurferDir, 'subjects');
    setenv('SUBJECTS_DIR', SUBJECTDIR);
    
    cwd = pwd;
    
    cd(SUBJECTDIR)
    try
        % LH
        cmd_lh = ['"' mri_vol2surf_bin '"', ...
              ' --src "' input_nii '"' ...
              ' --mni152reg' ...
              ' --projfrac ' projfrac ...
              ' --trgsubject ' fsaverage, ...
              ' --interp nearest', ...
              ' --hemi lh', ...
              ' --o "' out_lh '"'];
        [exitCode, outmsg] = system(cmd_lh);
        if exitCode ~= 0
            error('mri_vol2surf:LH', ...
                'Error running mri_vol2surf for LH:\n%s', outmsg);
        end
    
        % RH
        cmd_rh = ['"' mri_vol2surf_bin '"', ...
             ' --src "' input_nii '"' ...
             ' --mni152reg' ...
             ' --projfrac ' projfrac ...
             ' --trgsubject ' fsaverage, ...
             ' --interp nearest', ...
             ' --hemi rh', ...
             ' --o "' out_rh '"'];
        [exitCode, outmsg] = system(cmd_rh);
        if exitCode ~= 0
            error('mri_vol2surf:RH', ...
                'Error running mri_vol2surf for RH:\n%s', outmsg);
        end
    
        % read data
        lh_vertex_val = squeeze(load_mgh(out_lh));
        rh_vertex_val = squeeze(load_mgh(out_rh));
    
    catch ME
        % clean up
        if exist(out_lh, 'file'), delete(out_lh); end
        if exist(out_rh, 'file'), delete(out_rh); end
        if exist(tmpDir, 'dir'), rmdir(tmpDir, 's'); end
        rethrow(ME);
    end
    
    % clean up temp files
    if exist(out_lh, 'file'), delete(out_lh); end
    if exist(out_rh, 'file'), delete(out_rh); end
    if exist(tmpDir, 'dir'), rmdir(tmpDir, 's'); end
    
    cd(cwd);
end
