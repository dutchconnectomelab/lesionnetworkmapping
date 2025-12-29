function S = load_nifti(filename, varargin)
    % Call your existing loader (Jimmy Shen's or FreeSurfer's)
    nii = load_nii(filename, varargin{:});  % or load_nii if that's the function name

    % Normalize field names
    S = nii;
    if isfield(nii, 'img') && ~isfield(nii, 'vol')
        S.vol = double(nii.img);
    end
end
