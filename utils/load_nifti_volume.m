function vol = load_nifti_volume(fileName)
    % LOAD_NIFTI_VOLUME Load a NIfTI image volume with a toolbox-aware fallback.
    %
    %   vol = load_nifti_volume(fileName)
    %
    %   Loads a NIfTI file using MATLAB's niftiread (if the Image Processing
    %   Toolbox is available). If niftiread is not available, falls back to
    %   FreeSurfer/External IO via load_nifti and returns the .vol field.
    %
    % Input:
    %   fileName : Path to a NIfTI file (.nii or .nii.gz).
    %
    % Output:
    %   vol      : 3-D (or 4-D) numeric array containing the NIfTI volume data.
    %
    % Dependencies:
    %   • Image Processing Toolbox (niftiread) OR
    %   • load_nifti (e.g., FreeSurfer-compatible NIfTI reader) on the MATLAB path.
    %
    % Example
    %   vol = load_nifti_volume(fullfile(projectDir,'data','brain','degree.nii.gz'));

    if license('test', 'Image_Toolbox') && exist('niftiread', 'file') == 2
        vol = niftiread(fileName);
    else
        niftiStruct = load_nifti(fileName); % Assuming load_nifti is in the path
        vol = niftiStruct.vol;
    end
end