function [repoData, repoDir] = setup_project(projectDir)
%% =========================================================================
%                        SETUP_PROJECT
%  Detect LNM repository root, add all required paths, and load atlas data
% =========================================================================

%% ------------------------- User-defined Paths ----------------------------

% if not set, MATLAB will search for the required paths in the projectDir
% as an alternative in the folder utils_leaddbs

% Define paths to required toolboxes
% leadDBS_path = 'PATH_TO_LEADBS/leaddbs/'; % (Lead-DBS toolbox: https://www.lead-dbs.org)

% spm_path   = 'PATH_TO_SPM/spm'; % (SPM toolbox required by Lead-DBS)

%% --------------------------- Add to MATLAB Path --------------------------

warning off

try
    if ~exist('leadDBS_path','var')
        leadDBS_path = (fullfile(projectDir,'utils_leaddbs','leaddbs'));
    end
    
    if ~exist('spm_path','var')
        spm_path = (fullfile(projectDir,'utils_leaddbs','spm'));
    end
catch
    error('Please set paths to LeadDBS and SPM12 (requirement for LeadDBS) in setup_project.m')
end

% OS detection
if ispc, osname='windows'; elseif ismac, osname='mac'; else, osname='linux'; end

warning off;

% Add LeadDBS and SPM
if isfolder(leadDBS_path) && ~contains(path, [leadDBS_path pathsep])
    addpath(genpath(leadDBS_path));
end
if isfolder(spm_path) && ~contains(path, [spm_path pathsep])
    addpath(genpath(spm_path));
end

% Add general utility folders
folderPath=fullfile(projectDir, 'utils');
if isfolder(folderPath) && ~contains(path, [folderPath pathsep])
    addpath((folderPath));
end

folderPath=fullfile(projectDir, 'utils', 'toolbox_functions');
if isfolder(folderPath) && ~contains(path, [folderPath pathsep])
    addpath(genpath(folderPath));
end

if strcmp(osname,'windows') %special case for windows
    folderPath=fullfile(projectDir, 'utils', 'toolbox_functions_windows');
    addpath((folderPath));
    rmpath(fullfile(projectDir, 'utils', 'toolbox_funtions'));
end

% % Add project subdirectories
% subdirs = {
%     'normative_connectome', ...
%     'lesions', ...
%     'leadDBS_output', ...
%     'utils', ...
%     'parcellations', ...
%     'FSL_brain', ...
%     'freesurfer_files', ...
%     'data', ...
%     'leaddbs', ...
%     'spm12'
% };
% 
% for i = 1:numel(subdirs)
%     pth = fullfile(projectDir, subdirs{i});
%     if isfolder(pth)
%         addpath(genpath(pth));
%     end
% end

%% =========================================================================
%                        Load Atlas and Metadata
% =========================================================================

% Atlas NIfTI
atlas_file = fullfile(projectDir, 'utils', 'parcellations', ...
    'schaefer1000-yeo7+melbourne_MNI152_2mm.nii.gz'); 
atlas_vol = load_nifti_volume(atlas_file);

% CTAB file (annotation table)
ctab_file = fullfile(projectDir, 'utils', 'parcellations', ...
    'schaefer1000-yeo7.annot.ctab');

opts = detectImportOptions(ctab_file, 'FileType', 'text');
opts.Delimiter     = {' ', '\t'};
opts.VariableNames = {'ID', 'ROI', 'R', 'G', 'B', 'A'};
opts.VariableTypes = {'double', 'string', 'double', 'double', 'double', 'double'};

ctab = readtable(ctab_file, opts);

%% =========================================================================
%                   Load Custom Colormaps (for visualization)
% =========================================================================

[custom_cmap, cmap_brain, cmap_reds, cmap_hot] = get_custom_colormaps();

%% =========================================================================
%                  Define Atlas ROI Ranges 
% =========================================================================

% NOTE: This is fixed for the Melbourne54 and Schaefer1000 atlas used 
atlas_ROIs_subcortical = 1:54;         % Melbourne subcortical regions; these have a fixed position in the data repository
atlas_ROIs_cortical    = 55:1054;      % Schaefer cortical regions; these have a fixed position in the data repository

%% =========================================================================
%                  Package parameters and settings
% =========================================================================

repoData.atlas_ROIs_cortical = atlas_ROIs_cortical;
repoData.atlas_ROIs_subcortical = atlas_ROIs_subcortical;
repoData.atlas_vol = atlas_vol;
repoData.atlas_file = atlas_file;

repoData.ctab = ctab;
repoData.ctab_file = ctab_file;

repoData.colormaps.cmap_hot = cmap_hot;
repoData.colormaps.cmap_brain = cmap_brain;
repoData.colormaps.custom_cmap = custom_cmap;
repoData.colormaps.cmap_reds = cmap_reds;

repoDir.leadDBS_path=leadDBS_path;
repoDir.projectDir = projectDir;



