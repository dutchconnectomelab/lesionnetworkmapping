%% ========================================================================
%  LNM compared to LNM on random lesions 
%  ------------------------------------------------------------------------
%  This script compares Lesion Network Mapping (LNM) on patient lesions associated
%  with a specific disorder to an LNM map derived from lesions across a
%  random mix of disorders.
%
%  Requirements - SET PATHS IN setup_project.m:
%   - leadDBS toolbox: https://github.com/netstim/leaddbs
%   - SPM12 toolbox: needed for leadDBS https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%   - GSP1000 dataset: needed for leadDBS https://dataverse.harvard.edu/dataverse/cohenlab
% ========================================================================

close all; clc; clear;

%% =========================================================================
%%                         USER INPUT : projectDir 
%% =========================================================================

% setup_project
% set LeadDBS and SPM12 path in setup_project.m

% Define project directory and toolbox paths
projectDir   = ''; %fill in projectDir here where the scripts are located

%% =========================================================================
%                     Load  Environment and Define Paths
% =========================================================================
if isempty(projectDir) %if still not set, use current directory; 
    projectDir = fileparts(mfilename ('fullpath'));
end

[repoData, repoDir] = setup_project(projectDir);
atlas_vol=repoData.atlas_vol;
ctab = repoData.ctab;

%% =========================================================================
%                       SELECT LESIONSET
% =========================================================================

% lesionSet = 'SyntheticLesions'

% Alternatives:
lesionSet = 'ADDICTION'; 
% lesionSet = 'AGENCY';
% lesionSet = 'APHASIA_recovery'
% lesionSet = 'EPILEPSY';
% lesionSet = 'LESYMAP';
% lesionSet = 'MIGRAINE'
% lesionSet = 'PARKINSONISM';
% lesionSet = 'AMUSIA';
% lesionSet = 'APRAXIA';
% lesionSet = 'CREATIVITY';
% lesionSet = 'DELUSION';
% lesionSet = 'MANIA';
% lesionSet = 'NEGLECT';
% lesionSet = 'OCD';
% lesionSet = 'STUTTERING';

lesionDir   = fullfile(projectDir, 'data', 'lesion_masks', lesionSet);
lesionDirAll   = fullfile(projectDir, 'data', 'lesion_masks', '*');

%% =========================================================================
%   lead-DBS LNM Procedure on lesions associated to the specific condition
% =========================================================================

% Step 1: Select lesions in MNI152 space
lesionFiles = dir(fullfile(lesionDir, '*.nii.gz'));
nLesions    = length(lesionFiles);
fprintf('%d lesions from %s selected\n', nLesions, lesionSet);

leadDBS_outputfolder = fullfile(projectDir, 'data', 'leadDBS_output', lesionSet);

% Step 2 : For each individual lesion run leadDBS pipeline (WARNING slow - this step takes around 10-20 minutes per lesion)
use_precompute = true; % Optional: set false to compute from scratch
if ~use_precompute   
    for ii = 1:nLesions
        lesionToProcess{1} = fullfile(lesionFiles(ii).folder, lesionFiles(ii).name);
        fprintf('Running leadDBS connectome mapper for lesion %d of %d\n', ii, nLesions)
        writeoutsinglefiles = 0; 
        cs_fmri_conseed_seed_tc(fullfile(leadDBS_path, 'connectomes'), ...
            'GSP1000', lesionToProcess, 'seed', writeoutsinglefiles, ...
            leadDBS_outputfolder, []);
    end
end

% Step 2b : Collect all leadDBS output t-maps
leadDBS_outputFiles = dir(fullfile(leadDBS_outputfolder, '*T_funcmap.nii.gz'));
% Alternative: avgR maps
%leadDBS_outputFiles = dir(fullfile(leadDBS_outputfolder, '*avgR_funcmap.nii.gz'));
LNM_leadDBS_voxel = zeros(nLesions,91,109,91); % 2mm MNI
for ii = 1:nLesions
    fprintf('Loading pre-computed leadDBS output for lesion %d of %d\n', ii, nLesions)
    vol = load_nifti_volume(fullfile(leadDBS_outputFiles(ii).folder, ...
        leadDBS_outputFiles(ii).name));
    LNM_leadDBS_voxel(ii,:,:,:) = vol;
end

% Step 3 Group-level analysis : Combine leadDBS outputs into group LNM sensitivity map
LNM_leadDBS_voxel_map = squeeze(nanmean(LNM_leadDBS_voxel,1));
LNM_leadDBS_atlas_map = map_voxel2schaeferMelbourne(LNM_leadDBS_voxel_map, ...
    'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab);

%% =========================================================================
%    Lead-DBS LNM Procedure on lesions associated with no-specific condition
% =========================================================================

% nLesionsMixed=250;
%Alternative:
nLesionsMixed = nLesions; % similar as original set size

% Step 1: Select lesions in MNI152 space
fprintf('%d lesions from %s selected\n', nLesionsMixed, lesionSet);

% select random lesions as the set of lesions across multiple conditions
leadDBS_outputfolderAll = fullfile(projectDir, 'data', 'leadDBS_output', '*');
leadDBS_outputFilesMixed = dir(fullfile(leadDBS_outputfolderAll, '*T_funcmap.nii.gz'));
% % Alternative:
% leadDBS_outputFilesMixed = dir(fullfile(leadDBS_outputfolderAll, '*avgR_funcmap.nii.gz'));
%pick random selection of number of lesions from the total set of clinically informed lesions, but mixed conditions
leadDBS_outputFilesMixed = leadDBS_outputFilesMixed(randperm(size(leadDBS_outputFilesMixed,1),nLesionsMixed)); 

% Step 2b : Collect all leadDBS output t-maps
LNM_leadDBS_mixed_voxel = zeros(nLesionsMixed,91,109,91); % 2mm MNI
for ii = 1:nLesionsMixed
    fprintf('Loading pre-computed leadDBS output for mixed lesion %d of %d\n', ii, nLesionsMixed)

    vol = load_nifti_volume(fullfile(leadDBS_outputFilesMixed(ii).folder, ...
        leadDBS_outputFilesMixed(ii).name)); 
    LNM_leadDBS_mixed_voxel(ii,:,:,:) = vol;
end

% Step 3 Group-level analysis : Combine leadDBS outputs into group LNM sensitivity map
LNM_leadDBS_mixed_voxel_map = squeeze(nanmean(LNM_leadDBS_mixed_voxel,1));
LNM_leadDBS_mixed_atlas_map = map_voxel2schaeferMelbourne(LNM_leadDBS_mixed_voxel_map, ...
    'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab);

%% =========================================================================
%      Compare original patient LNM with randomized spin lesions
% =========================================================================

%% visualize LNM on brain surface
clear zscore*
zscore_LNM_orig = zscore(LNM_leadDBS_atlas_map); 
zscore_LNM_mixed = zscore(LNM_leadDBS_mixed_atlas_map); %plot against the mean of all random spins

plot_rh_surface_twotiles(projectDir,lesionSet, ...
    zscore_LNM_orig, zscore_LNM_mixed,... %data of panel 1 and 2
    'Patient vs mixed lesions LNM', 'original LNM map', 'LNM based on mixed lesions',...  % main title, and titles of panels
    'original LNM map (z-score)', 'LNM on mixed lesions (z-score)') %x and y label of third panel












