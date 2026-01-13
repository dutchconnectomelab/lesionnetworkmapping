%% ========================================================================
%  LNM Linear Form Step 1, 2 and group-level analysis Step 3
%  ------------------------------------------------------------------------
%  This script demonstrates the linear form of the Lesion Network Mapping 
%  (LNM) comparing leadDBS output and the streamlined LNM = sum(M*C) approach.
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

% lesionSet = 'SingleLesion';

% Alternatives:
lesionSet = 'SyntheticLesions'
% lesionSet = 'ADDICTION'; 
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

%% =========================================================================
%               SECTION 1: lead-DBS LNM Procedure 
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
        cs_fmri_conseed_seed_tc(fullfile(repoDir.leadDBS_path, 'connectomes'), ...
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

% Step 3 Group-level analysis: Combine leadDBS outputs into group LNM sensitivity map
LNM_leadDBS_voxel_map = squeeze(nanmean(LNM_leadDBS_voxel,1));

% Optional:
% group_T_threshold=7;
% tmask_pos = squeeze(mean(LNM_leadDBS_voxel>group_T_threshold));
% tmask_neg = -squeeze(mean(LNM_leadDBS_voxel<-group_T_threshold));
% LNM_leadDBS_voxel_map = tmask_pos + tmask_neg;

%% =========================================================================
%              SECTION 2: linear form LNM = sum(M*C)
% =========================================================================

% Step 1: Map lesions to atlas space
all_lesions_atlas = zeros(nLesions,1054);
for ii = 1:nLesions
    fprintf('Mapping lesion %d of %d to atlas space\n', ii, nLesions)
    vol = load_nifti_volume(fullfile(lesionFiles(ii).folder, lesionFiles(ii).name));
    lesionT=0; if nnz(vol) > 1000; lesionT=0.5;end %if the lesion is large 
    lesionVec = map_voxel2schaeferMelbourne(vol, ...
        'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab) > lesionT;
    all_lesions_atlas(ii,:) = lesionVec;
end

% Step 2a Load group-level connectome matrix (atlas space)
C = load(fullfile(projectDir, 'data', 'normative_connectome', ...
    'GSP1000_FC_Schaefer1000Melbourne54.mat')).FC; % functional connectivity matrix

% Optional: Fisher r-to-z
% C = load(fullfile(projectDir, 'data', 'normative_connectome', ...
%     'GSP1000_FCz_Schaefer1000Melbourne54.mat')).FC; % functional connectivity matrix

% Step 2b For each individual lesion: compute lesion FC map by averaging rows of C corresponding to lesion
% % % LNM_streamlined_atlas = nan(nLesions,1054);
for ii = 1:nLesions
    rows_matching_lesions = find(all_lesions_atlas(ii,:));
    LNM_streamlined_atlas(ii,:) = mean(C(rows_matching_lesions,:),1);
end

% Step 3 Group-level analysis : Combine outputs into group LNM sensitivity map
LNM_streamlined_atlas_map = nanmean(LNM_streamlined_atlas,1);

% Optional: linear matrix form LNM=sum(M*C)
M = all_lesions_atlas ./ sum(all_lesions_atlas,2);
LNM_equation3_atlas_map = nansum(M*C,1);

%% =========================================================================
%          SECTION 3: Compare the two implementations of LNM
% =========================================================================

%% Compare voxelwise leadDBS and streamlined versions
% Map voxelwise leadDBS LNM map to atlas for comparison
LNM_leadDBS_atlas_map = map_voxel2schaeferMelbourne(LNM_leadDBS_voxel_map, ...
    'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab);

%% =========================================================================
%          %% visualize on brain surface
% =========================================================================

zscore_LNM_voxel = zscore(LNM_leadDBS_atlas_map); 
zscore_LNM_eq3 = zscore(LNM_equation3_atlas_map);

plot_rh_surface_twotiles(projectDir,lesionSet, ...
    zscore_LNM_voxel, zscore_LNM_eq3,... %data of panel 1 and 2
    'Voxel-wise vs Atlas-based LNM', 'LNM leadDBS', 'LNM = sum(M*C)',...  % main title, and titles of panels
    'LNM leadDBS (z-score)', 'LNM = sum(M*C) (z-score)') %x and y label of third panel

% plot_lh_surface_twotiles(projectDir,lesionSet, ...
%     zscore_LNM_voxel, zscore_LNM_eq3,... %data of panel 1 and 2
%     'Voxel-wise vs Atlas-based LNM', 'LNM leadDBS', 'LNM = sum(M*C)',...  % main title, and titles of panels
%     'LNM leadDBS (z-score)', 'LNM = sum(M*C) (z-score)') %x and y label of third panel
