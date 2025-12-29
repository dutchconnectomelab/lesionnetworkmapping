%% ========================================================================
%  sLNM Linear Form Step 1, 2 and group-level analysis Step 3
%  ------------------------------------------------------------------------
%  This script demonstrates the linear form of the Symptom Lesion Network Mapping 
%  (sLNM) comparing leadDBS output and the atlas-based sLNM = sv * M * C approach.
%
%  Requirements:
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
%                     Load Environment and Define Paths
% =========================================================================
if isempty(projectDir) %if still not set, use current directory; 
    projectDir = fileparts(mfilename ('fullpath'));
end
[repoData, repoDir] = setup_project(projectDir);
atlas_vol=repoData.atlas_vol;
ctab = repoData.ctab;

%% =========================================================================
%                 SELECT DISORDER and CLINICAL DATA Sv
% =========================================================================

% lesionSet = 'SyntheticLesions'; %in the case of syntheticLesions sv = [1:100]

%Alternatives:
% lesionSet = 'APHASIA_RECOVERY'; %sv = patient aphasia score https://openneuro.org/datasets/ds004884/versions/1.0.1
lesionSet = 'MDDDLPFC_TMS'; % sv = symptom improvement after TMS treatment 
% lesionSet = 'LESYMAP'; % sv = symptom scorese are synthetic made by LESYMAP

% load in Clinical Values
lesionDir   = fullfile(projectDir, 'data', 'lesion_masks', lesionSet);
ClinicalFile = dir(fullfile(lesionDir, '*ClinVal.txt'));
sv = readtable(fullfile(ClinicalFile(1).folder, ClinicalFile(1).name),'ReadVariableNames', false).Var1;

%% =========================================================================
%                         SECTION 1: lead-DBS sLNM Procedure 
% =========================================================================

% Step 1: Select lesions in MNI152 space
lesionFiles = dir(fullfile(lesionDir, '*.nii.gz'));
nLesions    = length(lesionFiles);
fprintf('%d lesions from %s selected\n', nLesions, lesionSet);

leadDBS_outputfolder = fullfile(projectDir, 'data', 'leadDBS_output', lesionSet);

% Step 2 : For each individual lesion run leadDBS pipeline (WARNING slow - this step takes around 10-20 per lesion)
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
% leadDBS_outputFiles = dir(fullfile(leadDBS_outputfolder, '*avgR_funcmap.nii.gz'));
LNM_leadDBS_voxel = zeros(nLesions,91,109,91); % 2mm MNI
for ii = 1:nLesions
    fprintf('Loading pre-computed leadDBS output for lesion %d of %d\n', ii, nLesions)
    vol = load_nifti_volume(fullfile(leadDBS_outputFiles(ii).folder, ...
        leadDBS_outputFiles(ii).name));
    vol(isnan(vol))=0; %clear out some possible NaNs
    LNM_leadDBS_voxel(ii,:,:,:) = vol;
end

% Step 2c : read in clinical values 
% this is done above

% Step 3 Group-level analysis : compute r-map with clinical values
sLNM_rmap_leadDBS_voxel=zeros(size(LNM_leadDBS_voxel,2:4));
for xx = 1:size(LNM_leadDBS_voxel,2)
    for yy = 1:size(LNM_leadDBS_voxel,3)
        for zz = 1:size(LNM_leadDBS_voxel,4)
            sLNM_rmap_leadDBS_voxel(xx,yy,zz)=corr(LNM_leadDBS_voxel(:,xx,yy,zz), sv);
        end
    end
end
% clean up NaN
sLNM_rmap_leadDBS_voxel(isnan(sLNM_rmap_leadDBS_voxel))=0;

%% =========================================================================
%             SECTION 2: Streamlined linear form sLNM = (sv*M*C) 
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
    'GSP1000_FC_Schaefer1000Melbourne54.mat')).FC;

% Optional: Fisher r-to-z
% C = load(fullfile(projectDir, 'data', 'normative_connectome', ...
%     'GSP1000_FCz_Schaefer1000Melbourne54.mat')).FC; % functional connectivity matrix

% Step 2c : read in clinical values
% this is done above

% Step 3 Group-level analysis : compute r-map with clinical values
M = all_lesions_atlas ./ sum(all_lesions_atlas,2); %make lesion matrix M (each leasion has total sum of 1)
M(isnan(M))=0;
sv = (sv-mean(sv))./std(sv); % 0-center sv (e.g. take z-score)
sLNM_equation4_atlas_map = sv'*(M*C);

%% =========================================================================
%         SECTION 3: %% Compare voxelwise leadDBS and streamlined version
% =========================================================================

% Map voxelwise leadDBS LNM map to atlas for comparison
sLNM_rmap_leadDBS_atlas_map = map_voxel2schaeferMelbourne(sLNM_rmap_leadDBS_voxel, ...
    'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab);

%% visualize LNM on brain surface
zscore_sLNM_voxel = zscore(sLNM_rmap_leadDBS_atlas_map); 
zscore_sLNM_eq4 = zscore(sLNM_equation4_atlas_map);

plot_rh_surface_twotiles(projectDir,lesionSet, ...
    zscore_sLNM_voxel, zscore_sLNM_eq4,... %data of panel 1 and 2
    'Voxel-wise sLNM vs sv*M*C', 'sLNM leadDBS', 'sLNM sv*M*C',... % main title, and titles of panels
    'sLNM leadDBS (z-score)', 'sLNM sv*M*C (z-score)') %x and y label of third panel

