%% ========================================================================
%  LNM comparison to LNM with randomized Lesions
%  ------------------------------------------------------------------------
%  This script uses the linear form of the Lesion Network Mapping 
%  (LNM) comparing LNM output between patient and randomized (spin) lesions 
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

lesionSet = 'SyntheticLesions'

% Alternatives:
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
%              LNM linear form LNM = sum(M*C)
% =========================================================================

% Step 1a: Select lesions in MNI152 space
lesionDir   = fullfile(projectDir, 'data', 'lesion_masks', lesionSet);
lesionFiles = dir(fullfile(lesionDir, '*.nii.gz'));
nLesions    = length(lesionFiles);
fprintf('%d lesions from %s selected\n', nLesions, lesionSet);

% Step 1b: Map lesions to atlas space
all_lesions_atlas = zeros(nLesions,1054);
for ii = 1:nLesions
    fprintf('Mapping lesion %d of %d to atlas space\n', ii, nLesions)
    vol = load_nifti_volume(fullfile(lesionFiles(ii).folder, lesionFiles(ii).name));
    lesionT=0; if nnz(vol) > 1000; lesionT=0.5;end %if the lesion is large 
    lesionVec = map_voxel2schaeferMelbourne(vol, ...
        'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab) > lesionT;
    all_lesions_atlas(ii,:) = lesionVec;
end

% Step 2 : Take the group matrix C 
C = load(fullfile(projectDir, 'data', 'normative_connectome', ...
    'GSP1000_FC_Schaefer1000Melbourne54.mat')).FC;

% Optional: Fisher r-to-z
% C = load(fullfile(projectDir, 'data', 'normative_connectome', ...
%     'GSP1000_FCz_Schaefer1000Melbourne54.mat')).FC; % functional connectivity matrix

% for spin we only can use cortical 
% (this may lead to small lesion sets, if lesions are mostly subcortical)
atlas_ROIs_cortical=repoData.atlas_ROIs_cortical;
C=C(atlas_ROIs_cortical,atlas_ROIs_cortical);
all_lesions_atlas = all_lesions_atlas(:,atlas_ROIs_cortical);

% Step 3 Group-level analysis : Combine outputs into group LNM sensitivity map
M = all_lesions_atlas ./ sum(all_lesions_atlas,2); % make lesion matrix M (each leasion has total sum of 1)
LNM_equation3_atlas_map = nansum(M*C,1);

%% =========================================================================
%                 Compute LNM on randomized lesions (spin)
% =========================================================================

% Spin-model cortical-only randomization
M = all_lesions_atlas;

nPerm=100; %set to 10,000 for full run;
% spin the lesions randomly, creating nPerm random sets
randomSeed=538; %use fixed seed for replication
% Optional: use random seed for fresh runs
% randomSeed=abs(round(randn()*10000)); % or
randomSeed=NaN;

% run all lesions with a new spin, more liberal model
% M_spin_all = spin_model(M, 'projectDir',projectDir, 'random_number_generator', randomSeed, 'nPerm', nPerm, 'method', 'spin_independent', 'verbose', true);

% Optional: run all lesions with the same set, more conservative model
M_spin_all = spin_model(M, 'projectDir',projectDir, 'random_number_generator', randomSeed, 'nPerm', nPerm, 'method', 'spin_same_angle', 'verbose', true);

% Optional: run all lesions with parcel centroid rotation as presented by vasa et al.
% M_spin_all = spin_model_vasa(M, 'projectDir',projectDir, 'random_number_generator', randomSeed, 'nPerm', nPerm, 'method', 'spin_same_angle', 'verbose', true);

% Use the sets of spin lesions to compute LNM
LNM_spin_all=zeros(nPerm,length(atlas_ROIs_cortical));
corrs_spin = zeros(nPerm, 1);

for ii=1:nPerm
    % get one of the spins
    M_spin = squeeze(M_spin_all(ii,:,:));

    % remove lesions that end up fully in medial wall (no parcels left in
    % rotated lesion)
    M_spin = M_spin(sum(M_spin, 2) > 0, :);

    % binarize
    M_spin = M_spin > 0;

    % normalize spin permuted lesion vectors
    M_spin_norm = M_spin ./ sum(M_spin, 2);

    % compute LNM
    LNM_spin_all(ii,:) = sum(M_spin_norm * C, 1);

    % correlation with original LNM
    corrs_spin(ii) = corr(LNM_spin_all(ii,:)', LNM_equation3_atlas_map');
end


%% =========================================================================
%      Compare original patient LNM with randomized spin lesions
% =========================================================================

%% visualize LNM on brain surface
clear zscore*
zscore_LNM_orig(:,atlas_ROIs_cortical) = zscore(LNM_equation3_atlas_map); 
zscore_LNM_spin(:,atlas_ROIs_cortical) = zscore(mean(LNM_spin_all,1)); %plot against the mean of all random spins

plot_rh_surface_twotiles(projectDir,lesionSet, ...
    zscore_LNM_orig, zscore_LNM_spin,... %data of panel 1 and 2
    'Patient vs spin lesions LNM', 'original LNM', 'LNM based on randomized lesions',...  % main title, and titles of panels
    'Patient vs spin lesions LNM  (z-score)', 'original LNM  (z-score)') %x and y label of third panel
