%% ========================================================================
%  sLNM using random Clinical Values and random Lesions can result in similar
%  maps as original sLNM maps with clinically informed lesions and sv 
%  ------------------------------------------------------------------------
%  This script demonstrates that sLNM using randomized Lesions and random Clinival Values can lead to 
%  systematic traces in published sLNM maps
%
%  Requirements:
%   - leadDBS toolbox: https://github.com/netstim/leaddbs
%   - SPM12 toolbox: needed for leadDBS https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ 
%   - GSP1000 dataset: needed for leadDBS https://dataverse.harvard.edu/dataverse/cohenlab
% ========================================================================

close all; clc; clear;
warning('off','all')

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
%                 SELECT ORIGINAL PUBLISHED sLNM map
% =========================================================================

% voxel-wise LNM maps
original_sLNMmap = 'SIDDIQI_MSDEPRESSION'

% Alternatives:
% original_sLNMmap = 'FOX_MDDDLPFCTMS';
% original_sLNMmap = 'JOUTSA_TREMORDBS';
% original_sLNMmap = 'SIDDIQI_ANXIETYi'
% original_sLNMmap = 'JOUTSA_ADDICTIONCONTRASTAB';
% original_sLNMmap = 'SIDDIQI_PTSD'
% original_sLNMmap = 'LI_OCDDBS'

% load in original sLNM map
rootDir = fullfile(projectDir, 'data', 'LNM_networks');
files = dir(fullfile(rootDir, '**', [original_sLNMmap, '_2mm.nii.gz']));
if isempty(files)
    error('File not found.'); end
map_path = fullfile(files(1).folder, files(1).name);

orig_sLNMmap = load_nifti_volume(map_path); 
orig_sLNMmap_atlas = map_voxel2schaeferMelbourne(orig_sLNMmap, ...
        'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab);

%% =========================================================================
%               TAKE RANDOM LESIONS AND RANDOM CLIN VALS
% =========================================================================

randomLesions = 'syntheticLesions';
lesionDir   = fullfile(projectDir, 'data', 'lesion_masks', randomLesions);
lesionFiles = dir(fullfile(lesionDir, '*.nii.gz'));
sv=randn(size(lesionFiles,1),1); %pick a random sv to start with.

%% =========================================================================
%                         sLNM = sv * M * C , with random M and random sv
% =========================================================================

nLesions=size(lesionFiles,1);

% Step 1: Map lesions to atlas space
all_lesions_atlas = zeros(nLesions,1054);
for ii = 1:nLesions
    fprintf('Mapping lesion %d of %d to atlas space\n', ii, nLesions)
    vol = load_nifti_volume(fullfile(lesionFiles(ii).folder, lesionFiles(ii).name));
    lesionT = 0; if nnz(vol) > 1000; lesionT=0.5;end %if the lesion is large 
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

%% ========== sLNM map with random lesions and random sv ==========
% Step 3 Group-level analysis : compute r-map with clinical values
% we can run the entire sLNM procedure with randomized clinical values; 
% as predicted from Equation 4, the outcomes are not random, but driven by
% the factors of C; we now use random lesions and random sv
M = all_lesions_atlas ./ sum(all_lesions_atlas,2); % make lesion matrix M (each leasion has total sum of 1)
M(isnan(M))=0;

%for replication
rng(333)

nPerm=1000; %number of random runs
nRandomLesions=100;
sLNM_equation4_atlas_maps_random=zeros(nPerm,size(C,1));

for iPerm=1:nPerm
    Mrandom = all_lesions_atlas(randperm(size(all_lesions_atlas,1), nRandomLesions),:);
    Mrandom = Mrandom ./ sum(Mrandom);
    Mrandom(isnan(Mrandom))=0;
    sv_random = sv(randperm(length(sv))); %randomize the clinival values in sv
    sv_random = sv_random(1:nRandomLesions);
    sLNM_equation4_atlas_maps_random(iPerm,:)=sv_random'*(Mrandom*C);
end

% while run with a randomized sv_random, all random maps still share a strong underlying latent factor(s) 
% we can find such latent factor(s):

% Option 1 : take the mean or alternative find the positive and negative dominant
% cluster (see Supplementary Information)
% c=modularity_und(corr(sLNM_equation4_atlas_maps_random')); % 2 dominant clusters which are each other opposite. 
% sLNM_equation4_atlas_maps_random_randgroup=sLNM_equation4_atlas_maps_random(c==randi(2),:); %take randomly 1 of the groups
% latentfactor = mean(sLNM_equation4_atlas_maps_random_randgroup);

% Alternative - we can also take a PCA and examine the latent factors in
% the random runs; if there would be no random run, explained var
% should be (1+sqrt(2))^2) = ~5.8, but it is much higher indicating a
% strong structure in the random runs
[coeff_vec,s,~,~,explained_variance]=pca(sLNM_equation4_atlas_maps_random);
latentfactor=coeff_vec(:,1);

% Optional: the sign of the pca1 is arbitrary - for visualisation we can make it equal
if corr(latentfactor(:), orig_sLNMmap_atlas(:)) < 0
    latentfactor=-latentfactor;
end

%% =========================================================================
%       %%  COMPARE original published sLNM and sLNM with random lesions and random sv
% =========================================================================

plot_rh_surface_twotiles(projectDir,original_sLNMmap, ...
    zscore(latentfactor), zscore(orig_sLNMmap_atlas), ... %data of panel 1 and 2
    'Comparison random vs original sLNMP (sign is arbitrary)', 'average of sLNM on random lesions and random sv','original sLNM map',... % main title, and titles of panels
    'sLNM random lesions random sv (z-score)', 'original sLNM map (z-score)') %x and y label of third panel















