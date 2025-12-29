%% ========================================================================
%  LNM Linear Form Step 1, 2 and group-level analysis Step 3
%  ------------------------------------------------------------------------
%  This script demonstrates how published LNM maps are formed by basic
%  elementary properties of the connectivity matrix C: 9 factors of degree, subcortical
%  degree, broad resting-state modules and primary gradients of FC
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
%                     Load  Environment and Define Paths
% =========================================================================
if isempty(projectDir) %if still not set, use current directory; 
    projectDir = fileparts(mfilename ('fullpath'));
end

[repoData, repoDir] = setup_project(projectDir);
atlas_vol=repoData.atlas_vol;
ctab = repoData.ctab;

%% =========================================================================
%               SELECT ORIGINAL PUBLISHED LNM/sLNM MAP
% =========================================================================

% Example map
original_LNMmap = 'JOUTSA_ADDICTIONGROUPA';

% Alternatives:
% voxel-wise LNM maps
% original_LNMmap = 'PINES_PSYCHOSIS';
% original_LNMmap = 'SIDDIQI_PTSD';
% original_LNMmap = 'DARBY_CRIMINALITY';
% original_LNMmap = 'GANOS_TICS'
% original_LNMmap = 'COTOVIO_OCD';
% original_LNMmap = 'KUTSCHE_APHANTASIA';
% original_LNMmap = 'SEGAL_OCD'; 
% original_LNMmap = 'SEGAL_SCZ'; 
% original_LNMmap = 'ZARIFKAR_APHASIA'
% original_LNMmap = 'COTOVIO_MANIASETA'
% original_LNMmap = 'CRISTOFORI_TBI'
% voxel-wise sLNM maps
% original_LNMmap = 'SIDDIQI_MSDEPRESSION';
% original_LNMmap = 'JOUTSA_ADDICTIONCONTRASTAB';
% original_LNMmap = 'FOX_MDDDLPFCTMS';
% original_LNMmap = 'SIDDIQI_PTSD'
% original_LNMmap = 'SIDDIQI_ANXIETYi'
% original_LNMmap = 'JOUTSA_TREMORDBS';
% original_LNMmap = 'LI_OCDDBS'

%% ========================================================================
%    Find and load in map from computed or downloaded LNM maps
% =========================================================================
% Load in order from various potential directories
rootDir = fullfile(projectDir, 'data', 'LNM_networks');
files = dir(fullfile(rootDir, '**', [original_LNMmap, '_2mm.nii.gz']));
if isempty(files)
    error('File not found.'); end
map_path = fullfile(files(1).folder, files(1).name);

% Convert load in voxel map and convert to Schaefer-Melbourne atlas space
orig_LNMmap = load_nifti_volume(map_path); 
orig_LNMmap_atlas = map_voxel2schaeferMelbourne(orig_LNMmap, ...
        'projectDir', projectDir, 'atlas_vol', atlas_vol, 'ctab',ctab);

%% =========================================================================
%                    Load Connectivity Matrix and Metadata
% =========================================================================
C = load(fullfile(projectDir, 'data', 'normative_connectome', ...
    'GSP1000_FC_Schaefer1000Melbourne54.mat')).FC; % functional connectivity matrix

% Optional: Fisher r-to-z
% C = load(fullfile(projectDir, 'data', 'normative_connectome', ...
%     'GSP1000_FCz_Schaefer1000Melbourne54.mat')).FC; % functional connectivity matrix

%% =========================================================================
%                     Compute Explanatory Variables
% =========================================================================

rng(555); % for replication
fprintf('Regression model on degree, subcortical degree, basic modules, and gradients (9 factors), \n');
% 1. Global Degree
deg_path = fullfile(projectDir, 'data', 'degree', 'GSP1000_degree_fisher_z.nii.gz');
global_degree = load_nifti_volume(deg_path); 
global_degree(global_degree == global_degree(1)) = 0;  % Remove small noise values
global_degree_atlas = map_voxel2schaeferMelbourne(global_degree, 'projectDir', projectDir, ...
    'atlas_vol', atlas_vol, 'ctab', ctab);

% 2. Subcortical Degree
subcortical_degree = mean(C(repoData.atlas_ROIs_subcortical, :), 1);

% 3. RSN Degree (by module)
rsn_modules = modularity_und(C);
unique_modules = unique(rsn_modules);
rsn_degree = zeros(length(unique_modules), size(C,2));

for ii = 1:length(unique_modules)
    rsn_degree(ii, :) = mean(C(rsn_modules == unique_modules(ii), :), 1);
end

% 4. PCA Gradients
[c, ~, ~, ~, explained_variance] = pca(C);
k = 3;  % Number of principal components to use
gradients = c(:, 1:k)';


%% =========================================================================
%                        Run Regression Model
% =========================================================================
% Linear regression to explain the LNM map with network features
X = [global_degree_atlas; subcortical_degree; rsn_degree; gradients]';
s = regstats(orig_LNMmap_atlas, X);

% Report R²
R2 = s.rsquare;
fprintf('Explained variance (R²) for original LNM : %s %.4f\n', original_LNMmap, R2);

% plot
X_with_intercept = [ones(size(X,1),1), X];  % add intercept
Y_pred = X_with_intercept * s.beta;          % predicted values

plot_rh_surface_twotiles(projectDir, original_LNMmap, ...
    zscore(Y_pred), zscore(orig_LNMmap_atlas), ... %data of panel 1 and 2
    'Regression model explains signal in (s)LNM map', 'predicted map based on simple regression model', 'original LNM map', ... % main title, and titles of panels
    'predicted map based on simple regression model (z-score)','original LNM map (z-score)') %x and y label of third panel

%plot_lh_surface_twotiles(projectDir, original_LNMmap, ...
%    zscore(Y_pred), zscore(orig_LNMmap_atlas), ... %data of panel 1 and 2
%    'Regression model explains signal in (s)LNM map', 'predicted map based on simple regression model', 'original LNM map', ... % main title, and titles of panels
%    'predicted map based on simple regression model (z-score)','original LNM map (z-score)') %x and y label of third panel








