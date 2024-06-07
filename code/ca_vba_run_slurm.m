function genfi_cvr_analysis_ca_slurm(f_data, Model,rootDir,num_perm, specific_seed,column_index,f_mask)
% Code showing the use of Voxel-based Analysis (VBA) with fitlm and 
% commonality function using SLURM. 
% Call this script with input variables using other environment (e.g. bash)
%  - f_data       - xlsx file containing subject data with variables
%                   specified in Model
%                   
%  - Model        - Equation of the model (see below)
%                  
%  - rootdir      - Root output director where Model results are saved
%  - f_permMatrix - txt file with predefined permuted matrix (sub x perm) for a set of 
%                      subjects and a number of permutations. First column
%                      should include original order (i.e. 1:num sub)
%  - column_index - indicating which RandOrder column to run on a worker.
% 
%
% Running the analysis requires two variables in the workspace, T and
% Model, similar to the usage of 'fitlm'
%
%% T       - 
%  Table containing all variables specified in Model (see below).
%
%% Model   - 
% A string specifying the linear model formula using Wilkinson notation.
% (https://uk.mathworks.com/help/stats/wilkinson-notation.html) 
% e.g. Equation:
% 'y ~ terms', where 
% 'y'       - defines the name of the response variable
% 'terms'   - defines the model using the predictor variable names and the operators.
% 
% Images (e.g. anatomical T1w scans) can be used as response and/or predictor variables. 
%
% Model = 'f_T1w ~ age + sex'; 
%
% Specifies analysis to predict T1w intensity values using age and sex as predictors.
% 
% VBA can dissociate between brainmaps and univariate variables based on 
% the prefix 'f_' of the variable in Model, e.g. f_T1w indicates that
% T.f_T1w contains filepaths to individual's T1-weighted processed images. 
%
% For large models with many covariates of no interest/confounders, it is 
% adviseable to indicate these variables with 'c_'. That way VBA will not
% generate outputs for these variables making it more efficient.
%
% e.g. 
% Model = 'f_T1w ~ Age + c_Sex'
% VBA will not save maps for unique effects of Sex and shared effects of
% Age and Sex.
%
%
% These dependencies are at https://github.com/kamentsvetanov/external
% .../mat/spm12
% and need to be loaded in Matlab's path
% 

%addpath(genpath('/rds/user/kat35/hpc-work/projects/public-code/CommonalityAnalysis/code'));
%addpath(genpath('/rds/user/kat35/hpc-work/projects/external/mat/spm12'));

% -----------------------
% Load permutation matrix
% -----------------------

% permMatrix = readmatrix(f_permMatrix);
% permMatrix = permMatrix';

T = readtable(f_data);
column_index
num_perm
specific_seed
Model
f_mask
% -------------------------------------------------------------------------
% Assemble cfg structure needed to run the analysis for a specific
% permutation order in predefRandOrder
% -------------------------------------------------------------------------
cfg                 = [];
cfg.model           = Model;
cfg.rootDir         = rootDir;%'/rds/project/rds-tbdABxTZvic/genfi_cvr/analysis/vba/temp'; 
cfg.f_mask          = f_mask;%fullfile(datadir,'mask.nii');
cfg.doCommonality   = 0;
% cfg.doRobust        = 1;
cfg.numPerm         = str2num(num_perm);
% cfg.predefRandOrder = permMatrix;
cfg.whichRandOrder  = str2num(column_index);
cfg.specificSeed    = str2num(specific_seed);
cfg
cfg                 = ca_vba_glm_fitlm(T,cfg);

%%

% f_data  = '/rds/user/kat35/hpc-work/projects/public-code/CommonalityAnalysis/data/rsfa/subject_info.xlsx';
% Model   = 'f_rsfa ~ Age + c_Gender';
% Model   = 'f_rsfa ~ Age*GeneStatus + c_Gender';
% Model   = 'f_rsfa ~ Age*GeneStatus + f_gm + c_Gender';
% rootDir = '/rds/project/rds-tbdABxTZvic/genfi_cvr/analysis/vba/temp';
% f_permMatrix = '/rds/project/rds-tbdABxTZvic/genfi_cvr/analysis/vba/temp/permutation_matrix.txt';
% column_index = 1;
% f_mask = '/rds/user/kat35/hpc-work/projects/public-code/CommonalityAnalysis/data/rsfa/mask.nii';
% f_mask = '/rds/project/rds-tbdABxTZvic/genfi_cvr/analysis/vba/ca/SPM12_GM_mask70tpm_91x109x91.nii';

% numSub          = 20;
% numPerm         = 100;
% permutedMatrix  = arrayfun(@(x) randperm(numSub), 1:numPerm*1.1, 'UniformOutput', false);
% permutedMatrix  = cell2mat(permutedMatrix')';
% % Remove columns that cointain the original order
% idx = corr(permutedMatrix,[1:numSub]')==1;
% permutedMatrix(:,idx)=[];
% permutedMatrix(:,1) = [1:numSub]'; % Set first column to original order
% % get the right size of pertmuted Matrix
% randOrder = permutedMatrix(:,1:numPerm);

% T.f_rsfa = regexprep(T.f_rsfa,'/home/kt03/Projects/public-code/CommonalityAnalysis',rootdir);
% T.f_rsfa = string(T.f_rsfa);
% T.SubID = string(T.SubID);
% Ensure CommonalityAnalysis is in Matlab's path
% rootdir = fileparts(fileparts(which('ca_demo_vba.m')));
% datadir = [rootdir '/data/rsfa'];
% T       = readtable(fullfile(datadir,'subject_info.xlsx'));
% T.f_rsfa = regexprep(T.f_rsfa,'/home/kt03/Projects/public-code/CommonalityAnalysis',rootdir);
% T.f_rsfa = string(T.f_rsfa);
% T.SubID = string(T.SubID);

% writetable(T,'subject_info.csv','QuoteStrings',true)
