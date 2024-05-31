function ca_demo_vba_slurm(predefRandOrder,whichRandOrder)
% Code showing the use of Voxel-based Analysis (VBA) with fitlm and 
% commonality function using SLURM. 
% Call this script with input variables using other environment (e.g. bash)
%  - predefRandOrder - predefined permuted matrix (sub x perm) for a set of 
%                      subjects and a number of permutations. First column
%                      should include original order (i.e. 1:num sub)
%  - whichRandOrder - indicating which RandOrder column to run on a worker.
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

 
% Ensure CommonalityAnalysis is in Matlab's path
rootdir = fileparts(fileparts(which('ca_demo_vba.m')));
datadir = [rootdir '/data/rsfa'];
load(fullfile(datadir,'subject_info.mat'));


T.f_rsfa = regexprep(T.f_rsfa,'/home/kt03/Projects/public-code/CommonalityAnalysis',rootdir);
Model = 'f_rsfa ~ Age + c_Sex';

% -------------------------------------------------------------------------
% Assemble cfg structure needed to run the analysis for a specific
% permutation order in predefRandOrder
% -------------------------------------------------------------------------
cfg                 = [];
cfg.model           = Model;
cfg.rootDir         = '/imaging/camcan/sandbox/kt03/temp/'; 
cfg.f_mask          = fullfile(datadir,'mask.nii');
cfg.doCommonality   = 1;
cfg.predefRandOrder = predefRandOrder;
cfg.whichRandOrder  = whichRandOrder;
cfg                 = ca_vba_glm_fitlm(T,cfg);