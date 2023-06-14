% Code showing the use of Voxel-based Analysis (VBA) with fitlm and commonality function.
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
% Depenendencies (download from https://github.com/kamentsvetanov/external):
% SPM 12 
% palm_quickperms
%
% These dependencies are at https://github.com/kamentsvetanov/external
% .../mat/spm12
% .../mat/palm
% and need to be loaded in Matlab's path
% 


clear 
rootdir = '/home/kt03/Projects/public-code/CommonalityAnalysis/data/rsfa/';
load(fullfile(rootdir,'subject_info.mat'));

Model = 'f_rsfa ~ Age + c_Sex';

% -------------------------------------------------
% Assemble cfg structure needed to run the analysis
% -------------------------------------------------
cfg                 = [];
cfg.model           = Model;
cfg.rootDir         = '/imaging/camcan/sandbox/kt03/temp/'; 
cfg.f_mask          = fullfile(rootdir,'mask.nii');
cfg.numPerm         = 100;
cfg.doCommonality   = 1;
cfg                 = ca_vba_glm_fitlm(T,cfg);

% ---------------------------
% Perform TFCE thresholding
% ---------------------------
cfg.tfce.path2data  = cfg.outDir;
cfg.tfce.typeStats  = 'tval'; 
cfg.tfce.Ns         = size(cfg.tbl,1);
cfg.tfce.Np         = size(cfg.tbl,2)-1;
cfg.tfce.th         = 1.5;
cfg                 = ca_vba_tfce_threshold(cfg);

%% ------------------------------------------------------------------------
% Extract information for significant TFCE clusters in a results Table 
% -------------------------------------------------------------------------
prefix  = 'tfce150';
cfg     = ca_vba_tfce_resultsTable(cfg,prefix);
fout    = fullfile(cfg.outDir,sprintf('resultsTable_%s.xlsx',prefix));
writetable(cfg.tfce.tableConcat,fout);
save(regexprep(fout,'xlsx','mat'),'cfg','T');

%% ------------------------------------------------------------------------ 
% Clean up some directories
% -------------------------------------------------------------------------
fn = dir(cfg.outDir);fn = fn([fn.isdir]);fn = {fn.name}';
fn(ismember(fn,{'..','.'}))=[];
for idir = 1:numel(fn)
    rmdir(fullfile(cfg.outDir,fn{idir}),'s');
end

%% ROI extraction (based on Clusters for a prespecified Contrast/Coefficient)
% Select clusters based on TFCE resutls for contasts of interest
cfg1         = cfg;
cfg1.conName = {'Age'};% {'Age','Sex'} % Select contrast of interest
cfg1.doMask  = 1; % whether or not to save masks of the clusters
N = ca_vba_tfce_extractROI(cfg1,T);
    
    