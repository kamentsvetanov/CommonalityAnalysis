% Code showing the use of Voxel-based Analysis with fitlm and commonality function.
%
% Running the analysis requires two variables in the workspace, T and
% Model, similar to the usage of 'fitlm'
% 
% T       - Table containing all variables specified in Model. Similar to
%           'fitlm', variables in T can be continuous or categorical
%           variables. For voxel-based analysis is required to specify
%           filepaths to individuals' images (e.g. anatomical T1w scans) in
%           a variable prefixed by 'f_', e.g. f_t1w for a structure
%           variable in T containing filepaths to anatomical scans.
%
% Model   - A string specifying the linear model formula using Wilkinson notation.
%           (https://uk.mathworks.com/help/stats/wilkinson-notation.html)
%           y ~ terms, where y is the name of the response variable, and
%           terms defines the model using the predictor variable names and
%           the operators.
%           Images (e.g. anatomical T1w scans) can be used as response and/
%           or predictor variables. 
%           Model = 'f_t1w ~ age + sex'; Specifies analysis to predict T1w
%           intensity using age and sex as predictors.
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

Model = 'f_rsfa ~ Age + Sex';

% -------------------------------------------------
% Assemble cfg structure needed to run the analysis
% -------------------------------------------------
cfg                 = [];
cfg.model           = Model;
cfg.rootDir         = '/imaging/camcan/sandbox/kt03/temp/'; 
cfg.f_mask          = fullfile(rootdir,'mask.nii');
cfg.numPerm         = 4;
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

%% ------------------------------------------------------------------------ 
% Clean up some directories
% -------------------------------------------------------------------------
fn = dir(cfg.outDir);fn = fn([fn.isdir]);fn = {fn.name}';
fn(ismember(fn,{'..','.'}))=[];
for idir = 1:numel(fn)
    rmdir(fullfile(cfg.outDir,fn{idir}),'s');
end
