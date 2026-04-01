% Code showing the use of Voxel-based Analysis (VBA) with fitlm and
% commonality function using SLURM


% -------------------------------------------------------------------------
% 1. Create a spreadsheet file (csv format) containign subject information
% with variable names specified for the statistical model. For example, for
% a model
%     Model = 'f_t1w ~ age + sex'; 
% The following variable names are expected f_t1w, age and sex

% -------------------------------------------------------------------------
% 2. Set up shell script to perform commonality analysis for a specific model
% using SLURM
% See vba_demo_slurm_hpc.sh as an example
% Usage:
% sbatch vba_demo_slurm_hpc.sh
% -------------------------------------------------------------------------

% ---------------------------
% 3. Perform TFCE thresholding
% ---------------------------
addpath(genpath('/rds/user/kat35/hpc-work/projects/public-code/CommonalityAnalysis/code'));
addpath(genpath('/rds/user/kat35/hpc-work/projects/external/mat/spm12'));
cfg.tfce.path2data  = cfg.outDir;
cfg.tfce.typeStats  = 'tval'; 
cfg.tfce.Ns         = size(cfg.tbl,1);
cfg.tfce.Np         = size(cfg.tbl,2)-1;
cfg.tfce.th         = 1.97;
cfg.tfce.doRunInSerial = 1;
cfg                 = vba_tfce_threshold(cfg);


%% ------------------------------------------------------------------------
% 4. Extract information for significant TFCE clusters in a results Table 
% -------------------------------------------------------------------------
prefix      = 'tfce197'; % Number defined in cfg.tfce.th
clustersize = 8;
cfg         = vba_tfce_resultsTable(cfg,prefix,clustersize);
fout        = fullfile(cfg.outDir,sprintf('resultsTable_%s.xlsx',prefix));
writetable(cfg.tfce.tableConcat,fout);
save(regexprep(fout,'xlsx','mat'),'cfg');

%% ------------------------------------------------------------------------
% 5. Extract subject-specific data for peaks from significant TFCE clusters 
% -------------------------------------------------------------------------
[N] = vba_tfce_extractROI(cfg);
