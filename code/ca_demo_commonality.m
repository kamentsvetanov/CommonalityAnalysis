% Code showing the use of commonality function for a single model. To run
% voxel-level commonality analysis refer to demo_vba.m

% Load data from the Holzinger and Swineford (1939) study, available in the MBESS R package


clear 
% Ensure CommonalityAnalysis is in Matlab's path
rootdir = fileparts(fileparts(which('ca_demo_vba.m')));
dat = readtable([rootdir '/data/hsdata.csv']);

nameDV  = 'paragrap';
nameIV  = {'general'  'sentence' 'wordc'    'wordm'};
DAT     = dat(:,{nameDV nameIV{:}});
DAT     = normalize(DAT);

% Specify Linear Model using Wilkinson notation
model = 'paragrap ~ general + sentence + wordc + wordm';

% Run MLR using Formula where the terms are in Wilkinson notation
mlr = fitlm(DAT,model);


% Run Commonality Analysis using fitlm output
tic
CA = ca_stats_commonality(mlr);
toc

tic
CA1 = ca_stats_commonality(mlr);
toc

% Run Commonality Analysis using Permutations
cfg             = [];
cfg.mlr         = mlr;
cfg.doPerm      = 1;
cfg.numPerm     = 1000;
cfg.runParfor   = 0;
tic
CA1 = ca_stats_commonality(cfg);
toc

% Run Commonality Analysis using Permutations nad Robust regression
cfg             = [];
cfg.mlr         = mlr;
cfg.doPerm      = 1;
cfg.numPerm     = 100;
cfg.runParfor   = 0;
cfg.doRobust    = 1;
tic
CA1 = ca_stats_commonality(cfg);
toc

