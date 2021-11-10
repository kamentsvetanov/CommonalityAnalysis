% Code showing the use of commonality function for a single model. To run
% voxel-level commonality analysis refer to demo_vba.m

% Load data from the Holzinger and Swineford (1939) study, available in the MBESS R package
dat = readtable('/home/kt03/Projects/public-code/CommonalityAnalysis/data/hsdata.csv');

nameDV  = 'paragrap';
nameIV  = {'general'  'sentence' 'wordc'    'wordm'};
DAT     = dat(:,{nameDV nameIV{:}});
DAT     = normalize(DAT);

% Specify Linear Model using Wilkinson notation
model = 'paragrap ~ general + sentence + wordc + wordm';

% Run MLR using Formula where the terms are in Wilkinson notation
mlr = fitlm(DAT,model);


% Run Commonality Analysis using fitlm output
CA = commonality(mlr);

% Run Commonality Analysis using Permutations
cfg             = [];
cfg.mlr         = mlr;
cfg.doPerm      = 1;
cfg.numPerm     = 100;
cfg.runParfor   = 0;
tic
CA1 = commonality(cfg);
toc
