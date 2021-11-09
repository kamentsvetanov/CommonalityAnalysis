% Code showing the use of commonality function.

% Load data from the Holzinger and Swineford (1939) study, available in the MBESS R package
dat = readtable('/home/kt03/Projects/public-code/CommonalityAnalysis/data/hsdata.csv');

nameDV  = 'paragrap';
nameIV  = {'general'  'sentence' 'wordc'    'wordm'};
DAT     = dat(:,{nameDV nameIV{:}});
DAT     = normalize(DAT);

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

%% 
% 
% tbl = mlr.Variables;
% modelR = 'Cattellcont ~ 1 + G + E + Age + M + H';
% mlrR =  fitlm(tbl,modelR);
% 
% 
% cfg     = [];
% cfg.mlr = mlrR;
% % cfg.dat = dat;
% % cfg.model = model;
% % cfg.nameDV = nameDV;
% cfg.nameIV = nameIV;
% cfg.doPerm = 0;
% cfg.numPerm = 100;
% tic
% CAr = kat_stats_commonality_perm(cfg);
% toc
% idxAgeU = ismember(CAr.Properties.RowNames,'Age');
% idxAgeU = ismember(CAr.Properties.RowNames,'Age');
% normValue = CAr.Coefficient(idxAgeU);
% 
% 
