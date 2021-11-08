% data from the Holzinger and Swineford (1939) study, readily
% available in the MBESS R package
dat = readtable('/imaging/camcan/sandbox/kt03/projects/private-code/kamen/statistical/commonality/hsdata.csv');

nameDV = 'paragrap';
nameIV = {'general'  'sentence' 'wordc'    'wordm'};

DAT = dat(:,{nameDV nameIV{:}});
DAT = normalize(DAT);
model = 'paragrap ~ general + sentence + wordc + wordm';

% Run MLR using Formula where the terms are in Wilkinson notation
mlr = fitlm(DAT,model);

% Run CA using fitlm output
CA = kat_stats_commonality(mlr);

% Run CA using data matrix and Formula where the terms are in Wilkinson
% notation
CA = kat_stats_commonality(dat,model);

% Run CA using a data matrix with variable names of dependent variable (DV) 
% and multipe independent variables (IV)
CA = kat_stats_commonality(dat,nameDV,nameIV);


cfg     = [];
cfg.mlr = mlr;
% cfg.dat = dat;
% cfg.model = model;
% cfg.nameDV = nameDV;
% cfg.nameIV = nameIV;
cfg.doPerm = 0;
cfg.numPerm = 100;
% cfg.normValue = normValue
cfg.runInSerial = 0;
tic
CA1 = kat_stats_commonality_perm(cfg);
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
