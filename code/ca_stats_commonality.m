 function [commonalityMatrix] = ca_stats_commonality(cfg)
% Function performing Commonality Analysis, which partitions R2 explained
% by all predictors in multiple linear regression into variance unique to
% each predictor and variance shared between each combination of predictors
% See Nimon et al 2013 OrgResMethods and Nimon et al 2008 BehavResMethods
%
% Based on R functions 'aps' and 'commonality', part of calc.yhat package
% https://cran.r-project.org/web/packages/yhat/yhat.pdf
% 
% Usage:
% [CA] = commonality_perm(cfg)
% 
% Parameters (fields in cfg):
% cfg.mlr        - Linear Model output from fitlm
%
% Optional parameters:
% cfg.doPerm      - Whether or not to run permutations (0 | 1; default,0)
% cfg.numPerm     - Number of permuations (default,1)
% cfg.runParfor   - Run usuing parpool 
% cfg.normValue   - Whether to renormalise Variance explained relative to a variable of interest (to normalization value)
% cfg.doRobust    - whether or not to do robust regression
%
% Script edit log (kat)
% -----------------------------
% 21/04/2022 - omit reporting covariates of no interest by prefixing with 'c_'
% 30/05/2024 - Efficient estimation of R2 for each submodel is computed by hand, not using fitlm.
%              Fitlm could be useful if robust regression is needed. Thus,
%              default is the fast implementation.

% --------------------------------------------------------------
% Required packages (https://github.com/kamentsvetanov/external)
% --------------------------------------------------------------
% kat_import('palm'); https://github.com/andersonwinkler/PALM

% Check whether cfg is structure or a linear model
strClass = class(cfg);
switch strClass
    case 'struct'
        mlr = cfg.mlr;
    case 'LinearModel'
        mlr = cfg;
    otherwise
        error('Input should be structure or Linear Model');
end

try parforArg   = cfg.runParfor;    catch parforArg     = 0; end % Run using parpool or not [0 | inf]
try numPerm     = cfg.numPerm;      catch numPerm       = 1; end % Number of permtuations 
try doPerm      = cfg.doPerm;       catch doPerm        = 0; end % Whether to Run permuations [0 | 1]
try normValue   = cfg.normValue;    catch normValue     = 1; end % Whether to renormalise Variance explained relative to a variable of interest (to normalization value)
try doRobust    = cfg.doRobust;     catch doRobust      = 0; end % Whether or not to perform robust regression 

% Set parforArg to inf if runParfor was 1
if parforArg == 1
    parforArg = Inf;
end

% For standard commonality (i.e. no permuations), reset numPerm to 1 and
% parpool argument
if doPerm==0
    numPerm     = 1;
    parforArg   = 0;
end


% Extract information from mlr
dv          = mlr.ResponseName;     % Name of dependent variable
ivlist      = mlr.PredictorNames;   % Predictor names
dataMatrix  = mlr.Variables;        % Data matrix
Ns          = size(dataMatrix,1);   % Number of samples
k           = length(ivlist);       % Determine the number of independent variables (n).
numcc       = 2^k - 1;              % Determine the number of commonality coefficients (2^n-1).
ilist       = find(contains(dataMatrix.Properties.VariableNames,ivlist));
ilist       = mlr.Variables.Properties.VariableNames(ilist);

dataMatrix  = dataMatrix(:,[dv ilist]);

% ----------------------------
% Check for NaN in dataMatrix
% ----------------------------
if any(any(ismissing(dataMatrix)))
    warning('NaNs in the data.');
end

% --------------------------------------------------------
% Generate an ID for each independent variable to 2^(n-1).
% --------------------------------------------------------
ivID = table([2.^(0:(k - 1))]'); 
ivID.Properties.RowNames = ivlist;
nvar = numel(ivID);
if (nvar < 2)
    error('Commonality analysis not conducted. Insufficient number of regressors specified.');
end

% -------------------------------------------------------------------------
% Generate a bitmap matrix containing the bit representation of each 
% commonality coefficient.
% -------------------------------------------------------------------------
PredBitMap = zeros(k, numcc);
for i = 1:numcc
    n = i;
    j = 1;
    while n ~= 0
        if (mod(n,2) == 1)
            PredBitMap(j, i) = 1;
        end
        n = fix(n/2);
        j = j + 1;  
    end
end

[~,apsBitMap] = sort(sum(PredBitMap,1));
PredBitMap = array2table(PredBitMap);
PredBitMap.Properties.RowNames = ivlist;

effectBitMap = PredBitMap;
[~,orderApsBitMap] = sort(apsBitMap);

% -------------------------------------------------------------------------
% Estimate All-Possible-Subsets (APS) Regression
% Use the bitmap matrix to compute the R2 value for each combination of
% independent variables.
% -------------------------------------------------------------------------
APSMatrix = nan(numcc, 2);
Rownames = strings(size(PredBitMap,2),1);

for i = 1:numcc
    model = [dv ' ~ 1 '];
    rownames = [];
     pred  = {};
    for j = 1:k
        bit = PredBitMap{j, i};
        if bit == 1
            model = [model ' + ' ivlist{j}];
            rownames{end+1} = ivlist{j};
            pred  = [pred ivlist{j}];
        end
    end
    Model{i} = model; % Sub-Models for each combination of predictors
    ModelPred{i} = pred; % Predictors in each sub-model
    Rownames{i} = strjoin(rownames,',');
end

% -----------------------
% Create Commonality List
% Use the bitmap matrix to generate the list of R2 values needed for each
% commonality coefficient.
% -----------------------
commonalityList = cell(numcc,1);
for i = 1:numcc
    bit = effectBitMap{1, i};
    if (bit == 1)
        ilist = [0 -ivID{1,1}];
    else
        ilist = ivID{1,1};
    end
    for j = 2:nvar
        bit = effectBitMap{j, i};
        if bit == 1
            alist = ilist;
            blist = commonality_genlist(ilist, -ivID{j,1});
            ilist = [alist blist];
        else
            ilist = commonality_genlist(ilist, ivID{j,1});
        end
    end
    ilist = ilist .* -1;
    commonalityList{i} = ilist;
end

% Predefine order for outter permutations
% % randOrder = palm_quickperms(Ns,[],numPerm);

coeffPerm = nan(numcc,numPerm);
APSr2     = nan(numcc,numPerm);

% TEST: remove covariates of not interest to speed up
% Covariates are indicated by 'c_'
% idxCovs                     = contains(Rownames,'c_');
% Rownames(idxCovs)           = [];
% Model(idxCovs)              = [];
% coeffPerm(idxCovs)          = [];
% APSr2(idxCovs)              = [];
% commonalityList(idxCovs)    = [];
% apsBitMap(idxCovs)          = []; 
% numcc                       = sum(~idxCovs);
% [~,orderApsBitMap]          = sort(apsBitMap);
% apsBitMap(orderApsBitMap)   = 1:numel(apsBitMap); % remap apsBitMap after removing entries for Covariates

parfor (iPerm = 1:numPerm, parforArg) 
% for iPerm = 1:numPerm % Permutation of subject labels. First iteration uses the correct labels
    
    % First iteration uses the correct labels
    if iPerm == 1
        tempOrder = 1:Ns;
    else
        tempOrder = randperm(Ns);
    end
%     tempOrder = randOrder(:,iPerm);
    dataTemp  = dataMatrix;
    dataTemp.(dv) = dataMatrix.(dv)(tempOrder); % Shuffle order of DV
    
    % -------------------------------------------------------------------------
    % Estimate All-Possible-Subsets (APS) Regression
    % Use the bitmap matrix to compute the R2 value for each combination of
    % independent variables.
    % -------------------------------------------------------------------------
%     APSMatrix = array2table(nan(numcc, 2),'RowNames',cellstr(Rownames),'VariableNames',{'k','R2'});
    if doRobust == 1
        APSMatrix = estimateAPSMatrix_mlr(dataTemp,Model,numcc,doRobust);    
    else
        y           = dataMatrix.(dv);
        mu_y        = mean(y);
        APSMatrix = estimateAPSMatrix_fast(dataTemp,ModelPred,numcc,Ns,y,mu_y);
    end


    % Reorder APSMartix accoring to apsBitMap
    % Store the R2 value based on an index that is computed by apsBitMap the IDs of the related IV.
    APSMatrix = APSMatrix(apsBitMap,:);


    % -----------------------------------
    % Now start with Commonality Analysis
    % -----------------------------------
    R2 = APSMatrix(orderApsBitMap);
    commonalityCoeff = nan(numcc,1);
    % -------------------------------------------------------------------------
    % Use the list of R2 values to compute each commonality coefficient and
    % store it in CommonalityMatrix.
    % ------------------------------------------------------------------
    for i = 1:numcc
        r2list = commonalityList{i};
        numlist = length(r2list);
        ccsum = 0;
        for j = 1:numlist
            indexs = r2list(j);
            indexu = abs(indexs);
            if indexu ~= 0
                ccvalue = R2(indexu); % commonalityMatrix{indexu, 'R2'};
                if indexs < 0
                    ccvalue = ccvalue * -1;
                end
                ccsum = ccsum + ccvalue;

            end
        end
        commonalityCoeff(i) = ccsum;
    end
 
    coeffPerm(:,iPerm) = commonalityCoeff; %commonalityMatrix{:,'Coefficient'};
    APSr2(:,iPerm) = R2;
end

varnames = {'Coefficient','PercentTotal','tVal','pVal'};
if  doPerm
    varnames = [varnames ,'pPerm','pPareto','tPerm','tPareto'];
end

%-Estimate T-scores and p-values
%-------------------------------
R2       = abs(coeffPerm(:,1));
r        = real(sqrt(R2));
DF       = mlr.DFE;
tR2      = r ./ sqrt((1-R2)./DF);
pVal     = (1-spm_Tcdf(tR2,DF))*2; % multiply by 2 for two-sided?

% tsquared = DF.*(R2) ./ (1 - R2);
% tsquared = DF.*(r) ./ (1 - r);
% tR = sqrt(tsquared);

commonalityMatrix = array2table(nan(numcc, numel(varnames)),'RowNames',cellstr(Rownames),'VariableNames',varnames);
commonalityMatrix.('Coefficient') = coeffPerm(:,1); % Coefficients reflected variance explained (i.e. R2, not r-value as in GLM)
commonalityMatrix{:,'PercentTotal'} = commonalityMatrix{:,'Coefficient'}./APSr2(numcc,1);
commonalityMatrix.tVal  = tR2;
commonalityMatrix.pVal  = pVal;
commonalityMatrix.APSr2 = APSr2(:,1);

if doPerm
    % Flip sign of distribution if coefficient for true ordering is
    % negative
    signflip = sign(coeffPerm(:,1));
    origdata = coeffPerm(:,1) .* signflip;
    nulldata = coeffPerm .* signflip;
%     origdata = abs(coeffPerm(:,1));
%     nulldata = abs(coeffPerm);
    pPerm    = sum(nulldata' > repmat(origdata',numPerm,1))/numPerm;
    pPareto  = palm_pareto(origdata',nulldata',0,0,0);
    tPareto  = abs(tinv(pPareto/2,Ns));
    tPerm    = abs(tinv(pPerm/2,Ns));
    commonalityMatrix{:,{'pPerm','pPareto','tPerm','tPareto'}} = [pPerm' pPareto' tPerm' tPareto'];
end


commonalityMatrix = commonalityMatrix(apsBitMap, :);
calist = commonalityMatrix.Properties.RowNames;

% Estimate Total effects
commonalityMatrix{end+1,:} = sum(commonalityMatrix{:,:});
commonalityMatrix.Properties.RowNames{end} = 'Total';
commonalityMatrix.APSr2(end) = nan;


% TO DO -Renormalise Variance explained relative to a variable of interest (to
% normaliation value) %

commonalityMatrix.R2Normalised = commonalityMatrix.Coefficient./normValue; 

% ---------------------------------------------
% Assign a direction T-stats of common effects
% ---------------------------------------------
idxcommon = contains(commonalityMatrix.Properties.RowNames,',');
signvar = sign(commonalityMatrix.Coefficient(idxcommon));
idxvars = regexp(commonalityMatrix.Properties.VariableNames,'^t');
idxvars = ~cellfun(@isempty,idxvars);
commonalityMatrix{idxcommon,idxvars} = commonalityMatrix{idxcommon,idxvars}.*repmat(signvar,1,sum(idxvars));

% ---------------------------------------------
% Assign a direction of the main/simple effects
% ---------------------------------------------
for ivar = 1:nvar
    idx=ismember(calist,ivlist(ivar));
    signvar = sign(mlr.Coefficients.tStat(ivar+1));
    commonalityMatrix{idx,:} = commonalityMatrix{idx,:}.*signvar;
end

% -----------------------------------------------------------------
% Exclude variables prefixed by 'c_' (i.e.covariates of no interest)
% -----------------------------------------------------------------
idxRemove = contains(commonalityMatrix.Properties.RowNames,'c_');
commonalityMatrix(idxRemove,:) = [];

%%
function [newlist] = commonality_genlist(ivlist,value)

numlist = length(ivlist);
newlist = ivlist;
newlist = 0;
for i = 1:numlist
    newlist(i) = abs(ivlist(i)) + abs(value);
    if (((ivlist(i) < 0) && (value >= 0)) || ((ivlist(i) >=  0) && (value < 0))) 
        newlist(i) = newlist(i) * -1;
    end
end

%%

function APSMatrix = estimateAPSMatrix_mlr(dataTemp,Model,numcc,doRobust);

APSMatrix = nan(numcc,1);
for i = 1:numcc
%         tic
    mlrTemp  = fitlm(dataTemp,Model{i},'RobustOpts',doRobust);
%         toc
    APSMatrix(i)    = mlrTemp.Rsquared.Ordinary;
%         APSMatrix{i, 'R2'} = mlrTemp.Rsquared.Ordinary;
%         APSMatrix{i, 'k'} = sum(PredBitMap{:,i});
    
    %-(not implemented) Possibly a more efficient way to estimate Total variance explained by the
    % model. Note that it does not work for interactions or squared
    % terms defined in Wilkinson annotation. Instead these should be
    % modelled by the user.
    
%         y = dataTemp.(dv);
%         X = [ones(Ns,1) dataTemp{:,ModelPred{i}}];
%         mu_y    = mean(y);
%         if doRobust
%             [betas stats]= robustfit(X,y,[],[],'off');
%             ss_res  = sum(stats.resid.^2);
%         else
% %             [b,bint,r,rint,stats] = regress(dataTemp.(dv),dataTemp{:,ModelPred{i}});
%             betas   = X\y; % add constant term to make identical to robust fit and fitlm
%             ss_res  = sum((X*betas - y).^2);
%         end
%         ss_tot  = sum((y - mu_y).^2);
%         R2      = 1 - ss_res/ss_tot;
%         APSMatrix(i) = 1 - ss_res/ss_tot;
%         tic
%         y = dataTemp.(dv);
%         X = dataTemp{:,ModelPred{i}};
%         mu_y    = mean(y);
%         betas   = x\y;
%         ss_tot  = sum((y - mu_y).^2);
%         ss_res  = sum((X*betas - y).^2);
%         R2      = 1 - ss_res/ss_tot;
%         toc
    
    
end

%%

function  APSMatrix = estimateAPSMatrix_fast(dataTemp,ModelPred,numcc,Ns,y,mu_y)

APSMatrix = nan(numcc,1);
for i = 1:numcc
    % Calculate R2 by hand
    X = [ones(Ns,1) dataTemp{:,ModelPred{i}}];
    betas   = X\y;
    ss_tot  = sum((y - mu_y).^2);
    ss_res  = sum((X*betas - y).^2);
    APSMatrix(i) = 1 - ss_res/ss_tot;
end
