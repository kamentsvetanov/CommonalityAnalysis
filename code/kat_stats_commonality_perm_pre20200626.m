function [commonalityMatrix] = kat_stats_commonality_perm(cfg)
% Function performing Commonality Analysis, which partitions R2 explained
% by all predictors in multiple linear regression into variance unique to
% each predictor and variance shared between each combination of predictors
% See Nimon et al 2013 OrgResMethods and Nimon et al 2008 BehavResMethods
%
% Canibalized from R functions 'aps' and 'commonality', part of calc.yhat package 
% 
% Three types of inputs:
% 1 - Linear Model such as an output from fitlm (1 input needed)
% 2 - Data Matrix with formula string, input 1 and 2, respectively 
% 3 - Data matrix with names of Dependent and Independnt variables (inputs 2 and 3, respectively) 

try mlr = cfg.mlr; end
try 
    try 
        mlr = fitlm(cfg.dat,cfg.model);
    catch
        mlr = fitlm(cfg.dat,'ResponseVar',cfg.nameDV,'PredictorVars',cfg.nameIV);
    end
end


try runInSerial = cfg.runInSerial; catch runInSerial = 1;      end
try numPerm = cfg.numPerm;         catch numPerm = 1;          end
try doPerm = cfg.doPerm;           catch doPerm =0;            end

if runInSerial
    parforArg = 0;
else
    parforArg = Inf;
end

if ~doPerm 
    numPerm = 1;
    parforArg = 0;
end

    
dv = mlr.ResponseName;
ivlist = mlr.PredictorNames;
dataMatrix = mlr.Variables;    
Ns    = mlr.NumObservations;
k     = length(ivlist); % Determine the number of independent variables (n).
numcc = 2^k - 1; % Determine the number of commonality coefficients (2^n-1).
ilist = find(contains(dataMatrix.Properties.VariableNames,ivlist));
ilist = mlr.Variables.Properties.VariableNames(ilist);

dataMatrix = dataMatrix(:,[dv ilist]);
ivID = table([2.^(0:(k - 1))]'); % Generate an ID for each independent variable to 2^(n-1).
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
    for j = 1:k
        bit = PredBitMap{j, i};
        if bit == 1
            model = [model ' + ' ivlist{j}];
            rownames{end+1} = ivlist{j};
        end
    end
    Model{i} = model;
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
            blist = kat_stats_commonality_genlist(ilist, -ivID{j,1});
            ilist = [alist blist];
        else
            ilist = kat_stats_commonality_genlist(ilist, ivID{j,1});
        end
    end
    ilist = ilist .* -1;
    commonalityList{i} = ilist;
end



%
% Predefine order for outter permutations
randOrder = palm_quickperms(Ns,[],numPerm);

coeffPerm = nan(numcc,numPerm);
APSr2     = nan(numcc,numPerm);

parfor (iPerm = 1:numPerm, parforArg)
%     for iPermOut = 1:numPerm % Permutation of subject labels. First iteration uses the correct labels

    tempOrder = randOrder(:,iPerm);
    dataTemp  = dataMatrix;
    dataTemp.(dv) = dataMatrix.(dv)(tempOrder);
    
    % -------------------------------------------------------------------------
    % Estimate All-Possible-Subsets (APS) Regression
    % Use the bitmap matrix to compute the R2 value for each combination of
    % independent variables.
    % -------------------------------------------------------------------------
    
%     APSMatrix = array2table(nan(numcc, 2),'RowNames',cellstr(Rownames),'VariableNames',{'k','R2'});
    APSMatrix = nan(numcc,1);
    for i = 1:numcc
        mlrTemp  = fitlm(dataTemp,Model{i});
        APSMatrix(i)    = mlrTemp.Rsquared.Ordinary;
%         APSMatrix{i, 'R2'} = mlrTemp.Rsquared.Ordinary;
%         APSMatrix{i, 'k'} = sum(PredBitMap{:,i});
    end

    % Reorder APSMartix accoring to apsBitMap
    % Store the R2 value based on an index that is computed by apsBitMap the IDs of the related IV.
    APSMatrix = APSMatrix(apsBitMap,:);


    % -----------------------------------
    % Now start with Commonality Analysis
    % -----------------------------------
%     commonalityMatrix = array2table(nan(numcc, 3),'RowNames',cellstr(Rownames),'VariableNames',{'R2','Coefficient','PercentTotal'});
%     commonalityMatrix(:,'R2') = APSMatrix(orderApsBitMap,'R2');
%     R2 = APSMatrix{orderApsBitMap,'R2'};
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
                ccvalue = R2(indexu);%commonalityMatrix{indexu, 'R2'};
                if indexs < 0
                    ccvalue = ccvalue * -1;
                end
                ccsum = ccsum + ccvalue;

            end
        end
        
%         commonalityMatrix{i, 'Coefficient'} = ccsum;
        commonalityCoeff(i) = ccsum;
    end
%     commonalityMatrix(:,'R2') = [];
%     commonalityMatrix = commonalityMatrix(apsBitMap, :);
%     commonalityMatrix{:,'PercentTotal'} = commonalityMatrix{:,'Coefficient'}./APSMatrix{numcc,2};
%     
    coeffPerm(:,iPerm) = commonalityCoeff; %commonalityMatrix{:,'Coefficient'};
    APSr2(:,iPerm) = R2;
%     dataTemp = [];
end

varnames = {'Coefficient','PercentTotal','tR2'};
if doPerm
    varnames = [varnames ,'pPerm','pPareto','tPerm','tPareto'];
end

R2       = coeffPerm(:,1);
r        = sqrt(R2);
DF       = mlr.DFE;
tR2       = r ./ sqrt((1-R2)./DF);

commonalityMatrix = array2table(nan(numcc, numel(varnames)),'RowNames',cellstr(Rownames),'VariableNames',varnames);
commonalityMatrix.('Coefficient') = coeffPerm(:,1);
commonalityMatrix{:,'PercentTotal'} = commonalityMatrix{:,'Coefficient'}./APSr2(numcc,1);
commonalityMatrix.tR2 = tR2;

if doPerm
    nulldata = coeffPerm(:,2:end);
    origdata = coeffPerm(:,1);
    pPerm    = sum(nulldata' > repmat(origdata',numPerm-1,1))/numPerm;
    pPareto  = palm_pareto(origdata',nulldata',0,0,0);
    tPareto  = abs(tinv(pPareto/2,DF));
    tPerm    = abs(tinv(pPerm/2,DF));
    commonalityMatrix{:,{'pPerm','pPareto','tPerm','tPareto','tR2'}} = [pPerm' pPareto' tPerm' tPareto' tR2];
end
commonalityMatrix = commonalityMatrix(apsBitMap, :);



            
