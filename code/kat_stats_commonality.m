function [commonalityMatrix] = kat_stats_commonality(varargin)
% Function performing Commonality Analysis, which partitions R2 explained
% by all predictors in multiple linear regression into variance unique to
% each predictor and variance shared between each combination of predictors
% See Nimon et al 2013 OrgResMethods and Nimon et al 2008 BehavResMethods
%
% Canibalized from R functions 'aps' and 'commonality', part of calc.yhat package 
% 
% Three types of inputs:
% 1 - Linear Model such as an output from fitlm (1 input needed)
%     mlr = fitlm(dat,model);
%     CA = kat_stats_commonality(mlr);
%
% 2 - Data Matrix with formula string, input 1 and 2, respectively 
%     CA = kat_stats_commonality(dat,model);
%
% 3 - Data matrix with names of Dependent and Independnt variables (inputs 2 and 3, respectively) 
%     CA = kat_stats_commonality(dat,nameDV,nameIV);

strClass = class(varargin{1});

switch strClass
    case 'LinearModel'
        mlr = varargin{1};
    case 'table'
        if contains(varargin{2},'~')
            mlr = fitlm(varargin{1},varargin{2});
        else
            mlr = fitlm(varargin{1},'ResponseVar',varargin{2},'PredictorVars',varargin{3});
        end
end
   

dv = mlr.ResponseName;
ivlist = mlr.PredictorNames;
dataMatrix = mlr.Variables;        
k     = length(ivlist); % Determine the number of independent variables (n).
numcc = 2^k - 1; % Determine the number of commonality coefficients (2^n-1).
ilist = find(contains(dataMatrix.Properties.VariableNames,ivlist));
ilist = mlr.Variables.Properties.VariableNames(ilist);

dataMatrix = dataMatrix(:,[dv ilist]);
ivID = table([2.^(0:(k - 1))]'); % Generate an ID for each independent variable to 2^(n-1).
ivID.Properties.RowNames = ivlist;


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

% -------------------------------------------------------------------------
% Estimate All-Possible-Subsets (APS) Regression
% Use the bitmap matrix to compute the R2 value for each combination of 
% independent variables.
% -------------------------------------------------------------------------
APSMatrix = nan(numcc, 2);
rownames = strings(size(PredBitMap,2),1);

for i = 1:numcc
    model = [dv ' ~ 1 '];
    for j = 1:k
        bit = PredBitMap{j, i};
        if bit == 1
            model = [model ' + ' ivlist{j}];
        end
    end
    mlrTemp  = fitlm(dataMatrix,model);
    APSMatrix(i, 2) = mlrTemp.Rsquared.Ordinary;
    APSMatrix(i, 1) = sum(PredBitMap{:,i});
    rownames{i} = strjoin(mlrTemp.PredictorNames,',');
end


APSMatrix = array2table(APSMatrix,'RowNames',cellstr(rownames),'VariableNames',{'k','R2'});
% Reorder APSMartix accoring to apsBitMap
% Store the R2 value based on an index that is computed by apsBitMap the IDs of the related IV.
APSMatrix = APSMatrix(apsBitMap,:);

% -----------------------------------
% Now start with Commonality Analysis
% -----------------------------------
nvar = numel(ivID);

if (nvar < 2) 
    error('Commonality analysis not conducted. Insufficient number of regressors specified.');
end

effectBitMap = PredBitMap;
[~,orderApsBitMap] = sort(apsBitMap);
commonalityMatrix = nan(numcc, 3);
commonalityMatrix(:,1) = APSMatrix{orderApsBitMap,2};

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
            ccvalue = commonalityMatrix(indexu, 1);
            if indexs < 0
                ccvalue = ccvalue * -1;
            end
            ccsum = ccsum + ccvalue;
            
        end
    end
    commonalityMatrix(i, 2) = ccsum;
end

commonalityMatrix(:,1) = [];
commonalityMatrix(:,1) = commonalityMatrix(apsBitMap, 1);
commonalityMatrix(:,2) = commonalityMatrix(:,1)./APSMatrix{numcc,2};
commonalityMatrix = array2table(commonalityMatrix,'RowNames',cellstr(rownames(apsBitMap)),'VariableNames',{'Coefficient','PercentTotal'});
  