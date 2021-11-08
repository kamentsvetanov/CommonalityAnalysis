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

try mlr = cfg.mlr; catch mlr = cfg; end
try 
    try 
        mlr = fitlm(cfg.dat,cfg.model);
    catch
        mlr = fitlm(cfg.dat,'ResponseVar',cfg.nameDV,'PredictorVars',cfg.nameIV);
    end
end

try runInSerial = cfg.runInSerial; catch runInSerial = 1;      end % Run using parpool or not
try numPerm = cfg.numPerm;         catch numPerm = 1;          end % Number of permtuations
try doPerm = cfg.doPerm;           catch doPerm = 0;           end % Whether to Run permuations
try normValue = cfg.normValue;     catch normValue = 1;        end % Whether to renormalise Variance explained relative to a variable of interest (to normaliation value)
 

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
Ns    = size(dataMatrix,1);
k     = length(ivlist); % Determine the number of independent variables (n).
numcc = 2^k - 1; % Determine the number of commonality coefficients (2^n-1).
ilist = find(contains(dataMatrix.Properties.VariableNames,ivlist));
ilist = mlr.Variables.Properties.VariableNames(ilist);

dataMatrix = dataMatrix(:,[dv ilist]);

% ----------------------------
% Check for NaN in dataMatrix
% ----------------------------
if any(any(ismissing(dataMatrix)))
    warning('NaNs in the data.');
end


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
% % randOrder = palm_quickperms(Ns,[],numPerm);

coeffPerm = nan(numcc,numPerm);
APSr2     = nan(numcc,numPerm);

parfor (iPerm = 1:numPerm, parforArg)
% for iPerm = 1:numPerm % Permutation of subject labels. First iteration uses the correct labels

    if iPerm == 1
        tempOrder = 1:Ns;
    else
        tempOrder = randperm(Ns);
    end
%     tempOrder = randOrder(:,iPerm);
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
                ccvalue = R2(indexu); % commonalityMatrix{indexu, 'R2'};
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
if  doPerm
    varnames = [varnames ,'pPerm','pPareto','tPerm','tPareto'];
end

R2       = abs(coeffPerm(:,1));
r        = real(sqrt(R2));
DF       = mlr.DFE;
tR2       = r ./ sqrt((1-R2)./DF);

% tsquared = DF.*(R2) ./ (1 - R2);
% tsquared = DF.*(r) ./ (1 - r);
% tR = sqrt(tsquared);

commonalityMatrix = array2table(nan(numcc, numel(varnames)),'RowNames',cellstr(Rownames),'VariableNames',varnames);
commonalityMatrix.('Coefficient') = coeffPerm(:,1);
commonalityMatrix{:,'PercentTotal'} = commonalityMatrix{:,'Coefficient'}./APSr2(numcc,1);
commonalityMatrix.tR2 = tR2;
commonalityMatrix.APSr2 = APSr2;

if doPerm
    % Flip sign of distribution if coefficient for true ordering is
    % negative
    signflip = sign(coeffPerm(:,1));
    origdata = coeffPerm(:,1) .* signflip;
    nulldata = coeffPerm .* signflip;
    pPerm    = sum(nulldata' > repmat(origdata',numPerm,1))/numPerm;
    pPareto  = palm_pareto(origdata',nulldata',0,0,0);
    tPareto  = abs(tinv(pPareto/2,Ns));
    tPerm    = abs(tinv(pPerm/2,Ns));
    commonalityMatrix{:,{'pPerm','pPareto','tPerm','tPareto','tR2'}} = [pPerm' pPareto' tPerm' tPareto' tR2];
end


commonalityMatrix = commonalityMatrix(apsBitMap, :);
calist = commonalityMatrix.Properties.RowNames;

% Estimate Total effects
commonalityMatrix{end+1,:} = sum(commonalityMatrix{:,:});
commonalityMatrix.Properties.RowNames{end} = 'Total';
commonalityMatrix.APSr2(end) = nan;


% Renormalise Variance explained relative to a variable of interest (to
% normaliation value)

commonalityMatrix.R2Normalised = commonalityMatrix.Coefficient./normValue;

%  Total variance explaine by Age Unique and Shared
 nameNormValue = 'Age';
 idx=contains(calist,nameNormValue);
 commonalityMatrix{end+1,:} = sum(commonalityMatrix{idx,:});
 commonalityMatrix.Properties.RowNames{end} = 'AgeUS';
 

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


function [newlist] = kat_stats_commonality_genlist(ivlist,value)

numlist = length(ivlist);
newlist = ivlist;
newlist = 0;
for i = 1:numlist
    newlist(i) = abs(ivlist(i)) + abs(value);
    if (((ivlist(i) < 0) && (value >= 0)) || ((ivlist(i) >=  0) && (value < 0))) 
        newlist(i) = newlist(i) * -1;
    end
end
