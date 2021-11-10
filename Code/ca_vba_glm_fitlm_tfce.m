function [cfg,outDir] = ca_vba_glm_fitlm_tfce(T,cfg) 
%
% A function for group-level multi-modal voxel-wise brain image analysis.
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
% This GLM-like approach uses fitlm to save nii image of coeff and p-values 
% for each variable and residuals (unexplained effects by predictors).
% This could be useful in instances with voxel-specific covariates, e.g. to
% estimate variance explained or residuals in RSFA maps (across subjects)
% after controlling for the effects of ASL maps (voxel-specific), HRV and
% other effects as in 
% Tsvetanov et al 2020 Psyhophysiology (https://doi.org/10.1111/psyp.13714)
% and
% Wu et al 2021 (preprint)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Script edit log (kat)
% ---------------------
% 02/10/2017 - function initiation
% 12/10/2021 - streamlining with WIlkinson notations
% 
% Author : Kamen Tsvetanov, Ph.D., Neurocognitive Ageing
% Affil. : Department of Psychology, University of Cambridge
% Email  : kamen.tsvetanov@gmail.com  
% Website: http://www.kamentsvetanov.com
%__________________________________________________________________________
% Copyright (C) Kamen Tsvetanov 2018


try rootDir         = cfg.rootDir;          catch rootDir       = pwd;  end % Root directory to output analysis results
try doZscore        = cfg.doZscore;         catch doZscore      = 1;    end % Whether or not to z-score data (Default,1)
try doCommonality   = cfg.doCommonality;    catch doCommonality = 0;    end % Whehter or not to run commonality analysis (Default, 0)
try doRobust        = cfg.doRobust;         catch doRobust      = 0;    end % Whether or not to run robust regression (Default, 0)
try numPerm         = cfg.numPerm;          catch numPerm       = 1000; end % Number of permuations
try f_mask          = cfg.f_mask;           catch error('Error. \nPlease provide brain mask.'); end
try Model           = cfg.model;            catch error('Error. \nPlease specify the model using Wilkinson notaions.'); end 


idxNumeric      = vartype('numeric'); 
idxCategorical  = vartype('categorical'); 

% -------------------------
% Get Variable names
% -------------------------
VarNames = strtrim(split(Model,["~","*","+","^"]));
ResponseVar = VarNames(1); 
PredictorVars = VarNames(2:end);

%-------------------------------------
% Identify subjects with missing data
% -----------------------------------
idxSub   = all(~ismissing(T(:,VarNames)),2);
tbl = T(idxSub,VarNames);

%----------------------------
% Import and apply Brain mask 
% ---------------------------
Ymask   = logical(spm_read_vols(spm_vol(f_mask)));
dimMask = size(Ymask);
idxMask = find(Ymask);

% --------------------------------------------------------------
% Import Variables/Modalities with Brain images (e.g. RSFA maps
% ---------------------------------------------------------------
idxMaps         = find(contains(VarNames,'f_'));
VarNamesMaps    = VarNames(idxMaps);
VarNamesGlobal  = VarNames(~contains(VarNames,'f_'));

numVox   = size(idxMask,1);
numSub      = size(tbl,1);
numMap      = numel(idxMaps);

Y = nan(numSub,numVox,numMap);

for iMap = 1:numMap
    nameMaps = VarNames{idxMaps(iMap)};
    V   = spm_vol(char(tbl.(nameMaps)));
    % Check that Mask and this modality dimensions match
    if ~isequal(dimMask,V(1).dim)
        strError = sprintf('Error. \nDimensions of Mask and %s do not match.',nameMaps);
        error(strError);
    end       
    y   = spm_read_vols(V);
    y   = permute(y,[4 1 2 3]);
    Y(:,:,iMap) = y(:,Ymask);
end


if doZscore
    tbl(:,idxNumeric) = normalize(tbl(:,idxNumeric));
    Y = normalize(Y,1);
end

 

% ------------------------------------------------------------------------
% Prepare the table for fitlm with one voxel for dependent and predictor
% maps and set a template to extract data after parpool
% ------------------------------------------------------------------------
Temp = tbl(:,VarNamesGlobal);
Temp{:,VarNamesMaps{:}}   = Y(:,ceil(numVox/2),:);
tbl = Temp;

mlr_temp     = fitlm(tbl,Model);
nameCoef     = mlr_temp.CoefficientNames;
numCoef      = numel(nameCoef);

   
  

% Clean up workspace to imporve parpool performance
clear T X X3D Xvec Y3D Yvec; 



%% make output folders for each  permutation

strModel = regexprep(Model,'f_','');
strModel = regexprep(strModel,'(\<[a-z])','${upper($1)}');
strModel = regexprep(strModel,{'+',' '},'');
strModel = regexprep(strModel,'~','_iv');
strModel = regexprep(strModel,'*','X');
strModel = ['dv' strModel '_nS' num2str(numSub) '_nPerm' num2str(numPerm)];

outDir = fullfile(rootDir,strModel);
mkdir(outDir);

 % make folder for tValue
for i=1:numCoef
    nameOutput = nameCoef{i};
    tValfldDir = fullfile(outDir,['tVal_',nameOutput]); 
    mkdir(tValfldDir);
    bValfldDir = fullfile(outDir,['bVal_',nameOutput]); 
    mkdir(bValfldDir);
 end
% ----------------------------------------------
% Predefine Commonality Analysis output maps
% ----------------------------------------------
if doCommonality
    cfg.mlr     = mlr_temp;
    cfg.doPerm  = 0;
    CA_temp     = ca_stats_commonality(cfg);
    nameVarCA   = CA_temp.Properties.RowNames;
    nameVarCA   = regexprep(nameVarCA,',','');
    numVarCA    = numel(nameVarCA); 
end

randOrder   = palm_quickperms(numSub,[],numPerm);

for iperm = 1:numPerm
    iperm
    
    tempOrder  = randOrder(:,iperm);
    datY         = Y;
    datY(:,:,1)  = datY(tempOrder,:,1); % Shuffle dependent variable
    coef         = zeros(numVox,numCoef);
    pvals        = ones(numVox,numCoef);
    tvals        = zeros(numVox,numCoef);
    Rvals        = zeros(numVox,1);
    if doCommonality
        tvalCA     = zeros(numVox,numVarCA);
        rvalCA     = zeros(numVox,numVarCA);
    end
    
    parfor iVox = 1:numVox
        Yvox = squeeze(datY(:,iVox,:));
        Yvox(:,1) = Yvox(tempOrder,:);
%         if ~isnan(YarrayZ(1,iVox)) && ~YarrayZ(1,iVox)==0 & (~(sum(table2array(Xvox))==0) | ~doPredMaps) 
            tbl_temp = tbl;    
            tbl_temp{:,VarNamesMaps{:}} = Yvox;
            
            mlr             = [];  
            mlr             = fitlm(tbl_temp,Model,'RobustOpts',doRobust);
            coef(iVox,:)    = mlr.Coefficients.Estimate;
            pvals(iVox,:)   = mlr.Coefficients.pValue;
            tvals(iVox,:)   = mlr.Coefficients.tStat;
            Rvals(iVox,:)   = mlr.Rsquared.Adjusted;

            if doCommonality
                cfgtemp     = cfg; cfgtemp.mlr = mlr;
                CA = ca_stats_commonality(cfgtemp);
                tvalCA(iVox,:) = CA{:,'tR2'};
                rvalCA(iVox,:) = CA{:,'Coefficient'};

            end
    %         disp(['Cool, the ' num2str(iVox) 'th voxel was completed!']);
%          end
    end
    % --------------------------------------------------------
    % Write results (coeff, pvals and residuals) to nii images
    % --------------------------------------------------------
    % For every permutation Write results (tvalCA) to nii images in separate folder,
    % so for 1000 permutations there should be 1000 folders
    % (perm_00001,perm_00002 etc), each one containing the tvalCA maps for
    % every effect (common and shared) 
    Vtemp = V(1);
    Vtemp.pinfo(1)=1;
    tcoefMap = zeros(Vtemp.dim);
    bcoefMap = zeros(Vtemp.dim);

    for i=1:numCoef
        nameOutput = nameCoef{i};

        % make folder for tValue
        rootDir = outDir;

        tValfldDir = fullfile(rootDir,['tVal_',nameOutput]); 
%         mkdir(tValfldDir);

        tcoefMap(idxMask) = tvals(:,i);  
        Vtemp.fname = fullfile(tValfldDir,sprintf('glmT_tVal_%s_perm_%.5d.nii',nameOutput,iperm));
        spm_write_vol(Vtemp,tcoefMap);

        % make folder for rValue
        bValfldDir = fullfile(rootDir,['bVal_',nameOutput]); 
%         mkdir(bValfldDir);
        bcoefMap(idxMask) = coef(:,i);  
        Vtemp.fname = fullfile(bValfldDir,sprintf('glmT_bVal_%s_perm_%.5d.nii',nameOutput,iperm));
        spm_write_vol(Vtemp,bcoefMap);
    end 

    if doCommonality
        CAVtemp = V(1);
        CAVtemp.pinfo(1)=1;
        CAtcoefMap = zeros(CAVtemp.dim);
        CArcoefMap = zeros(CAVtemp.dim);
        for i=1:numVarCA 
            nameOutput = nameVarCA{i};

            % make folder for tValue
            rootDir = outDir;
            tValfldDir = fullfile(rootDir,['tValCA_',nameOutput]); 
            mkdir(tValfldDir);
            CAtcoefMap(idxMask) = tvalCA(:,i);  
            CAVtemp.fname = fullfile(tValfldDir,sprintf('glmT_tValCA_%s_perm_%.5d.nii',nameOutput,iperm));
            spm_write_vol(CAVtemp,CAtcoefMap);

            % make folder for rValue
            rValfldDir = fullfile(rootDir,['bValCA_',nameOutput]); 
            mkdir(rValfldDir);
            CArcoefMap(idxMask) = rvalCA(:,i);  
            CAVtemp.fname = fullfile(rValfldDir,sprintf('glmT_rValCA_%s_perm_%.5d.nii',nameOutput,iperm));
            spm_write_vol(CAVtemp,CArcoefMap);
        end 
    end     
end
   
% Assemble and save output structure
cfg.tbl         = tbl;
cfg.randOrder   = randOrder;
cfg.outDir      = outDir;
cfg.mlr         = mlr_temp;
fout = fullfile(outDir,'analysis_cfg.mat')
save(fout,'cfg');

%% Perform TFCE thresholding
% cfg = [];
cfg.tfce.path2data = outDir;
cfg.tfce.typeStats = 'tVal'; 
cfg.tfce.Ns = size(cfg.tbl,1);
cfg.tfce.Np = size(cfg.tbl,2)-1;
ca_vba_tfce_threshold(cfg);

if doCommonality
    cfg.typeStats = 'tValCA';
    ca_vba_tfce_threshold(cfg);
end
