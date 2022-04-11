function [cfg,outDir] = ca_vba_glm_fitlm(T,cfg) 
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
% Wu et al 2021 (https://www.biorxiv.org/content/10.1101/2021.11.10.468042v1)
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
% 10/01/2022 - estimate main effects as 'One-Sample' by specifying 'y ~ 1'
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


modeltype       = 'regression'; % 'one-sample' also possible if 'y~1'

idxNumeric      = vartype('numeric'); 
idxCategorical  = vartype('categorical'); 

% -------------------------
% Get Variable names
% -------------------------
VarNames = strtrim(split(Model,["~","*","+","^"]));
ResponseVar = VarNames(1); 
PredictorVars = VarNames(2:end);

% -------------------------
% '1' in PredictorVars requires estimation of main effect, e.g. do not
% zscore data
% --------------------------
idx1 = ismember(VarNames,'1'); 
if any(idx1)
    VarNames(idx1)  = [];
    doZscore        = 0;
    modeltype       = 'one-sample'; 
end
    

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
Temp{:,VarNamesMaps} = squeeze(Y(:,ceil(numVox/2),:));
tbl = Temp;

mlr_temp     = fitlm(tbl,Model);
nameCoef     = mlr_temp.CoefficientNames;
numCoef      = numel(nameCoef);

% Clean up workspace to imporve parpool performance
clear T X X3D Xvec Y3D Yvec; 

%% make output folders for each  permutation

strModel = regexprep(Model,'f_','');
strModel = regexprep(strModel,'c_','');
strModel = regexprep(strModel,'(\<[a-z])','${upper($1)}');
strModel = regexprep(strModel,{'+',' '},'');
strModel = regexprep(strModel,'~','_iv');
strModel = regexprep(strModel,'*','X');
strModel = ['dv' strModel '_n' num2str(numSub) '_nPerm' num2str(numPerm)];

% ----------------------------------------------
% Predefine Commonality Analysis or GLM output maps
% ----------------------------------------------
if doCommonality
    outDir      = fullfile(rootDir,['ca_' strModel]);
    cfg.mlr     = mlr_temp;
    cfg.doPerm  = 0;
    CA_temp     = ca_stats_commonality(cfg);
    nameVarCA   = CA_temp.Properties.RowNames;
    nameVarCA   = regexprep(nameVarCA,',','');
    numVarCA    = numel(nameVarCA);    
    nameCoef    = nameVarCA;
    numCoef     = numVarCA;
else
    outDir = fullfile(rootDir,['glm_' strModel]);
end

% -----------------
% Make folders
% -----------------
for i=1:numCoef
    if ~contains(nameCoef{i},'c_') 
        nameOutput = regexprep(nameCoef{i},'f_','');
        mkdir(fullfile(outDir,['tval_',nameOutput]));
    %     mkdir(fullfile(outDir,['bval_',nameOutput]));
    end
end

%% ------------------------------------------------------------------------
% Permutations 
randOrder   = palm_quickperms(numSub,[],numPerm);
switch modeltype
    case 'one-sample' % INCOMPLETE. Estimate Main/Average/Group effect using one-sample t-test
        cfg_temp = cfg;
        ca_vba_util_onesample(cfg_temp);
       
    case 'regression'
        
        for iperm = 1:numPerm
%             iperm

            tempOrder   = randOrder(:,iperm);
            datY        = Y;
            tvals       = zeros(numVox,numCoef);

            if doCommonality
                % Loop through all voxels
                parfor iVox = 1:numVox
                    Yvox = squeeze(datY(:,iVox,:));
                    % Ensure all subjects have non-nan values
                    if sum(sum(isnan(Yvox)))==0
                        tbl_temp = tbl;    
                        tbl_temp{:,VarNamesMaps} = Yvox;
                        tbl_temp(:,ResponseVar) = tbl_temp(tempOrder,ResponseVar);% Shuffle dependent variable

                        mlr             = [];  
                        mlr             = fitlm(tbl_temp,Model,'RobustOpts',doRobust);
                        cfgtemp         = cfg; 
                        cfgtemp.mlr     = mlr;
                        CA              = ca_stats_commonality(cfgtemp);
                        tvals(iVox,:)   = CA{:,'tR2'};
        %                 bvals(iVox,:)   = CA{:,'Coefficient'};
                    end           
                end

            % ---------------------------------------    
            % Perform and save standard GLM analysis
            % ---------------------------------------
            else
               parfor iVox = 1:numVox
                    Yvox = squeeze(datY(:,iVox,:));
        %             Yvox = Yvox(tempOrder,:);

                    % Ensure all subjects have non-nan values
                    if sum(sum(isnan(Yvox)))==0
                        tbl_temp = tbl;    
                        tbl_temp{:,VarNamesMaps} = Yvox;
                        tbl_temp(:,ResponseVar) = tbl_temp(tempOrder,ResponseVar);% Shuffle dependent variable

                        mlr             = [];  
                        mlr             = fitlm(tbl_temp,Model,'RobustOpts',doRobust);
        %                 bvals(iVox,:)   = mlr.Coefficients.Estimate;
                        tvals(iVox,:)   = mlr.Coefficients.tStat;
                    end
                end
            end

            % --------------------------------------------------------
            % Write results (coeff, pvals and residuals) to nii images
            % --------------------------------------------------------
            % For every permutation Write results (tvalCA) to nii images in separate folder,
            % so for 1000 permutations there should be 1000 folders
            % (perm_00001,perm_00002 etc), each one containing the tvalCA maps for
            % every effect (common and shared) 

            for i=1:numCoef
                % Output coefficients not containing covariates 
                % (i.e. predictors prefixed by 'c_')
                if ~contains(nameCoef{i},'c_') 
                    nameOutput      = regexprep(nameCoef{i},'f_','');
                    Vtemp           = V(1);
                    Vtemp.pinfo(1)  =1;
                    tmap            = zeros(Vtemp.dim);
            %         bmap            = zeros(Vtemp.dim);

                    % Save t-stats
                    tdir            = fullfile(outDir,['tval_',nameOutput]); 
                    tmap(idxMask)   = tvals(:,i);  
                    Vtemp.fname     = fullfile(tdir,sprintf('results_null_%.5d.nii',iperm));
                    spm_write_vol(Vtemp,tmap);

                    % Save beta coefficients (r-values if data was z-scoredd)
            %         bdir            = fullfile(outDir,['bval_',nameOutput]); 
            %         bmap(idxMask)   = bvals(:,i);  
            %         Vtemp.fname     = fullfile(bdir,sprintf('results_null_%.5d.nii',iperm));
            %         spm_write_vol(Vtemp,bmap);
                end
            end  
        end
end

% ----------------------------------
% Assemble and save output structure
% ----------------------------------
cfg.tbl         = tbl;
cfg.randOrder   = randOrder;
cfg.outDir      = outDir;
cfg.mlr         = mlr_temp;
fout = fullfile(outDir,'analysis_cfg.mat');
save(fout,'cfg');
