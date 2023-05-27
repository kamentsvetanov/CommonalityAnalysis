function [N] = ca_vba_tfce_extractROI(cfg,T)
% Extract Cluster values for each participant and brainmap type for a given
% coordinate as specified in TFCE results table
%
% Inputs:
% ----------
% cfg         - cfg variable with varios fields defined in ca_vba_tfce_resultsTable. 
%
% cfg.doMasks - Save masks for each cluster (Optional, Default - no save)
%
% cfg.conName - Contrast names of interest used to define Clusters
%               (optional)
%
% T           - Subject-specific T table
%
% Outputs:
% -----------
%   cfg   - same as input inluding additional subfiedl statTable.
%   
%   N     - Same as table T where cluster values are appended as additional
%           variables with the following naming convention
%           Modality_Contrast_Hemisphere_RegionName
%
% Author : Kamen Tsvetanov, Ph.D., Neurocognitive Ageing
% Affil. : Department of Clinical Neurosciences, University of Cambridge
% Email  : kamen.tsvetanov@gmail.com  
% Website: http://www.kamentsvetanov.com
% Date   : 27 May 2023; Last revision: 
%__________________________________________________________________________
% Copyright (C) 2023, Kamen Tsvetanov

try doMasks = cfg.doMasks; catch doMasks = 0; end
try contr   = cfg.conName; catch contr   = {}; end

% Make a new table across all coefficients and clusters
% -----------------------------------------------------
namecoeff = fieldnames(cfg.tfce.results);
R = table();
for icoeff = 1:numel(namecoeff)
    coeff  = string(namecoeff{icoeff});
    tbl = cfg.tfce.results.(coeff);
    
    if ~ismissing(tbl.cluster_name) % Check if table is not empty
        tbl.coeff_name = repmat(coeff,size(tbl,1),1);
        R      = [R;tbl];
    end
end
R.id = [1:size(R,1)]';

% Extract subject average values for each cluster in each imaging data set
% ------------------------------------------------------------------------
varnameMaps = regexp(cfg.model,('f_\w*'),'match');

% Create Cluster names based on coefficient name, cluster number and node name
nodenames = strcat(R.coeff_name,'_',R.cluster_name);
nodenames = regexprep(nodenames,'\.','_');
nodenames = regexprep(nodenames,{' ','-'},'');

% (optional) Select a subsest of contrasts, as some are not of interest
% ---------------------------------------------------------------------
% contr     = {'Age'}; % Here we focus only clusters with effects in Ageing
% contr     = {'Age','Sex'}; % Here we focus only clusters with effects in Ageing and Sex
if ~isempty(contr)
    idx       = contains(nodenames,strcat(contr,'_'));
    nodenames = nodenames(idx);
    R         = R(idx,:);
end

roi             = [];
roi.outDir      = [cfg.outDir '/ROI']; mkdir(roi.outDir);
roi.space       = 'cluster';

N = T;
for imap = 1:numel(varnameMaps)
    varname         = varnameMaps{imap};
    newvarname      = strcat(varname(3:end),'_',nodenames);%[varname(3:end) '_' nodename];
    roi.VOI         = {newvarname{:};R.cluster_voxels{:}}';
    roi.datafiles   = T.(varname); % List file paths to subjects' first level contrast images

    % Option to save Cluster mask in MNI space
    if imap == 1
        roi.doSaveMask  = doMasks;
    else
        roi.doSaveMask  = 0;
    end

    tblNode = ca_vba_util_extractROI(roi);
    N       = [N tblNode];
end
