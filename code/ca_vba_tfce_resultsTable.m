function [cfg] = ca_vba_tfce_resultsTable(cfg,prefix)
% Post-process TFCE result maps to generate results table, similar to SPM
% Results tables
% 
% Usage:
% [cfg] = ca_vba_tfce_resultsTable(cfg,prefix)
% 
% 
% Inputs:
%   cfg   - cfg variable with varios fields define in ca_vba_glm_fitlm.m and 
%           ca_vba_tfce_threshold.m. 
%   prefix- string defining the naming pattern used to select statistical
%           maps of clusters surviving TFCE, i.e. outputs from ca_vba_tfce_threshold.m
%           These are likely to be 'tfce***', where *** was the
%           defined threshold level. Note that 'prefix' is further defined
%           by adding type of statistics used for TFCE below.
%
% Outputs:
%   cfg   - same as input inluding additional subfiedl statTable.
%   
%   results - Results Table, separate for each contrast image including 
%           t-scores, pvalues and coordinates for clusters and peaks
%           surviving TFCE. It also includes number and index of voxels
%           corresponding to a given cluster.
%
%
% Author : Kamen Tsvetanov, Ph.D., Neurocognitive Ageing
% Affil. : Department of Clinical Neurosciences, University of Cambridge
% Email  : kamen.tsvetanov@gmail.com  
% Website: http://www.kamentsvetanov.com
% Date   : 16 April 2022; Last revision: 
%__________________________________________________________________________
% Copyright (C) 2022, Kamen Tsvetanov
%
% Revision Log:
%
%
% ----------------------- BEGIN CODE ------------------------

voxthresh   = .1; % generally T-score, e.g. T=1.97 is equivalent to p-value =0.05
clustersize =  9; % cluster size in number of voxels
npeaks      = 10; % maximum number of peaks 


%-Extract some information 
%-------------------------
typeStats = cfg.tfce.typeStats;

try prefix = cfg.prefix; catch prefix = ['tfce' num2str(cfg.tfce.th*100)]; end % Default takes the threshold level used for TFCE

%-Select files based on type of statistics used for tfce
%-------------------------------------------------------
prefix  = [prefix '_' typeStats '_'];

fn   = rdir(fullfile(cfg.outDir,[prefix '*']));fn = {fn.name}';

Vaal = spm_vol('/imaging/camcan/sandbox/kt03/projects/external/mat/peak_nii/aal_MNI_V4.img'); 
Yaal = spm_read_vols(Vaal);
aalInv = inv(Vaal(1).mat);
Naal = load('/imaging/camcan/sandbox/kt03/projects/external/mat/peak_nii/aal_MNI_V4_List.mat');
Naal = Naal.ROI;
D    = 8; % size of ROI to draw around peak 



%-Predefine template Results Table 
%----------------------------------
sz = [1,6];
vartypes = ["string","double","double","double","double","double"];
varnames = ["cluster_name","cluster_nvox","cluster_tval","cluster_pval","peak_tval","peak_pval"];
Ttemp    = table('Size',sz,'VariableTypes',vartypes,'VariableNames',varnames);
Results  = [];

for icoef = 1:numel(fn)
    fname = fn{icoef};
    [d,f] = fileparts(fname);
    V     = spm_vol(fname);
    [Y, volXYZmm] = spm_read_vols(V);
    
    
    %-load P-value image
    fname_pval = regexprep(fname,typeStats,'pval'); % Select pval maps

    

    [resultTbl] = ca_vba_get_cluster_maxima(fname,fname_pval, voxthresh, clustersize, npeaks);
    
    %-Get inverse transform:
    %-----------------------
    Vinv  = inv(V(1).mat);
    
    %-Transform from Volume voxel space to AAL voxel space:
    %------------------------------------------------------
    aalXYZvox  = aalInv(1:3,1:3)*volXYZmm + repmat(aalInv(1:3,4),1,size(volXYZmm,2));
    
    %-Transform from voxel to mm space in AAL {mm}
    %--------------------------------------------------------------------------
    aalXYZmm      = Vaal(1).mat(1:3,:)*[aalXYZvox; ones(1, size(aalXYZvox,2))];
    Q          = ones(1,size(     aalXYZmm,2));%ones(1,size(xSPM.XYZmm,2));
    O          = ones(1,size(     aalXYZmm,2));
    
   
    P = spm_read_vols(spm_vol(regexprep(fname,typeStats,'pval'))); % Select pval maps
    CC = bwconncomp(Y,cfg.tfce.C);
    
   
    Tclust = Ttemp;
    for iC = 1:CC.NumObjects
        idx             = CC.PixelIdxList{iC};
        tval            = Y(idx);
        pval            = P(idx);
       
        
%          % get maxima
%         [N,Z,M,A,~] = spm_max(abs(Y(idx)),volXYZvox(:,idx));
%         [N,Z,M,A,~] = spm_max(Y1(:),volXYZvox(:,ind));
%         
        
        Tclust.peak_pval(iC)    = min(pval);
        Tclust.peak_tval(iC)    = max(tval); 
        [~,I]                   = max(tval); 
        Tclust.peak_xyzMNI(iC,:)= volXYZmm(:,idx(I))';       % Peak MNI coordinates
        peak_xyzAal             = ceil(aalXYZvox(:,idx(I)));  % Peak Vox coordinates in AAL resolution
        peak_xyzAalmm           = ceil(aalXYZmm(:,idx(I)));  % Peak Vox coordinates in AAL resolution
        idxAal                  = Yaal(peak_xyzAal(1),peak_xyzAal(2),peak_xyzAal(3));
        
        idx_neighbours          = find(sum((     aalXYZmm - peak_xyzAalmm*Q).^2) <= D^2);
        
        %-get AAL membership for neighbouring voxels to the Peak
        %----------------------------------------------------------
        idxAal = [];
        for ivox = 1:numel(idx_neighbours)
            temp_xyzvox     = ceil(aalXYZvox(:,idx_neighbours(ivox)));
            idxAal(ivox) = Yaal(temp_xyzvox(1),temp_xyzvox(2),temp_xyzvox(3));
        end
        idxAal = mode(nonzeros(idxAal)); % Pick the most represented region
        if isnan(idxAal)
            Tclust.peak_name(iC)  = 'unknown';
        else
            Tclust.peak_name(iC)  = string(Naal(idxAal).Nom_C);
        end
        Tclust.peak_voxels{iC}= idx_neighbours;
        
        %-get AAL membership for neighbouring voxels to the Peak
        %-------------------------------------------------------
        idxAal = [];
        for ivox = 1:numel(idx)
            temp_xyzvox     = ceil(aalXYZvox(:,idx(ivox)));
            idxAal(ivox) = Yaal(temp_xyzvox(1),temp_xyzvox(2),temp_xyzvox(3));
        end
        idxAal = mode(nonzeros(idxAal)); % Pick the most represented region
         if isnan(idxAal)
            Tclust.cluster_name(iC)  = 'unknown';
        else
            Tclust.cluster_name(iC)  = string(Naal(idxAal).Nom_C);
        end
%         Tclust.cluster_name(iC)  = string(Naal(idxAal).Nom_C);
        
        Tclust.cluster_nvox(iC)  = numel(idx);        
        Tclust.cluster_tval(iC)  = mean(tval);
        Tclust.cluster_pval(iC)  = mean(pval);
        Tclust.cluster_voxels{iC}= idx;
    end
    contrastname = (regexprep(f,prefix,''));
    Results.(contrastname) = Tclust;
    ResTbl.(contrastname) = resultTbl;
    resultTbl.contrast_name = strings(size(resultTbl,1),1);
    resultTbl.contrast_name(1) = contrastname;
    
    if icoef==1
        ResTblConcat = [resultTbl];
    else
        ResTblConcat = [ResTblConcat;resultTbl];
    end
end
ResTblConcat = movevars(ResTblConcat,'contrast_name','before','cluster_level[k]');
cfg.tfce.results    = Results;
cfg.tfce.table      = ResTbl;
cfg.tfce.tableConcat= ResTblConcat;

