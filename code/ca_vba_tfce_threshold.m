function [cfg] = ca_vba_tfce_threshold(cfg)
% Apply TFCE thresholding using matlab_tfce_transform on output from
% voxel-based analysis with permutations
% https://github.com/markallenthornton/MatlabTFCE

% --------------------------------------------------------------
% Required packages (https://github.com/kamentsvetanov/external)
% --------------------------------------------------------------
% addpath(genpath('/imaging/camcan/sandbox/kt03/projects/private-code/kamen'));
% kat_import('bwt');
% kat_import('load_nii');
% kat_import('tfce');
% kat_import('spm12');
% kat_import('matlabcentral');

% ------------------
% TCFE options
% ------------------
try tfce        = cfg.tfce;        catch           end
try H           = tfce.H;          catch H  = 2;   end % height exponent, default = 2
try E           = tfce.E;          catch E  = 0.5; end % extent exponent, default = 0.5
try C           = tfce.C;          catch C  = 26;  end % connectivity, default = 26 (6 = surface, 18 = edge, 26 = corner)
try dh          = tfce.dh;         catch dh = 0.1; end %  step size for cluster formation, default = .1
try th          = tfce.th;         catch th = 1.97;end % Threshold level to apply correction
try clustSize   = tfce.clustersize;catch clustSize=5;end % Minimum number of voxels forming a cluster
try Ns          = tfce.Ns;         catch Ns = size(cfg.tbl,1); fprintf('Selecting sample size based on subjects in tbl.\n'); end % Number of subject
try Np          = tfce.Np;         catch Np = size(cfg.tbl,2)-1; fprintf('Selecting number predictors based on variables in tbl.\n'); end % Number of predictors
try path2data   = tfce.path2data;  catch error('Specify path to permutations.\n'); end
try typeStats   = tfce.typeStats;  catch error('Specify folder prefix, based on statistic type. [tval | rval | bval].\n'); end % This is the pattern used to select directories, e.g. bVal, tVal glmTCA

try showIntercept = tfce.showIntercept; catch showIntercept = 0; end % Whether to show results for (Intercept). Default = 0, no 

nameCoefficients    = dir(fullfile(path2data,[typeStats '*']));
nameCoefficients    = nameCoefficients([nameCoefficients.isdir]);
nameCoefficients    = {nameCoefficients.name}';
numCoeff            = numel(nameCoefficients);
df                  = Ns - Np; % Degrees of freedom

% -------------------------------------------------
% Exclude the Intercept, as it's rarely of interest
% -------------------------------------------------
if ~showIntercept == 1
    idx_intercept = contains(nameCoefficients,'(Intercept)');
    nameCoefficients(idx_intercept) = [];
end
  

% ------------------------------
% Loop through each coefficient
% ------------------------------
for iCoeff = 1:numel(nameCoefficients)
    
    namecoeff   = nameCoefficients{iCoeff};
    fname       = ca_rdir(fullfile(path2data,namecoeff,'results_null_*.nii')); fname = {fname.name}';
    nperm       = size(fname,1);

    % Get t-stats for observed data (1st image)
    V = spm_vol(fname{1});
    Yreal = spm_read_vols(V);
    Ypos = Yreal;Ypos(Yreal<0)=0;
    Yneg = Yreal;Yneg(Yreal>0)=0;
    [trealPos]  = ca_matlab_tfce_transform(Ypos,H,E,C,dh,th);
    [trealNeg]  = ca_matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);
    treal       = trealPos-trealNeg;    
    
    % cycle through permutations
    nvox = numel(Yreal);
    exceedancesPos = zeros(size(Yreal));
    exceedancesNeg = zeros(size(Yreal));
%     parfor(p = 1:nperm,parworkers)
    parfor p = 1:nperm
        V           = spm_vol(fname{p});
        Y           = spm_read_vols(V);
        Ypos = Y; Ypos(Y<0)=0;
        Yneg = Y; Yneg(Y>0)=0;
        [tnullPos]  = ca_matlab_tfce_transform(Ypos,H,E,C,dh,th);
        [tnullNeg]  = ca_matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);
        % compare maxima to t-values and increment as appropriate
        curexceeds      = max(tnullPos(:)) >= trealPos;
        exceedancesPos  = exceedancesPos + curexceeds;
        curexceeds      = max(tnullNeg(:)) >= trealNeg;
        exceedancesNeg  = exceedancesNeg + curexceeds;
        maxNeg(p) = max(tnullNeg(:));
    end

    correctedPos = exceedancesPos./(nperm);
    correctedNeg = exceedancesNeg./(nperm);

    corrected4d = []; corrected = [];
    corrected4d(:,:,:,1)    = correctedPos;
    corrected4d(:,:,:,2)    = correctedNeg;
    corrected               = min(corrected4d,[],4);
    Yreal(corrected>0.05)   = 0;

    % ---------------------------------------------------------
    % Remove clusters smaller than certain number of voxels
    % ---------------------------------------------------------
    cc = bwconncomp(abs(Yreal)>=th,C);
    for iclust = 1:cc.NumObjects
        if numel(cc.PixelIdxList{iclust})<clustSize
            Yreal(cc.PixelIdxList{iclust}) = 0;
        end
    end
    
    % ----------------------------
    % Write output image
    % ----------------------------
    fout = fullfile(path2data,sprintf('tfce%3.f_%s.nii',th*100,namecoeff));
    V.fname = fout;
    V.private = [];
    spm_write_vol(V,Yreal);
    
    fout = regexprep(fout,typeStats,'pval');
    V.fname = fout;
    V.private = [];
    spm_write_vol(V,corrected);
end

%% Prepare output
tfce.H  = H;
tfce.E  = E;
tfce.C  = C;
tfce.dh = dh;
tfce.th = th;
tfce.numPerm = nperm;
cfg.tfce= tfce;

save(fullfile(cfg.outDir,'analysis_cfg.mat'),'cfg','tfce');

% function [tPos tNeg] = tfce_transform(fname,H,E,C,dh,th)
% 
% V = spm_vol(fname);
% Y = spm_read_vols(V);
% 
% Ypos = Y;Ypos(Y<0)=0;
% Yneg = Y;Yneg(Y>0)=0;
% [tPos] = matlab_tfce_transform(Ypos,H,E,C,dh,th);
% [tNeg] = matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);