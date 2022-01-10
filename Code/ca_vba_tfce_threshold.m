function kat_vba_tfce_threshold(cfg)
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
try H           = tfce.H;          catch H = 2;    end % height exponent, default = 2
try E           = tfce.E;          catch E = 0.5;  end % extent exponent, default = 0.5
try C           = tfce.C;          catch C = 26;   end % connectivity, default = 26 (6 = surface, 18 = edge, 26 = corner)
try dh          = tfce.dh;         catch dh = 0.2; end %  step size for cluster formation, default = .1
try th          = tfce.th;         catch th = 1.97;end % Threshold level to apply correction
try Ns          = tfce.Ns;         catch error('Sample size is needed.'); end % Number of subject
try Np          = tfce.Np;         catch error('Number of predictors is needed.'); end % Number of predictors
try path2data   = tfce.path2data;  catch error('Specify path to permutations.'); end
try typeStats   = tfce.typeStats;  catch error('Specify folder prefix, based on statistic type. [tVal | rVal | bVal]'); end % This is the pattern used to select directories, e.g. bVal, tVal glmTCA

nameCoefficients    = dir(fullfile(path2data,[typeStats '*']));
nameCoefficients    = nameCoefficients([nameCoefficients.isdir]);
nameCoefficients    = {nameCoefficients.name}';
numCoeff            = numel(nameCoefficients);
df                  = Ns - Np; % Degrees of freedom

% ------------------------------
% Loop through each coefficient
% ------------------------------
for iCoeff = 1:numel(nameCoefficients)
    
    namecoeff   = nameCoefficients{iCoeff};
    fname       = rdir(fullfile(path2data,namecoeff,'results_null_*.nii')); fname = {fname.name}';
    nperm       = size(fname,1);

    % Get t-stats for observed data (1st image)
    V = spm_vol(fname{1});
    Y = spm_read_vols(V);
    Ypos = Y;Ypos(Y<0)=0;
    Yneg = Y;Yneg(Y>0)=0;
    [trealPos]  = matlab_tfce_transform(Ypos,H,E,C,dh,th);
    [trealNeg]  = matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);
    treal       = trealPos-trealNeg;
    
    % cycle through permutations
    nvox = numel(Y);
    exceedancesPos = zeros(size(Y));
    exceedancesNeg = zeros(size(Y));
%     parfor(p = 1:nperm,parworkers)
    parfor p = 1:nperm
        V = spm_vol(fname{p});
        Y = spm_read_vols(V);
        Ypos = Y;Ypos(Y<0)=0;
        Yneg = Y;Yneg(Y>0)=0;
        [tnullPos] = matlab_tfce_transform(Ypos,H,E,C,dh,th);
        [tnullNeg] = matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);
        % compare maxima to t-values and increment as appropriate
        curexceeds      = max(tnullPos(:)) >= trealPos;
        exceedancesPos  = exceedancesPos + curexceeds;
        curexceeds      = min(-tnullNeg(:)) <= -trealNeg;
        exceedancesNeg  = exceedancesNeg + curexceeds;
    end

    correctedPos = exceedancesPos./(nperm);
    correctedNeg = exceedancesNeg./(nperm);

    corrected4d = []; corrected = [];
    corrected4d(:,:,:,1) = correctedPos;
    corrected4d(:,:,:,2) = correctedNeg;
    corrected = min(corrected4d,[],4);
%     corrected = min([correctedPos; correctedNeg],4);
    treal(corrected>0.05) = 0;

    % Write output image
    fout = fullfile(cfg.outDir,sprintf('%s_tfce%d.nii',namecoeff,th*100));
    V.fname = fout;
    spm_write_vol(V,treal);
end


% function [tPos tNeg] = tfce_transform(fname,H,E,C,dh,th)
% 
% V = spm_vol(fname);
% Y = spm_read_vols(V);
% 
% Ypos = Y;Ypos(Y<0)=0;
% Yneg = Y;Yneg(Y>0)=0;
% [tPos] = matlab_tfce_transform(Ypos,H,E,C,dh,th);
% [tNeg] = matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);