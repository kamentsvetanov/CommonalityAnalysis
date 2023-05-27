function [VOIout] = kat_extract_VOI(cfg)
% Function to extract signals for given coordinates or a mask across list
% of images

% Author : Kamen Tsvetanov, Ph.D., Neurocognitive Ageing
% Affil. : Department of Clinical Neurosciences, University of Cambridge
% Email  : kamen.tsvetanov@gmail.com  
% Website: http://www.kamentsvetanov.com
% Date   : 27 May 2023; Last revision: 
%__________________________________________________________________________
% Copyright (C) 2023, Kamen Tsvetanov

% ------------
% Upack cfg
% ------------
try doSaveMask  = cfg.doSaveMask;   catch doSaveMask = 0;   end
try space       = cfg.space;        catch space = 'sphere'; end % 'sphere' | 'box' | 'mask'
try voisize     = cfg.size;         catch voisize = 6;      end % Size of VOI in mm space


Datafiles   = cfg.datafiles; % List of images
VOI         = cfg.VOI;%[30.79 38.15 29.85]';% VOI coordinates in mm space


numNode     = size(VOI,1);
nodeinfo    = VOI(:,2);
nodename    = VOI(:,1); % Name of the nodes


% --------------------------------------------------------------------------
% read and vectorize all subjects images
% --------------------------------------------------------------------------
V = spm_vol(Datafiles);
V = [V{:}];
Y4d = spm_read_vols(V);
Y = permute(Y4d,[4 1 2 3]);
Y = Y(:,:);


% --------------------------------------------
% Get XYZz space by loading on subject's image
% --------------------------------------------
[Z,XYZz] = spm_read_vols(V(1));

% Get inverse transform:
Yinv  = inv(V(1).mat);

% Transform ROI XYZ (in voxel space):
XYZ   = Yinv(1:3,1:3)*XYZz + repmat(Yinv(1:3,4),1,size(XYZz,2));

%-Voxels in entire search volume {mm}
%--------------------------------------------------------------------------
XYZmm      = V(1).mat(1:3,:)*[XYZ; ones(1, size(XYZ,2))];
Q          = ones(1,size(     XYZmm,2));%ones(1,size(xSPM.XYZmm,2));
O          = ones(1,size(     XYZmm,2));
%     FWHM       = xSPM.FWHM;

VOIout = nan(size(Datafiles,1),size(VOI,1));
for iVoi = 1:numNode
    
   
    switch space
        case 'sphere' %-Sphere based on a signle voxel
        %----------------------------------------------------------------------
            xyzmm = nodeinfo{iVoi}';
            if ~isfield(cfg,'size')
                D  = spm_input('radius of VOI {mm}',-2);
            else
                D  = voisize;
            end
            str    = sprintf('%0.1fmm sphere',D);
            if D == 1
                j      = find(sum((     XYZmm - xyzmm*Q).^2) <= 2);
                k      = find(sum((     XYZmm - xyzmm*O).^2) <= 2);
            else
                j      = find(sum((     XYZmm - xyzmm*Q).^2) <= D^2);
                k      = find(sum((     XYZmm - xyzmm*O).^2) <= D^2);
            end
        %         D      = D./xSPM.VOX;
        
            % Define name
            % -----------
            strCoord = [];
            for i =1:3
                if sign(xyzmm(i))>0
                    strCoord = [strCoord '+'];
                else
                    strCoord = [strCoord '-'];
                end
                strCoord = [strCoord num2str(abs(xyzmm(i)))];
            end
            fname = fullfile(cfg.outDir,sprintf('mask_%s_xyz%s_sphere%d.nii',nodename{iVoi},strCoord,voisize));
        
        case 'cluster' %- Cluster defined by many voxels, e.g. irregular stats ROI
        % ---------------------------------------------------------------------
            xyzmm = nodeinfo{iVoi}';
            j = xyzmm;
            
            %-Define name
            % -----------
             fname = fullfile(cfg.outDir,sprintf('mask_%s_cluster%d.nii',nodename{iVoi},numel(xyzmm)));
            

        case 'box' %-Box
        %----------------------------------------------------------------------
            xyzmm = nodeinfo{iVoi}';
            if ~isfield(cfg,'voisize')
                D  = spm_input('box dimensions [k l m] {mm}',-2);
            else
                D  = cfg.voisize;
            end
            if length(D)~=3, D = ones(1,3)*D(1); end
            str    = sprintf('%0.1f x %0.1f x %0.1f mm box',D(1),D(2),D(3));
            j      = find(all(abs(xSPM.XYZmm - xyzmm*Q) <= D(:)*Q/2));
            k      = find(all(abs(     XYZmm - xyzmm*O) <= D(:)*O/2));
            D      = D./xSPM.VOX;
            
            % Define name
            % -----------
            strCoord = [];
            for i =1:3
                if sign(xyzmm(i))>0
                    strCoord = [strCoord '+'];
                else
                    strCoord = [strCoord '-'];
                end
                strCoord = [strCoord num2str(abs(xyzmm(i)))];
            end
            fname = fullfile(cfg.outDir,sprintf('mask_%s_xyz%s_box%d.nii',nodename{iVoi},strCoord,voisize));

        case 'mask' %-Mask Image - TO DO
        %----------------------------------------------------------------------
            
            if isempty(nodeinfo{iVoi})
                [VM,sts] = spm_select([1 Inf],'image','Image defining search volume');
                if ~sts, TabDat = []; xSVC = []; return; end
            else
                VM = nodeinfo{iVoi};
            end
            D      = spm_vol(VM);
%             if numel(D) > 1
%                 fprintf('Computing union of all masks.\n');
%                 spm_check_orientations(D);
%                 D2 = struct(...
%                     'fname',   ['virtual_SVC_mask' spm_file_ext],...
%                     'dim',     D(1).dim,...
%                     'dt',      [spm_type('uint8') spm_platform('bigend')],...
%                     'mat',     D(1).mat,...
%                     'n',       1,...
%                     'pinfo',   [1 0 0]',...
%                     'descrip', 'SVC mask');
%                 D2.dat     = false(D2.dim);
%                 for i=1:numel(D)
%                     D2.dat = D2.dat | spm_read_vols(D(i));
%                 end
%                 D2.dat     = uint8(D2.dat);
%                 D  = D2;
%             end
%             str    = spm_file(D.fname,'short30');
%             str    = regexprep(str, {'\\' '\^' '_' '{' '}'}, ...
%                 {'\\\\' '\\^' '\\_' '\\{' '\\}'}); % Escape TeX special characters
%             str    = sprintf('image mask: %s',str); 
            
%             VOX    = sqrt(sum(D.mat(1:3,1:3).^2));
%             FWHM   = FWHM.*(xSPM.VOX./VOX);
%             XYZ    = D.mat \ [xSPM.XYZmm; ones(1, size(xSPM.XYZmm, 2))];
%             j      = find(spm_sample_vol(D, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
%             XYZ    = D.mat \ [     XYZmm; ones(1, size(     XYZmm, 2))];
%             k      = find(spm_sample_vol(D, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
            
            Ym = spm_read_vols(D);
            j = logical(Ym(:));
            fname = fullfile(cfg.outDir,sprintf('mask_%s_predefinedMask.nii',nodename{iVoi}));
    end

    Ymask = zeros(V(1).dim);
    Ymask(j) = 1;

    % create a new mask(optional)
    if doSaveMask
       
        Vmask = V(1);
        Vmask.fname = fname;
        spm_write_vol(Vmask,Ymask);
    end
    
    VOIout(:,iVoi) = nanmean(Y(:,j),2);
end

VOIout = array2table(VOIout,'VariableNames',nodename);


