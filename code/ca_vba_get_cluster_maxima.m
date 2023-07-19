function [resultTbl] = ca_vba_get_cluster_maxima(fname,fname_pval, voxthresh, clustersize, npeaks)
%function [resTbl] = get_cluster_maxima(map, voxthresh, sizethresh, npeaks)
% Returns table with up to npeaks maxima for clusters retrieved from map
% and defined by voxthresh and sizethresh.
% 
% Inputs:
%   fname       - full filepath to nii image
%   voxthresh   - generally T-score
%   clustersize - cluster size in number of voxels
%   npeaks      - maximum number of peaks 
%
%   e.g.:
%   fname       = 'my_file.nii';
%   voxthresh   = 2.5;
%   sizethresh  = 10;
%   npeaks      = 3;
%
% Adapted from: https://www.researchgate.net/post/How_to_extract_peak_coordinates_from_a_given_cluster_with_MATLAB_as_in_the_SPM_results_viewer
%
% Author : Kamen Tsvetanov, Ph.D., Neurocognitive Ageing
% Affil. : Department of Clinical Neurosciences, University of Cambridge
% Email  : kamen.tsvetanov@gmail.com  
% Website: http://www.kamentsvetanov.com
%__________________________________________________________________________
% Copyright (C) Kamen Tsvetanov 2022

%-Load AAL atlas information from peak_nii toolbox
%--------------------------------------------------
p = mfilename('fullpath');
p = regexp(p,'^.*/CommonalityAnalysis/code/','match');
p = fullfile(p{1},'external/aal');
Vaal    = spm_vol(fullfile(p,'aal_MNI_V4.img'));
Yaal    = spm_read_vols(Vaal);
aalInv  = inv(Vaal(1).mat);
InfoAAL = load(fullfile(p,'aal_MNI_V4_List.mat'));
InfoAAL = InfoAAL.ROI;

%-load image 
%-----------
V        = spm_vol(fname);
[Y, XYZ] = spm_read_vols(V);

%-load image 
%-----------
Vpval = spm_vol(fname_pval);
[Ypval, XYZpval]= spm_read_vols(Vpval);

%-voxel space
%-------------
dim     = size(Y);
[R,C,P] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
RCP     = [R(:)';C(:)';P(:)'];

%-Transform from Volume voxel space to AAL voxel space:
%------------------------------------------------------
aalXYZvox  = ceil(aalInv(1:3,1:3)*XYZ + repmat(aalInv(1:3,4),1,size(XYZ,2)));
aalXYZvoxAll = aalXYZvox;

%-Transform from voxel to mm space in AAL {mm}
%--------------------------------------------------------------------------
aalXYZmm      = Vaal(1).mat(1:3,:)*[aalXYZvox; ones(1, size(aalXYZvox,2))];
aalXYZmmAll   = aalXYZmm;
Q          = ones(1,size(     aalXYZmm,2));%ones(1,size(xSPM.XYZmm,2));
O          = ones(1,size(     aalXYZmm,2));
D    = 8; % size of ROI to draw around peak               

%-thresholding
%-------------
ind     = abs(Y) > voxthresh;
Y       = Y(ind);
Ypval   = Ypval(ind);
RCP     = RCP(:, ind);
XYZ     = XYZ(:, ind);
aalXYZvox = aalXYZvox(:,ind);
aalXYZmm = aalXYZmm(:,ind);

%-get maxima for positive effects
%------------
[N,Z,M,A,~] = spm_max(Y(:), RCP);

%-get maxima for Negative effects
%------------
[Nn,Zn,Mn,An,~] = spm_max(-Y(:), RCP);

% Check that regions from positive and negative contrast have unique
% numbers
% idClustPos = unique(A);
% idClustNeg = unique(An);
% if sum(ismember(idClustPos,idClustNeg))
% %     idxRepeat = ismember(idClustPos,idClustNeg);
% %     idClustNeg = idClustPos(idxRepeat)
%     error('Negative and positive contrast have same ID cluster. This has to be avoided!');
% end

% Concatenate Positive and Negative effects
if ~isempty(N) & ~isempty(Nn)
    N = [N; Nn];
    Z = [Z; -Zn];
    M = [M Mn];
    A = [A; An+max(A)];%Add max number in A to An to avoid merging of clusters with identical IDs

elseif isempty(N) & ~isempty(Nn) % No positive effects, only negative effects
    N = Nn;
    Z = -Zn;
    M = Mn;
    A = An;%Add max number in A to An to avoid merging of clusters with identical IDs
end
    
%-get MNI coordinates for all clusters
%--------------------------------------
if ~isempty(M)
    [~, idxCoord] = ismember(M.', RCP.', 'rows');
    XYZpeak     = XYZ(1:3, idxCoord);
    aalXYZvoxPeak = aalXYZvox(:,idxCoord);
    aalXYZmmPeak= aalXYZmm(:,idxCoord);
    % Get associated pvalues
%     [~, idxCoordPval] = ismember(M.', RCP.', 'rows');
    Pval = Ypval(idxCoord); 
end

%-create results table
%----------------------
resTbl = {'cluster-level','voxel-level',[],'coordinates (mm)',[],[]; ...
'k','tVal','pVal','x','y','z'};

%-Predefine template Results Table 
%----------------------------------
sz = [1,7];
vartypes = ["string","double","double","double","double","double","double"];
varnames = ["cluster_name","cluster_level[k]","peak_tval","peak_pval","x_coord","y_coord","z_coord"];
resultTbl= table('Size',sz,'VariableTypes',vartypes,'VariableNames',varnames);
Results  = [];


curRow = 1;
for iClust = 1:max(A)
    % get associated voxels
    idxClust = A == iClust;
    
    % get associated coordinates
    xyzClust = XYZpeak(:,idxClust);
    xyzClustAALvox = aalXYZvoxPeak(:,idxClust);
    xyzClustAALmm = aalXYZmmPeak(:,idxClust);
    
    % sort by value
    B = Z(idxClust);
    [~,I]  = sort(abs(B),'descend');
    B      = B(I);  
    pClust = Pval(idxClust);
    pClust = pClust(I); % Reorder according to B
    
    % cluster extent
    n = N(find(idxClust,1));
    
    % cluster size threshold for printing
    if n > clustersize
%         resTbl{curRow,1} = n;
        resultTbl{curRow,"cluster_level[k]"} = n;
        % iterate over peaks
        for iPeak = 1:min(length(I), npeaks)
            % value
%             resTbl{curRow,2} = B(iPeak);
            resultTbl{curRow,"peak_tval"} = B(iPeak);
            resultTbl{curRow,"peak_pval"} = pClust(iPeak);
            % coordinates
            xyz = xyzClust(:, I(iPeak));
%             resTbl{curRow,3} = xyz(1);
%             resTbl{curRow,4} = xyz(2);
%             resTbl{curRow,5} = xyz(3);
            resultTbl{curRow,"x_coord"} = xyz(1);
            resultTbl{curRow,"y_coord"} = xyz(2);
            resultTbl{curRow,"z_coord"} = xyz(3);

            %-Get AAL label
            %---------------------------------------------------------------
            xyzAAL = xyzClustAALvox(:,I(iPeak));
            idxAal = Yaal(xyzAAL(1),xyzAAL(2),xyzAAL(3));
            xyzAALmm = xyzClustAALmm(:,I(iPeak));
            idx_neighbours          = find(sum((     aalXYZmmAll - xyzAALmm*Q).^2) <= D^2);
            
            %-get AAL membership for neighbouring voxels to the Peak
            %----------------------------------------------------------
            idxAal = [];
            for ivox = 1:numel(idx_neighbours)
                temp_xyzvox  = ceil(aalXYZvoxAll(:,idx_neighbours(ivox)));
                idxAal(ivox) = Yaal(temp_xyzvox(1),temp_xyzvox(2),temp_xyzvox(3));
            end
            idxAal = mode(nonzeros(idxAal)); % Pick the most represented region
            if isnan(idxAal)
                resultTbl{curRow,"cluster_name"}   = 'unknown';
            else
                resultTbl{curRow,"cluster_name"} = string(InfoAAL(idxAal).Nom_C);
            end
            
            curRow = curRow + 1;
        end
    end
end

% Nan all zero values for cluster level[k]
resultTbl{resultTbl.("cluster_level[k]")==0,"cluster_level[k]"}=nan;

resultTbl = movevars(resultTbl,'cluster_level[k]','before','cluster_name');
% display table
% disp(resTbl);
% save csv
% [p,n,~] = fileparts(map);
% writetable(table(resTbl), fullfile(p, [n '_peaks.csv']), 'WriteVariableNames', 0);
end