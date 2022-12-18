function ca_vba_util_onesample(cfg)
% Untily function in Commonality toolbox to compute main/average effect
% using one-sample test to estimation
% % One-sample T permutation (flip signs not re-label conditions)
% N.B. incomplete!!!!

Y(isnan(Y)) = 0;
H = 2;
E = 0.5;
C = 26;
dh = 0.2;
th = 1.97;
% Get original correlations
numSub = size(Y,1);
numVar = size(Y,2);
meanX  = nanmean(Y, 1);
stdX   = nanstd(Y,0, 1);
t_orig = sqrt(numSub).*meanX./stdX;
t_orig(abs(t_orig)==inf) = 0; % Reset Inf to 0;
t_orig(isnan(t_orig)) = 0; % Reset NanNs to 0;
dat0 = zeros(size(Y));

Ypos = t_orig;Ypos(t_orig<0)=0;
Yneg = t_orig;Yneg(t_orig>0)=0;
[trealPos]  = matlab_tfce_transform(Ypos,H,E,C,dh,th);
[trealNeg]  = matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);
treal       = trealPos-trealNeg;

% cycle through permutations
nvox = numel(t_orig);
exceedancesPos = zeros(size(t_orig));
exceedancesNeg = zeros(size(t_orig));

%randomly select 1's and -1's
randOrder(randOrder<numSub/2) = -1;
randOrder(randOrder>=numSub/2)= 1;

parfor iperm=1:numPerm
%             iperm
    %randomly select 1's and -1's
    random_flip     = randOrder(:,iperm);% randsample([1, -1],size(Y,1),true);%([1, -1],nsub) | ([1 -1],nsub,true)
    mat_flipped     = repmat(random_flip,1,nvox).*Y;

    %calculate 15 t stats for permuted data
    mn_loop = nanmean(mat_flipped, 1);
    sd_loop = nanstd(mat_flipped,0, 1);
    t_loop = sqrt(numSub).*mn_loop./sd_loop;
    t_loop(abs(t_loop)==inf) = 0; % Reset Inf to 0;
    t_loop(isnan(t_loop)) = 0; % Reset NanNs to 0;

    Ypos = t_loop;Ypos(t_loop<0)=0;
    Yneg = t_loop;Yneg(t_loop>0)=0;
    [tnullPos] = matlab_tfce_transform(Ypos,H,E,C,dh,th);
    [tnullNeg] = matlab_tfce_transform(abs(Yneg),H,E,C,dh,th);
    % compare maxima to t-values and increment as appropriate
    curexceeds      = max(tnullPos(:)) >= trealPos;
    exceedancesPos  = exceedancesPos + curexceeds;
    curexceeds      = min(-tnullNeg(:)) <= -trealNeg;
    exceedancesNeg  = exceedancesNeg + curexceeds;

end
correctedPos = exceedancesPos./(numPerm);
correctedNeg = exceedancesNeg./(numPerm);

corrected4d = []; corrected = [];
corrected4d(:,:,:,1) = correctedPos;
corrected4d(:,:,:,2) = correctedNeg;
corrected = min(corrected4d,[],4);
%     corrected = min([correctedPos; correctedNeg],4);
treal(corrected>0.05) = 0;
% Write output image

Vtemp           = V(1);
Vtemp.pinfo(1)  =1;
Vtemp.fname     = fullfile(outDir,sprintf('%s_tfce%d.nii',nameOutput,th*100));
tmap            = zeros(Vtemp.dim);
tmap(idxMask)   = treal;  
spm_write_vol(Vtemp,tmap);

%          % Save t-stats
%         tdir            = fullfile(outDir,['tval_(Intercept)',]); 
%         mkdir(tdir);
%         Vtemp           = V(1);
%         Vtemp.pinfo(1)  =1;
%         Vtemp.fname     = fullfile(tdir,sprintf('results_null_%.5d.nii',1));
%         tmap            = zeros(Vtemp.dim);
%         tmap(idxMask)   = t_orig;  
%         spm_write_vol(Vtemp,tmap);