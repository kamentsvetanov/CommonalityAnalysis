function [cfg] = kat_stats_perm_maxT(cfg)
% DESCRIPTION
% Permutaion based adjustment of p- and T-values based on the following
% statistic tests: one-sample, two-sample and regression
% 
% INPUTs:
% ---------
% cfg: a structure with the following fields
%   .nperm - number of permutations
%   .type  - type of permutation test: 
%           'one-sample' 
%           'two-sample-unpaired'
%           'two-sample-paired'
%           'regression'
%   .data  - array of data organised depending on type
%           'one-sample': n x m array, where n are observations and m is
%               defferent samples/tests
%           'two-sample-upaired': n x m, where n are observations and m = 2
%               for two unpaired/unmatched samples, e.g. young and older group of
%               same Ns
%           'two-sample-paired': same as upaired test, but here the
%               samples come from two independent measurements within a 
%               subject, e.g. conditions
%           'regression': n x m array, where n are obsevations and m are 
%               variables specified in 'model' (see below)
%   .model - this option applies to "regression" and two-sample-unpaired test. In the former, the model
%               specifies the variables. In the latter, the model specifies
%               specified indices for semparation of the two groups to be compared.
%
%           
%
%
%
%
%
% Adusted from Jannete Mumford
% Assume you have 15 ROIs and for each ROI you want to look at the
% hypothesis test for the overall mean (over 20 measures).  Although we could use a bonferroni
% test to correct for all 15 tests, we might be better off using a
% permuataion based test.
%
% Steps of the max-statistic based permutation test
% Step 1:  Calculate test statistics for all 15 tests without use of
% permutations
%
% Step 2:  Randomly select some observations and multiply their values by
% negative 1.  Note, if the null were true in this case, then it shouldn't
% matter if we change the sign of our values or not.
%
% Step 3: Calculate the 15 test statistics and store the value of the
% maximum of the 15 statistics
%
% Step 4: Repeat 2 & 3 many, many times (let's do 5000)
% Step 5:  Compare the original test statistics from step 1 to the cutoff
% corresponding to the 95th percentile of the max tstat distribution
%
% load in data.  Note that each column holds the 20 mesures for a single ROI
% roi_vals=load('/imaging/camcan/sandbox/kt03/toolboxes/permutations/roi_data.mat', '-ascii');

try nperm = cfg.nperm; catch nperm = 5000; end


dat      = full(cfg.data); 
numObs   = size(dat,1);
numVar   = size(dat,2);
tmax     = zeros(1,nperm);
testtype = cfg.type;



%=================================
%Steps 2 & 3 &4
switch testtype
    case 'one-sample' % One-sample T permutation (flip signs not re-label conditions)
        % Get original correlations
        meanX  = mean(dat, 1);
        stdX   = std(dat,0, 1);
        t_orig = sqrt(numObs).*meanX./stdX;
        dat0 = zeros(size(dat));
        [H,P,CI,stat]=ttest(dat,dat0,'tail','both');
        t_orig2 = stat.tstat;
        p_orig = 1-tcdf(abs(t_orig), numVar-1); 

        parfor ii=1:nperm
        %     ii
            %randomly select 1's and -1's
            random_flip     = randsample([1, -1],numObs,true);%([1, -1],nsub) | ([1 -1],nsub,true)
            random_flip_mat = random_flip'*ones(1,numVar);
            mat_flipped     = random_flip_mat.*dat;

            %calculate 15 t stats for permuted data
            mn_loop = mean(mat_flipped, 1);
            sd_loop = std(mat_flipped,0, 1);
            t_loop = sqrt(numObs).*mn_loop./sd_loop;


            %Cacluate the max t stat for this iteration
            tmax(ii) = max(abs(t_loop));
        end
    case 'two-sample-unpaired' % re-label randomly across rows
        % In this case we use Wilcoxon Rank-Sum test. Note that number of
        % observations in both samples could be different
        idx         = logical(cfg.model);
        dat1        = dat(idx,:); % Select dataset 1
        dat2        = dat(~idx,:); % Select dataset 2
        [h,p,ci,st] = ttest2(dat1,dat2);%[p,h,stats] = ranksum(dat1,dat2);
        p_orig      = p;
        t_orig      = st.tstat;
      
        parfor ii = 1:nperm
            %randomly select 1's and -1's
            newdat      = Shuffle(Shuffle(dat));
            dat1        = newdat(idx,:); % Select dataset 1
            dat2        = newdat(~idx,:); % Select dataset 2
            [h,p,ci,st] = ttest2(dat1,dat2);%[p,h,stats] = ranksum(dat1,dat2);
            tmax(ii)    = max(abs(st.tstat));
%             random_flip_mat = 
        end
    case 'two-sample-paired' % re-label within a row
        % in this case we use Wilcoxo Sign-Rank Test
        numSample = size(dat,1);
        dat1 = dat(:,1);
        dat2 = dat(:,2);
        [p,h,stats] = signrank(dat1,dat2);
        t_orig      = abs(tinv(p/2,numSample)) * sign(stats.zval);
        p_orig      = p;
        parfor ii = 1:nperm
            newDat = Shuffle(Shuffle(dat,2));
            dat1 = newDat(:,1);
            dat2 = newDat(:,2)
            [p,h,stats] = signrank(dat1,dat2);
            tmax(ii)    = max(abs(tinv(p/2,numSample)));
        end
        
        numVar = 1;
    case 'regression' % regression
        if isfield(cfg,'model')
            M_orig = cfg.model;
        else
            M_orig = dat(:,1);
            dat      = dat(:,2:end);
            numVar   = size(dat,2);
        end
        contr       = [0 1 zeros(1,size(M_orig,2)-1)]'; % Set the contrast
        I      = ones(size(M_orig,1),1);
        M      = [I  M_orig];
        [t_orig, p_orig] = kat_stats_spm_ancova(M,[],dat,contr);%
       
        parfor ii=1:nperm
            M = [I Shuffle(Shuffle(M_orig(:,1))) M_orig(:,2:end)];% Randomize Y or order of subjects  M_orig(randperm(size(M_orig,1)),1)
%             t_loop1 = corr(M,dat)
            [t_loop] = kat_stats_spm_ancova(M,[],dat,contr);
          
            tmax(ii) = max(abs(t_loop)); %Cacluate the max t stat for this iteration
           
        end
end


% Step 5:  Calculate the percentiles for each statistic based on
% permutation-based null dist
p_perm = zeros(numVar,1);
for ii = 1:numVar
    p_perm(ii)=sum(tmax>=abs(t_orig(ii)))/nperm;
end

% stat.p_orig = p_orig;
cfg.p_perm = p_perm';
cfg.tstat  = t_orig;
cfg.p_orig = p_orig;
try cfg.tstat2  = t_orig2; end
% 
% [p_orig', p_perm]
% 
% sig_pcer = p_orig<0.05;
% sig_bon  = p_orig<0.05/nobs;
% sig_perm = p_perm<0.05;
% 
% sig_vals = [sig_pcer; sig_bon; sig_perm']'
% sum(sig_vals)

