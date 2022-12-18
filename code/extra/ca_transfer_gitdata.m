

% Load S and T from Tsvtetanov et al 2020
load('/imaging/camcan/sandbox/kt03/archived/2020TsvetanovPsyP/data/settings_02_proc_data_R1_20200617.mat');

% set up new table with a subset of individuals
tbl        = [];
tbl.f_rsfa = T.f_rsfaAromaW1_norm_wcm2f11_10_128_s12;
tbl.Age    = T.Age;
tbl.Sex    = T.GenderNum;
% tbl.SubID  = T.SubCCIDc;
tbl        = struct2table(tbl);

% Remove subjects with missing rsfa data
idx = cellfun(@isempty,tbl.f_rsfa,'UniformOutput',false);
idx = [idx{:}]';
tbl(idx,:) = [];

% Take n young and n old individuals
n = 10;
ns = size(tbl,1);
tbl = sortrows(tbl,'Age');
tbl = tbl([1:n ns-n+1:ns],:);

% Copy rsfa maps to commonality repo
outdir = '/home/kt03/Projects/public-code/CommonalityAnalysis/data/rsfa';
tblOut = tbl;
for isub = 1:size(tbl,1)
    f_in = tbl.f_rsfa{isub};
    f_in = regexprep(f_in,'projects/rsfa_cv2','archived/2020TsvetanovPsyP');
    f_out = fullfile(outdir,sprintf('sub%04.f.nii',isub));
    copyfile(f_in,f_out);

    tblOut.f_rsfa{isub} = f_out;
    tblOut.SubID{isub}  = sprintf('sub%04.f.nii',isub);
end
T = tblOut;
save(fullfile(outdir,'subject_info.mat'),'T');