#!/bin/bash

model='f_rsfa ~ Age + c_Sex'
rootDir=/imaging/camcan/sandbox/kt03/temp/
f_mask=/imaging/camcan/sandbox/kt03/projects/public-code/CommonalityAnalysis/data/rsfa/mask.nii
f_table=/imaging/camcan/sandbox/kt03/projects/public-code/CommonalityAnalysis/data/rsfa/subject_info.mat
numPerm=100
doCommonality=1


