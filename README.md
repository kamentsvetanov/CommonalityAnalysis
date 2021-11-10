## Commonality Analysis for Neuroimaging

Commonality analysis for neuroimaging. Code supporting Wu et al, "Cerebral blood flow predicts multiple demand network activity and fluid intelligence across the lifespan".

### Background
This voxel-wise GLM-like approach uses fitlm to save nii image of coeff and p-values for each variable and residuals (unexplained effects by predictors). This could be useful in instances with voxel-specific covariates, e.g. in Tsvetanov et al 2020 Psyhophysiology (https://doi.org/10.1111/psyp.13714) we estimate variance explained and residuals in RSFA maps (across subjects) after controlling for the effects of ASL maps (voxel-specific), T1w maps and other effects.

We extended this voxel-wise approach to commonality analysis in Wu et al 2021 (preprint link to follow).


![image](./Figures/Figure_1.png)
