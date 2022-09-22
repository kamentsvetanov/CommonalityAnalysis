## Commonality Analysis for Neuroimaging

Commonality analysis for neuroimaging. Code supporting Wu et al, "Cerebral blood flow predicts multiple demand network activity and fluid intelligence across the adult lifespan".

### Background
This voxel-wise GLM-like approach uses MATLAB's fitlm to save nii image of coefficients and p-values for each variable and residuals. This could be useful in instances with voxel-specific covariates. For example, in [Tsvetanov et al 2020](https://doi.org/10.1111/psyp.13714), we estimated variance explained and residuals in RSFA maps (across subjects) after controlling for regionally-speciffic effects of ASL maps, T1w maps in addition to other systemic effects.

We extended this voxel-wise approach to commonality analysis in [Wu et al 2022](https://www-sciencedirect-com.ezp.lib.cam.ac.uk/science/article/pii/S0197458022002044#sec0022).


![image](./Figures/Figure_1.png)
