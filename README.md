# Including Near-Infrared Spectroscopy Covariates via Sparse Bayesian Infinite Factor Models

This repository contains the R scripts used for the short paper *"Including Near-Infrared Spectroscopy Covariates via Sparse Bayesian Infinite Factor Models"*, co-authored with Alessandro Lanteri and Guillaume Kon Kam King, and submitted to SIS 2026.

The code is organized into six scripts covering preprocessing, benchmark models, the SBIFM implementation, factor interpretation, and the figures reported in the paper.

## Repository structure

### `01_data_preprocess.R`

Preprocesses the biscuit NIRS data. Observation 23 is removed as an outlier, following the literature (Brown et al., 2001). The script also includes the optional wavelength subsampling step, allowing the analysis to be run either on the full spectrum ($ p \approx 700 \u007f$) or on the reduced 256-wavelength setting. The response variable is selected by choosing one of the four constituents: fat, flour, sugar, or water.

### `02_pls.R`

Fits Partial Least Squares (PLS) models using the `plsR` package after preprocessing. The script assumes that the wavelength setting and response variable have already been selected.

### `03_glmnet_models.R`

Implements multi-output Ridge and LASSO regression. LASSO corresponds to `alpha = 1`, while Ridge corresponds to `alpha = 0`. The script includes the full preprocessing pipeline internally, so it can be run directly without first executing `01_data_preprocess.R`. It returns the four MSPE values, one for each response variable. To work on the full dataset with `p = 700`, remove the subsampling block; to switch from LASSO to Ridge, change `alpha` from 1 to 0.

### `04_sbifm_model.R`

Contains the Sparse Bayesian Infinite Factor Model (SBIFM) code of Bhattacharya and Dunson (2011), developed by Ewan Poworoznek and adapted with the Durante specification. Additional code has been added to compute regression coefficients `beta`, MSPE, the log-likelihood trace, and the evolution of the number of active factors for diagnostic purposes. This script should be run after choosing the response and the data setting in `01_data_preprocess.R`.

### `05_mostPredFactor.R`

Computes the most predictive factor using Ewan Poworoznek’s `infinitefactor` package and the rotational-alignment procedure from *Efficiently Resolving Rotational Ambiguity in Bayesian Matrix Sampling with Matching* (Poworoznek et al.). The script also produces the final plot used to identify the index of the most predictive factor.

### `06_plots.R`

Produces the figures included in the paper: the MSPE dotplot comparing the two dimensional settings across response variables, the beta plots from SBIFM and Ridge, and the most predictive factor visualization.

## Notes

* The data are publicly available in the original R packages and are not redistributed here.
* The repository is kept intentionally minimal and focused on reproducibility of the paper results.
* SBIFM is computationally intensive; the code is provided for full reproducibility, while benchmark models can be run in a much shorter time.

## Reference

If you use this code, please cite the corresponding short paper and the original methodological papers cited within the scripts.
