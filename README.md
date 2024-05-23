# Robust statistical methods for M/EEG multidimensional data

- GLM including OLS, IRLS, WLS optimization for time (ERP) and frequency (power spectra) domains. GLM deals better accounts for within-subject variance.

- nonparametric statistics, including bootstrap and permutation approaches

- robust corrections for type 1 error (maximum likelihood, spatiotemporal clustering, threshold-free cluster enhancement (TFCE)

- 95% confidence intervals (CIs) and 95% Bayesian high-density intervals (HDIs)

- data visualization (mass-univariate, scalp topography, course plot)


## Comparing GLM/no-GLM results on ERP data (N = 78), p<0.001

Vectorized approach to significantly increase computation cost. Vectorized code performs matrix operations directly, which are generally more efficient and less error-prone than iterating through individual elements. Matrix multiplication and other linear algebra operations in MATLAB are highly optimized. The core computations, such as the design matrix multiplication (XTWX, XTWY), fitting the GLM (pinv(XTWX) * XTWY), and computing residuals (data_reshaped - Yhat), are the same in both approaches.


Optimization options: OLS, IRLS, WLS

WLS options
- Pernet's principal components projection (PCP): https://www.humanbrainmapping.org/files/Aperture%20Neuro/Accepted%20Works%20PDF/7_51_Cyril_ElectroEncephaloGraphy_robust_statistical_linear.pdf
- Hubert M-estimator: Combines squared error for small residuals and absolute error for large residuals, offering a balance between efficiency and robustness.
- Tukey's Biweight: Reduces the influence of outliers more aggressively by setting large residuals to zero weight.

# Raw ERP (no GLM)
<img width="40%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_RAW_corrected.png"> <img width="30%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_RAW_corrected_peak-channel">

Cluster 1: 228 to 948 ms. Peak effect: channel P3 at 436 ms (t = 14) 

Cluster 2: 132 to 216 ms. Peak effect: channel P8 at 156 ms (t = 5.7) 

### GLM with OLS optimization
<img width="40%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_OLS_corrected.png"> <img width="30%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_OLS_corrected_peak-channel.png">

Cluster 1: 228 to 948 ms. Peak effect: channel P4 at 460 ms (t = 13.4) 

Cluster 2: 132 to 224 ms. Peak effect: channel P8 at 156 ms (t = 5.9) 

### GLM with IRLS optimization

<img width="40%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_IRLS_corrected.png"> <img width="30%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_IRLS_corrected_peak-channel.png">

Cluster 1: 228 to 948 ms. Peak effect: channel P3 at 448 ms (t = 13.7) 

Cluster 2: 132 to 224 ms. Peak effect: channel O1 at 176 ms (t = 6.2) 

### WLS - Tukey biweight function

<img width="40%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_WLS-Tukey_corrected.png"> <img width="30%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_WLS-Tukey_corrected_peak-channel.png">

Cluster 1: 232 to 948 ms. Peak effect: channel P4 at 460 ms (t = 13.9) 

Cluster 2: 132 to 220 ms. Peak effect: channel O1 at 172 ms (t = 5.7) 

### WLS - Hubert M-estimator

<img width="40%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_WLS-Hubert_corrected.png"> <img width="30%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_WLS-Hubert_corrected_peak-channel.png">

Cluster 1: 228 to 948 ms. Peak effect: channel P3 at 448 ms (t = 13.7) 

Cluster 2: 132 to 224 ms. Peak effect: channel O1 at 176 ms (t = 6) 

### WLS - Pernet's PCP

<img width="40%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_WLS-PCP_corrected.png"> <img width="30%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_WLS-PCP_corrected_peak-channel.png">


Cluster 1: 224 to 948 ms. Peak effect: channel P3 at 396 ms (t = 14) 
Cluster 2: 152 to 212 ms. Peak effect: channel PO3 at 176 ms (t = 6) 
