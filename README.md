# Robust statistical methods for M/EEG multidimensional data

- GLM including OLS, IRLS, WLS optimization for time (ERP) and frequency (power spectra) domains. GLM deals better accounts for within-subject variance.

- nonparametric statistics, including bootstrap and permutation approaches

- robust corrections for type 1 error (maximum likelihood, spatiotemporal clustering, threshold-free cluster enhancement (TFCE)

- 95% confidence intervals (CIs) and 95% Bayesian high-density intervals (HDIs)

- data visualization (mass-univariate, scalp topography, course plot)


## Comparing GLM/no-GLM results on ERP data (N = 78), p<0.001

# Raw ERP (no GLM)
<img width="50%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_RAW_corrected.png">
<img width="30%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_RAW_corrected_peak-channel">

Cluster 1: 228 to 948 ms. Peak effect: channel P3 at 436 ms (t = 14) 
Cluster 2: 132 to 216 ms. Peak effect: channel P8 at 156 ms (t = 5.7) 

# GLM with OLS optimization
<img width="50%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_OLS_corrected.png">
<img width="50%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_OLS_corrected_peak-channel.png">

Cluster 1: 228 to 948 ms. Peak effect: channel P4 at 460 ms (t = 13.4) 
Cluster 2: 132 to 224 ms. Peak effect: channel P8 at 156 ms (t = 5.9) 

# GLM with IRLS optimization

<img width="50%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_IRLS_corrected.png">
<img width="50%" src="https://github.com/amisepa/eeg_robust_statistics/blob/main/outputs/result_unpleasant-neutral_GLM_IRLS_corrected_peak-channel.png">

