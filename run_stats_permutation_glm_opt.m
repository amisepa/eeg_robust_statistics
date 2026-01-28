function [betas_obs, tvals_obs, tvals_H0, w_obs, pvals_obs, pvals_H0] = run_stats_permutation_glm_opt( ...
    Y_all, X, condition_col, nSub, nPerm, varargin)

% RUN_STATS_PERMUTATION_GLM_OPT - Permutation-based GLM with full inference
%
% Performs trial-level GLM with permutation testing for non-parametric
% inference. Returns observed statistics, null distributions, and 
% permutation-based p-values for multiple comparison correction.
%
% SYNTAX:
%   [betas_obs, tvals_obs, tvals_H0, w_obs, pvals_obs, pvals_H0] = ...
%       run_stats_permutation_glm_opt(Y_all, X, condition_col, nSub, nPerm)
%
%   [...] = run_stats_permutation_glm_opt(..., 'Name', Value)
%
% INPUTS:
%   Y_all         - [nChan x nTime x nObs] Trial-level data
%                   nChan: number of channels/electrodes
%                   nTime: number of time points
%                   nObs:  total number of trials (across all subjects)
%
%   X             - [nObs x nP] Design matrix
%                   nP: number of predictors (typically intercept + condition)
%                   Example: [ones(nObs,1), condition_col]
%
%   condition_col - [nObs x 1] Binary condition labels (0/1)
%                   Used for within-subject permutation shuffling
%
%   nSub          - Scalar. Number of subjects
%
%   nPerm         - Scalar. Number of permutations (e.g., 5000)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'Subjects'    - [nObs x 1] Subject index for each trial
%                   Required if trials are unbalanced across subjects
%                   Example: [1 1 1 2 2 2 3 3 3...]'
%                   Default: assumes equal trials per subject
%
%   'Method'      - Estimation method. Options:
%                   'OLS'  - Ordinary Least Squares (default)
%                   'IRLS' - Iteratively Reweighted Least Squares (robust)
%                   'WLS'  - Weighted Least Squares (robust)
%
%   'WeightType'  - Weight function for WLS method. Options:
%                   'PCP'   - Principal Component Projection (default)
%                   'Huber' - Huber weights
%                   'Tukey' - Tukey bisquare weights
%
%   'Progress'    - Logical. Display progress messages
%                   Default: true
%
% OUTPUTS:
%   betas_obs     - [nChan x nTime x nP] Observed regression coefficients
%                   Beta weights for each predictor at each channel/time point
%
%   tvals_obs     - [nChan x nTime x nP] Observed t-statistics
%                   T-values for testing each predictor
%
%   tvals_H0      - [nChan x nTime x nP x nPerm] Null distribution of t-values
%                   Permutation-based null distribution for each predictor
%
%   w_obs         - [nObs x 1] Observation weights
%                   For OLS: all ones
%                   For IRLS/WLS: computed robust weights
%
%   pvals_obs     - [nChan x nTime x nP] Observed permutation p-values
%                   Two-tailed p-values computed from permutation distribution
%                   p = proportion of |t_H0| >= |t_obs|
%
%   pvals_H0      - [nChan x nTime x nP x nPerm] Null distribution of p-values
%                   P-values for each permutation (for MCC procedures)
%                   Each perm tested against all other perms
%
% ALGORITHM:
%   1. Fits GLM to observed data using specified method (OLS/IRLS/WLS)
%   2. Performs within-subject permutation of condition labels
%   3. Refits GLM for each permutation to build null distribution
%   4. Computes permutation p-values (two-tailed):
%      - pvals_obs: observed statistic vs null distribution
%      - pvals_H0:  each permutation vs other permutations
%
% PERMUTATION SCHEME:
%   Within-subject label shuffling: condition labels are randomly shuffled
%   within each subject to preserve within-subject correlation structure
%   while breaking the condition-response association under H0.
%
% EXAMPLE 1: Basic two-condition comparison
%   % Data: crash [12 x 251 x 14], nocrash [12 x 251 x 14]
%   Y_all = cat(3, crash, nocrash);
%   condition_col = [ones(14,1); zeros(14,1)];
%   X = [ones(28,1), condition_col];
%   subj_idx = [1:14, 1:14]';
%   
%   [betas, tvals, tvals_H0, w, pvals, pvals_H0] = ...
%       run_stats_permutation_glm_opt(Y_all, X, condition_col, 14, 5000, ...
%       'Subjects', subj_idx, 'Method', 'OLS');
%   
%   % Apply multiple comparison correction
%   mask = compute_mcc(tvals(:,:,2), pvals(:,:,2), ...
%                      tvals_H0(:,:,2,:), pvals_H0(:,:,2,:), ...
%                      'cluster', 0.05, chanlocs);
%
% EXAMPLE 2: Robust regression with unbalanced trials
%   [betas, tvals, tvals_H0, w, pvals, pvals_H0] = ...
%       run_stats_permutation_glm_opt(Y_all, X, condition_col, 20, 10000, ...
%       'Subjects', subj_idx, 'Method', 'IRLS', 'Progress', true);
%
% NOTES:
%   - Assumes condition predictor is in column 2 of design matrix X
%   - For two-tailed tests, uses absolute value of t-statistics
%   - Permutation p-values are exact (not parametric approximations)
%   - pvals_H0 is needed for cluster-based or max-based corrections
%   - Larger nPerm (5000-10000) recommended for stable p-value estimates
%
% DEPENDENCIES:
%   - irls_glm.m (if using 'IRLS' method)
%   - compute_wls_weights.m (if using 'WLS' method)
%
% SEE ALSO:
%   compute_mcc, permutation_test, cluster_correction
%
% REFERENCES:
%   Pernet, C. R., et al. (2015). LIMO EEG: A toolbox for hierarchical 
%   linear modeling of electroencephalographic data. Computational 
%   Intelligence and Neuroscience.
%
%   Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing 
%   of EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190.
% 
% Cedric Cannard, Jan 2026

% ---------------- options ----------------
p = inputParser;
addParameter(p, 'Subjects', []);
addParameter(p, 'Method', 'OLS');
addParameter(p, 'WeightType', 'PCP');
addParameter(p, 'Progress', true);
parse(p, varargin{:});
subj_idx   = p.Results.Subjects;
method     = upper(p.Results.Method);
weightType = upper(p.Results.WeightType);
show_prog  = p.Results.Progress;

[nChan, nTime, nObs] = size(Y_all);
nP = size(X,2);

% ---------------- subject indexing ----------------
if isempty(subj_idx)
    trials_per_sub = nObs / nSub;
    if mod(trials_per_sub,1) ~= 0
        error('Provide Subjects index if trials are unbalanced.');
    end
    subj_idx = repelem((1:nSub)', trials_per_sub);
else
    subj_idx = subj_idx(:);
end

% ---------------- vectorize data ----------------
Y = reshape(permute(Y_all, [3 1 2]), nObs, nChan*nTime);
rankX = rank(X);
dof   = max(nObs - rankX, 1);

% ---------------- estimation: observed GLM ----------------
if show_prog
    fprintf('Observed GLM (%s)\n', method);
end

switch method
    case 'OLS'
        B = pinv(X' * X) * (X' * Y);
        Res = Y - X * B;
        w_obs = ones(nObs,1);
    case 'IRLS'
        [B, Res, w_obs] = irls_glm(X, Y);
    case 'WLS'
        w_obs = compute_wls_weights(X, Y, weightType);
        W = diag(w_obs);
        B = pinv(X' * W * X) * (X' * W * Y);
        Res = Y - X * B;
    otherwise
        error('Unknown method.');
end

mse = sum(bsxfun(@times, Res.^2, w_obs), 1) / dof;
SE  = sqrt(bsxfun(@times, diag(pinv(X' * diag(w_obs) * X)), mse));
Tobs = B ./ SE;

betas_obs = permute(reshape(B', nChan, nTime, nP), [1 2 3]);
tvals_obs = permute(reshape(Tobs', nChan, nTime, nP), [1 2 3]);

% ---------------- permutations ----------------
if show_prog
    fprintf('Running %d permutations...\n', nPerm);
end

tvals_H0 = nan(nChan, nTime, nP, nPerm);
sub_rows = arrayfun(@(s) find(subj_idx==s), 1:nSub, 'uni', false);

parfor iPerm = 1:nPerm
    cond_perm = condition_col;
    % within-subject label shuffle
    for s = 1:nSub
        rows = sub_rows{s};
        cond_perm(rows) = cond_perm(rows(randperm(numel(rows))));
    end
    
    Xp = X;
    Xp(:,2) = cond_perm;   % assumes condition is column 2
    
    % Initialize variables to avoid parfor warnings
    Bp = [];
    Resp = [];
    
    switch method
        case 'OLS'
            Bp = pinv(Xp' * Xp) * (Xp' * Y);
            Resp = Y - Xp * Bp;
        case 'IRLS'
            [Bp, Resp] = irls_glm(Xp, Y, w_obs);
        case 'WLS'
            W = diag(w_obs);
            Bp = pinv(Xp' * W * Xp) * (Xp' * W * Y);
            Resp = Y - Xp * Bp;
    end
    
    msep = sum(bsxfun(@times, Resp.^2, w_obs), 1) / dof;
    SEp  = sqrt(bsxfun(@times, diag(pinv(Xp' * diag(w_obs) * Xp)), msep));
    Tp   = Bp ./ SEp;
    tvals_H0(:,:,:,iPerm) = permute(reshape(Tp', nChan, nTime, nP), [1 2 3]);
end

if show_prog
    disp('Permutation GLM complete.');
end

% ============ COMPUTE PERMUTATION P-VALUES ============
if show_prog
    fprintf('Computing permutation p-values...\n');
end

% Compute p-values from permutation distribution (two-tailed)
pvals_obs = zeros(nChan, nTime, nP);
pvals_H0 = zeros(nChan, nTime, nP, nPerm);

% Progress tracking
total_iterations = nP * nChan * nTime;
completed = 0;
last_percent = 0;

if show_prog
    fprintf('Progress: 0%%');
end

for iP = 1:nP
    for iChan = 1:nChan
        for iTime = 1:nTime
            t_obs = tvals_obs(iChan, iTime, iP);
            t_null = squeeze(tvals_H0(iChan, iTime, iP, :));
            
            % Two-tailed: proportion of |t_null| >= |t_obs|
            pvals_obs(iChan, iTime, iP) = mean(abs(t_null) >= abs(t_obs));
            
            % For H0 distribution: each perm vs all other perms
            for iPerm = 1:nPerm
                t_this = t_null(iPerm);
                t_others = t_null([1:iPerm-1, iPerm+1:end]);
                pvals_H0(iChan, iTime, iP, iPerm) = mean(abs(t_others) >= abs(t_this));
            end
            
            % Update progress
            if show_prog
                completed = completed + 1;
                current_percent = floor(100 * completed / total_iterations);
                if current_percent > last_percent && mod(current_percent, 10) == 0
                    fprintf('\b\b\b\b%3d%%', current_percent);
                    last_percent = current_percent;
                end
            end
        end
    end
end

if show_prog
    fprintf('\b\b\b\b100%%\n');
    fprintf('P-value computation complete.\n');
    fprintf('Summary: %d channels, %d time points, %d predictors\n', nChan, nTime, nP);
    fprintf('Min p-value (observed): %.6f\n', min(pvals_obs(:)));
    fprintf('Max p-value (observed): %.6f\n', max(pvals_obs(:)));
end

end
