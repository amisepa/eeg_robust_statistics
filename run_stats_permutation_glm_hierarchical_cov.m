function [betas_obs, tvals_obs, tvals_H0, w_obs, pvals_obs, pvals_H0, dz_map, beta_subj_obs_3d] = ...
    run_stats_permutation_glm_hierarchical_cov(Y_all, X, condition_col, nSub, nPerm, varargin)

% RUN_STATS_PERMUTATION_GLM_HIERARCHICAL_COVARIATE
%
% Hierarchical (two-level) GLM with permutation testing, extended to support
% subject-level covariates at Level 2.
%
% Level 1: Within-subject trial-level GLM → subject-specific condition effects
% Level 2: Between-subject regression of condition effects on covariates
%
% SYNTAX:
%   [betas_obs, tvals_obs, tvals_H0, w_obs, pvals_obs, pvals_H0, dz_map, beta_subj] = ...
%       run_stats_permutation_glm_hierarchical_covariate(Y_all, X, condition_col, nSub, nPerm, ...
%       'Subjects', subj_idx, 'Covariate', cov_vals, ...)
%
% INPUTS (same as hierarchical version, plus):
%   Y_all, X, condition_col, nSub, nPerm - see run_stats_permutation_glm_hierarchical
%
% OPTIONAL NAME-VALUE PAIRS:
%   'Subjects'      - [nTrials x 1] Subject index for each trial (REQUIRED)
%   'Method'        - 'OLS' | 'IRLS' | 'WLS' (default: 'OLS')
%   'WeightType'    - 'PCP' | 'Huber' | 'Tukey' (for WLS)
%   'Progress'      - Logical (default: true)
%   'Parallel'      - Logical (default: true)
%
%   'Covariate'     - [nSub x 1] Continuous covariate for Level 2.
%                     If empty, defaults to one-sample t-test (intercept only).
%                     The covariate is z-scored internally.
%
%   'GroupSplit'     - 'none' | 'median' | 'custom'
%                     'median' = median-split covariate into two groups
%                     'custom' = use GroupLabels directly
%                     (default: 'none', i.e., use continuous covariate)
%
%   'GroupLabels'    - [nSub x 1] Binary group labels (0/1) for custom split.
%                     Used only when GroupSplit = 'custom'.
%
% OUTPUTS:
%   betas_obs       - [nChan x nTime x nP2] Level 2 betas
%                     nP2 = number of Level 2 predictors
%                     Col 1 = intercept (mean condition effect)
%                     Col 2 = covariate effect (if Covariate provided)
%
%   tvals_obs       - [nChan x nTime x nP2] Level 2 t-statistics
%
%   tvals_H0        - [nChan x nTime x nP2 x nPerm] Null distribution
%
%   w_obs           - [nTrials x 1] Level 1 trial weights
%
%   pvals_obs       - [nChan x nTime x nP2] Permutation p-values
%
%   pvals_H0        - [nChan x nTime x nP2 x nPerm] Null p-value distributions
%
%   dz_map          - [nChan x nTime] Cohen's d_z (mean/std of subject betas)
%
%   beta_subj_obs_3d - [nSub x nChan x nTime] Subject-specific condition betas
%                      Useful for post-hoc plotting and diagnostics
%
% PERMUTATION LOGIC:
%   - Level 1 betas are computed ONCE from the real data (not re-permuted)
%   - Permutations occur at Level 2 only, with separate schemes per predictor:
%     * Intercept (Predictor 1): Sign-flip subject betas (tests mean ≠ 0)
%     * Covariate (Predictor 2+): Shuffle covariate across subjects (tests slope ≠ 0)
%   - Sign-flipping is the exchangeability-valid approach for one-sample tests
%   - Label-shuffling is valid for testing between-subject predictors
%
% EXAMPLE - Continuous covariate (meditation experience):
%   meditation_scores = [8 4 1 2 3 6 0 7 4 6 6 0 2 6 2 6]'; % 16 subjects
%   [betas, tvals, tvals_H0, w, pvals, pvals_H0, dz, beta_sub] = ...
%       run_stats_permutation_glm_hierarchical_covariate(Y_all, X, condition_col, 16, 5000, ...
%       'Subjects', subj_idx, 'Method', 'WLS', 'WeightType', 'Huber', ...
%       'Covariate', meditation_scores);
%   % tvals(:,:,1) = intercept (main condition effect)
%   % tvals(:,:,2) = covariate effect (meditation modulates condition?)
%
% EXAMPLE - Median split:
%   [betas, tvals, tvals_H0, w, pvals, pvals_H0] = ...
%       run_stats_permutation_glm_hierarchical_covariate(Y_all, X, condition_col, 16, 5000, ...
%       'Subjects', subj_idx, 'Method', 'WLS', 'WeightType', 'Huber', ...
%       'Covariate', meditation_scores, 'GroupSplit', 'median');
%
% Cedric Cannard, 2026

% ---------------- Parse options ----------------
p = inputParser;
addParameter(p, 'Subjects', [], @(x) ~isempty(x));
addParameter(p, 'Method', 'OLS');
addParameter(p, 'WeightType', 'PCP');
addParameter(p, 'Progress', true);
addParameter(p, 'Parallel', true);
addParameter(p, 'Covariate', []);
addParameter(p, 'GroupSplit', 'none');
addParameter(p, 'GroupLabels', []);
parse(p, varargin{:});

subj_idx     = p.Results.Subjects;
method       = upper(p.Results.Method);
weightType   = upper(p.Results.WeightType);
show_prog    = p.Results.Progress;
use_parallel = p.Results.Parallel;
% covariate    = p.Results.Covariate(:);
covariate    = p.Results.Covariate;
if isvector(covariate), covariate = covariate(:); end   % column-ize only if vector
groupSplit   = lower(p.Results.GroupSplit);
groupLabels  = p.Results.GroupLabels(:);

if nargin < 5 || isempty(nPerm), nPerm = 1000; end
if isempty(subj_idx)
    error('Subjects index is required for hierarchical modeling.');
end

[nChan, nTime, nTrials] = size(Y_all);
subj_idx = subj_idx(:);

% Validate inputs
assert(length(subj_idx) == nTrials, 'Subjects index must match number of trials.');
assert(length(condition_col) == nTrials, 'condition_col must match number of trials.');

% Subject indexing
unique_subs = unique(subj_idx);
if length(unique_subs) ~= nSub
    warning('Number of unique subjects (%d) differs from nSub (%d). Using unique count.', ...
            length(unique_subs), nSub);
    nSub = length(unique_subs);
end
sub_rows = arrayfun(@(s) find(subj_idx == s), unique_subs, 'uni', false);
trials_per_sub = cellfun(@length, sub_rows);

% ============ BUILD LEVEL 2 DESIGN MATRIX ============
if ~isempty(covariate)
    % assert(length(covariate) == nSub, 'Covariate must have nSub elements.');
    assert(size(covariate, 1) == nSub, 'Covariate must have nSub rows. Got %d rows for %d subjects.', size(covariate,1), nSub);

    if size(covariate, 2) > 1 && ~strcmpi(groupSplit, 'none')
        error('GroupSplit only supported for single covariates. Got %d columns.', size(covariate,2));
    end
    
    switch groupSplit
        case 'median'
            med_val = median(covariate);
            grp = double(covariate > med_val);
            % Handle ties at median: assign to lower group
            grp(covariate == med_val) = 0;
            if show_prog
                fprintf('Median split at %.2f: Group 0 (n=%d), Group 1 (n=%d)\n', ...
                        med_val, sum(grp==0), sum(grp==1));
            end
            % Center the binary variable
            X2 = [ones(nSub,1), grp - mean(grp)];
            cov_label = 'Group (median split)';
            
        case 'custom'
            assert(~isempty(groupLabels) && length(groupLabels)==nSub, ...
                   'GroupLabels must have nSub elements for custom split.');
            grp = groupLabels;
            X2 = [ones(nSub,1), grp - mean(grp)];
            cov_label = 'Group (custom)';
            
        % otherwise  % 'none' — continuous covariate
        %     % Z-score covariate for interpretability
        %     cov_z = (covariate - mean(covariate)) / std(covariate);
        %     X2 = [ones(nSub,1), cov_z];
        %     cov_label = 'Covariate (z-scored)';
        otherwise  % 'none' — continuous covariate(s)
            nCov = size(covariate, 2);
            cov_z = (covariate - mean(covariate, 1)) ./ std(covariate, 0, 1);  % z-score each column
            X2 = [ones(nSub, 1), cov_z];   % [nSub x (1 + nCov)]
            if nCov == 1
                cov_label = 'Covariate (z-scored)';
            else
                cov_label = sprintf('%d covariates (z-scored)', nCov);
            end
    end
else
    % No covariate: intercept-only → one-sample t-test
    X2 = ones(nSub, 1);
    cov_label = 'Intercept only (one-sample t-test)';
end

nP2 = size(X2, 2);

if show_prog
    fprintf('\n=== HIERARCHICAL GLM WITH COVARIATE ===\n');
    fprintf('Level 2 design: %s (%d predictors)\n', cov_label, nP2);
end

% ============ PARALLEL SETUP ============
if use_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        if show_prog, fprintf('Starting parallel pool...\n'); end
        try
            pool = parpool;
            if show_prog, fprintf('Parallel pool: %d workers.\n', pool.NumWorkers); end
        catch
            if show_prog, warning('Could not start parallel pool. Running serially.'); end
            use_parallel = false;
        end
    else
        if show_prog, fprintf('Parallel pool: %d workers.\n', pool.NumWorkers); end
    end
end

% ============ VECTORIZE DATA ============
Y = reshape(permute(Y_all, [3 1 2]), nTrials, nChan*nTime);

% ============ LEVEL 1: WITHIN-SUBJECT MODELING ============
if show_prog
    fprintf('\nLevel 1: Within-subject modeling (%s)\n', method);
    fprintf('  Subjects: %d | Trials: %d | Channels: %d | Time: %d\n', ...
            nSub, nTrials, nChan, nTime);
    fprintf('  Trials/subject: min=%d, max=%d, mean=%.1f\n', ...
            min(trials_per_sub), max(trials_per_sub), mean(trials_per_sub));
    level1_start = tic;
end

beta_subj_obs = zeros(nSub, nChan*nTime);
w_obs = ones(nTrials, 1);

for iSub = 1:nSub
    rows = sub_rows{iSub};
    Y_sub = Y(rows, :);
    X_sub = X(rows, :);
    
    switch method
        case 'OLS'
            B_sub = pinv(X_sub' * X_sub) * (X_sub' * Y_sub);
        case 'IRLS'
            [B_sub, ~, w_sub] = irls_glm(X_sub, Y_sub);
            w_obs(rows) = w_sub;
        case 'WLS'
            w_sub = compute_wls_weights(X_sub, Y_sub, weightType);
            W_sub = diag(w_sub);
            B_sub = pinv(X_sub' * W_sub * X_sub) * (X_sub' * W_sub * Y_sub);
            w_obs(rows) = w_sub;
    end
    
    % Extract condition effect (column 2 of design matrix)
    beta_subj_obs(iSub, :) = B_sub(2, :);
    
    if show_prog && mod(iSub, max(1, ceil(nSub/10))) == 0
        fprintf('  Subject %d/%d\n', iSub, nSub);
    end
end

if show_prog
    fprintf('Level 1 complete in %.2f seconds\n', toc(level1_start));
end

% Reshape to [nSub x nChan x nTime]
beta_subj_obs_3d = permute(reshape(beta_subj_obs', nChan, nTime, nSub), [3 1 2]);

% ============ LEVEL 2: GROUP-LEVEL WITH COVARIATE (OBSERVED) ============
if show_prog
    fprintf('\nLevel 2: Group-level analysis\n');
    level2_start = tic;
end

% Fit Level 2 GLM: beta_subj = X2 * B2 + error
% beta_subj is [nSub x nVox], X2 is [nSub x nP2]
nVox = nChan * nTime;
B2 = pinv(X2' * X2) * (X2' * beta_subj_obs);  % [nP2 x nVox]
Res2 = beta_subj_obs - X2 * B2;
dof2 = max(nSub - rank(X2), 1);
mse2 = sum(Res2.^2, 1) / dof2;                 % [1 x nVox]
cov_XtX_inv = diag(pinv(X2' * X2));            % [nP2 x 1]
SE2 = sqrt(cov_XtX_inv * mse2);                % [nP2 x nVox]
T2_obs = B2 ./ SE2;                            % [nP2 x nVox]

% Reshape outputs
betas_obs = zeros(nChan, nTime, nP2);
tvals_obs = zeros(nChan, nTime, nP2);
for iP = 1:nP2
    betas_obs(:,:,iP) = reshape(B2(iP,:), nChan, nTime);
    tvals_obs(:,:,iP) = reshape(T2_obs(iP,:), nChan, nTime);
end

% Cohen's d_z (from raw subject betas, independent of covariate)
betas_mean = squeeze(mean(beta_subj_obs_3d, 1));
betas_std  = squeeze(std(beta_subj_obs_3d, 0, 1));
dz_map = betas_mean ./ betas_std;

if show_prog
    fprintf('Level 2 complete in %.2f seconds\n', toc(level2_start));
end

% ============ PERMUTATIONS AT LEVEL 2 ============
% Two permutation schemes run simultaneously:
%   - Predictor 1 (intercept): Sign-flip subject betas (tests mean ≠ 0)
%   - Predictor 2+ (covariate): Shuffle covariate assignment (tests slope ≠ 0)
% These are independent and both run within the same loop for efficiency.

if show_prog
    fprintf('\n=== PERMUTATION TESTING (Level 2) ===\n');
    if nP2 > 1
        fprintf('Scheme: sign-flip (intercept) + label-shuffle (covariate)\n');
    else
        fprintf('Scheme: sign-flip (intercept only)\n');
    end
    fprintf('Permutations: %d\n', nPerm);
    perm_start = tic;
end

tvals_H0 = zeros(nChan, nTime, nP2, nPerm);

% Pre-extract Level 2 data for permutation
beta_subj_for_perm = beta_subj_obs;  % [nSub x nVox]

if use_parallel
    % Parallel permutations
    parfor iPerm = 1:nPerm
        tmp = zeros(nChan, nTime, nP2);
        
        % --- Intercept null: sign-flip ---
        signs = 2*(randi(2, nSub, 1) - 1) - 1;  % random +1/-1
        Y2_signed = bsxfun(@times, beta_subj_for_perm, signs);
        B2_int = pinv(X2' * X2) * (X2' * Y2_signed);
        Res2_int = Y2_signed - X2 * B2_int;
        mse2_int = sum(Res2_int.^2, 1) / dof2;
        SE2_int = sqrt(cov_XtX_inv * mse2_int);
        T2_int = B2_int ./ SE2_int;
        tmp(:,:,1) = reshape(T2_int(1,:), nChan, nTime);
        
        % --- Covariate null: label-shuffle ---
        if nP2 > 1
            perm_idx = randperm(nSub);
            X2p = X2(perm_idx, :);
            B2_cov = pinv(X2p' * X2p) * (X2p' * beta_subj_for_perm);
            Res2_cov = beta_subj_for_perm - X2p * B2_cov;
            mse2_cov = sum(Res2_cov.^2, 1) / dof2;
            cov_inv_p = diag(pinv(X2p' * X2p));
            SE2_cov = sqrt(cov_inv_p * mse2_cov);
            T2_cov = B2_cov ./ SE2_cov;
            for iP = 2:nP2
                tmp(:,:,iP) = reshape(T2_cov(iP,:), nChan, nTime);
            end
        end
        
        tvals_H0(:,:,:,iPerm) = tmp;
    end
else
    % Serial permutations
    if show_prog
        h_wait = waitbar(0, 'Starting permutations...', 'Name', 'Permutation Progress');
        cleanup_waitbar = onCleanup(@() close(h_wait));
    end
    
    for iPerm = 1:nPerm
        % --- Intercept null: sign-flip ---
        signs = 2*(randi(2, nSub, 1) - 1) - 1;
        Y2_signed = bsxfun(@times, beta_subj_for_perm, signs);
        B2_int = pinv(X2' * X2) * (X2' * Y2_signed);
        Res2_int = Y2_signed - X2 * B2_int;
        mse2_int = sum(Res2_int.^2, 1) / dof2;
        SE2_int = sqrt(cov_XtX_inv * mse2_int);
        T2_int = B2_int ./ SE2_int;
        tvals_H0(:,:,1,iPerm) = reshape(T2_int(1,:), nChan, nTime);
        
        % --- Covariate null: label-shuffle ---
        if nP2 > 1
            perm_idx = randperm(nSub);
            X2p = X2(perm_idx, :);
            B2_cov = pinv(X2p' * X2p) * (X2p' * beta_subj_for_perm);
            Res2_cov = beta_subj_for_perm - X2p * B2_cov;
            mse2_cov = sum(Res2_cov.^2, 1) / dof2;
            cov_inv_p = diag(pinv(X2p' * X2p));
            SE2_cov = sqrt(cov_inv_p * mse2_cov);
            T2_cov = B2_cov ./ SE2_cov;
            for iP = 2:nP2
                tvals_H0(:,:,iP,iPerm) = reshape(T2_cov(iP,:), nChan, nTime);
            end
        end
        
        if show_prog && mod(iPerm, max(1, floor(nPerm/100))) == 0
            progress = iPerm / nPerm;
            elapsed = toc(perm_start);
            remaining = elapsed / progress - elapsed;
            waitbar(progress, h_wait, sprintf('Perm %d/%d | ~%.1fm remaining', ...
                    iPerm, nPerm, remaining/60));
        end
    end
end

if show_prog
    perm_time = toc(perm_start);
    fprintf('Permutations complete in %.2f minutes (%.4f sec/perm)\n', ...
            perm_time/60, perm_time/nPerm);
end

% ============ COMPUTE PERMUTATION P-VALUES ============
if show_prog
    fprintf('\nComputing permutation p-values...\n');
    pval_start = tic;
end

pvals_obs = zeros(nChan, nTime, nP2);
pvals_H0 = zeros(nChan, nTime, nP2, nPerm);

for iP = 1:nP2
    for iChan = 1:nChan
        for iTime = 1:nTime
            t_obs = tvals_obs(iChan, iTime, iP);
            t_null = squeeze(tvals_H0(iChan, iTime, iP, :));
            
            % Two-tailed p-value
            pvals_obs(iChan, iTime, iP) = mean(abs(t_null) >= abs(t_obs));
            
            % H0 distribution of p-values (for MCC)
            for iPerm = 1:nPerm
                t_this = t_null(iPerm);
                t_others = t_null([1:iPerm-1, iPerm+1:end]);
                pvals_H0(iChan, iTime, iP, iPerm) = mean(abs(t_others) >= abs(t_this));
            end
        end
    end
end

if show_prog
    pval_time = toc(pval_start);
    fprintf('P-values complete in %.2f seconds\n', pval_time);
    
    fprintf('\n=== SUMMARY ===\n');
    fprintf('Data: %d channels × %d time points × %d subjects\n', nChan, nTime, nSub);
    fprintf('Level 2: %s\n', cov_label);
    fprintf('Permutations: %d (sign-flip for intercept', nPerm);
    if nP2 > 1
        fprintf(', label-shuffle for covariate)\n');
    else
        fprintf(')\n');
    end
    for iP = 1:nP2
        pv = pvals_obs(:,:,iP);
        fprintf('\nPredictor %d:\n', iP);
        fprintf('  Min p: %.6f | Median p: %.4f\n', min(pv(:)), median(pv(:)));
        fprintf('  < 0.05: %d (%.1f%%) | < 0.01: %d (%.1f%%)\n', ...
                sum(pv(:)<0.05), 100*mean(pv(:)<0.05), ...
                sum(pv(:)<0.01), 100*mean(pv(:)<0.01));
    end
end

end