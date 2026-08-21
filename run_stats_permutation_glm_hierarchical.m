function [betas_obs, tvals_obs, tvals_H0, w_obs, pvals_obs, pvals_H0, dz_map] = ...
    run_stats_permutation_glm_hierarchical(Y_all, X, condition_col, nSub, nPerm, varargin)

% RUN_STATS_PERMUTATION_GLM_HIERARCHICAL - Hierarchical GLM with permutation testing
%
% Performs hierarchical (two-level) GLM with permutation testing:
%   Level 1: Within-subject trial-level modeling (subject-specific effects)
%   Level 2: Between-subject group-level modeling (random effects)
%
% This approach is more powerful when you have multiple trials per subject
% and want to model both within-subject and between-subject variability.
%
% SYNTAX:
%   [betas_obs, tvals_obs, tvals_H0, w_obs, pvals_obs, pvals_H0] = ...
%       run_stats_permutation_glm_hierarchical(Y_all, X, condition_col, nSub, nPerm)
%
%   [...] = run_stats_permutation_glm_hierarchical(..., 'Name', Value)
%
% INPUTS:
%   Y_all         - [nChan x nTime x nTrials] Trial-level data
%                   nTrials = total trials across all subjects
%
%   X             - [nTrials x nP] Design matrix
%                   Typically: [ones(nTrials,1), condition_col, subject_dummies...]
%
%   condition_col - [nTrials x 1] Binary condition labels (0/1)
%                   Used for within-subject permutation shuffling
%
%   nSub          - Scalar. Number of subjects
%
%   nPerm         - Scalar. Number of permutations (default: 1000)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'Subjects'    - [nTrials x 1] Subject index for each trial (REQUIRED)
%
%   'Method'      - Estimation method for Level 1. Options:
%                   'OLS'  - Ordinary Least Squares (default)
%                   'IRLS' - Iteratively Reweighted Least Squares (robust)
%                   'WLS'  - Weighted Least Squares (robust)
%
%   'WeightType'  - Weight function for WLS method. Options:
%                   'PCP'   - Principal Component Projection (default)
%                   'Huber' - Huber weights
%                   'Tukey' - Tukey bisquare weights
%
%   'Progress'    - Logical. Display progress messages (default: true)
%
%   'Parallel'    - Logical. Use parallel computing for permutations (default: true)
%                   Requires Parallel Computing Toolbox
%
% OUTPUTS:
%   betas_obs     - [nChan x nTime] Group-level condition effect
%                   Average condition effect across subjects
%
%   tvals_obs     - [nChan x nTime] Observed t-statistics
%                   T-values testing group-level condition effect
%
%   tvals_H0      - [nChan x nTime x nPerm] Null distribution of t-values
%                   Permutation-based null distribution
%
%   w_obs         - [nTrials x 1] Trial-level observation weights
%                   For OLS: all ones
%                   For IRLS/WLS: computed robust weights
%
%   pvals_obs     - [nChan x nTime] Observed permutation p-values
%                   Two-tailed p-values from permutation distribution
%
%   pvals_H0      - [nChan x nTime x nPerm] Null distribution of p-values
%                   P-values for each permutation (for MCC procedures)
%
% ALGORITHM:
%   Level 1 (Within-subject):
%     For each subject s:
%       1. Fit GLM: Y_s = X_s * beta_s + error_s
%       2. Extract subject-specific condition effect: beta_s(condition)
%       3. Optionally use robust estimation (IRLS/WLS)
%
%   Level 2 (Between-subject):
%     1. Stack subject-specific effects: [beta_1, beta_2, ..., beta_S]
%     2. Test if mean(beta_s) ≠ 0 using one-sample t-test
%     3. Permutation: shuffle condition labels within each subject
%
% PERMUTATION SCHEME:
%   Within-subject label shuffling at Level 1, then recompute Level 2
%   statistics. Preserves within-subject correlation structure.
%
% EXAMPLE:
%   % Trial-level data: 30 trials/subject/condition, 14 subjects
%   % crash_trials: [12 x 251 x 30 x 14]
%   % nocrash_trials: [12 x 251 x 30 x 14]
%
%   crash_all = reshape(crash_trials, 12, 251, []);     % [12 x 251 x 420]
%   nocrash_all = reshape(nocrash_trials, 12, 251, []); % [12 x 251 x 420]
%   Y_all = cat(3, crash_all, nocrash_all);             % [12 x 251 x 840]
%
%   condition_col = [ones(420,1); zeros(420,1)];
%   X = [ones(840,1), condition_col];
%   subj_idx = [repmat((1:14)', 30, 1); repmat((1:14)', 30, 1)];
%
%   [betas, tvals, tvals_H0, w, pvals, pvals_H0] = ...
%       run_stats_permutation_glm_hierarchical(Y_all, X, condition_col, 14, 1000, ...
%       'Subjects', subj_idx, 'Method', 'IRLS');
%
% NOTES:
%   - More powerful than single-level when trials >> subjects
%   - Accounts for both within- and between-subject variability
%   - Subject-specific effects allow for individual differences
%   - Robust methods at Level 1 protect against outlier trials
%
% SEE ALSO:
%   run_stats_permutation_glm_opt, compute_robust_erp

% ---------------- options ----------------
p = inputParser;
addParameter(p, 'Subjects', [], @(x) ~isempty(x));
addParameter(p, 'Method', 'OLS');
addParameter(p, 'WeightType', 'PCP');
addParameter(p, 'Progress', true);
addParameter(p, 'Parallel', true);
parse(p, varargin{:});

subj_idx   = p.Results.Subjects;
method     = upper(p.Results.Method);
weightType = upper(p.Results.WeightType);
show_prog  = p.Results.Progress;
use_parallel = p.Results.Parallel;

% Set default nPerm if not provided
if nargin < 5 || isempty(nPerm)
    nPerm = 1000;
end

if isempty(subj_idx)
    error('Subjects index is required for hierarchical modeling.');
end

[nChan, nTime, nTrials] = size(Y_all);
subj_idx = subj_idx(:);

% ---------------- validate inputs ----------------
if length(subj_idx) ~= nTrials
    error('Length of Subjects index must match number of trials.');
end

if length(condition_col) ~= nTrials
    error('Length of condition_col must match number of trials.');
end

% ---------------- parallel setup ----------------
if use_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        if show_prog
            fprintf('Starting parallel pool...\n');
        end
        try
            pool = parpool;
            if show_prog
                fprintf('Parallel pool started with %d workers.\n', pool.NumWorkers);
            end
        catch
            if show_prog
                warning('Could not start parallel pool. Running serially.');
            end
            use_parallel = false;
        end
    else
        if show_prog
            fprintf('Using existing parallel pool with %d workers.\n', pool.NumWorkers);
        end
    end
else
    if show_prog
        fprintf('Parallel computing disabled. Running serially.\n');
    end
end

% ---------------- subject indexing ----------------
unique_subs = unique(subj_idx);
if length(unique_subs) ~= nSub
    warning('Number of unique subjects (%d) differs from nSub (%d). Using unique count.', ...
        length(unique_subs), nSub);
    nSub = length(unique_subs);
end

% Get trial indices for each subject
sub_rows = arrayfun(@(s) find(subj_idx == s), unique_subs, 'uni', false);
trials_per_sub = cellfun(@length, sub_rows);

% ---------------- vectorize data ----------------
Y = reshape(permute(Y_all, [3 1 2]), nTrials, nChan*nTime);

% ---------------- Level 1: Within-subject modeling (OBSERVED) ----------------
if show_prog
    fprintf('\n=== HIERARCHICAL GLM ===\n');
    fprintf('Level 1: Within-subject modeling (%s)\n', method);
    fprintf('  Subjects: %d\n', nSub);
    fprintf('  Total trials: %d\n', nTrials);
    fprintf('  Trials/subject: min=%d, max=%d, mean=%.1f\n', ...
        min(trials_per_sub), max(trials_per_sub), mean(trials_per_sub));
    fprintf('  Channels: %d, Time points: %d\n', nChan, nTime);
    level1_start = tic;
end

% Store subject-specific condition effects
beta_subj_obs = zeros(nSub, nChan*nTime);  % [nSub x nVox]
w_obs = ones(nTrials, 1);

for iSub = 1:nSub
    rows = sub_rows{iSub};
    Y_sub = Y(rows, :);
    X_sub = X(rows, :);

    % Fit within-subject GLM
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

    % Extract condition effect (assumes column 2)
    beta_subj_obs(iSub, :) = B_sub(2, :);

    if show_prog && mod(iSub, max(1, ceil(nSub/10))) == 0
        fprintf('  Subject %d/%d (%.1f%%)\n', iSub, nSub, 100*iSub/nSub);
    end
end

if show_prog
    fprintf('Level 1 complete in %.2f seconds\n', toc(level1_start));
end

% ---------------- Level 2: Group-level analysis (OBSERVED) ----------------
if show_prog
    fprintf('\nLevel 2: Group-level analysis\n');
    level2_start = tic;
end

% Reshape to [nSub x nChan x nTime]
beta_subj_obs_3d = permute(reshape(beta_subj_obs', nChan, nTime, nSub), [3 1 2]);

% Group mean and standard error
betas_obs_mean = squeeze(mean(beta_subj_obs_3d, 1));  % [nChan x nTime]
betas_obs_std = squeeze(std(beta_subj_obs_3d, 0, 1)); % [nChan x nTime]
betas_obs_se = betas_obs_std / sqrt(nSub);

% Compute Cohen's d_z directly from subject-specific effects
dz_map = betas_obs_mean ./ betas_obs_std;  % This is d_z!

% One-sample t-test against zero
tvals_obs = betas_obs_mean ./ betas_obs_se;  % [nChan x nTime]
betas_obs = betas_obs_mean;

if show_prog
    fprintf('Level 2 complete in %.2f seconds\n', toc(level2_start));
end

% ---------------- Permutations ----------------
if show_prog
    fprintf('\n=== PERMUTATION TESTING ===\n');
    fprintf('Running %d permutations%s\n', nPerm, ...
        ternary(use_parallel, ' (parallel)', ' (serial)'));
    fprintf('  Each permutation fits %d subject-level GLMs\n', nSub);
    perm_start = tic;
end

tvals_H0 = nan(nChan, nTime, nPerm);

% Progress tracking variables
if show_prog && ~use_parallel
    % Only use waitbar for serial execution
    h_wait = waitbar(0, 'Starting permutations...', 'Name', 'Permutation Progress');
    cleanup_waitbar = onCleanup(@() close(h_wait));
end

% Milestone tracking for parallel
completed_milestones = false(10, 1);  % Track 10%, 20%, ..., 100%

if use_parallel
    % Parallel execution
    parfor iPerm = 1:nPerm
        % Permute condition labels within each subject
        cond_perm = condition_col;
        for s = 1:nSub
            rows = sub_rows{s};
            cond_perm(rows) = cond_perm(rows(randperm(numel(rows))));
        end

        % Update design matrix
        Xp = X;
        Xp(:, 2) = cond_perm;

        % Level 1: Within-subject modeling (permuted)
        beta_subj_perm = zeros(nSub, nChan*nTime);

        for iSub = 1:nSub
            rows = sub_rows{iSub};
            Y_sub = Y(rows, :);
            Xp_sub = Xp(rows, :);

            % Initialize to avoid parfor warnings
            B_sub = [];

            switch method
                case 'OLS'
                    B_sub = pinv(Xp_sub' * Xp_sub) * (Xp_sub' * Y_sub);

                case 'IRLS'
                    [B_sub, ~] = irls_glm(Xp_sub, Y_sub, w_obs(rows));

                case 'WLS'
                    W_sub = diag(w_obs(rows));
                    B_sub = pinv(Xp_sub' * W_sub * Xp_sub) * (Xp_sub' * W_sub * Y_sub);
            end

            beta_subj_perm(iSub, :) = B_sub(2, :);
        end

        % Level 2: Group-level analysis (permuted)
        beta_subj_perm_3d = permute(reshape(beta_subj_perm', nChan, nTime, nSub), [3 1 2]);
        betas_perm_mean = squeeze(mean(beta_subj_perm_3d, 1));
        betas_perm_std = squeeze(std(beta_subj_perm_3d, 0, 1));
        betas_perm_se = betas_perm_std / sqrt(nSub);

        tvals_H0(:, :, iPerm) = betas_perm_mean ./ betas_perm_se;
    end

    % After parfor completes, report completion
    if show_prog
        perm_time = toc(perm_start);
        fprintf('\nAll permutations complete!\n');
        fprintf('  Total time: %.2f minutes (%.2f sec/perm)\n', ...
            perm_time/60, perm_time/nPerm);
    end

else
    % Serial execution with waitbar
    for iPerm = 1:nPerm
        % Permute condition labels within each subject
        cond_perm = condition_col;
        for s = 1:nSub
            rows = sub_rows{s};
            cond_perm(rows) = cond_perm(rows(randperm(numel(rows))));
        end

        % Update design matrix
        Xp = X;
        Xp(:, 2) = cond_perm;

        % Level 1: Within-subject modeling (permuted)
        beta_subj_perm = zeros(nSub, nChan*nTime);

        for iSub = 1:nSub
            rows = sub_rows{iSub};
            Y_sub = Y(rows, :);
            Xp_sub = Xp(rows, :);

            B_sub = [];

            switch method
                case 'OLS'
                    B_sub = pinv(Xp_sub' * Xp_sub) * (Xp_sub' * Y_sub);

                case 'IRLS'
                    [B_sub, ~] = irls_glm(Xp_sub, Y_sub, w_obs(rows));

                case 'WLS'
                    W_sub = diag(w_obs(rows));
                    B_sub = pinv(Xp_sub' * W_sub * Xp_sub) * (Xp_sub' * W_sub * Y_sub);
            end

            beta_subj_perm(iSub, :) = B_sub(2, :);
        end

        % Level 2: Group-level analysis (permuted)
        beta_subj_perm_3d = permute(reshape(beta_subj_perm', nChan, nTime, nSub), [3 1 2]);
        betas_perm_mean = squeeze(mean(beta_subj_perm_3d, 1));
        betas_perm_std = squeeze(std(beta_subj_perm_3d, 0, 1));
        betas_perm_se = betas_perm_std / sqrt(nSub);

        tvals_H0(:, :, iPerm) = betas_perm_mean ./ betas_perm_se;

        % Update waitbar
        if show_prog && mod(iPerm, max(1, floor(nPerm/100))) == 0
            progress = iPerm / nPerm;
            elapsed = toc(perm_start);
            remaining = elapsed / progress - elapsed;
            waitbar(progress, h_wait, ...
                sprintf('Permutation %d/%d (%.1f%%) | Elapsed: %.1fm | Remaining: ~%.1fm', ...
                iPerm, nPerm, 100*progress, elapsed/60, remaining/60));
        end
    end

    if show_prog
        perm_time = toc(perm_start);
        fprintf('\nPermutation complete in %.2f minutes (%.2f sec/perm)\n', ...
            perm_time/60, perm_time/nPerm);
    end
end

% ============ COMPUTE PERMUTATION P-VALUES ============
if show_prog
    fprintf('\n=== COMPUTING P-VALUES ===\n');
    pval_start = tic;
    h_wait_pval = waitbar(0, 'Computing p-values: 0%', 'Name', 'P-value Progress');
    cleanup_pval = onCleanup(@() close_if_valid(h_wait_pval));
end

tvals_obs = reshape(tvals_obs, nChan, nTime);
tvals_H0  = reshape(tvals_H0,  nChan, nTime, nPerm);

pvals_obs = zeros(nChan, nTime);
pvals_H0  = zeros(nChan, nTime, nPerm);

total_iter = nChan * nTime;
completed  = 0;
last_pct   = 0;

for iChan = 1:nChan
    for iTime = 1:nTime
        t_obs  = tvals_obs(iChan, iTime);
        t_null = squeeze(tvals_H0(iChan, iTime, :));   % [nPerm x 1]

        % Observed p-value. The observed statistic is included in the null
        % (Phipson & Smyth 2010): mean(...) alone can return exactly 0, which
        % is not a possible p-value for a finite permutation test.
        pvals_obs(iChan, iTime) = (1 + sum(abs(t_null) >= abs(t_obs))) / (1 + nPerm);

        % Vectorized H0 p-values: each perm vs all others
        abs_null = abs(t_null);                          % [nPerm x 1]
        total_ge = sum(abs_null >= abs_null', 1);        % [1 x nPerm]
        pvals_H0(iChan, iTime, :) = (total_ge - 1) / (nPerm - 1);

        % Progress
        completed   = completed + 1;
        current_pct = floor(100 * completed / total_iter);
        if show_prog && current_pct > last_pct && mod(current_pct, 10) == 0
            elapsed   = toc(pval_start);
            remaining = elapsed / completed * (total_iter - completed);
            fprintf('  P-values: %d%% | Elapsed: %.1fs | Remaining: ~%.1fs\n', ...
                current_pct, elapsed, remaining);
            waitbar(completed/total_iter, h_wait_pval, ...
                sprintf('P-values: %d%% | ~%.1fs remaining', current_pct, remaining));
            last_pct = current_pct;
        end
    end
end

if show_prog
    pval_time = toc(pval_start);
    fprintf('P-values complete in %.2f seconds\n', pval_time);
    close_if_valid(h_wait_pval)

    fprintf('\n=== SUMMARY ===\n');
    fprintf('Data dimensions: %d channels × %d time points\n', nChan, nTime);
    fprintf('Subjects: %d (avg %.1f trials/subject)\n', nSub, mean(trials_per_sub));
    fprintf('Permutations: %d\n', nPerm);
    fprintf('Total processing time: %.2f minutes\n', (toc(perm_start) + pval_time)/60);
    fprintf('\nObserved p-values:\n');
    fprintf('  Min: %.6f\n', min(pvals_obs(:)));
    fprintf('  Max: %.6f\n', max(pvals_obs(:)));
    fprintf('  Median: %.6f\n', median(pvals_obs(:)));
    fprintf('  < 0.05: %d voxels (%.2f%%)\n', ...
        sum(pvals_obs(:) < 0.05), 100*mean(pvals_obs(:) < 0.05));
    fprintf('  < 0.01: %d voxels (%.2f%%)\n', ...
        sum(pvals_obs(:) < 0.01), 100*mean(pvals_obs(:) < 0.01));

end

end

% Helper function for conditional expression
function out = ternary(condition, true_val, false_val)
if condition
    out = true_val;
else
    out = false_val;
end
end

% helper
function close_if_valid(h)
    if ishandle(h)
        close(h)
    end
end
