function [betas_obs, tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation_glm( ...
    Y_all, X, phase_col, nSub, nPerm, varargin)

% Parse optional inputs
p = inputParser;
addParameter(p, 'Duration', []);
addParameter(p, 'Group', []);
parse(p, varargin{:});
duration_col = p.Results.Duration;
group_col    = p.Results.Group;

% Z-score continuous covariates for stability
if ~isempty(duration_col)
    duration_col = zscore(duration_col);
end
if ~isempty(group_col) && ~islogical(group_col) && ~all(ismember(group_col,[0 1]))
    group_col = zscore(group_col);
end

nChan = size(Y_all,1);
nPredictors = size(X,2);
nObs = size(Y_all,3);

% Auto-detect trials per subject
trials_per_sub = nObs / nSub;
if mod(trials_per_sub,1) ~= 0
    error('Number of observations is not evenly divisible by number of subjects');
end

% Precompute for observed X
XtX_inv = inv(X'*X);

% ---------------------------
% Observed GLM
% ---------------------------
fprintf('Running statistics on observed Betas...\n');
betas_obs = nan(nChan, size(Y_all,2), nPredictors);
tvals     = nan(nChan, size(Y_all,2), nPredictors);
pvals     = nan(nChan, size(Y_all,2), nPredictors);

parfor iChan = 1:nChan
    fprintf(' - channel %g/%g\n', iChan, nChan);

    betas_ch = nan(size(Y_all,2), nPredictors);
    tvals_ch = nan(size(Y_all,2), nPredictors);
    pvals_ch = nan(size(Y_all,2), nPredictors);

    for iFreq = 1:size(Y_all,2)
        y = squeeze(Y_all(iChan,iFreq,:));
        b = X \ y;  % OLS solution
        betas_ch(iFreq,:) = b;

        resid = y - X*b;
        dof = nObs - rank(X);
        mse = sum(resid.^2) / dof;
        se  = sqrt(diag(XtX_inv * mse));
        tvals_ch(iFreq,:) = b ./ se;
        pvals_ch(iFreq,:) = 2 * (1 - tcdf(abs(tvals_ch(iFreq,:)), dof));
    end

    betas_obs(iChan,:,:) = betas_ch;
    tvals(iChan,:,:)     = tvals_ch;
    pvals(iChan,:,:)     = pvals_ch;
end

% ---------------------------
% Permutation GLM
% ---------------------------
fprintf('Running permutation statistics (H0)...\n');
tvals_H0 = nan(nChan, size(Y_all,2), nPredictors, nPerm);
pvals_H0 = nan(nChan, size(Y_all,2), nPredictors, nPerm);

parfor iPerm = 1:nPerm
    fprintf('   - permutation %g/%g\n', iPerm, nPerm);

    % Permute within-subject phase labels
    phase_perm = phase_col;
    for s = 1:nSub
        idx = ((s-1)*trials_per_sub + 1):(s*trials_per_sub);
        phase_perm(idx) = phase_perm(idx(randperm(trials_per_sub)));
    end

    % Build permuted X
    X_perm = build_permuted_X_simple(phase_perm, duration_col, group_col);
    XtX_inv_perm = inv(X_perm'*X_perm);

    tvals_perm_i = nan(nChan, size(Y_all,2), nPredictors);
    pvals_perm_i = nan(nChan, size(Y_all,2), nPredictors);

    for iChan = 1:nChan
        for iFreq = 1:size(Y_all,2)
            y = squeeze(Y_all(iChan,iFreq,:));
            b = X_perm \ y;
            resid = y - X_perm*b;
            dof = nObs - rank(X_perm);
            mse = sum(resid.^2) / dof;
            se  = sqrt(diag(XtX_inv_perm * mse));
            tvals_perm_i(iChan,iFreq,:) = b ./ se;
            pvals_perm_i(iChan,iFreq,:) = 2 * (1 - tcdf(abs(tvals_perm_i(iChan,iFreq,:)), dof));
        end
    end

    tvals_H0(:,:,:,iPerm) = tvals_perm_i;
    pvals_H0(:,:,:,iPerm) = pvals_perm_i;
end

disp('Done with GLM permutation statistics.');
end

%% Subfunction: build permuted X
function Xp = build_permuted_X_simple(phase_perm, duration_col, group_col)
% Phase contrast coding
recov_vs_ret = zeros(size(phase_perm));
recov_vs_ret(phase_perm == 3) = 1;
recov_vs_ret(phase_perm == 2) = -1;

% Base predictors
Xp = [ ...
    ones(size(phase_perm)), ...
    double(phase_perm == 2), ...
    recov_vs_ret ...
];

% Optional covariates
if ~isempty(duration_col)
    Xp(:,end+1) = duration_col;
end
if ~isempty(group_col)
    Xp(:,end+1) = group_col;
    Xp(:,end+1) = group_col .* double(phase_perm == 2);
    Xp(:,end+1) = group_col .* recov_vs_ret;
end
end
