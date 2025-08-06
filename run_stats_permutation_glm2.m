function [betas_obs, tvals_obs, pvals_obs, tvals_H0, pvals_H0] = run_stats_permutation_glm2( ...
    data_cell, X, phase_col, nPerm, varargin)
% RUN_GLM_VECTORIZED_PERM
%   Vectorized GLM (OLS) across channels × freqs × trials with permutations.
%   Input data is a cell array {nSub}, each entry [nChan × nFreq × nTrials_subj].
%
% INPUTS:
%   data_cell  - cell array {nSub} with [nChan × nFreq × nTrials] per subject
%   X          - design matrix [nObs × nPredictors]
%   phase_col  - phase coding vector [nObs × 1]
%   nPerm      - number of permutations
%
% OPTIONAL:
%   'Duration' - duration covariate (optional, z-scored)
%   'Group'    - group covariate (optional, z-scored or binary)
%
% OUTPUTS:
%   betas_obs  - [nChan × nFreq × nPredictors]
%   tvals_obs  - [nChan × nFreq × nPredictors]
%   pvals_obs  - same dims
%   tvals_H0   - [nChan × nFreq × nPredictors × nPerm]
%   pvals_H0   - same dims

% Parse optional inputs
p = inputParser;
addParameter(p, 'Duration', []);
addParameter(p, 'Group', []);
parse(p, varargin{:});
duration_col = p.Results.Duration;
group_col    = p.Results.Group;

% Z-score continuous covariates
if ~isempty(duration_col)
    duration_col = zscore(duration_col);
end
if ~isempty(group_col) && ~islogical(group_col) && ~all(ismember(group_col,[0 1]))
    group_col = zscore(group_col);
end

nSub = numel(data_cell);
nChan = size(data_cell{1}, 1);
nFreq = size(data_cell{1}, 2);
nPredictors = size(X, 2);
nObs = size(X, 1);

% Detect number of trials per subject automatically
trials_per_sub = nObs / nSub;
if mod(trials_per_sub,1) ~= 0
    error('Number of trials per subject is not an integer. Check phase_col or data alignment.');
end

% Concatenate all subjects' trials along 3rd dim (trials axis)
data_all = cat(3, data_cell{:}); % [nChan × nFreq × nObs]

% Reshape to 2D for vectorized GLM: (trials × features)
Y = reshape(permute(data_all, [3 1 2]), nObs, nChan*nFreq);

% ------------------------------
% Observed GLM (OLS)
% ------------------------------
fprintf('Running observed GLM (vectorized OLS)...\n');
XtX_inv = inv(X' * X);
B = XtX_inv * (X' * Y); % [nPredictors × features]
Yhat = X * B;
resid = Y - Yhat;
dof = nObs - rank(X);
mse = sum(resid.^2, 1) / dof; % 1 × features
SE  = sqrt(bsxfun(@times, diag(XtX_inv), mse)); % [nPredictors × features]
T   = B ./ SE; % [nPredictors × features]
P   = 2 * (1 - tcdf(abs(T), dof));

% Reshape outputs
betas_obs = permute(reshape(B', nChan, nFreq, nPredictors), [1 2 3]);
tvals_obs = permute(reshape(T', nChan, nFreq, nPredictors), [1 2 3]);
pvals_obs = permute(reshape(P', nChan, nFreq, nPredictors), [1 2 3]);

% ------------------------------
% Permutations
% ------------------------------
fprintf('Running %d permutations...\n', nPerm);
tvals_H0 = nan(nChan, nFreq, nPredictors, nPerm);
pvals_H0 = nan(nChan, nFreq, nPredictors, nPerm);

parfor iPerm = 1:nPerm
    if mod(iPerm, 50) == 0
        fprintf(' Permutation %d/%d\n', iPerm, nPerm);
    end

    % Permute within subject
    phase_perm = phase_col;
    for s = 1:nSub
        idx = ((s-1)*trials_per_sub + 1):(s*trials_per_sub);
        phase_perm(idx) = phase_perm(idx(randperm(trials_per_sub)));
    end

    % Build permuted design
    X_perm = build_permuted_X_simple(phase_perm, duration_col, group_col);

    % GLM on permuted data (vectorized)
    XtX_inv_perm = inv(X_perm' * X_perm);
    Bp = XtX_inv_perm * (X_perm' * Y);
    Yhatp = X_perm * Bp;
    residp = Y - Yhatp;
    dofp = nObs - rank(X_perm);
    msep = sum(residp.^2, 1) / dofp;
    SEp  = sqrt(bsxfun(@times, diag(XtX_inv_perm), msep));
    Tp   = Bp ./ SEp;
    Pp   = 2 * (1 - tcdf(abs(Tp), dofp));

    % Store
    tvals_H0(:,:,:,iPerm) = permute(reshape(Tp', nChan, nFreq, nPredictors), [1 2 3]);
    pvals_H0(:,:,:,iPerm) = permute(reshape(Pp', nChan, nFreq, nPredictors), [1 2 3]);
end

disp('Done with vectorized GLM + permutations.');
end

%% Helper

function Xp = build_permuted_X_simple(phase_perm, duration_col, group_col)

recov_vs_ret = zeros(size(phase_perm));
recov_vs_ret(phase_perm == 3) = 1;
recov_vs_ret(phase_perm == 2) = -1;

Xp = [ ...
    ones(size(phase_perm)), ...
    double(phase_perm == 2), ...
    recov_vs_ret ...
];

if ~isempty(duration_col)
    Xp(:,end+1) = duration_col;
end
if ~isempty(group_col)
    Xp(:,end+1) = group_col;
    Xp(:,end+1) = group_col .* double(phase_perm == 2);
    Xp(:,end+1) = group_col .* recov_vs_ret;
end
end
