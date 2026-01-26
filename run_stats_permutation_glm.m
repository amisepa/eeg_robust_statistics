function [betas_obs, tvals_obs, pvals_obs, tvals_H0, pvals_H0] = run_stats_permutation_glm( ...
    data_cell, X, phase_col, nPerm, varargin)
% RUN_STATS_PERMUTATION_GLM
% Vectorized OLS GLM across channels × features with within-subject permutations.
%
% Outputs:
% - tvals_H0 is the full null distribution of t-values (correct for later corrections).
% - pvals_obs are empirical permutation p-values (two-tailed, uncorrected).
% - pvals_H0 is returned empty (tcdf-based p-values are omitted by design).
%
% INPUTS:
%   data_cell  - cell array {nSub}, each [nChan × nFeat × nTrials_subj]
%   X          - observed design matrix [nObs × nPredictors]
%   phase_col  - coding vector [nObs × 1] permuted within subject
%   nPerm      - number of permutations
%
% OPTIONAL:
%   'Duration' - duration covariate (optional, z-scored)
%   'Group'    - group covariate (optional, z-scored unless binary/logic)
%   'Progress' - true/false (default true)
%
% OUTPUTS:
%   betas_obs  - [nChan × nFeat × nPredictors]
%   tvals_obs  - [nChan × nFeat × nPredictors]
%   pvals_obs  - [nChan × nFeat × nPredictors] empirical permutation p-values
%   tvals_H0   - [nChan × nFeat × nPredictors × nPerm]
%   pvals_H0   - [] (omitted)

% Parse optional inputs
p = inputParser;
addParameter(p, 'Duration', []);
addParameter(p, 'Group', []);
addParameter(p, 'Progress', true);
parse(p, varargin{:});

duration_col = p.Results.Duration;
group_col    = p.Results.Group;
show_prog    = p.Results.Progress;

% Z-score continuous covariates
if ~isempty(duration_col)
    duration_col = zscore(duration_col(:));
end
if ~isempty(group_col)
    group_col = group_col(:);
    if ~islogical(group_col) && ~all(ismember(unique(group_col), [0 1]))
        group_col = zscore(group_col);
    end
end

% Dimensions
nSub = numel(data_cell);
nChan = size(data_cell{1}, 1);
nFeat = size(data_cell{1}, 2);
nPredictors = size(X, 2);
nObs = size(X, 1);

if numel(phase_col) ~= nObs
    error('phase_col must be length nObs');
end
phase_col = phase_col(:);

% Detect trials per subject
trials_per_sub = nObs / nSub;
if mod(trials_per_sub, 1) ~= 0
    error('nObs/nSub is not an integer. Check phase_col or data alignment.');
end
trials_per_sub = round(trials_per_sub);

% Concatenate all subjects: [nChan × nFeat × nObs]
data_all = cat(3, data_cell{:});
if size(data_all,3) ~= nObs
    error('Concatenated data has %d obs but X expects %d obs', size(data_all,3), nObs);
end

% Vectorize features: Y [nObs × (nChan*nFeat)]
Y = reshape(permute(data_all, [3 1 2]), nObs, nChan*nFeat);
nF = size(Y,2);

% ------------------------------
% Observed GLM (OLS)
% ------------------------------
if show_prog
    fprintf('Running observed GLM (vectorized OLS)...\n');
end

XtX_inv = pinv(X' * X);
rankX   = rank(X);
dof     = max(nObs - rankX, 1);

B    = XtX_inv * (X' * Y);       % [nPred × nF]
Res  = Y - X * B;                % [nObs × nF]
mse  = sum(Res.^2, 1) / dof;     % [1 × nF]
SE   = sqrt(bsxfun(@times, diag(XtX_inv), mse)); % [nPred × nF]
Tobs = B ./ SE;                  % [nPred × nF]

% reshape observed outputs
betas_obs = permute(reshape(B',    nChan, nFeat, nPredictors), [1 2 3]);
tvals_obs = permute(reshape(Tobs', nChan, nFeat, nPredictors), [1 2 3]);

% ------------------------------
% Permutations: build full H0 of t-values
% ------------------------------
if show_prog
    fprintf('Running %d permutations...\n', nPerm);
end

tvals_H0 = nan(nChan, nFeat, nPredictors, nPerm, 'like', data_all);

parfor iPerm = 1:nPerm

    % 1) permute phase labels within subject
    phase_perm = phase_col;
    for s = 1:nSub
        idx = ((s-1)*trials_per_sub + 1):(s*trials_per_sub);
        phase_perm(idx) = phase_perm(idx(randperm(trials_per_sub)));
    end

    % 2) build permuted design
    X_perm = build_permuted_X_simple(phase_perm, duration_col, group_col);

    % 3) permuted OLS fit
    XtX_inv_perm = pinv(X_perm' * X_perm);
    dofp         = max(nObs - rank(X_perm), 1);

    Bp   = XtX_inv_perm * (X_perm' * Y);
    Resp = Y - X_perm * Bp;

    msep = sum(Resp.^2, 1) / dofp;
    SEp  = sqrt(bsxfun(@times, diag(XtX_inv_perm), msep));
    Tp   = Bp ./ SEp;

    % store null t-values
    tvals_H0(:,:,:,iPerm) = permute(reshape(Tp', nChan, nFeat, nPredictors), [1 2 3]);

    if show_prog && mod(iPerm, 50) == 0
        fprintf('  permutation %d/%d\n', iPerm, nPerm);
    end
end

% ------------------------------
% Empirical permutation p-values for observed effects (two-tailed)
% ------------------------------
% p = (1 + sum(|Tperm| >= |Tobs|)) / (nPerm + 1)
if show_prog
    fprintf('Computing empirical permutation p-values (two-tailed, uncorrected)...\n');
end

Tobs_abs = abs(Tobs); % [nPred x nF]

Th0 = reshape(permute(tvals_H0, [3 1 2 4]), nPredictors, nF, nPerm); % [nPred x nF x nPerm]
Th0_abs = abs(Th0);

p_emp = nan(nPredictors, nF, 'like', data_all);
for j = 1:nPredictors
    Hj = squeeze(Th0_abs(j,:,:));   % [nF x nPerm]
    oj = Tobs_abs(j,:).';           % [nF x 1]
    cnt = sum(Hj >= oj, 2);         % [nF x 1]
    p_emp(j,:) = (1 + cnt(:)') / (nPerm + 1);
end

pvals_obs = permute(reshape(p_emp', nChan, nFeat, nPredictors), [1 2 3]);

% keep signature
pvals_H0 = [];

if show_prog
    disp('Done: permutation GLM (full H0 t-values retained; no tcdf outputs).');
end
end

%% Helper
function Xp = build_permuted_X_simple(phase_perm, duration_col, group_col)
phase_perm = phase_perm(:);

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
