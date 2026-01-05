function [betas_obs, tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation_glm( ...
    Y_all, X, condition_col, nSub, nPerm, varargin)
% RUN_STATS_PERMUTATION_GLM
% Vectorized OLS GLM with within-subject permutations and maxT (t-max) FWE correction.
%
% Fixes versus original:
% 1) Permutation inference is empirical (no tcdf-based p-values inside permutations)
% 2) Multiple-comparisons control via maxT across all (chan x feat) for each predictor
% 3) Supports >=2 observations per subject via within-subject shuffling of condition labels
%
% INPUTS
%   Y_all         [nChan x nFeat x nObs]
%   X             [nObs x nP] canonical: [1, condition, group, condition.*group, duration(optional), ...]
%   condition_col [nObs x 1] 0/1 coding for condition (must match X(:,2) in observed design)
%   nSub          number of subjects
%   nPerm         number of permutations
%
% Name-value:
%   'Group'       group_col vector length nObs (required if nP>=3)
%   'Subjects'    subj_idx vector length nObs (optional; inferred if evenly balanced)
%   'Duration'    duration_col vector length nObs (required if nP>=5)
%   'Progress'    true/false (default true)
%   'StoreH0'     true/false store full tvals_H0 (default true, can be memory heavy)
%
% OUTPUTS
%   betas_obs  [nChan x nFeat x nP]
%   tvals      [nChan x nFeat x nP]
%   pvals      [nChan x nFeat x nP] FWE-corrected p-values via maxT per predictor
%   tvals_H0   [nChan x nFeat x nP x nPerm] (optional; empty if StoreH0=false)
%   pvals_H0   returned as [] (kept for backward compatibility; empirical p-values are in pvals)

% ---------------- options ----------------
p = inputParser;
addParameter(p, 'Group', []);
addParameter(p, 'Subjects', []);
addParameter(p, 'Duration', []);
addParameter(p, 'Progress', true);
addParameter(p, 'StoreH0', true);
parse(p, varargin{:});

group_col    = p.Results.Group;
subj_idx     = p.Results.Subjects;
duration_col = p.Results.Duration;
show_prog    = p.Results.Progress;
storeH0      = p.Results.StoreH0;

% ---------------- checks ----------------
[nChan, nFeat, nObs] = size(Y_all);
nP = size(X,2);

if numel(condition_col) ~= nObs
    error('condition_col must be length nObs');
end
condition_col = condition_col(:);

if ~all(ismember([0 1], unique(condition_col)))
    error('condition_col must contain both 0 and 1');
end

if nP >= 3 && isempty(group_col)
    error('X has a group column but Group was not provided');
end
if nP >= 5 && isempty(duration_col)
    error('X has a duration column but Duration was not provided');
end

if ~isempty(group_col),    group_col    = group_col(:);    end
if ~isempty(duration_col), duration_col = duration_col(:); end

% subject index
if isempty(subj_idx)
    trials_per_sub = nObs / nSub;
    if mod(trials_per_sub,1) ~= 0
        error('Observations are not evenly divisible by nSub. Provide Subjects index.');
    end
    subj_idx = repelem((1:nSub)', trials_per_sub);
else
    subj_idx = subj_idx(:);
    if numel(subj_idx) ~= nObs
        error('Subjects vector must be length nObs');
    end
    u = unique(subj_idx);
    if numel(u) ~= nSub
        error('Subjects index does not contain nSub unique subjects');
    end
    counts = arrayfun(@(s) sum(subj_idx==s), u);
    if any(counts ~= counts(1))
        error('Each subject must have the same number of observations for within-subject permutations');
    end
    trials_per_sub = counts(1);
end

if trials_per_sub < 2
    error('Need at least 2 observations per subject for within-subject contrast');
end

% precompute subject row lists
sub_rows = arrayfun(@(s) find(subj_idx==s), (1:nSub)', 'uni', false);

% ---------------- vectorize data ----------------
% Y: [nObs x (nChan*nFeat)]
Y = reshape(permute(Y_all, [3 1 2]), nObs, nChan*nFeat);
nF = size(Y,2);

% ---------------- observed OLS fit ----------------
if show_prog, fprintf('Observed fit: vectorized OLS\n'); end

XtX     = X' * X;
XtX_inv = pinv(XtX);
rankX   = rank(X);
dof     = max(nObs - rankX, 1);

B    = XtX_inv * (X' * Y);     % [nP x nF]
Yhat = X * B;
Res  = Y - Yhat;
mse  = sum(Res.^2, 1) / dof;   % [1 x nF]
SE   = sqrt(bsxfun(@times, diag(XtX_inv), mse)); % [nP x nF]
T    = B ./ SE;                % [nP x nF]

% reshape observed outputs
betas_obs = permute(reshape(B', nChan, nFeat, nP), [1 2 3]);
tvals     = permute(reshape(T', nChan, nFeat, nP), [1 2 3]);

% ---------------- permutation null + maxT ----------------
if show_prog, fprintf('Running %d within-subject permutations...\n', nPerm); end

if storeH0
    tvals_H0 = nan(nChan, nFeat, nP, nPerm, 'like', Y_all);
else
    tvals_H0 = [];
end

% maxT per predictor and permutation
maxT_H0 = nan(nP, nPerm);

% columns 6+ are copied verbatim; if you include condition-dependent terms there,
% you must rebuild them yourself inside the permutation loop.
X_tail = [];
if nP > 5
    X_tail = X(:,6:end);
end

parfor iPerm = 1:nPerm
    % 1) permute condition labels within subject (preserves counts)
    cond_perm = condition_col;
    for s = 1:nSub
        rows = sub_rows{s};
        cond_perm(rows) = cond_perm(rows(randperm(numel(rows))));
    end

    % 2) rebuild permuted design matrix in canonical order
    Xp = zeros(nObs, nP);
    Xp(:,1) = 1;
    if nP >= 2, Xp(:,2) = cond_perm; end
    if nP >= 3, Xp(:,3) = double(group_col); end
    if nP >= 4, Xp(:,4) = double(cond_perm) .* double(group_col); end
    if nP >= 5, Xp(:,5) = duration_col; end
    if ~isempty(X_tail), Xp(:,6:end) = X_tail; end

    % 3) vectorized OLS on permuted design
    XtXp     = Xp' * Xp;
    XtXp_inv = pinv(XtXp);
    dof_p    = max(nObs - rank(Xp), 1);

    Bp   = XtXp_inv * (Xp' * Y);         % [nP x nF]
    Resp = Y - Xp * Bp;
    msep = sum(Resp.^2, 1) / dof_p;      % [1 x nF]
    SEp  = sqrt(bsxfun(@times, diag(XtXp_inv), msep));
    Tp   = Bp ./ SEp;                    % [nP x nF]

    if storeH0
        tvals_H0(:,:,:,iPerm) = permute(reshape(Tp', nChan, nFeat, nP), [1 2 3]);
    end

    % maxT across all features for each predictor
    maxT_H0(:, iPerm) = max(abs(Tp), [], 2);
end

% ---------------- empirical maxT FWE p-values ----------------
% pvals(chan,feat,p) = mean( maxT_H0(p,:) >= abs(T_obs(p,chan,feat)) )
pvals = nan(nChan, nFeat, nP, 'like', Y_all);

Tobs = reshape(permute(tvals, [3 1 2]), nP, nF); % [nP x nF]
for j = 1:nP
    thr = maxT_H0(j, :); % [1 x nPerm]
    p_j = mean(thr(:) >= abs(Tobs(j,:)), 1); % [1 x nF]
    pvals(:,:,j) = reshape(p_j, nChan, nFeat);
end

% keep output for backward compatibility
pvals_H0 = [];

if show_prog, disp('Done: permutation GLM with maxT FWE correction.'); end
end
