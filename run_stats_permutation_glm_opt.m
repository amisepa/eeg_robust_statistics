function [betas_obs, tvals_obs, tvals_H0, w_obs] = run_stats_permutation_glm_opt_limo_arch( ...
    Y_all, X, condition_col, nSub, nPerm, varargin)
% RUN_STATS_PERMUTATION_GLM_OPT_LIMO_ARCH
%
% LIMO-style architecture:
%   - One observed GLM fit (OLS / IRLS / WLS)
%   - Permutation-based inference (within-subject label shuffling)
%   - No parametric tcdf/fcdf p-values
%   - Full H0 distributions returned for downstream correction
%
% INPUTS
%   Y_all         [nChan x nTime x nObs] trial-level data
%   X             [nObs x nP] design matrix
%   condition_col [nObs x 1] binary condition labels
%   nSub          number of subjects
%   nPerm         number of permutations
%
% OPTIONS
%   'Subjects'    subj_idx [nObs x 1]
%   'Method'      'OLS' | 'IRLS' | 'WLS'
%   'WeightType'  'PCP' | 'Huber' | 'Tukey' (for WLS)
%   'Progress'    true/false

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
end
