function [betas_obs, tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation_glm_opt( ...
    Y_all, X, condition_col, nSub, nPerm, varargin)
% RUN_STATS_PERMUTATION_GLM_OPT
% Vectorized OLS with optional robustness for a paired two-condition GLM and OLS permutations.
%
% INPUTS
%   Y_all         [nChan x nFeat x nObs] observations stacked as [A; B] per subject
%   X             [nObs x nP] canonical order: [1, condition, group, condition.*group, duration(optional), ...]
%   condition_col length nObs, 0 for condition A, 1 for condition B, pattern [0;1] per subject
%   nSub          number of subjects
%   nPerm         number of permutations
%
% Name-value options
%   'Group'        group_col vector (0 novice, 1 expert)
%   'Subjects'     subj_idx vector length nObs (optional; inferred if evenly balanced)
%   'Duration'     duration_col vector length nObs (optional; required if X has that column)
%   'Progress'     true/false to print progress (default true)
%   'RobustMode'   'none' | 'wls' | 'irls'  (default 'none')
%
% Notes
%   1. WLS computes one weight per subject from cross-feature residual energy, then applies it to both rows.
%   2. IRLS refines each feature with Huber weights on the observed fit only. Permutations stay OLS for speed and a clean null.
%   3. Outputs match your original function signature.

% 1) Parse options
p = inputParser;
addParameter(p, 'Group', []);
addParameter(p, 'Subjects', []);
addParameter(p, 'Duration', []);
addParameter(p, 'Progress', true);
addParameter(p, 'RobustMode', 'none');   % 'none' | 'wls' | 'irls'
parse(p, varargin{:});
group_col    = p.Results.Group;
subj_idx     = p.Results.Subjects;
duration_col = p.Results.Duration;
show_prog    = p.Results.Progress;
mode         = lower(p.Results.RobustMode);

% 2) Basic checks and shapes
[nChan, nFreq, nObs] = size(Y_all);
nP = size(X,2);

if numel(condition_col) ~= nObs
    error('condition_col must be length nObs');
end
condition_col = condition_col(:);

if nP >= 3 && isempty(group_col)
    error('X has a group column but Group was not provided');
end
if nP >= 5 && isempty(duration_col)
    error('X has a duration column but Duration was not provided');
end

if ~isempty(group_col),    group_col    = group_col(:);    end
if ~isempty(duration_col), duration_col = duration_col(:); end

% Subjects indexing and pairing checks
if isempty(subj_idx)
    trials_per_sub = nObs / nSub;
    if mod(trials_per_sub,1) ~= 0
        error('Observations are not evenly divisible by nSub. Provide Subjects index.');
    end
    subj_idx = repelem((1:nSub)', trials_per_sub);
else
    if numel(subj_idx) ~= nObs
        error('Subjects must be length nObs');
    end
    u = unique(subj_idx(:));
    counts = arrayfun(@(s) sum(subj_idx==s), u);
    if any(counts ~= counts(1))
        error('Each subject must have the same number of observations');
    end
    trials_per_sub = counts(1);
end

if trials_per_sub < 2
    error('Need at least 2 observations per subject for within-subject designs');
end
if ~all(ismember([0 1], unique(condition_col)))
    error('condition_col must contain 0 and 1');
end

% 3) Vectorize features
nFeat = nChan * nFreq;
Y = reshape(permute(Y_all, [3 1 2]), nObs, nFeat);  % [obs x features]

% Precompute useful quantities
XtX     = X' * X;
XtX_inv = pinv(XtX);
rankX   = rank(X);
dof_obs = max(nObs - rankX, 1);

% 4) Observed fit: OLS or robust
switch mode
    case 'none'
        if show_prog, fprintf('Observed fit: OLS\n'); end
        % OLS
        B    = XtX_inv * (X' * Y);      % [nP x nFeat]
        Yhat = X * B;
        Res  = Y - Yhat;
        mse  = sum(Res.^2, 1) / dof_obs;
        SE   = sqrt(bsxfun(@times, diag(XtX_inv), mse));
        T    = B ./ SE;
        P    = 2 * (1 - tcdf(abs(T), dof_obs));

    case 'wls'
        if show_prog, fprintf('Observed fit: subject-paired WLS\n'); end
        % 4.1 Quick OLS to estimate residual energy
        B0 = XtX_inv * (X' * Y);
        R0 = Y - X * B0;                               % [nObs x nFeat]
        rE = sqrt(median(R0.^2, 2));                   % robust energy per obs
        % Collapse to one weight per subject, assign back to both rows
        rE_sub = accumarray(subj_idx, rE, [], @median);   % [nSub x 1]
        s      = 1.4826 * mad(rE_sub, 1) + eps;
        z      = (rE_sub - median(rE_sub)) ./ s;
        c      = 1.345;                                % Huber constant
        w_sub  = 1 ./ max(1, abs(z)/c);
        w_sub  = max(w_sub, 1e-6);
        W_vec  = w_sub(subj_idx);                      % [nObs x 1], same weight for paired rows

        % 4.2 Weighted OLS
        sqrtW = sqrt(W_vec);
        Xw = bsxfun(@times, X, sqrtW);
        Yw = bsxfun(@times, Y, sqrtW);
        XtXw     = Xw' * Xw;
        XtXw_inv = pinv(XtXw);
        B    = XtXw_inv * (Xw' * Yw);
        Yhat = X * B;
        Res  = Y - Yhat;
        mse  = sum(bsxfun(@times, Res.^2, W_vec), 1) / dof_obs;
        SE   = sqrt(bsxfun(@times, diag(XtXw_inv), mse));
        T    = B ./ SE;
        P    = 2 * (1 - tcdf(abs(T), dof_obs));

    case 'irls'
        if show_prog, fprintf('Observed fit: IRLS with Huber weights\n'); end
        % 4.1 Initialize with OLS
        B = XtX_inv * (X' * Y);
        % 4.2 Feature-wise IRLS
        maxIt = 25; tol = 1e-4; c = 1.345;
        parfor f = 1:nFeat
            b = B(:,f);
            y = Y(:,f);
            for it = 1:maxIt
                r = y - X*b;
                s = 1.4826 * mad(r, 1) + eps;
                u = r ./ s;
                w = 1 ./ max(1, abs(u)/c);          % Huber
                w = max(w, 1e-6);
                % Weighted least squares step using sqrt weights
                sw = sqrt(w);
                Xwf = bsxfun(@times, X, sw);
                ywf = y .* sw;
                bw  = pinv(Xwf' * Xwf) * (Xwf' * ywf);
                if norm(bw - b, 2) < tol * (1 + norm(b,2)), b = bw; break; end
                b = bw;
            end
            B(:,f) = b;
        end
        Yhat = X * B;
        Res  = Y - Yhat;
        mse  = sum(Res.^2, 1) / dof_obs;
        SE   = sqrt(bsxfun(@times, diag(XtX_inv), mse));  % plug-in with unweighted XtX_inv
        T    = B ./ SE;
        P    = 2 * (1 - tcdf(abs(T), dof_obs));

    otherwise
        error('RobustMode must be ''none'', ''wls'', or ''irls''');
end

% 5) Reshape observed outputs
betas_obs = permute(reshape(B', nChan, nFreq, nP), [1 2 3]);
tvals     = permute(reshape(T', nChan, nFreq, nP), [1 2 3]);
pvals     = permute(reshape(P', nChan, nFreq, nP), [1 2 3]);

% 6) Permutations under H0 (OLS, vectorized)
if show_prog, fprintf('Running permutation statistics (H0)...\n'); end
tvals_H0 = nan(nChan, nFreq, nP, nPerm, 'like', Y_all);
pvals_H0 = nan(nChan, nFreq, nP, nPerm, 'like', Y_all);

% Precompute subject row indices for paired swap
sub_rows = arrayfun(@(s) find(subj_idx==s), 1:nSub, 'uni', false);

parfor iPerm = 1:nPerm
    % 6.1 Paired swap of condition within subject
    cond_perm = condition_col;
    for s = 1:nSub
        rows = sub_rows{s};
        if rand > 0.5
            cond_perm(rows) = cond_perm(rows(end:-1:1));
        end
    end

    % 6.2 Rebuild permuted X in canonical order
    Xp = zeros(size(X));
    Xp(:,1) = 1;
    if nP >= 2, Xp(:,2) = cond_perm; end
    if nP >= 3, Xp(:,3) = double(group_col); end
    if nP >= 4, Xp(:,4) = double(cond_perm) .* double(group_col); end
    if nP >= 5, Xp(:,5) = duration_col; end
    if nP > 5,  Xp(:,6:end) = X(:,6:end); end

    % 6.3 Vectorized OLS on permuted design
    XtXp     = Xp' * Xp;
    XtXp_inv = pinv(XtXp);
    dof_p    = max(nObs - rank(Xp), 1);

    Bp   = XtXp_inv * (Xp' * Y);
    Resp = Y - Xp * Bp;
    msep = sum(Resp.^2, 1) / dof_p;
    SEp  = sqrt(bsxfun(@times, diag(XtXp_inv), msep));
    Tp   = Bp ./ SEp;
    Pp   = 2 * (1 - tcdf(abs(Tp), dof_p));

    tvals_H0(:,:,:,iPerm) = permute(reshape(Tp', nChan, nFreq, nP), [1 2 3]);
    pvals_H0(:,:,:,iPerm) = permute(reshape(Pp', nChan, nFreq, nP), [1 2 3]);

    if show_prog
        step = max(1, floor(nPerm/20));
        if mod(iPerm, step) == 0
            fprintf('   - permutation %d/%d\n', iPerm, nPerm);
        end
    end
end

if show_prog, disp('Done with GLM permutation statistics.'); end
end
