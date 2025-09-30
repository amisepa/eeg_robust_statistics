function [betas_obs, tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation_glm( ...
    Y_all, X, condition_col, nSub, nPerm, varargin)
% RUN_STATS_PERMUTATION_GLM
% Y_all         [nChan x nFeat x nObs]
% X             [nObs x nP] in canonical order:
%               [1, condition, group, condition.*group, duration(optional)]
% condition_col length nObs, 0/1, pattern [0;1] per subject
% nSub          number of subjects
% nPerm         number of permutations
%
% Name-value:
%   'Group'       group_col vector (0 novice, 1 expert)
%   'Subjects'    subj_idx vector length nObs (optional; inferred)
%   'Duration'    duration_col vector length nObs (optional)
%   'Progress'    true/false print progress (default true)

% --------------- options ---------------
p = inputParser;
addParameter(p, 'Group', []);
addParameter(p, 'Subjects', []);
addParameter(p, 'Duration', []);
addParameter(p, 'Progress', true);
parse(p, varargin{:});
group_col    = p.Results.Group;
subj_idx     = p.Results.Subjects;
duration_col = p.Results.Duration;
show_prog    = p.Results.Progress;

% --------------- checks ---------------
[nChan, nFreq, nObs] = size(Y_all);
nP = size(X,2);

if numel(condition_col) ~= nObs
    error('condition_col must be length nObs');
end
condition_col = condition_col(:);

% check canonical layout requirements
if nP >= 3 && isempty(group_col)
    error('X has a group column but Group was not provided');
end
if nP >= 5 && isempty(duration_col)
    error('X has a duration column but Duration was not provided');
end

if ~isempty(group_col),    group_col    = group_col(:);    end
if ~isempty(duration_col), duration_col = duration_col(:); end

if isempty(subj_idx)
    trials_per_sub = nObs / nSub;
    if mod(trials_per_sub,1) ~= 0
        error('Observations are not evenly divisible by nSub. Provide Subjects index.');
    end
    subj_idx = repelem((1:nSub)', trials_per_sub);
else
    if numel(subj_idx) ~= nObs
        error('Subjects vector must be length nObs');
    end
    u = unique(subj_idx(:));
    counts = arrayfun(@(s) sum(subj_idx==s), u);
    if any(counts ~= counts(1))
        error('Each subject must have the same number of observations for paired permutations');
    end
    trials_per_sub = counts(1);
end
if trials_per_sub < 2
    error('Need at least 2 observations per subject for within-subject contrast');
end
if ~all(ismember([0 1], unique(condition_col)))
    error('condition_col must contain both 0 and 1');
end

% --------------- observed fit ---------------
if show_prog, fprintf('Running statistics on observed Betas...\n'); end
betas_obs = nan(nChan, nFreq, nP);
tvals     = nan(nChan, nFreq, nP);
pvals     = nan(nChan, nFreq, nP);

XtX_inv = pinv(X' * X);
rankX   = rank(X);
dof_obs = max(nObs - rankX, 1);

parfor iChan = 1:nChan
    betas_ch = nan(nFreq, nP);
    tvals_ch = nan(nFreq, nP);
    pvals_ch = nan(nFreq, nP);
    for iFreq = 1:nFreq
        y = squeeze(Y_all(iChan, iFreq, :));
        b = X \ y;
        resid = y - X*b;

        mse = sum(resid.^2) / dof_obs;
        se  = sqrt(diag(XtX_inv * mse));
        se(se==0) = Inf;

        t_here = (b ./ se);
        p_here = 2 * (1 - tcdf(abs(t_here), dof_obs));

        betas_ch(iFreq,:) = b.';
        tvals_ch(iFreq,:) = t_here.';
        pvals_ch(iFreq,:) = p_here.';
    end
    betas_obs(iChan,:,:) = betas_ch;
    tvals(iChan,:,:)     = tvals_ch;
    pvals(iChan,:,:)     = pvals_ch;
end

% --------------- permutations ---------------
if show_prog, fprintf('Running permutation statistics (H0)...\n'); end
tvals_H0 = nan(nChan, nFreq, nP, nPerm);
pvals_H0 = nan(nChan, nFreq, nP, nPerm);

% precompute subject row lists
sub_rows = arrayfun(@(s) find(subj_idx==s), 1:nSub, 'uni', false);

parfor iPerm = 1:nPerm
    % swap condition within subject
    cond_perm = condition_col;
    for s = 1:nSub
        rows = sub_rows{s};
        if rand > 0.5
            cond_perm(rows) = cond_perm(rows(end:-1:1));
        end
    end

    % rebuild X_perm deterministically in canonical order
    X_perm = zeros(size(X));
    X_perm(:,1) = 1;
    if nP >= 2
        X_perm(:,2) = cond_perm;
    end
    if nP >= 3
        X_perm(:,3) = double(group_col);
    end
    if nP >= 4
        X_perm(:,4) = double(cond_perm) .* double(group_col);
    end
    if nP >= 5
        X_perm(:,5) = duration_col;
    end
    % if user added more columns beyond 5, copy them verbatim
    if nP > 5
        X_perm(:,6:end) = X(:,6:end);
    end

    XtX_inv_perm = pinv(X_perm' * X_perm);
    dof_perm     = max(nObs - rank(X_perm), 1);

    tvals_perm_i = nan(nChan, nFreq, nP);
    pvals_perm_i = nan(nChan, nFreq, nP);

    for iChan2 = 1:nChan
        for iFeat2 = 1:nFreq
            y = squeeze(Y_all(iChan2, iFeat2, :));
            b = X_perm \ y;
            resid = y - X_perm*b;

            mse = sum(resid.^2) / dof_perm;
            se  = sqrt(diag(XtX_inv_perm * mse));
            se(se==0) = Inf;

            t_here = (b ./ se);
            p_here = 2 * (1 - tcdf(abs(t_here), dof_perm));

            tvals_perm_i(iChan2,iFeat2,:) = t_here;
            pvals_perm_i(iChan2,iFeat2,:) = p_here;
        end
    end

    tvals_H0(:,:,:,iPerm) = tvals_perm_i;
    pvals_H0(:,:,:,iPerm) = pvals_perm_i;

    if show_prog
        step = max(1, floor(nPerm/20));
        if mod(iPerm, step) == 0
            fprintf('   - permutation %d/%d\n', iPerm, nPerm);
        end
    end
end

if show_prog, disp('Done with GLM permutation statistics.'); end
end
