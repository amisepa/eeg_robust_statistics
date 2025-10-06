function [tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation(data1, data2, nPerm, method, dpt)
% RUN_STATS_PERMUTATION
% Pairwise permutation testing with synchronized permutations across channels.
% Direction is data1 minus data2. Robust to NaNs and zero variance.
%
% Inputs
%   data1   [nChan x nTimes x nSub] for dpt, or [nChan x nTimes x n1] for idpt
%   data2   same shape as data1 for dpt, or [nChan x nTimes x n2] for idpt
%   nPerm   number of permutations. Default 1000
%   method  'mean' or 'trimmed mean'. Default 'trimmed mean' (20 percent)
%   dpt     'dpt' for paired design or 'idpt' for independent groups
%
% Outputs
%   tvals       [nChan x nTimes]
%   pvals       [nChan x nTimes]
%   tvals_H0    [nChan x nTimes x nPerm]
%   pvals_H0    [nChan x nTimes x nPerm]

% Defaults
if nargin < 3 || isempty(nPerm), nPerm = 1000; end
if nargin < 4 || isempty(method), method = 'trimmed mean'; end
if nargin < 5 || isempty(dpt), error("Specify 'dpt' or 'idpt'"); end
trim = 0.20;

% Validate dimensions and safely promote 2D to 3D
if ndims(data1) == 2 && ndims(data2) == 2
    % Interpret as [nChan x nSub]  ->  [nChan x 1 x nSub]
    [nChan1, nSub1] = size(data1);
    [nChan2, nSub2] = size(data2);
    if nChan1 ~= nChan2
        error('data1 and data2 must have the same number of channels.')
    end
    data1 = reshape(data1, nChan1, 1, nSub1);
    data2 = reshape(data2, nChan2, 1, nSub2);
elseif ndims(data1) == 3 && ndims(data2) == 3
    % ok as is
else
    error('Inputs must be 2D [nChan x nSub] or 3D [nChan x nTimes x nSub].')
end

[nChan, nTimes, nSub] = size(data1);
if size(data2,1) ~= nChan || size(data2,2) ~= nTimes
    error('data1 and data2 must match in [nChan x nTimes].')
end

tvals    = nan(nChan, nTimes);
pvals    = nan(nChan, nTimes);
tvals_H0 = nan(nChan, nTimes, nPerm);
pvals_H0 = nan(nChan, nTimes, nPerm);

% Sizes
switch lower(dpt)
    case 'dpt'
        if size(data2,3) ~= nSub
            error('Paired case requires the same number of subjects in data1 and data2.')
        end
    case 'idpt'
        n1 = size(data1,3);
        n2 = size(data2,3);
    otherwise
        error("Unknown dependency type: use 'dpt' or 'idpt'")
end

% Observed stats
method_is_mean = strcmpi(method,'mean');
is_dpt = strcmpi(dpt,'dpt');
[hObs, qObs] = make_parfor_waitbar(nChan, 'Observed stats');
if method_is_mean
    parfor iChan = 1:nChan
        % Build [n x nTimes] once per channel, no squeeze in inner loop
        X1 = permute(data1(iChan,:,:), [3 2 1]);  % [n x nTimes]
        X2 = permute(data2(iChan,:,:), [3 2 1]);  % [n x nTimes]
        [tv, pv] = t_mean_matrix_fast(X1, X2, dpt);
        tvals(iChan,:) = tv;
        pvals(iChan,:) = pv;
        send(qObs, 1);
    end
else
    parfor iChan = 1:nChan
        for t = 1:nTimes
            % Always initialize temporaries for parfor analysis
            x1c = []; x2c = [];
            tv  = 0;  pv  = 1;

            x1 = data1(iChan, t, :); x1 = x1(:);
            x2 = data2(iChan, t, :); x2 = x2(:);

            if is_dpt
                keep = ~(isnan(x1) | isnan(x2));
                if any(keep)
                    x1c = x1(keep); x2c = x2(keep);
                end
            else
                if any(~isnan(x1)), x1c = x1(~isnan(x1)); end
                if any(~isnan(x2)), x2c = x2(~isnan(x2)); end
            end

            if ~isempty(x1c) && ~isempty(x2c)
                [tv, pv] = t_yuen_safe(x1c, x2c, dpt, trim);
            end

            tvals(iChan, t) = tv;
            pvals(iChan, t) = pv;
        end
        send(qObs, 1);
    end
end
close(hObs);

% Synchronized permutations
swap_mask = [];
idx_perm  = [];
if is_dpt
    swap_mask = rand(size(data1,3), nPerm) > 0.5;   % [nSub x nPerm]
else
    pooledN  = size(data1,3) + size(data2,3);
    idx_perm = zeros(nPerm, pooledN, 'uint32');
    for perm = 1:nPerm
        idx_perm(perm,:) = uint32(randperm(pooledN));
    end
end

% ================= Null distribution =================
NpermTot = nChan * nPerm;
[hPerm, qPerm] = make_parfor_waitbar(NpermTot, 'Permutation nulls');

% Precompute sizes once
n1 = size(data1, 3);
n2 = size(data2, 3);
pooledN = n1 + n2;

% Build permutation helpers
if is_dpt
    % [n1 x nPerm] logical swap mask
    swap_mask = rand(n1, nPerm) > 0.5;
else
    % [nPerm x pooledN] permutation indices as doubles
    idx_perm = zeros(nPerm, pooledN);
    for k = 1:nPerm
        idx_perm(k,:) = randperm(pooledN);
    end
end


% parfor over channels
parfor iChan = 1:nChan
    A = data1(iChan,:,:);  % [1 x nTimes x n1]
    B = data2(iChan,:,:);  % [1 x nTimes x n2]

    for perm = 1:nPerm
        if is_dpt
            % within-subject swaps
            permA = A; permB = B;
            sw = swap_mask(:,perm);            % [n1 x 1] logical
            if any(sw)
                permA(:,:,sw) = B(:,:,sw);
                permB(:,:,sw) = A(:,:,sw);
            end
        else
            % independent groups: reindex pooled
            X = cat(3, A, B);                  % [1 x nTimes x pooledN]
            idx = idx_perm(perm,:);            % 1 x pooledN
            permA = X(:,:,idx(1:n1));          % [1 x nTimes x n1]
            permB = X(:,:,idx(n1+1:pooledN));  % [1 x nTimes x n2]
        end

        if method_is_mean
            X1 = permute(permA, [3 2 1]);      % [n x nTimes]
            X2 = permute(permB, [3 2 1]);      % [n x nTimes]
            [tv, pv] = t_mean_matrix_fast(X1, X2, dpt);
        else
            tv = nan(1, nTimes);
            pv = nan(1, nTimes);
            for t = 1:nTimes
                xx1 = permA(1,t,:); xx1 = xx1(:);
                xx2 = permB(1,t,:); xx2 = xx2(:);
                if is_dpt
                    keep = ~(isnan(xx1) | isnan(xx2));
                    xx1c = xx1(keep); xx2c = xx2(keep);
                else
                    xx1c = xx1(~isnan(xx1));
                    xx2c = xx2(~isnan(xx2));
                end
                [tv(t), pv(t)] = t_yuen_safe(xx1c, xx2c, dpt, trim);
            end
        end

        tvals_H0(iChan,:,perm) = tv;
        pvals_H0(iChan,:,perm) = pv;

        send(qPerm, 1);                        % one tick per inner iter
    end
end

close(hPerm);

end

% ===================== helpers =====================

function [tvec, pvec] = t_mean_matrix_fast(X1, X2, dpt)
% X1, X2 are [n x nTimes], NaN allowed
eps_se = 1e-12;
nTimes = size(X1,2);
if strcmpi(dpt,'dpt')
    mask = ~isnan(X1) & ~isnan(X2);       % [n x nTimes]
    D = X1 - X2; D(~mask) = NaN;
    n  = sum(mask,1);                      % 1 x nTimes
    md = mean(D,1,'omitmissing');
    sd = std(D,0,1,'omitmissing');
    se = sd ./ sqrt(max(n,1));

    % guards
    tiny = ~isfinite(se) | (se < eps_se);
    t   = md ./ se;
    % substitute Inf when mean not tiny but se tiny
    fix = tiny & isfinite(md) & (abs(md) >= eps_se);
    t(fix) = sign(md(fix)).*Inf;
    t(tiny & ~fix) = 0;

    df = max(n - 1, 1);
    p = 2 * (1 - tcdf(abs(t), df));

    tvec = t; pvec = p;
else
    nA = sum(~isnan(X1),1); nB = sum(~isnan(X2),1);
    mA = mean(X1,1,'omitmissing');     
    mB = mean(X2,1,'omitmissing');
    vA = var(X1,0,1,'omitmissing');    
    vB = var(X2,0,1,'omitmissing');

    se2 = vA./max(nA,1) + vB./max(nB,1);
    se  = sqrt(se2);
    md  = mA - mB;

    tiny = ~isfinite(se) | (se < eps_se);
    t    = md ./ se;
    fix  = tiny & isfinite(md) & (abs(md) >= eps_se);
    t(fix) = sign(md(fix)).*Inf;
    t(tiny & ~fix) = 0;

    num = se2.^2;
    den = (vA.^2) ./ (max(nA.^2 .* max(nA-1,1), 1)) + ...
        (vB.^2) ./ (max(nB.^2 .* max(nB-1,1), 1));
    df  = num ./ max(den, eps_se);
    df  = max(df, 1);

    p = 2 * (1 - tcdf(abs(t), df));

    tvec = t; pvec = p;
end
end


function [tstat, p] = t_yuen_safe(x1, x2, dpt, gamma)
% Two sided Yuen test with safety for small N and NaNs, data1 minus data2
eps_se = 1e-12;
if strcmpi(dpt,'dpt')
    d = x1 - x2;
    [tm, wv, h] = trim_winsor_stats(d, gamma);
    if h <= 1 || ~isfinite(wv)
        tstat = 0; p = 1; return
    end
    se = sqrt( (numel(d) * wv) / (h * (h - 1)) );
    if ~isfinite(se) || se < eps_se
        if abs(tm) < eps_se, tstat = 0; p = 1; else, tstat = sign(tm)*Inf; p = 1/numel(d); end
        return
    end
    tstat = tm ./ se;
    df = max(h - 1, 1);
    p = 2 * (1 - tcdf(abs(tstat), df));
else
    [tmx, wvx, hx] = trim_winsor_stats(x1, gamma);
    [tmy, wvy, hy] = trim_winsor_stats(x2, gamma);
    if hx <= 1 || hy <= 1 || ~isfinite(wvx) || ~isfinite(wvy)
        tstat = 0; p = 1; return
    end
    vx  = (numel(x1) * wvx) / (hx * (hx - 1));
    vy  = (numel(x2) * wvy) / (hy * (hy - 1));
    se2 = vx + vy;
    if ~isfinite(se2) || se2 < eps_se
        md = tmx - tmy;
        if abs(md) < eps_se, tstat = 0; p = 1; else, tstat = sign(md)*Inf; p = 1/(numel(x1)+numel(x2)); end
        return
    end
    tstat = (tmx - tmy) ./ sqrt(se2);
    df = (se2.^2) / ( (vx.^2)/(hx - 1) + (vy.^2)/(hy - 1) );
    df = max(df, 1);
    p = 2 * (1 - tcdf(abs(tstat), df));
end
end

function [tm, wv, h] = trim_winsor_stats(v, gamma)
% Trim to central (1-2*gamma) and compute winsorized variance
v = v(~isnan(v));
n = numel(v);
if n == 0
    tm = NaN; wv = NaN; h = 0; return
end
v = sort(v);
g = floor(gamma * n);
h = n - 2*g;
if h <= 0
    tm = mean(v); wv = var(v,0); h = 0; return
end
core = v(g+1:n-g);
tm = mean(core);
if g > 0
    w = v;
    w(1:g)     = v(g+1);
    w(n-g+1:n) = v(n-g);
else
    w = v;
end
wv = sum( (w - mean(w)).^2 ) / max(n - 1, 1);
if isnan(wv), wv = 0; end
end

% ===== put this helper as a SUBFUNCTION at the end of the file =====
function [h, q] = make_parfor_waitbar(N, titleStr)
% Waitbar that updates from parfor workers via a DataQueue
h = waitbar(0, sprintf('%s 0%%', titleStr));
q = parallel.pool.DataQueue;
cnt = 0;  % lives on client
afterEach(q, @update);
    function update(~)
        cnt = cnt + 1;
        frac = min(cnt / N, 1);
        if isvalid(h)
            waitbar(frac, h, sprintf('%s %2.0f%%', titleStr, 100*frac));
        end
    end
end