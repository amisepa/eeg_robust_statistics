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

[nChan, nTimes, ~] = size(data1);
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
        nSub = size(data1,3);
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
[hObs, qObs] = make_parfor_waitbar(nChan, 'Observed stats');
parfor iChan = 1:nChan
    for t = 1:nTimes
        x1 = squeeze(data1(iChan, t, :));
        x2 = squeeze(data2(iChan, t, :));

        switch lower(dpt)
            case 'dpt'
                % remove subjects with NaN in either vector
                keep = ~(isnan(x1) | isnan(x2));
                x1c = x1(keep); x2c = x2(keep);
            case 'idpt'
                x1c = x1(~isnan(x1));
                x2c = x2(~isnan(x2));
        end

        if strcmpi(method,'mean')
            [tv, pv] = t_mean_safe(x1c, x2c, dpt);
        elseif strcmpi(method,'trimmed mean')
            [tv, pv] = t_yuen_safe(x1c, x2c, dpt, trim);
        else
            error("Unknown method: choose 'mean' or 'trimmed mean'")
        end

        tvals(iChan, t) = tv;
        pvals(iChan, t) = pv;
    end
    send(qObs, 1);                           % notify progress
end
close(hObs);

% Synchronized permutations
swap_mask = [];
idx_perm  = [];
if strcmpi(dpt, 'dpt')
    swap_mask = rand(size(data1,3), nPerm) > 0.5;   % [nSub x nPerm]
else
    pooledN  = size(data1,3) + size(data2,3);
    idx_perm = zeros(nPerm, pooledN, 'uint32');
    for perm = 1:nPerm
        idx_perm(perm,:) = uint32(randperm(pooledN));
    end
end

% Null distribution
NpermTot = nChan * nPerm;                    % total iterations: channel × permutation
[hPerm, qPerm] = make_parfor_waitbar(NpermTot, 'Permutation nulls');
parfor iChan = 1:nChan
    A = data1(iChan,:,:);   % [1 x nTimes x n]
    B = data2(iChan,:,:);

    n1_local = size(A,3);
    n2_local = size(B,3);

    for perm = 1:nPerm
        if strcmpi(dpt, 'dpt')
            permA = A; permB = B;
            sw = swap_mask(:,perm);
            if any(sw)
                permA(:,:,sw) = B(:,:,sw);
                permB(:,:,sw) = A(:,:,sw);
            end
        else
            X = cat(3, A, B);
            idx = idx_perm(perm,:);
            permA = X(:,:,idx(1:n1_local));
            permB = X(:,:,idx(n1_local+1:n1_local+n2_local));
        end

        if strcmpi(method,'mean')
            % vectorized across time with NaN handling
            X1 = squeeze(permA(1,:,:));  % [n x nTimes] with possible NaN
            X2 = squeeze(permB(1,:,:));
            [tv, pv] = t_mean_matrix_safe(X1, X2, dpt);
        else
            % trimmed mean per time with NaN handling
            tv = nan(1, nTimes);
            pv = nan(1, nTimes);
            for t = 1:nTimes
                xx1 = squeeze(permA(1,t,:));
                xx2 = squeeze(permB(1,t,:));
                if strcmpi(dpt,'dpt')
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

        send(qPerm, 1);                      % notify progress
    end
end
close(hPerm);

end

% ===================== helpers =====================

function [tstat, p] = t_mean_safe(x1, x2, dpt)
% Two sided t test, data1 minus data2, robust to NaN and zero variance
eps_se = 1e-12;
if strcmpi(dpt,'dpt')
    n = numel(x1);
    if n < 2
        tstat = 0; p = 1; return
    end
    d  = x1 - x2;
    md = mean(d,'omitnan');
    sd = std(d,0,'omitnan');
    se = sd ./ sqrt(n);
    if ~isfinite(se) || se < eps_se
        if abs(md) < eps_se, tstat = 0; p = 1;
        else, tstat = sign(md)*Inf; p = 1/n; end
        return
    end
    tstat = md ./ se;
    df = max(n - 1, 1);
    p = 2 * (1 - tcdf(abs(tstat), df));
else
    nA = numel(x1); nB = numel(x2);
    if nA < 2 || nB < 2
        tstat = 0; p = 1; return
    end
    mA = mean(x1,'omitnan');  mB = mean(x2,'omitnan');
    vA = var(x1,0,'omitnan'); vB = var(x2,0,'omitnan');
    se2 = vA./nA + vB./nB; se = sqrt(se2);
    if ~isfinite(se) || se < eps_se
        md = mA - mB;
        if abs(md) < eps_se, tstat = 0; p = 1;
        else, tstat = sign(md)*Inf; p = 1/(nA+nB); end
        return
    end
    tstat = (mA - mB) ./ se;
    % Welch Satterthwaite df, guard zeros
    num = se2.^2;
    den = (vA.^2)/(nA^2*(nA-1) + (nA==1)) + (vB.^2)/(nB^2*(nB-1) + (nB==1));
    df = num./max(den, eps_se);
    df = max(df, 1);
    p = 2 * (1 - tcdf(abs(tstat), df));
end
end

function [tvec, pvec] = t_mean_matrix_safe(X1, X2, dpt)
% Columnwise t for matrices [n x nTimes], NaN safe, data1 minus data2
nTimes = size(X1,2);
tvec = nan(1,nTimes); pvec = nan(1,nTimes);
for t = 1:nTimes
    x1 = X1(:,t); x2 = X2(:,t);
    [tvec(t), pvec(t)] = t_mean_safe(x1(~isnan(x1)), x2(~isnan(x2)), dpt);
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