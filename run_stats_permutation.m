function [tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation(data1, data2, nPerm, method, dpt)
% RUN_STATS_PERMUTATION
% Pairwise permutation testing with synchronized permutations across channels.
% Direction is always data1 minus data2.
%
% Syntax
%   [tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation(data1, data2, nPerm, method, dpt)
%
% Description
%   Computes observed pointwise test statistics and permutation null distributions
%   for a pairwise comparison suitable for TFCE or t max familywise error control.
%   For dependent data it uses subject wise swaps. For independent data it uses
%   pooled label shuffles. At each permutation, the same randomization is applied
%   to all channels which preserves dependence required by TFCE and t max.
%
% Direction
%   All t values are for data1 minus data2.
%   p values are two sided.
%
% Inputs
%   data1   [nChan x nTimes x nSub] for dpt, or [nChan x nTimes x n1] for idpt
%   data2   same shape as data1 for dpt, or [nChan x nTimes x n2] for idpt
%   nPerm   number of permutations. Default 1000
%   method  'mean' or 'trimmed mean'. Default 'trimmed mean' (20 percent)
%   dpt     'dpt' for paired design or 'idpt' for independent groups
%
% Outputs
%   tvals       [nChan x nTimes] observed test statistic for data1 minus data2
%   pvals       [nChan x nTimes] observed pointwise two sided p value
%   tvals_H0    [nChan x nTimes x nPerm] null test statistics
%   pvals_H0    [nChan x nTimes x nPerm] null pointwise p values
%
% Examples
%   % Paired within subject comparison
%   [t,p,tH0,pH0] = run_stats_permutation(PSD_A, PSD_B, 2000, 'mean', 'dpt');
%   % Independent groups
%   [t,p,tH0,pH0] = run_stats_permutation(PSD_grp1, PSD_grp2, 5000, 'trimmed mean', 'idpt');
%
% Notes
%   1) Synchronized permutations make TFCE and t max valid over the field.
%   2) p values here are pointwise. Use your MCC routine for FWE control.
%   3) No dependency on external yuend or limo_ttest. Everything is computed here.
%
% Author: Cedric Cannard
% Version: 2025 10 05

% Defaults
if nargin < 3 || isempty(nPerm), nPerm = 1000; end
if nargin < 4 || isempty(method), method = 'trimmed mean'; end
if nargin < 5 || isempty(dpt), error("Specify 'dpt' or 'idpt'"); end
trim = 0.20; % for trimmed mean

% Accept 2D [nTimes x nSub] by promoting to [1 x nTimes x nSub]
if ndims(data1) == 2
    data1 = reshape(data1, 1, size(data1,1), size(data1,2));
    data2 = reshape(data2, 1, size(data2,1), size(data2,2));
end
if size(data1,3) == 1 && size(data1,2) > 1
    warning('Time dimension missing or collapsed, check input shape')
end

[nChan, nTimes, nSubOrN1] = size(data1);

tvals    = nan(nChan, nTimes);
pvals    = nan(nChan, nTimes);
tvals_H0 = nan(nChan, nTimes, nPerm);
pvals_H0 = nan(nChan, nTimes, nPerm);

% Sizes
switch lower(dpt)
    case 'dpt'
        nSub = nSubOrN1;
        if size(data2,3) ~= nSub
            error('Paired input requires same number of subjects in data1 and data2.')
        end
    case 'idpt'
        n1 = nSubOrN1;
        n2 = size(data2,3);
    otherwise
        error("Unknown dependency type: use 'dpt' or 'idpt'")
end

% Observed stats
disp('Computing observed test statistics...')
parfor iChan = 1:nChan
    for t = 1:nTimes
        x1 = squeeze(data1(iChan, t, :));
        x2 = squeeze(data2(iChan, t, :));

        switch lower(method)
            case 'mean'
                if strcmpi(dpt, 'dpt')
                    % paired mean, data1 minus data2
                    d  = x1 - x2;
                    n  = numel(d);
                    md = mean(d);
                    sd = std(d, 0);
                    se = sd ./ sqrt(n);
                    tv = md ./ se;
                    df = n - 1;
                    pv = 2 * (1 - tcdf(abs(tv), df));
                else
                    % independent mean, data1 minus data2
                    nA = numel(x1); nB = numel(x2);
                    mA = mean(x1);  mB = mean(x2);
                    vA = var(x1, 0); vB = var(x2, 0);
                    se = sqrt(vA./nA + vB./nB);
                    tv = (mA - mB) ./ se;
                    df = (vA./nA + vB./nB).^2 ./ ((vA.^2)./(nA^2*(nA-1)) + (vB.^2)./(nB^2*(nB-1)));
                    pv = 2 * (1 - tcdf(abs(tv), df));
                end

            case 'trimmed mean'
                if strcmpi(dpt, 'dpt')
                    [tv, df, pv] = t_yuen_paired_d1minusd2(x1, x2, trim);
                else
                    [tv, df, pv] = t_yuen_independent_d1minusd2(x1, x2, trim);
                end

            otherwise
                error("Unknown method: choose 'mean' or 'trimmed mean'")
        end

        tvals(iChan, t) = tv;
        pvals(iChan, t) = pv;
    end
end

% Synchronized permutations
disp('Preparing synchronized permutations...')
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
disp('Running permutation test...')
parfor iChan = 1:nChan
    A = data1(iChan,:,:);   % [1 x nTimes x n]
    B = data2(iChan,:,:);

    n1_local = size(A,3);
    n2_local = size(B,3);

    for perm = 1:nPerm
        if strcmpi(dpt, 'dpt')
            % subject wise swaps
            permA = A; permB = B;
            sw = swap_mask(:,perm);
            if any(sw)
                permA(:,:,sw) = B(:,:,sw);
                permB(:,:,sw) = A(:,:,sw);
            end
        else
            % pooled label shuffle
            X = cat(3, A, B);
            idx = idx_perm(perm,:);
            permA = X(:,:,idx(1:n1_local));
            permB = X(:,:,idx(n1_local+1:n1_local+n2_local));
        end

        if strcmpi(method, 'mean')
            X1 = squeeze(permA(1,:,:))';  % [n x nTimes]
            X2 = squeeze(permB(1,:,:))';
            if strcmpi(dpt, 'dpt')
                % paired, data1 minus data2
                D  = X1 - X2;
                n  = size(D,1);
                md = mean(D,1);
                sd = std(D,0,1);
                se = sd ./ sqrt(n);
                tv = md ./ se;
                df = n - 1;
                pv = 2 * (1 - tcdf(abs(tv), df));
            else
                % independent, data1 minus data2
                nA = size(X1,1); nB = size(X2,1);
                mA = mean(X1,1);  mB = mean(X2,1);
                vA = var(X1,0,1); vB = var(X2,0,1);
                se = sqrt(vA./nA + vB./nB);
                tv = (mA - mB) ./ se;
                df = (vA./nA + vB./nB).^2 ./ ((vA.^2)./(nA^2*(nA-1)) + (vB.^2)./(nB^2*(nB-1)));
                pv = 2 * (1 - tcdf(abs(tv), df));
            end
            tvals_H0(iChan,:,perm) = tv;
            pvals_H0(iChan,:,perm) = pv;

        else
            % Trimmed mean per time using data1 minus data2
            ttmp = nan(1, nTimes);
            ptmp = nan(1, nTimes);
            for t = 1:nTimes
                xx1 = squeeze(permA(1,t,:));
                xx2 = squeeze(permB(1,t,:));
                if strcmpi(dpt, 'dpt')
                    [tv, ~, pv] = t_yuen_paired_d1minusd2(xx1, xx2, trim);
                else
                    [tv, ~, pv] = t_yuen_independent_d1minusd2(xx1, xx2, trim);
                end
                ttmp(t) = tv; ptmp(t) = pv;
            end
            tvals_H0(iChan,:,perm) = ttmp;
            pvals_H0(iChan,:,perm) = ptmp;
        end
    end
end

disp('Permutation statistics completed.')
end

% ===================== SUBFUNCTIONS =====================

function [tstat, df, p] = t_yuen_independent_d1minusd2(x, y, gamma)
% Yuen two sample trimmed mean test for data1 minus data2
x = x(:); y = y(:);
[tmx, wvx, hx] = trim_winsor_stats(x, gamma);
[tmy, wvy, hy] = trim_winsor_stats(y, gamma);

vx  = (numel(x) * wvx) / (hx * (hx - 1));
vy  = (numel(y) * wvy) / (hy * (hy - 1));
se2 = vx + vy;

tstat = (tmx - tmy) ./ sqrt(se2);
df = (se2.^2) / ( (vx.^2)/(hx - 1) + (vy.^2)/(hy - 1) );
p = 2 * (1 - tcdf(abs(tstat), df));
end

function [tstat, df, p] = t_yuen_paired_d1minusd2(x, y, gamma)
% Paired Yuen test for data1 minus data2
d = x(:) - y(:);
[tm_d, wv_d, h_d] = trim_winsor_stats(d, gamma);
se = sqrt( (numel(d) * wv_d) / (h_d * (h_d - 1)) );
tstat = tm_d ./ se;
df = h_d - 1;
p = 2 * (1 - tcdf(abs(tstat), df));
end

function [tm, wv, h] = trim_winsor_stats(v, gamma)
% Returns trimmed mean tm, winsorized variance wv, and h = n - 2g
v = sort(v(~isnan(v)));
n = numel(v);
if n < 3
    tm = mean(v); if isnan(tm), tm = 0; end
    wv = var(v,0); if isnan(wv), wv = 0; end
    h  = max(n - 2*floor(gamma*n), 1);
    return
end
g = floor(gamma * n);
h = n - 2*g;
core = v(g+1:n-g);
tm = mean(core);

% Winsorize
if g > 0
    w = v;
    w(1:g)     = v(g+1);
    w(n-g+1:n) = v(n-g);
else
    w = v;
end
wv = sum( (w - mean(w)).^2 ) / (n - 1);
end
