function [mask, pvals, max_th] = correct_max(tvals, tvals_H0, pthresh)
if nargin < 3 || isempty(pthresh), pthresh = 0.05; end

% collect per-permutation maxima of TFCE scores
maxval = squeeze(max(max(tvals_H0,[],1),[],2));   % [nPerm x 1]
maxval = maxval(isfinite(maxval));
if isempty(maxval), error('correct_max: empty max distribution'); end

% threshold for FWER alpha via order statistic
sortmax = sort(maxval(:));
nboot   = numel(sortmax);
U       = max(1, min(round((1 - pthresh) * nboot), nboot));
max_th  = sortmax(U);

mask = tvals >= max_th;

% FWER p-values: P(max_H0 >= t_obs) with (k+1)/(n+1)
cmp   = maxval >= tvals(:).';                % [nPerm x nVox]
pflat = (sum(cmp,1) + 1) ./ (nboot + 1);     % 1 x nVox
pvals = reshape(pflat, size(tvals));
end
