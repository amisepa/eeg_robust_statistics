function [mask, pvals, max_th] = correct_max(tfce_obs, tfce_H0, pthresh)
% CORRECT_MAX  Familywise correction from TFCE null maxima
%
%   [mask, pvals, max_th] = correct_max(tfce_obs, tfce_H0, pthresh)
%
% Inputs
%   tfce_obs : observed TFCE map [X x Y] or [X x Y x Z]
%   tfce_H0  : TFCE null maps stacked along 3rd or 4th dim
%   pthresh  : familywise alpha (default = 0.05)
%
% Outputs
%   mask     : significant elements after FWER correction
%   pvals    : corrected p-values for all elements
%   max_th   : TFCE threshold at alpha
%
% Notes
%   Uses max-TFCE null distribution to ensure strong FWER control.
%
%   p(i) = ( #perms with max(TFCE_H0) >= tfce_obs(i) + 1 ) / (nPerm + 1)
%

if nargin < 3 || isempty(pthresh), pthresh = 0.05; end

% determine permutation dimension automatically
nDim = ndims(tfce_H0);
if nDim == 3
    permDim = 3;
elseif nDim == 4
    permDim = 4;
else
    error('tfce_H0 must be 3D or 4D');
end

% gather maxima across permutations
maxval = squeeze(max(tfce_H0,[],1));            % collapse first dimension
maxval = squeeze(max(maxval,[],1));             % collapse remaining spatial dims
maxval = maxval(isfinite(maxval));              % remove NaNs/Infs

if isempty(maxval)
    error('No finite values found in tfce_H0');
end

% sort and get cutoff
sortmax = sort(maxval(:));
n = numel(sortmax);
cutIndex = max(1, ceil((1 - pthresh) * n));
max_th = sortmax(cutIndex);

% FWER mask
mask = tfce_obs >= max_th;

% corrected p-values
nPerm = numel(maxval);
flat_obs = tfce_obs(:);
pflat = zeros(size(flat_obs));

for i = 1:numel(flat_obs)
    pflat(i) = (sum(maxval >= flat_obs(i)) + 1) / (nPerm + 1);
end

pvals = reshape(pflat, size(tfce_obs));
end
