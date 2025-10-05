function [mask, pvals, max_th] = correct_max(tvals, tvals_H0, pthresh)
% CORRECT_MAX
% One sided t max correction using the distribution of maxima under H0.
% p-values are FWER controlled: p(i,j) = mean(max_H0 >= tvals(i,j))

% inputs
if nargin < 3 || isempty(pthresh), pthresh = 0.05; end

% collect per-permutation maxima
[na, nb, nboot] = size(tvals_H0);
maxval = nan(nboot,1);
for k = 1:nboot
    data = tvals_H0(:,:,k);
    maxval(k) = max(data(:));
end

% drop Inf and NaN
maxval = maxval(isfinite(maxval));
if isempty(maxval)
    error('correct_max: empty max distribution after removing NaN/Inf')
end

% sort and get threshold
sortmaxM = sort(maxval(:));
nboot = numel(sortmaxM);
U = round((1 - pthresh) * nboot);
U = max(1, min(U, nboot));           % clamp
max_th = sortmaxM(U);
mask   = tvals >= max_th;

fprintf('Maximum threshold = %g\n', max_th)

% FWER p-values: P(max_H0 >= t_obs)
pvals = nan(na, nb);
for i = 1:na
    for j = 1:nb
        pvals(i,j) = mean(sortmaxM >= tvals(i,j));
        if pvals(i,j) == 0
            pvals(i,j) = 1 / nboot;  % smallest attainable
        end
    end
end

% figure for diagnostics when no detections
if ~any(mask(:))
    figure('Name','Results under H0 after max-correction')
    plot(sortmaxM, 'LineWidth', 3); grid on; hold on;

    % mark threshold
    idxU = U;
    plot(idxU, max_th, 'r*', 'LineWidth', 5)
    txt = ['Bootstrap threshold ' num2str(max_th) ' ->'];
    text(idxU, max_th, txt, 'FontSize', 12, 'HorizontalAlignment', 'right');

    % mark observed field maximum
    [tmax, loc] = max(tvals(:));
    plot(loc, tmax, 'r*', 'LineWidth', 5)
    txt = ['Maximum observed: ' num2str(tmax) ' ->'];
    text(loc, tmax, txt, 'FontSize', 12, 'HorizontalAlignment', 'right');

    title('Maxima under H0', 'FontSize', 12)
    xlabel('Sorted bootstrap iterations', 'FontSize', 12)
    ylabel('t max', 'FontSize', 12)
    box on; set(gca,'Layer','Top')
end
end
