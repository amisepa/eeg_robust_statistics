function [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_thresh, sep_thresh, polarity_split)
% PULL_CLUSTERS Extract significant clusters after MCC. 
% 
%   Identifies labeled clusters from a binary significance mask, evaluates
%   their spatial/temporal/frequency bounds, polarity, and peak effect.
%   Optionally merges adjacent or overlapping clusters that share consistent
%   polarity. Returns logical masks, cluster characteristics, ES, and a 
%   structured cluster summary.
%
% USAGE:
%   [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_gap_thresh)
%
% INPUTS:
%   mask             - Binary significance mask [channels x time/freq]
%   stats            - Observed t-values or test statistic map [same size as mask]
%   xaxis            - X-axis values (e.g., frequency in Hz or time in ms)
%   chanlocs         - EEGLAB-style channel location structure
%   datatype         - 'frequency' or 'time'
%   grp              - 'dpt' (within-group) or 'idpt' (between-group)
%   n                - Cell array with subject counts: {n1} or {n1, n2}
%   merge_thresh     - Maximum gap (in xaxis units) for merging clusters
%
% OUTPUTS:
%   mask             - Cell array of logical masks for each merged cluster
%   clust_summary    - Table summarizing each cluster (bounds, peaks, effect sizes)
%
% Author: Cedric Cannard, 2021–2025

% Early exit if nothing is significant
if sum(mask, 'all') == 0
    disp('No significant differences, nothing to plot.')
    mask = {};
    clust_summary = table;
    return
end

% sep_thresh must exist, or set a default
if ~exist('sep_thresh','var') || isempty(sep_thresh)
    sep_thresh = 2;   % Hz, tweak to taste
end

% For when some positive and negative clusters are overlapping and mixed up
% this separates them
% connectivity: minimal (no diagonals)
conn = conndef(ndims(mask), 'minimal');
if polarity_split
    warning("Re-labeling clusters by polarity (positive vs negative t-values)")
    pos_mask = (mask > 0) & (stats > 0);
    neg_mask = (mask > 0) & (stats < 0);
    pos_labels = bwlabeln(pos_mask, conn);
    neg_labels = bwlabeln(neg_mask, conn);
    if max(pos_labels(:)) > 0
        neg_labels(neg_labels > 0) = neg_labels(neg_labels > 0) + max(pos_labels(:));
    end
    label_mask = pos_labels + neg_labels;
else
    % if mask is logical, this preserves sign-less labeling
    label_mask = bwlabeln(mask, conn);
end

% Label each cluster
% label_mask = bwlabeln(mask);

% number of labeled blobs, not number of elements
n_cluster = max(label_mask(:));

% convert a gap in Hz to bins
dx = median(diff(xaxis));
min_gap_bins = max(1, round(sep_thresh / dx));

% --- Split each labeled blob by frequency gaps before bounds/peaks ---
new_masks = {};
for iClust = 1:n_cluster
    mask_i = (label_mask == iClust);

    % 1) 1D signature across frequency
    sig = any(mask_i, 1);                  % 1 x F logical

    % 2) Fill tiny gaps only if the SE has positive length
    if exist('imclose','file')
        se_len = max(0, min_gap_bins-1);
        if se_len > 0
            sig_filled = imclose(sig, ones(1, se_len));
        else
            sig_filled = sig;
        end
    else
        sig_filled = sig;
    end

    % 3) Label contiguous true runs along frequency
    if any(sig_filled)
        [L, numRuns] = bwlabel(sig_filled);   % 1D labeling
    else
        L = zeros(size(sig_filled));
        numRuns = 0;
    end

    % 4) Produce submasks per run
    if numRuns <= 1
        new_masks{end+1} = mask_i; %#ok<AGROW>
    else
        for r = 1:numRuns
            cols = find(L == r);
            if ~isempty(cols)
                sub = false(size(mask_i));
                sub(:, cols) = mask_i(:, cols);
                if any(sub(:))
                    new_masks{end+1} = sub; %#ok<AGROW>
                end
            end
        end
    end
end

% Replace with split set
mask_clusters = new_masks;
n_cluster = numel(mask_clusters);

% Guard for when min_gap_bins == 1 (structuring element becomes ones(1,0))
if exist('imclose','file') && min_gap_bins > 1
    se = ones(1, min_gap_bins-1);
    sig_filled = imclose(sig, se);
else
    sig_filled = sig;
end

% now compute bounds and peaks from mask_clusters as you already do
cluster_start = nan(1,n_cluster);
cluster_end   = nan(1,n_cluster);
clust_maxval  = nan(1,n_cluster);
clust_maxchan = nan(1,n_cluster);
clust_maxfreq = nan(1,n_cluster);
for iClust = 1:n_cluster
    mask_i = mask_clusters{iClust};
    tmp = stats .* mask_i;
    tmp(tmp == Inf | tmp == -Inf) = NaN;

    sigframes = any(mask_i, 1);
    cluster_start(iClust) = find(sigframes, 1, 'first');
    cluster_end(iClust)   = find(sigframes, 1, 'last');

    [V, type] = max([abs(min(tmp(:))) max(tmp(:))]);
    clust_maxval(iClust) = V * (-1)^(type==1);
    [clust_maxchan(iClust), clust_maxfreq(iClust)] = ind2sub(size(tmp), find(tmp == clust_maxval(iClust), 1));
end

% Sort clusters by peak xaxis value
[~, idx] = sort(xaxis(clust_maxfreq));
cluster_start   = cluster_start(idx);
cluster_end     = cluster_end(idx);
clust_maxchan   = clust_maxchan(idx);
clust_maxfreq   = clust_maxfreq(idx);
clust_maxval    = clust_maxval(idx);
mask_clusters   = mask_clusters(idx);
clust_bounds    = [cluster_start; cluster_end]';

% convert merge_thresh from x units to bin indices
dx = median(diff(xaxis));
gap = max(0, round(merge_thresh / dx));

% Merge clusters with consistent polarity if requested
if merge_thresh > 0
    N = size(clust_bounds, 1);

    % stable polarity from peak t
    clust_sign = sign(clust_maxval);

    % helper: overlap or adjacency in index space
    overlaps_or_adjacent = @(a,b) (b(1) <= a(2)+gap) && (b(2) >= a(1)-gap);
    contained = @(a,b) (b(1) >= a(1) && b(2) <= a(2)) || (a(1) >= b(1) && a(2) <= b(2));

    % build adjacency matrix for clusters that should be merged
    A = false(N,N);
    for i = 1:N-1
        for j = i+1:N
            if clust_sign(i) ~= clust_sign(j), continue, end
            a = clust_bounds(i,:); b = clust_bounds(j,:);
            if overlaps_or_adjacent(a,b) || contained(a,b)
                A(i,j) = true; A(j,i) = true;
            end
        end
    end

    % connected components = merge groups
    if any(A(:))
        G = graph(A);
        comp = conncomp(G);             % 1..K labels
    else
        comp = 1:N;                     % no merges
    end

    K = max(comp);
    merged_mask_clusters = cell(1,K);
    merged_bounds = nan(K,2);
    merged_maxval = nan(1,K);
    merged_maxchan = nan(1,K);
    merged_maxfreq = nan(1,K);

    for k = 1:K
        idx = find(comp == k);
        % merge masks
        m = false(size(mask_clusters{1}));
        for t = idx(:)'
            m = m | mask_clusters{t};
        end
        merged_mask_clusters{k} = m;

        % merged bounds
        merged_bounds(k,1) = min(clust_bounds(idx,1));
        merged_bounds(k,2) = max(clust_bounds(idx,2));

        % representative peak = largest |t|
        [~, ii] = max(abs(clust_maxval(idx)));
        pick = idx(ii);
        merged_maxval(k)  = clust_maxval(pick);
        merged_maxchan(k) = clust_maxchan(pick);
        merged_maxfreq(k) = clust_maxfreq(pick);
    end
else
    merged_mask_clusters = mask_clusters;
    merged_bounds        = clust_bounds;
    merged_maxval        = clust_maxval;
    merged_maxchan       = clust_maxchan;
    merged_maxfreq       = clust_maxfreq;
end

% Report and summarize
if merge_thresh > 0
    mode_str = 'Merged';
else
    mode_str = 'No merge';
end
fprintf('\nCluster Summary (%s mode):\n', mode_str);
effect_sizes = nan(1, numel(merged_mask_clusters));

for i = 1:numel(merged_mask_clusters)
    peak_val = xaxis(merged_maxfreq(i));
    nElectrodes = sum(any(merged_mask_clusters{i}, 2));
    t = merged_maxval(i);  

    % Compute effect size
    if strcmpi(grp, 'dpt')
        d = t / sqrt(n{1});
        if n{1} < 30
            % d = d * (1 - (3 / (4*n{1} - 5)));         % Hedges' g
            % % d = sqrt(t^2 / (t^2 + n{1} - 1));         % Use r² equivalent for within-subject single-trials
            d = t^2 / (t^2 + n{1} - 1);                 % or use η² (eta-squared) for within-subject single-trials 
        end
    elseif strcmpi(grp, 'idpt')
        d = t * sqrt(1/n{1} + 1/n{2});
        if n{1} < 30
            df = n{1} + n{2} - 2;
            d = d * (1 - (3 / (4*df - 1)));  % Hedges' g
        end
    else
        error("grp must be either 'dpt' or 'idpt'.")
    end

    effect_sizes(i) = d;

    if strcmpi(datatype, 'frequency')
        fprintf('Cluster %d: %g to %g Hz (%g channels). Peak: %s at %g Hz (t = %.2f; ES = %.2f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), ...
            nElectrodes, chanlocs(merged_maxchan(i)).labels, peak_val, merged_maxval(i), d);
    else
        fprintf('Cluster %d: %g to %g ms (%g channels). Peak: %s at %g ms (t = %.2f; ES = %.2f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), ...
            nElectrodes, chanlocs(merged_maxchan(i)).labels, peak_val, merged_maxval(i), d);
    end
end

% Construct summary table
nMerged = numel(merged_mask_clusters);
clust_summary = table;
clust_summary.Cluster = (1:nMerged)';
clust_summary.Start = xaxis(merged_bounds(:,1))';
clust_summary.End   = xaxis(merged_bounds(:,2))';
clust_summary.Peak  = xaxis(merged_maxfreq)';
clust_summary.Tvalue = merged_maxval';
clust_summary.Channel = string({chanlocs(merged_maxchan).labels})';
clust_summary.Polarity = repmat("", nMerged, 1);
clust_summary.Polarity(merged_maxval > 0) = "Positive";
clust_summary.Polarity(merged_maxval < 0) = "Negative";
clust_summary.NumElectrodes = cellfun(@(m) sum(any(m, 2)), merged_mask_clusters)';

if strcmpi(datatype, 'time')
    clust_summary.Properties.VariableNames(2:4) = {'Start_ms', 'End_ms', 'Peak_ms'};
elseif strcmpi(datatype, 'frequency')
    clust_summary.Properties.VariableNames(2:4) = {'Start_Hz', 'End_Hz', 'Peak_Hz'};
end

% Return final mask
mask = merged_mask_clusters;

