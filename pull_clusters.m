function [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_thresh, polarity_split)
% PULL_CLUSTERS Extract and optionally merge significant clusters after MCC.
%
%   This function identifies labeled clusters from a binary significance mask,
%   optionally separates them by polarity (positive vs negative effects),
%   evaluates their spatial/temporal/frequency properties, and merges adjacent
%   clusters if they are contiguous and have consistent polarity.
%
%   It returns logical masks for each cluster and a summary table with detailed
%   cluster properties including bounds, polarity, peak effect, and effect size.
%
% USAGE:
%   [mask, clust_summary] = pull_clusters( ...
%       mask_corr, tvals, xaxis, chanlocs, datatype, grp, n, merge_gap_thresh, polarity_split)
%
% INPUTS:
%   mask             - Binary significance mask [nChans x nFreqs or nTimes]
%   stats            - Test statistic map (e.g., t-values), same size as mask
%   xaxis            - Vector of frequency or time values (e.g., 1:64 Hz)
%   chanlocs         - EEGLAB-style structure with channel labels
%   datatype         - String: 'frequency' or 'time'
%   grp              - String: 'dpt' for within-subject or 'idpt' for between-subject
%   n                - Cell array of subject counts: {n1} or {n1, n2}
%   merge_gap_thresh - Maximum gap (in xaxis steps) allowed for merging clusters (set to 0 to disable)
%   polarity_split   - Boolean flag to enable polarity-based cluster labeling (e.g., required for TFCE)
%
% OUTPUTS:
%   mask             - Cell array of logical masks (1 per merged cluster)
%   clust_summary    - Table summarizing:
%                         - Cluster number
%                         - Start, end, and peak (Hz or ms)
%                         - T-value at peak
%                         - Effect size (Cohen's d or Hedges' g)
%                         - Peak channel
%                         - Cluster polarity ('Positive' or 'Negative')
%                         - Number of contributing electrodes
%
% NOTES:
%   - Use polarity_split = true when correction method (e.g., TFCE) returns a single mask that includes both positive and negative effects.
%   - Merged clusters must be both adjacent (within merge_gap_thresh) and share the same polarity.
%   - Peak effect size is based on the maximum absolute t-value within each cluster.
%
% EXAMPLE:
%   [mask, clust_tbl] = pull_clusters(mask_tfce, tvals, freqs, chanlocs, ...
%                                     'frequency', 'dpt', {30}, 3, true);
%
% Author: Cedric Cannard, 2021–2025

% Skip if empty
if sum(mask, 'all') == 0
    disp('No significant differences, nothing to plot.')
    mask = {};
    clust_summary = table;
    return
end

% --- Step 1: Separate positive and negative clusters ---
if polarity_split
    warning("Re-labeling clusters by polarity (positive vs negative t-values)")
    
    pos_mask = (mask > 0) & (stats > 0);
    neg_mask = (mask > 0) & (stats < 0);

    pos_labels = bwlabeln(pos_mask);
    neg_labels = bwlabeln(neg_mask);

    if max(pos_labels(:)) > 0
        neg_labels(neg_labels > 0) = neg_labels(neg_labels > 0) + max(pos_labels(:));
    end

    mask = pos_labels + neg_labels;
end

% --- Step 2: Extract initial clusters ---
n_cluster     = max(mask(:));
cluster_start = nan(1,n_cluster);
cluster_end   = nan(1,n_cluster);
clust_maxval  = nan(1,n_cluster);
clust_maxchan = nan(1,n_cluster);
clust_maxfreq = nan(1,n_cluster);
mask_clusters = cell(1, n_cluster);
for iClust = 1:n_cluster
    mask_i = (mask == iClust);
    mask_clusters{iClust} = mask_i;

    tmp = stats .* mask_i;
    tmp(tmp == Inf | tmp == -Inf) = NaN;

    sigframes = any(mask_i, 1);
    cluster_start(iClust) = find(sigframes, 1, 'first');
    cluster_end(iClust)   = find(sigframes, 1, 'last');

    [V, type] = max([abs(min(tmp(:))) max(tmp(:))]);
    clust_maxval(iClust) = V * (-1)^(type == 1);
    [clust_maxchan(iClust), clust_maxfreq(iClust)] = ind2sub(size(tmp), find(tmp == clust_maxval(iClust), 1));
end

% --- Step 3: Sort by frequency or time ---
[~, idx] = sort(xaxis(clust_maxfreq), 'ascend');
cluster_start   = cluster_start(idx);
cluster_end     = cluster_end(idx);
clust_maxchan   = clust_maxchan(idx);
clust_maxfreq   = clust_maxfreq(idx);
clust_maxval    = clust_maxval(idx);
mask_clusters   = mask_clusters(idx);
clust_bounds    = [cluster_start; cluster_end]';

% --- Step 4: Merge adjacent clusters if enabled ---
if merge_thresh > 0
    n_clusters = size(clust_bounds, 1);
    merged_mask_clusters = {};
    merged_bounds = [];
    merged_maxval = [];
    merged_maxchan = [];
    merged_maxfreq = [];
    used = false(1, n_clusters);
    merge_id = 1;

    for i = 1:n_clusters
        if used(i), continue, end
        curr_mask  = mask_clusters{i};
        curr_start = clust_bounds(i,1);
        curr_end   = clust_bounds(i,2);
        curr_val   = clust_maxval(i);
        curr_chan  = clust_maxchan(i);
        curr_freq  = clust_maxfreq(i);
        curr_vals  = stats(curr_mask);
        curr_sign  = sign(median(curr_vals(~isnan(curr_vals))));

        used(i) = true;

        for j = i+1:n_clusters
            if used(j), continue, end
            test_mask = mask_clusters{j};
            is_overlap = any(curr_mask & test_mask, 'all');
            is_adjacent = clust_bounds(j,1) <= (curr_end + merge_thresh);
            if ~(is_overlap || is_adjacent), continue, end

            test_vals = stats(test_mask);
            test_sign = sign(median(test_vals(~isnan(test_vals))));

            if test_sign == curr_sign
                curr_mask = curr_mask | test_mask;
                curr_start = min(curr_start, clust_bounds(j,1));
                curr_end   = max(curr_end, clust_bounds(j,2));
                if abs(clust_maxval(j)) > abs(curr_val)
                    curr_val  = clust_maxval(j);
                    curr_chan = clust_maxchan(j);
                    curr_freq = clust_maxfreq(j);
                end
                used(j) = true;
            end
        end

        merged_mask_clusters{merge_id} = curr_mask;
        merged_bounds(merge_id,:)     = [curr_start curr_end];
        merged_maxval(merge_id)       = curr_val;
        merged_maxchan(merge_id)      = curr_chan;
        merged_maxfreq(merge_id)      = curr_freq;
        merge_id = merge_id + 1;
    end
else
    merged_mask_clusters = mask_clusters;
    merged_bounds = clust_bounds;
    merged_maxval = clust_maxval;
    merged_maxchan = clust_maxchan;
    merged_maxfreq = clust_maxfreq;
end

% --- Step 5: Report and summarize ---
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
            % Optional: choose based on what metric you want
            % d = d * (1 - (3 / (4*n{1} - 5)));         % Hedges' g
            % d = sqrt(t^2 / (t^2 + n{1} - 1));         % Use r² equivalent for within-subject single-trials  (proportion of variance explained)
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

% Create summary table
nMerged = numel(merged_mask_clusters);
clust_summary = table;
clust_summary.Cluster = (1:nMerged)';
clust_summary.Start = xaxis(merged_bounds(:,1))';
clust_summary.End   = xaxis(merged_bounds(:,2))';
clust_summary.Peak  = xaxis(merged_maxfreq)';
clust_summary.Tvalue = merged_maxval';
clust_summary.EffectSize = effect_sizes';
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

% Final output
mask = merged_mask_clusters;
