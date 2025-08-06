function [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_thresh)
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

% Label each cluster separately
label_mask = bwlabeln(mask);
n_cluster = max(label_mask(:));

% Preallocate
cluster_start   = nan(1, n_cluster);
cluster_end     = nan(1, n_cluster);
clust_maxval    = nan(1, n_cluster);
clust_maxchan   = nan(1, n_cluster);
clust_maxfreq   = nan(1, n_cluster);
mask_clusters   = cell(1, n_cluster);
for iClust = 1:n_cluster
    mask_i = (label_mask == iClust);
    mask_clusters{iClust} = mask_i;

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

% Merge clusters with consistent polarity if requested
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
        curr_start = clust_bounds(i, 1);
        curr_end   = clust_bounds(i, 2);
        curr_val   = clust_maxval(i);
        curr_chan  = clust_maxchan(i);
        curr_freq  = clust_maxfreq(i);
        curr_sign  = sign(median(stats(curr_mask), 'omitnan'));

        used(i) = true;

        for j = i+1:n_clusters
            if used(j), continue, end
            test_mask = mask_clusters{j};
            is_overlap = any(curr_mask & test_mask, 'all');
            is_adjacent = clust_bounds(j, 1) <= (curr_end + merge_thresh);
            if ~(is_overlap || is_adjacent), continue, end

            new_mask = curr_mask | test_mask;
            new_sign = sign(median(stats(new_mask), 'omitnan'));

            if new_sign == curr_sign
                curr_mask  = new_mask;
                curr_start = min(curr_start, clust_bounds(j,1));
                curr_end   = max(curr_end, clust_bounds(j,2));
                if abs(clust_maxval(j)) > abs(curr_val)
                    curr_val  = clust_maxval(j);
                    curr_chan = clust_maxchan(j);
                    curr_freq = clust_maxfreq(j);
                end
                used(j) = true;
                curr_sign = new_sign;
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
    merged_bounds        = clust_bounds;
    merged_maxval        = clust_maxval;
    merged_maxchan       = clust_maxchan;
    merged_maxfreq       = clust_maxfreq;
end

% Output cluster summary
if merge_thresh > 0
    mode_str = 'Merged';
else
    mode_str = 'No merge';
end
fprintf('\nCluster Summary (%s mode):\n', mode_str);
for i = 1:numel(merged_mask_clusters)
    peak_val = xaxis(merged_maxfreq(i));
    nElectrodes = sum(any(merged_mask_clusters{i}, 2));
    % Compute effect size
    if strcmpi(grp, 'dpt')
        d = merged_maxval(i) / sqrt(n{1});
        if n{1} < 30
            d = d * (1 - (3 / (4*n{1} - 5)));
        end
    elseif strcmpi(grp, 'idpt')
        d = merged_maxval(i) * sqrt(1/n{1} + 1/n{2});
        if n{1} < 30
            df = n{1} + n{2} - 2;
            d = d * (1 - (3 / (4*df - 1)));
        end
    else
        error("grp must be either 'dpt' or 'idpt'.")
    end

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
