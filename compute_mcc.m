function [mask, clust_summary] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, pthresh, xaxis, datatype, chanlocs, grp, n)
% COMPUTE_MCC Perform multiple comparisons correction on t-values.
%
%   [mask, clust_summary] = compute_mcc(...)
%
% INPUTS:
%   tvals          - Observed t-values [channels x freq/time]
%   pvals          - Uncorrected p-values
%   tvals_H0       - Null t-values [chan x time x boot]
%   pvals_H0       - Null p-values [chan x time x boot]
%   mcctype        - Type of correction: 0=uncorrected, 1=max-t, 2=cluster, 3=TFCE
%   pthresh        - Significance threshold (e.g., 0.05)
%   neighbormatrix - Matrix of spatial adjacency
%   xaxis          - X-axis values (e.g., frequency or time)
%   datatype       - 'frequency' or 'time' (used for labeling)
%   chanlocs       - EEGLAB-style structure with channel labels
%   grp            - specify if group is 'dpt' or 'idpt'
%   n              - number of subjects
%
% OUTPUTS:
%   mask           - Binary mask of significant clusters
%   pcorr          - Corrected p-values
%   nClust         - Number of significant clusters
%   summary_table  - Table of merged cluster summaries
%   mask_clusters  - Cell array of binary masks (1 per merged cluster)

% Author: Cedric Cannard, modified July 2025

merge_gap_thresh = 5;

% Prep outputs
mask = [];
pcorr = [];
clust_bounds = [];
clust_maxchan = [];
clust_maxfreq = [];
clust_maxval = [];
nClust = [];
clust_summary = table;

nChan = size(tvals,1);
nTimes = size(tvals,2);
stats = tvals; % alias for clarity

% Ensure xaxis exists and has correct length
if nargin < 8 || isempty(xaxis)
    warning('xaxis not provided. Using default 1:nTimes.')
    xaxis = 1:nTimes;
end

% Find electrdode neighbors
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs);


% ------------------------- MCC Type -----------------------------
switch mcctype
    case 0 % Uncorrected
        pcorr = pvals;
        mask = pcorr <= pthresh;
        nClust = length(unique(bwlabeln(mask))) - 1;

    case 1 % Max-t
        [mask, pcorr] = correct_max(abs(tvals), abs(tvals_H0), pthresh);

    case 2 % Cluster correction
        [mask, pcorr] = correct_cluster(tvals.^2, pvals, tvals_H0.^2, pvals_H0, neighbormatrix, mcctype, pthresh);
        nClust = length(unique(mask)) - 1;

    case 3 % TFCE
        ndim = ndims(tvals);
        nBoot = size(tvals_H0, 3);
        tfce_H0_score = nan(nChan, nTimes, nBoot);
        fprintf("Performing threshold-free cluster enhancement on H0 data (%g nBoot permutations)...\n", nBoot)
        parfor b = 1:nBoot
            fprintf(' - boot %g\n', b)
            tfce_H0_score(:,:,b) = limo_tfce(ndim, squeeze(tvals_H0(:,:,b)), neighbormatrix, 0);
        end
        tfce_score = limo_tfce(ndim, tvals, neighbormatrix);
        [mask, pcorr] = correct_max(tfce_score, tfce_H0_score, pthresh);
        nClust = length(unique(mask)) - 1;
end

if sum(mask,'all') == 0
    disp('No significant differences, nothing to plot.')
    return
end
% if isempty(mask)
%     disp('Empty mask, computing one from corrected p-values.');
%     mask = pcorr < alpha;
% end

% % Fix polarity conflict: re-label clusters by polarity
% if ~isempty(mask) && any(mask(:))
%     warning("Overwriting mask to address cluster ploratiy conflict!")
% 
%     % 1. Get positive and negative t-values within the mask
%     pos_mask = (mask > 0) & (tvals > 0);
%     neg_mask = (mask > 0) & (tvals < 0);
% 
%     % 2. Label each separately
%     pos_labels = bwlabeln(pos_mask);
%     neg_labels = bwlabeln(neg_mask);
% 
%     % 3. Offset neg_labels to avoid overlap
%     if max(pos_labels(:)) > 0
%         neg_labels(neg_labels > 0) = neg_labels(neg_labels > 0) + max(pos_labels(:));
%     end
% 
%     % 4. Recombine into final labeled mask
%     mask = pos_labels + neg_labels;
% 
%     % 5. Update cluster count
%     nClust = length(unique(mask)) - 1;
% end

% Get clusters properties (start/end, peak, channel, frame/freq)
n_cluster     = max(mask(:));
cluster_start = nan(1,n_cluster);
cluster_end   = nan(1,n_cluster);
clust_maxval  = nan(1,n_cluster);
clust_maxchan = nan(1,n_cluster);
clust_maxfreq = nan(1,n_cluster);

for iClust = 1:n_cluster
    % Binary mask for this cluster
    mask_i = (mask == iClust);
    mask_clusters{iClust} = mask_i;

    % Restrict stats to this cluster only
    tmp = stats .* mask_i;
    tmp(tmp == Inf | tmp == -Inf) = NaN;

    % % Get bounds of the cluster
    % sigframes = sum(tmp,1);
    % cluster_start(iClust) = find(sigframes, 1, 'first');
    % cluster_end(iClust) = find(sigframes, 1, 'last');

    % Use mask (not weighted by t-values) to find cluster bounds
    sigmask = mask == iClust;
    sigframes = any(sigmask, 1);  % any significant channel at each frame
    cluster_start(iClust) = find(sigframes, 1, 'first');
    cluster_end(iClust)   = find(sigframes, 1, 'last');

    % Get max value and location
    [V, type] = max([abs(min(tmp(:))) max(tmp(:))]);
    if type == 2
        clust_maxval(iClust) = V(1);
    else
        V = -V;
        clust_maxval(iClust) = V(1);
    end
    [clust_maxchan(iClust), clust_maxfreq(iClust)] = ind2sub(size(tmp), find(tmp == V(1), 1));
end

% Convert frequency or latency values at cluster peaks
cluster_peak_xvals = xaxis(clust_maxfreq);

% Sort clusters by frequency or time
[~, idx] = sort(cluster_peak_xvals, 'ascend');
cluster_start     = cluster_start(idx);
cluster_end       = cluster_end(idx);
clust_maxchan   = clust_maxchan(idx);
clust_maxfreq   = clust_maxfreq(idx);
clust_maxval    = clust_maxval(idx);
mask_clusters     = mask_clusters(idx);  % sort masks accordingly

% Output cluster bounds
clust_bounds = [cluster_start; cluster_end]';

% AUTO-MERGE SIMILAR CLUSTERS (with polarity re-evaluation after each merge)
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

    % Start new merged cluster with cluster i
    curr_mask = mask_clusters{i};
    curr_start = clust_bounds(i,1);
    curr_end   = clust_bounds(i,2);
    curr_val   = clust_maxval(i);
    curr_chan  = clust_maxchan(i);
    curr_freq  = clust_maxfreq(i);

    % Evaluate current polarity
    curr_vals = stats(curr_mask);
    curr_vals = curr_vals(~isnan(curr_vals));
    curr_sign = sign(median(curr_vals));  % robust polarity

    used(i) = true;

    % Check each other cluster j for potential merge
    for j = i+1:n_clusters
        if used(j), continue, end

        % Check overlap or adjacency
        test_mask = mask_clusters{j};
        is_overlap = any(curr_mask & test_mask, 'all');
        is_adjacent = clust_bounds(j,1) <= (curr_end + merge_gap_thresh);

        if ~(is_overlap || is_adjacent), continue, end

        % Simulate merged mask
        new_mask = curr_mask | test_mask;
        new_vals = stats(new_mask);
        new_vals = new_vals(~isnan(new_vals));
        new_sign = sign(median(new_vals));

        % Only allow merge if merged polarity is consistent
        if new_sign == curr_sign
            curr_mask = new_mask;
            curr_start = min(curr_start, clust_bounds(j,1));
            curr_end   = max(curr_end, clust_bounds(j,2));

            % Update peak if stronger
            if abs(clust_maxval(j)) > abs(curr_val)
                curr_val  = clust_maxval(j);
                curr_chan = clust_maxchan(j);
                curr_freq = clust_maxfreq(j);
            end

            used(j) = true;

            % Update polarity to new cluster (important!)
            curr_sign = new_sign;
        end
    end

    % Store final merged cluster
    merged_mask_clusters{merge_id} = curr_mask;
    merged_bounds(merge_id,:)     = [curr_start curr_end];
    merged_maxval(merge_id)       = curr_val;
    merged_maxchan(merge_id)      = curr_chan;
    merged_maxfreq(merge_id)      = curr_freq;
    merge_id = merge_id + 1;
end


% Output merged clusters summary
fprintf('\nMerged Cluster Summary:\n');
for i = 1:numel(merged_mask_clusters)
    peak_freq_val = xaxis(merged_maxfreq(i));
    nElectrodes = sum(any(merged_mask_clusters{i}, 2)); % number of electrodes in cluster

    % Cohen's d
    if strcmpi(grp, 'dpt') % within-group comparisons 
        d = merged_maxval(i) / sqrt(n{:});   % Cohen's dz = t / √n
    elseif strcmpi(grp, 'idpt') % between group comparisons
        d = merged_maxval(i) * sqrt(1/n{1} + 1/n{2});    % Cohen's d = t × √(1/n₁ + 1/n₂)
    else
        error("grp must be either 'dpt' or 'idpt'.")
    end


    if strcmpi(datatype, 'frequency')
        fprintf('Cluster %d: %g to %g Hz (%g channels). Peak effect: %s at %g Hz (t = %.2f; Cohen''s d = %.2f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), nElectrodes, ...
            chanlocs(merged_maxchan(i)).labels, peak_freq_val, merged_maxval(i), d);
    else
        fprintf('Cluster %d: %g to %g ms (%g channels). Peak effect: %s at %g ms (t = %.2f; Cohen''s d = %.2f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), nElectrodes, ...
            chanlocs(merged_maxchan(i)).labels, peak_freq_val, merged_maxval(i), d);
    end
end

% Create summary table
clust_summary = table;
nMerged = numel(merged_mask_clusters);
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

% Rename columns based on datatype
if strcmpi(datatype, 'time')
    clust_summary.Properties.VariableNames(2:4) = {'Start_ms', 'End_ms', 'Peak_ms'};
elseif strcmpi(datatype, 'frequency')
    clust_summary.Properties.VariableNames(2:4) = {'Start_Hz', 'End_Hz', 'Peak_Hz'};
end

% Assign final mask_clusters to merged ones
mask = merged_mask_clusters;
clust_bounds = merged_bounds;
clust_maxchan = merged_maxchan;
clust_maxfreq = merged_maxfreq;
clust_maxval = merged_maxval;
