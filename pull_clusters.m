function [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_thresh, sep_thresh, polarity_split, es_metric)
% PULL_CLUSTERS Extract significant clusters after MCC. 
% 
%   Identifies labeled clusters from a binary significance mask, evaluates
%   their spatial/temporal/frequency bounds, polarity, and peak effect.
%   Optionally merges adjacent or overlapping clusters that share consistent
%   polarity. Returns logical masks, cluster characteristics, ES, and a 
%   structured cluster summary.
%
% USAGE:
%   [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_thresh, sep_thresh, polarity_split, es_metric)
%
% INPUTS:
%   mask             - Binary significance mask [channels x time/freq]
%   stats            - Observed t-values or test statistic map [same size as mask]
%   xaxis            - X-axis values (e.g., frequency in Hz or time in ms)
%   chanlocs         - EEGLAB-style channel location structure
%   datatype         - 'frequency' or 'time'
%   grp              - 'dpt' (within-group) or 'idpt' (between-group)
%   n                - Cell array with subject counts: {n1} or {n1, n2}
%   merge_thresh     - Maximum gap in xaxis units (ms or Hz) for merging 
%                       clusters (default = 0, no merging)
%   sep_thresh       - Minimum width in bins/frames to keep a cluster
%                       (default = 1)
%   polarity_split   - If true, makes polarity exclusive per x bin to get
%                       both positive and negative clusters as separate (default = false)
%   es_metric        - Effect size (ES) metric to compute: 'd' (cohen's d; default),
%                       'g (hedge g)', 'eta2', 'r'
%
% OUTPUTS:
%   mask             - Cell array of logical masks for each merged cluster
%   clust_summary    - Table summarizing each cluster (bounds, peaks, effect sizes)
%
% Author: Cedric Cannard, 2021–2025

% 1) Early exit if nothing is significant
if nnz(mask) == 0
    disp('No significant differences, nothing to plot.')
    mask = {};
    clust_summary = table;
    return
end

% 2) Defaults for thresholds
if ~exist('sep_thresh','var') || isempty(sep_thresh)
    sep_thresh = 1;   % in frames or freq bins
end
if ~exist('merge_thresh','var') || isempty(merge_thresh)
    merge_thresh = 0; % in xaxis units (ms or Hz); 0 = no merging
end
if ~exist('polarity_split','var') || isempty(polarity_split)
    polarity_split = false;
end

% Default ES metric and options
if ~exist('es_metric','var') || isempty(es_metric)
    es_metric = 'd';   % 'd', 'g', 'eta2', 'r'
end

if ~exist('es_opts','var') || isempty(es_opts)
    es_opts = struct();      % optional: es_opts.df, es_opts.welch, es_opts.s1, es_opts.s2
end
if strcmpi(es_metric, 'r')
    es_opts.df = n{:} - 2;
end

% Optional polarity exclusivity per x bin
if polarity_split
    [~, idx_row] = max(abs(stats), [], 1);
    lin_idx = sub2ind(size(stats), idx_row, 1:numel(idx_row));
    sign_at_absmax = sign(stats(lin_idx));
    sgn_pos = repmat(sign_at_absmax > 0, size(mask,1), 1);
    sgn_neg = repmat(sign_at_absmax < 0, size(mask,1), 1);
    pos_mask_all = mask & sgn_pos;
    neg_mask_all = mask & sgn_neg;

    conn = 4;
    pos_labels = bwlabel(pos_mask_all, conn);
    neg_labels = bwlabel(neg_mask_all, conn);
    if max(pos_labels(:)) > 0
        neg_labels(neg_labels > 0) = neg_labels(neg_labels > 0) + max(pos_labels(:));
    end
    label_img = pos_labels + neg_labels;
else
    conn = 4;
    label_img = bwlabel(mask, conn);
end

% Extract initial cluster masks
n_cluster = max(label_img(:));
mask_clusters = cell(1,0);
for iClust = 1:n_cluster
    m = (label_img == iClust);
    if any(m(:))
        mask_clusters{end+1} = m; %#ok<AGROW>
    end
end

if isempty(mask_clusters)
    disp('No clusters were labeled.')
    mask = {};
    clust_summary = table;
    return
end

% Consolidate clusters that overlap in time (bwlabel may fragment across
% non-adjacent channel indices that are part of the same temporal event)
changed = true;
while changed
    changed = false;
    for ii = 1:numel(mask_clusters)-1
        if ~any(mask_clusters{ii}(:)), continue, end
        sig_i = any(mask_clusters{ii}, 1);   % 1 x T
        for jj = ii+1:numel(mask_clusters)
            if ~any(mask_clusters{jj}(:)), continue, end
            sig_j = any(mask_clusters{jj}, 1);
            if any(sig_i & sig_j)  % temporal overlap
                mask_clusters{ii} = mask_clusters{ii} | mask_clusters{jj};
                mask_clusters{jj} = false(size(mask_clusters{jj}));
                sig_i = any(mask_clusters{ii}, 1);
                changed = true;
            end
        end
    end
end
mask_clusters = mask_clusters(cellfun(@(m) any(m(:)), mask_clusters));

% % The above code merges any two clusters that share even a single frequency bin in any channel.
% min_chan_overlap = max(2, round(0.1 * nChan));  % e.g. at least 10% of channels
% changed = true;
% while changed
%     changed = false;
%     for ii = 1:numel(mask_clusters)-1
%         if ~any(mask_clusters{ii}(:)), continue, end
%         sig_i = any(mask_clusters{ii}, 1);
%         for jj = ii+1:numel(mask_clusters)
%             if ~any(mask_clusters{jj}(:)), continue, end
%             sig_j = any(mask_clusters{jj}, 1);
%             overlap_bins = sig_i & sig_j;
%             % Count how many channels are active in both at overlapping bins
%             if any(overlap_bins)
%                 chan_overlap = sum(any(mask_clusters{ii}(:, overlap_bins), 2) & ...
%                                    any(mask_clusters{jj}(:, overlap_bins), 2));
%                 if chan_overlap >= min_chan_overlap
%                     mask_clusters{ii} = mask_clusters{ii} | mask_clusters{jj};
%                     mask_clusters{jj} = false(size(mask_clusters{jj}));
%                     sig_i = any(mask_clusters{ii}, 1);
%                     changed = true;
%                 end
%             end
%         end
%     end
% end

% Split clusters by internal gaps along x using sep_thresh
min_gap_bins = max(1, round(sep_thresh));
new_masks = {};
for iClust = 1:numel(mask_clusters)
    m = mask_clusters{iClust};
    sig = any(m, 1);

    if exist('imclose','file') == 2
        se_len = max(0, min_gap_bins-1);
        if se_len > 0
            sig_filled = imclose(sig, ones(1, se_len));
        else
            sig_filled = sig;
        end
    else
        sig_filled = sig;
    end

    if any(sig_filled)
        [L, numRuns] = bwlabel(sig_filled);
    else
        L = zeros(size(sig_filled));
        numRuns = 0;
    end

    if numRuns <= 1
        new_masks{end+1} = m; %#ok<AGROW>
    else
        for r = 1:numRuns
            cols = find(L == r);
            if ~isempty(cols)
                sub = false(size(m));
                sub(:, cols) = m(:, cols);
                if any(sub(:))
                    new_masks{end+1} = sub; %#ok<AGROW>
                end
            end
        end
    end
end
mask_clusters = new_masks;

% Prune clusters whose width is below sep_thresh
min_width_bins = max(1, round(sep_thresh));
keep = true(1, numel(mask_clusters));
for i = 1:numel(mask_clusters)
    sigcols = any(mask_clusters{i}, 1);
    if ~any(sigcols)
        keep(i) = false;
        continue
    end
    first_col = find(sigcols, 1, 'first');
    last_col  = find(sigcols, 1, 'last');
    width_bins = last_col - first_col + 1;
    if width_bins < min_width_bins
        keep(i) = false;
    end
end
mask_clusters = mask_clusters(keep);

if isempty(mask_clusters)
    disp('All clusters were pruned by sep_thresh.')
    mask = {};
    clust_summary = table;
    return
end

% Compute bounds and peaks
n_cluster = numel(mask_clusters);
cluster_start = nan(1,n_cluster);
cluster_end   = nan(1,n_cluster);
clust_maxval  = nan(1,n_cluster);
clust_maxchan = nan(1,n_cluster);
clust_maxfreq = nan(1,n_cluster);

for iClust = 1:n_cluster
    m = mask_clusters{iClust};
    tmp = stats .* m;
    tmp(~m) = NaN;

    sigframes = any(m, 1);
    cluster_start(iClust) = find(sigframes, 1, 'first');
    cluster_end(iClust)   = find(sigframes, 1, 'last');

    [mx_pos, idx_pos] = max(tmp(:), [], 'omitnan');
    [mx_neg, ~]       = min(tmp(:), [], 'omitnan');
    if abs(mx_neg) > abs(mx_pos)
        clust_maxval(iClust) = mx_neg;
        [clust_maxchan(iClust), clust_maxfreq(iClust)] = ind2sub(size(tmp), find(tmp == mx_neg, 1));
    else
        clust_maxval(iClust) = mx_pos;
        [clust_maxchan(iClust), clust_maxfreq(iClust)] = ind2sub(size(tmp), idx_pos);
    end
end

% Sort by peak x value
[~, idx] = sort(xaxis(clust_maxfreq));
cluster_start   = cluster_start(idx);
cluster_end     = cluster_end(idx);
clust_maxchan   = clust_maxchan(idx);
clust_maxfreq   = clust_maxfreq(idx);
clust_maxval    = clust_maxval(idx);
mask_clusters   = mask_clusters(idx);
clust_bounds    = [cluster_start; cluster_end]';

% ---- Merge clusters by gap in xaxis units (sequential greedy) ----
% Merge if:
%   (1) gap <= merge_thresh, AND
%   (2) same peak channel, AND
%   (3) the t-values at the peak channel do NOT cross zero in the gap
if merge_thresh > 0
    clust_sign = sign(clust_maxval);
    i = 1;
    while i < size(clust_bounds, 1)
        % Compute gap in xaxis units between end of cluster i and start of i+1
        gap_xunits = xaxis(clust_bounds(i+1, 1)) - xaxis(clust_bounds(i, 2));

        if gap_xunits <= merge_thresh && clust_maxchan(i) == clust_maxchan(i+1)
            % Check for zero-crossing in the gap at the shared peak channel
            gap_cols = (clust_bounds(i, 2)+1):(clust_bounds(i+1, 1)-1);

            do_merge = true;
            if ~isempty(gap_cols)
                ch = clust_maxchan(i);
                gap_vals = stats(ch, gap_cols);
                % Veto if t-values cross zero relative to cluster polarity
                if clust_sign(i) ~= 0 && any(sign(gap_vals) == -clust_sign(i))
                    do_merge = false;
                end
            end

            if do_merge
                % Merge cluster i+1 into i
                clust_bounds(i, 2) = max(clust_bounds(i, 2), clust_bounds(i+1, 2));
                mask_clusters{i} = mask_clusters{i} | mask_clusters{i+1};

                % Fill the gap columns in the mask for the shared channels
                if ~isempty(gap_cols)
                    shared_chans = any(mask_clusters{i}, 2);
                    mask_clusters{i}(shared_chans, gap_cols) = true;
                end

                % Update peak to the one with larger absolute t
                [~, pick] = max(abs([clust_maxval(i), clust_maxval(i+1)]));
                if pick == 2
                    clust_maxval(i)  = clust_maxval(i+1);
                    clust_maxchan(i) = clust_maxchan(i+1);
                    clust_maxfreq(i) = clust_maxfreq(i+1);
                end
                clust_sign(i) = sign(clust_maxval(i));

                % Remove cluster i+1
                clust_bounds(i+1, :)  = [];
                mask_clusters(i+1)    = [];
                clust_maxval(i+1)     = [];
                clust_maxchan(i+1)    = [];
                clust_maxfreq(i+1)    = [];
                clust_sign(i+1)       = [];
            else
                i = i + 1;
            end
        else
            i = i + 1;
        end
    end
end

% Assign final merged variables
merged_mask_clusters = mask_clusters;
merged_bounds        = clust_bounds;
merged_maxval        = clust_maxval;
merged_maxchan       = clust_maxchan;
merged_maxfreq       = clust_maxfreq;

% Report and summarize
if merge_thresh > 0, mode_str = 'Merged'; else, mode_str = 'No merge'; end
fprintf('\nCluster Summary (%s mode):\n', mode_str);

ES_full = nan(1, numel(merged_mask_clusters));
ES_name = strings(1, numel(merged_mask_clusters));
df_used = nan(1, numel(merged_mask_clusters));

for i = 1:numel(merged_mask_clusters)
    peak_val    = xaxis(merged_maxfreq(i));
    nElectrodes = sum(any(merged_mask_clusters{i}, 2));
    t           = merged_maxval(i);

    [ES_full(i), ES_name(i), df_used(i)] = effect_size_from_t(t, grp, n, es_metric, es_opts);

    if strcmpi(datatype, 'frequency')
        fprintf('Cluster %d: %g to %.0f Hz (%g channels). Peak: %s at %.0f Hz (t = %.2f; %s = %.2f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), ...
            nElectrodes, chanlocs(merged_maxchan(i)).labels, peak_val, t, ES_name(i), ES_full(i));
    elseif strcmpi(datatype, 'nonlinear')
        fprintf('Cluster %d: scales %.0f to %.0f (%g channels). Peak: %s at scale %.0f (t = %.2f; %s = %.2f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), ...
            nElectrodes, chanlocs(merged_maxchan(i)).labels, peak_val, t, ES_name(i), ES_full(i));
    else  % time as default
        fprintf('Cluster %d: %.0f to %.0f ms (%g channels). Peak: %s at %.0f ms (t = %.2f; %s = %.2f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), ...
            nElectrodes, chanlocs(merged_maxchan(i)).labels, peak_val, t, ES_name(i), ES_full(i));
        
    end
end

% Summary table
nMerged = numel(merged_mask_clusters);
clust_summary = table;
clust_summary.Cluster          = (1:nMerged).';
clust_summary.Start            = xaxis(merged_bounds(:,1)).';
clust_summary.End              = xaxis(merged_bounds(:,2)).';
clust_summary.Peak             = xaxis(merged_maxfreq).';

clust_summary.Tvalue_full      = merged_maxval(:);
clust_summary.Tvalue           = round(merged_maxval(:), 3);

clust_summary.ES_full          = ES_full(:);
clust_summary.ES               = round(ES_full(:), 3);
clust_summary.ES_metric        = repmat(string(es_metric), nMerged, 1);
clust_summary.df_used          = df_used(:);

clust_summary.Channel          = string({chanlocs(merged_maxchan).labels}).';
clust_summary.Polarity         = strings(nMerged,1);
clust_summary.Polarity(merged_maxval > 0) = "Positive";
clust_summary.Polarity(merged_maxval < 0) = "Negative";
clust_summary.NumElectrodes    = cellfun(@(m) sum(any(m, 2)), merged_mask_clusters).';

clust_summary.Properties.VariableNames(2:4) = {'Start','End','Peak'};

% return final mask
mask = merged_mask_clusters;

end



%% Helper

function [es, name, df_out] = effect_size_from_t(t, grp, n, metric, es_opts)
% EFFECT_SIZE_FROM_T - Compute effect sizes from t-statistics
%
% FORMULAS:
%   Paired (dpt):
%     d_z = t / sqrt(n)          [Cohen's d_z for repeated measures]
%     g_z = d_z * J              [Hedges' g_z with bias correction]
%     
%   Independent (idpt):
%     d = t * sqrt(1/n1 + 1/n2)  [Cohen's d pooled]
%     g = d * J                   [Hedges' g]
%
% NOTE: For paired designs, d_z = mean(diff) / sd(diff) = t / sqrt(n)

if nargin < 5 || isempty(es_opts), es_opts = struct(); end
if ~isfield(es_opts,'df'), es_opts.df = []; end
if ~isfield(es_opts,'welch'), es_opts.welch = false; end
if ~isfield(es_opts,'s1'), es_opts.s1 = []; end
if ~isfield(es_opts,'s2'), es_opts.s2 = []; end
metric = lower(metric);

switch lower(grp)
    case 'dpt'   % paired/dependent
        n1 = n{1};
        df = n1 - 1;
        if ~isempty(es_opts.df), df = es_opts.df; end
        
        dz = t / sqrt(n1);
        J  = 1 - (3 / (4*df - 1));
        
        switch metric
            case 'd'
                es = dz;
                name = "d_z";
            case 'g'
                es = dz * J;
                name = "g_z";
            case 'eta2'
                es = t^2 / (t^2 + df);
                name = "eta2";
            case 'r'
                eta = t^2 / (t^2 + df);
                es = sign(t) * sqrt(eta);
                name = "r";
            otherwise
                error('Unknown es_metric: %s', metric);
        end
        df_out = df;

    case 'idpt'  % independent/between-subjects
        n1 = n{1}; n2 = n{2};

        if es_opts.welch
            assert(~isempty(es_opts.s1) && ~isempty(es_opts.s2), ...
                'Welch conversion requires es_opts.s1 and es_opts.s2.');
            s1 = es_opts.s1; s2 = es_opts.s2;
            mdiff = t * sqrt(s1^2/n1 + s2^2/n2);
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2));
            d  = mdiff / sp;
            df_w = welch_df(s1, s2, n1, n2);
            if ~isempty(es_opts.df), df_w = es_opts.df; end
            J  = 1 - (3 / (4*df_w - 1));
            
            switch metric
                case 'd', es = d;      name = "d";
                case 'g', es = d * J;  name = "g";
                case 'eta2', es = t^2 / (t^2 + df_w); name = "eta2";
                case 'r',    eta = t^2 / (t^2 + df_w); es = sign(t)*sqrt(eta); name = "r";
                otherwise, error('Unknown es_metric: %s', metric);
            end
            df_out = df_w;

        else
            df = n1 + n2 - 2;
            if ~isempty(es_opts.df), df = es_opts.df; end
            d  = t * sqrt(1/n1 + 1/n2);
            J  = 1 - (3 / (4*df - 1));
            
            switch metric
                case 'd', es = d;           name = "d";
                case 'g', es = d * J;       name = "g";
                case 'eta2', es = t^2 / (t^2 + df); name = "eta2";
                case 'r',    eta = t^2 / (t^2 + df); es = sign(t)*sqrt(eta); name = "r";
                otherwise, error('Unknown es_metric: %s', metric);
            end
            df_out = df;
        end

    otherwise
        error("grp must be 'dpt' or 'idpt'.")
end
end

function dfw = welch_df(s1, s2, n1, n2)
v1 = s1^2 / n1; v2 = s2^2 / n2;
dfw = (v1 + v2)^2 / ( (v1^2 / (n1 - 1)) + (v2^2 / (n2 - 1)) );
end