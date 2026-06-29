function [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_thresh, sep_thresh, polarity_split, es_metric)
% PULL_CLUSTERS Extract significant clusters after MCC. 
% 
%   Identifies labeled clusters from a binary significance mask, evaluates
%   their spatial/temporal/frequency bounds, polarity, and peak effect.
%   Optionally merges adjacent or overlapping clusters that share consistent
%   polarity. Returns logical masks, cluster characteristics, ES, and a 
%   structured cluster summary.
%
%   When xaxis has only one value (uniscale data, e.g., a single entropy
%   measure across channels), clustering is performed in channel space only.
%   Clusters are then contiguous sets of electrodes (using 1D connectivity),
%   and the summary reports channel lists rather than x-axis bounds.
%
%   For 2D [channels x time/freq] data, polarity_split operates POST-HOC:
%   clusters are first extracted from the full binary mask, then each
%   cluster is assigned a polarity label based on the sign of its peak
%   t-value. This avoids fragmentation caused by pre-hoc column-wise
%   polarity gating, which fails when opposite-polarity effects coexist
%   at different scalp regions within the same frequency/time range.
%
% USAGE:
%   [mask, clust_summary] = pull_clusters(mask, stats, xaxis, chanlocs, datatype, grp, n, merge_thresh, sep_thresh, polarity_split, es_metric)
%
% INPUTS:
%   mask             - Binary significance mask [channels x time/freq]
%   stats            - Observed t-values or test statistic map [same size as mask]
%   xaxis            - X-axis values (e.g., frequency in Hz or time in ms)
%   chanlocs         - EEGLAB-style channel location structure
%   datatype         - 'frequency', 'time', 'nonlinear', or 'scalar'
%   grp              - 'dpt' (within-group) or 'idpt' (between-group)
%   n                - Cell array with subject counts: {n1} or {n1, n2}
%   merge_thresh     - Maximum gap in xaxis units (ms or Hz) for merging 
%                       clusters (default = 0, no merging)
%   sep_thresh       - Minimum width in bins/frames to keep a cluster
%                       (default = 0)
%   polarity_split   - For uniscale data: if true, splits positive and
%                       negative channel clusters before labeling.
%                       For 2D data: polarity is always assigned post-hoc
%                       from each cluster's peak t-value sign; this flag
%                       has no effect on the 2D path (default = true).
%   es_metric        - Effect size (ES) metric to compute: 'd' (cohen's d; default),
%                       'g (hedge g)', 'eta2', 'r'
%
% OUTPUTS:
%   mask             - Cell array of logical masks for each merged cluster
%   clust_summary    - Table summarizing each cluster (bounds, peaks, effect sizes)
%
% Author: Cedric Cannard, 2021-2025

% 1) Orient inputs: function expects [nChan x nScales] where nChan == numel(chanlocs).
%    For uniscale data this means [nChan x 1]. If the input arrives transposed
%    (e.g. [1 x nChan] from a row-vector t-test output), flip it here.
nChan = numel(chanlocs);
if size(stats,1) ~= nChan && size(stats,2) == nChan
    stats = stats';
    mask  = mask';
end

% Early exit if nothing is significant
if nnz(mask) == 0
    disp('No significant differences, nothing to plot.')
    mask = {};
    clust_summary = table;
    return
end

% 2) Defaults for thresholds
if ~exist('sep_thresh','var') || isempty(sep_thresh)
    sep_thresh = 0;
end
if ~exist('merge_thresh','var') || isempty(merge_thresh)
    merge_thresh = 0;
end
if ~exist('polarity_split','var') || isempty(polarity_split)
    polarity_split = true;
end
if ~exist('es_metric','var') || isempty(es_metric)
    es_metric = 'd';
end
if ~exist('es_opts','var') || isempty(es_opts)
    es_opts = struct();
end
if strcmpi(es_metric, 'r')
    es_opts.df = n{:} - 2;
end

% =========================================================================
% UNISCALE PATH: single value per channel [nChan x 1].
% Detected when xaxis is empty/scalar, datatype is 'scalar', or after
% orientation normalization above stats has only one column.
% Multiscale data [nChan x nScales] goes to the standard path.
% =========================================================================
isUniscale = size(stats, 2) == 1 || isempty(xaxis) || ...
             numel(xaxis) == 1   || strcmpi(datatype, 'scalar');

if isUniscale

    sig_vec = mask(:, 1);   % [nChan x 1] logical
    t_vec   = stats(:, 1);  % [nChan x 1] t-values

    if ~any(sig_vec)
        disp('No significant channels.')
        mask = {};
        clust_summary = table;
        return
    end

    % Optional polarity split along channel dimension (meaningful here
    % because there is no second spatial dimension to absorb the conflict)
    if polarity_split
        t_pos = sig_vec & (t_vec > 0);
        t_neg = sig_vec & (t_vec < 0);
        [lbl_pos, np] = bwlabel(t_pos, 4);
        [lbl_neg, nn] = bwlabel(t_neg, 4);
        lbl_neg(lbl_neg > 0) = lbl_neg(lbl_neg > 0) + np;
        lbl_vec = lbl_pos + lbl_neg;
        n_raw = np + nn;
    else
        [lbl_vec, n_raw] = bwlabel(sig_vec, 4);
    end

    % Build cell array of channel-index vectors per cluster
    chan_clusters = cell(1, n_raw);
    for k = 1:n_raw
        chan_clusters{k} = find(lbl_vec == k);
    end
    chan_clusters = chan_clusters(~cellfun(@isempty, chan_clusters));

    if isempty(chan_clusters)
        disp('No spatial clusters found.')
        mask = {};
        clust_summary = table;
        return
    end

    % Prune clusters smaller than sep_thresh channels
    keep = cellfun(@(c) numel(c) >= max(1, round(sep_thresh)), chan_clusters);
    chan_clusters = chan_clusters(keep);

    if isempty(chan_clusters)
        disp('All spatial clusters pruned by sep_thresh.')
        mask = {};
        clust_summary = table;
        return
    end

    % Find peak channel (largest |t|) per cluster
    n_clust = numel(chan_clusters);
    peak_chan   = nan(1, n_clust);
    peak_tval   = nan(1, n_clust);
    clust_masks = cell(1, n_clust);

    for k = 1:n_clust
        chans = chan_clusters{k};
        t_sub = t_vec(chans);
        [~, pk] = max(abs(t_sub));
        peak_chan(k) = chans(pk);
        peak_tval(k) = t_sub(pk);

        m = false(size(mask));       % [nChan x 1]
        m(chans, 1) = true;
        clust_masks{k} = m;
    end

    % Sort clusters by peak channel index
    [~, sidx] = sort(peak_chan);
    peak_chan     = peak_chan(sidx);
    peak_tval     = peak_tval(sidx);
    clust_masks   = clust_masks(sidx);
    chan_clusters = chan_clusters(sidx);

    % Report and build summary table
    fprintf('\nSpatial Cluster Summary (uniscale, %s):\n', es_metric);

    ES_full = nan(1, n_clust);
    ES_name = strings(1, n_clust);
    df_used = nan(1, n_clust);

    for k = 1:n_clust
        [ES_full(k), ES_name(k), df_used(k)] = effect_size_from_t(peak_tval(k), grp, n, es_metric, es_opts);
        fprintf('Cluster %d: %d channels. Peak: %s (t = %.2f; %s = %.2f)\n', ...
            k, numel(chan_clusters{k}), ...
            chanlocs(peak_chan(k)).labels, peak_tval(k), ES_name(k), ES_full(k));
    end

    % Build summary table
    clust_summary = table;
    clust_summary.Cluster       = (1:n_clust).';
    clust_summary.NumElectrodes = cellfun(@numel, chan_clusters).';
    clust_summary.PeakChannel   = string({chanlocs(peak_chan).labels}).';
    clust_summary.Tvalue_full   = peak_tval(:);
    clust_summary.Tvalue        = round(peak_tval(:), 3);
    clust_summary.ES_full       = ES_full(:);
    clust_summary.ES            = round(ES_full(:), 3);
    clust_summary.ES_metric     = repmat(string(es_metric), n_clust, 1);
    clust_summary.df_used       = df_used(:);
    clust_summary.Polarity      = strings(n_clust, 1);
    clust_summary.Polarity(peak_tval > 0) = "Positive";
    clust_summary.Polarity(peak_tval < 0) = "Negative";

    mask = clust_masks;
    return
end

% =========================================================================
% STANDARD PATH: 2D [channels x time/freq] data.
%
% Polarity is handled POST-HOC: bwlabel runs on the full binary mask so
% that spatially distinct opposite-polarity effects occupying the same
% frequency/time columns do not fragment into thin alternating slices.
% Each cluster's polarity is then assigned from the sign of its peak
% t-value after all merging and pruning is complete.
% =========================================================================

% Validate xaxis length against stats columns
nCols = size(stats, 2);
if numel(xaxis) ~= nCols
    warning('pull_clusters: xaxis length (%d) does not match stats columns (%d). Reindexing to bin indices.', numel(xaxis), nCols);
    xaxis = 1:nCols;
end

% Label connected components on the full binary mask (no pre-hoc polarity split)
conn      = 4;
label_img = bwlabel(mask, conn);

% Extract initial cluster masks
n_cluster    = max(label_img(:));
mask_clusters = cell(1, 0);
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

% Consolidate clusters that overlap in time/freq columns
changed = true;
while changed
    changed = false;
    for ii = 1:numel(mask_clusters)-1
        if ~any(mask_clusters{ii}(:)), continue, end
        sig_i = any(mask_clusters{ii}, 1);
        for jj = ii+1:numel(mask_clusters)
            if ~any(mask_clusters{jj}(:)), continue, end
            sig_j = any(mask_clusters{jj}, 1);
            if any(sig_i & sig_j)  % temporal overlap
                % Only consolidate if peak polarities match
                tmp_i = stats .* mask_clusters{ii}; tmp_i(~mask_clusters{ii}) = NaN;
                tmp_j = stats .* mask_clusters{jj}; tmp_j(~mask_clusters{jj}) = NaN;
                [mx_i_pos, ~] = max(tmp_i(:), [], 'omitnan');
                [mx_i_neg, ~] = min(tmp_i(:), [], 'omitnan');
                peak_i = (abs(mx_i_neg) > abs(mx_i_pos)) * mx_i_neg + ...
                         (abs(mx_i_pos) >= abs(mx_i_neg)) * mx_i_pos;
                [mx_j_pos, ~] = max(tmp_j(:), [], 'omitnan');
                [mx_j_neg, ~] = min(tmp_j(:), [], 'omitnan');
                peak_j = (abs(mx_j_neg) > abs(mx_j_pos)) * mx_j_neg + ...
                         (abs(mx_j_pos) >= abs(mx_j_neg)) * mx_j_pos;
                if sign(peak_i) ~= sign(peak_j)
                    continue  % opposite polarity, keep separate
                end
                mask_clusters{ii} = mask_clusters{ii} | mask_clusters{jj};
                mask_clusters{jj} = false(size(mask_clusters{jj}));
                sig_i = any(mask_clusters{ii}, 1);
                changed = true;
            end
        end
    end
end

% Split clusters by internal gaps along x using sep_thresh
min_gap_bins = max(1, round(sep_thresh));
new_masks = {};
for iClust = 1:numel(mask_clusters)
    m   = mask_clusters{iClust};
    sig = any(m, 1);

    if exist('imclose','file') == 2
        se_len = max(0, min_gap_bins - 1);
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
        L       = zeros(size(sig_filled));
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
    first_col  = find(sigcols, 1, 'first');
    last_col   = find(sigcols, 1, 'last');
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
n_cluster     = numel(mask_clusters);
cluster_start = nan(1, n_cluster);
cluster_end   = nan(1, n_cluster);
clust_maxval  = nan(1, n_cluster);
clust_maxchan = nan(1, n_cluster);
clust_maxfreq = nan(1, n_cluster);

for iClust = 1:n_cluster
    m   = mask_clusters{iClust};
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
[~, idx]      = sort(xaxis(clust_maxfreq));
cluster_start = cluster_start(idx);
cluster_end   = cluster_end(idx);
clust_maxchan = clust_maxchan(idx);
clust_maxfreq = clust_maxfreq(idx);
clust_maxval  = clust_maxval(idx);
mask_clusters = mask_clusters(idx);
clust_bounds  = [cluster_start; cluster_end]';

% ---- Merge clusters by gap in xaxis units (sequential greedy) ----
if merge_thresh > 0
    clust_sign = sign(clust_maxval);
    i = 1;
    while i < size(clust_bounds, 1)
        gap_xunits = xaxis(clust_bounds(i+1, 1)) - xaxis(clust_bounds(i, 2));

        if gap_xunits <= merge_thresh && clust_sign(i) == clust_sign(i+1)
            % Same polarity and within gap threshold: candidate for merge
            gap_cols = (clust_bounds(i, 2)+1):(clust_bounds(i+1, 1)-1);

            do_merge = true;
            if ~isempty(gap_cols)
                [~, pick] = max(abs([clust_maxval(i), clust_maxval(i+1)]));
                ch_ref   = clust_maxchan([i, i+1]);
                ch_ref   = ch_ref(pick);
                gap_vals = stats(ch_ref, gap_cols);
                if clust_sign(i) ~= 0 && any(sign(gap_vals) == -clust_sign(i))
                    do_merge = false;
                end
            end

            if do_merge
                chans_i   = any(mask_clusters{i},   2);
                chans_ip1 = any(mask_clusters{i+1}, 2);

                clust_bounds(i, 2)  = max(clust_bounds(i, 2), clust_bounds(i+1, 2));
                mask_clusters{i}    = mask_clusters{i} | mask_clusters{i+1};

                if ~isempty(gap_cols)
                    shared_chans = chans_i & chans_ip1;
                    mask_clusters{i}(shared_chans, gap_cols) = true;
                end

                [~, pick] = max(abs([clust_maxval(i), clust_maxval(i+1)]));
                if pick == 2
                    clust_maxval(i)  = clust_maxval(i+1);
                    clust_maxchan(i) = clust_maxchan(i+1);
                    clust_maxfreq(i) = clust_maxfreq(i+1);
                end
                clust_sign(i) = sign(clust_maxval(i));

                clust_bounds(i+1, :) = [];
                mask_clusters(i+1)   = [];
                clust_maxval(i+1)    = [];
                clust_maxchan(i+1)   = [];
                clust_maxfreq(i+1)   = [];
                clust_sign(i+1)      = [];
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
clust_summary.Cluster       = (1:nMerged).';
clust_summary.Start         = xaxis(merged_bounds(:,1)).';
clust_summary.End           = xaxis(merged_bounds(:,2)).';
clust_summary.Peak          = xaxis(merged_maxfreq).';

clust_summary.Tvalue_full   = merged_maxval(:);
clust_summary.Tvalue        = round(merged_maxval(:), 3);

clust_summary.ES_full       = ES_full(:);
clust_summary.ES            = round(ES_full(:), 3);
clust_summary.ES_metric     = repmat(string(es_metric), nMerged, 1);
clust_summary.df_used       = df_used(:);

clust_summary.Channel       = string({chanlocs(merged_maxchan).labels}).';
clust_summary.Polarity      = strings(nMerged, 1);
clust_summary.Polarity(merged_maxval > 0) = "Positive";
clust_summary.Polarity(merged_maxval < 0) = "Negative";
clust_summary.NumElectrodes = cellfun(@(m) sum(any(m, 2)), merged_mask_clusters).';

clust_summary.Properties.VariableNames(2:4) = {'Start', 'End', 'Peak'};

% Return final mask
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
if ~isfield(es_opts, 'df'),    es_opts.df    = []; end
if ~isfield(es_opts, 'welch'), es_opts.welch = false; end
if ~isfield(es_opts, 's1'),    es_opts.s1    = []; end
if ~isfield(es_opts, 's2'),    es_opts.s2    = []; end
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
                es   = dz;
                name = "d_z";
            case 'g'
                es   = dz * J;
                name = "g_z";
            case 'eta2'
                es   = t^2 / (t^2 + df);
                name = "eta2";
            case 'r'
                eta  = t^2 / (t^2 + df);
                es   = sign(t) * sqrt(eta);
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
            s1   = es_opts.s1; s2 = es_opts.s2;
            mdiff = t * sqrt(s1^2/n1 + s2^2/n2);
            sp    = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2));
            d     = mdiff / sp;
            df_w  = welch_df(s1, s2, n1, n2);
            if ~isempty(es_opts.df), df_w = es_opts.df; end
            J     = 1 - (3 / (4*df_w - 1));

            switch metric
                case 'd',    es = d;           name = "d";
                case 'g',    es = d * J;        name = "g";
                case 'eta2', es = t^2 / (t^2 + df_w); name = "eta2";
                case 'r',    eta = t^2 / (t^2 + df_w); es = sign(t)*sqrt(eta); name = "r";
                otherwise,   error('Unknown es_metric: %s', metric);
            end
            df_out = df_w;

        else
            df = n1 + n2 - 2;
            if ~isempty(es_opts.df), df = es_opts.df; end
            d  = t * sqrt(1/n1 + 1/n2);
            J  = 1 - (3 / (4*df - 1));

            switch metric
                case 'd',    es = d;           name = "d";
                case 'g',    es = d * J;        name = "g";
                case 'eta2', es = t^2 / (t^2 + df); name = "eta2";
                case 'r',    eta = t^2 / (t^2 + df); es = sign(t)*sqrt(eta); name = "r";
                otherwise,   error('Unknown es_metric: %s', metric);
            end
            df_out = df;
        end

    otherwise
        error("grp must be 'dpt' or 'idpt'.")
end
end

function dfw = welch_df(s1, s2, n1, n2)
v1  = s1^2 / n1; v2 = s2^2 / n2;
dfw = (v1 + v2)^2 / ( (v1^2 / (n1 - 1)) + (v2^2 / (n2 - 1)) );
end