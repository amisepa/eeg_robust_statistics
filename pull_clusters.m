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
%                       (default = 0.5)
%   sep_thresh       - Minimum bandwidth (in xaxis units) to keep a cluster
%                       (default = 0.5)
%   polarity_split   - If true, makes polarity exclusive per x bin
%                       (default = true)
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
    sep_thresh = .5;   % in xaxis units
end
if ~exist('merge_thresh','var') || isempty(merge_thresh)
    merge_thresh = .5; % in xaxis units
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
% es_opts.df     numeric, overrides df used for eta2 and r, and J in Hedges
% es_opts.welch  true for Welch t in independent groups
% es_opts.s1,s2  group SDs (only needed to convert Welch t to d/g)


% 3) Optional polarity exclusivity per x bin
%    Enforce that at each column only one polarity can contribute to clusters
if polarity_split
    [~, idx_row] = max(abs(stats), [], 1);
    lin_idx = sub2ind(size(stats), idx_row, 1:numel(idx_row));
    sign_at_absmax = sign(stats(lin_idx));         % 1 x F in {-1,0,+1}
    sgn_pos = repmat(sign_at_absmax > 0, size(mask,1), 1);
    sgn_neg = repmat(sign_at_absmax < 0, size(mask,1), 1);
    pos_mask_all = mask & sgn_pos;
    neg_mask_all = mask & sgn_neg;

    conn = 4;                                      % 2-D connectivity
    pos_labels = bwlabel(pos_mask_all, conn);
    neg_labels = bwlabel(neg_mask_all, conn);
    if max(pos_labels(:)) > 0
        neg_labels(neg_labels > 0) = neg_labels(neg_labels > 0) + max(pos_labels(:));
    end
    label_img = pos_labels + neg_labels;
else
    % Standard labeling without polarity separation
    conn = 4;
    label_img = bwlabel(mask, conn);
end

% 4) Extract initial cluster masks
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

% 5) Split clusters by internal gaps along x using sep_thresh
dx = median(diff(xaxis));
min_gap_bins = max(1, round(sep_thresh / dx));
new_masks = {};
for iClust = 1:numel(mask_clusters)
    m = mask_clusters{iClust};

    % 1D signature across columns
    sig = any(m, 1);                      % 1 x F logical

    % Fill gaps up to sep_thresh - 1 columns if imclose exists
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

    % Label runs
    if any(sig_filled)
        [L, numRuns] = bwlabel(sig_filled);
    else
        L = zeros(size(sig_filled));
        numRuns = 0;
    end

    % Create submasks per run
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

% 6) Prune clusters whose bandwidth is below sep_thresh
min_width_bins = max(1, round(sep_thresh / dx));
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

% 7) Compute bounds and peaks
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

    % Peak with sign preserved
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

% 8) Merge clusters by gap and polarity
gap_bins = max(0, round(merge_thresh / dx));
if merge_thresh > 0
    N = size(clust_bounds, 1);
    clust_sign = sign(clust_maxval);
    overlaps_or_adjacent = @(a,b) (b(1) <= a(2)+gap_bins) && (b(2) >= a(1)-gap_bins);
    contained = @(a,b) (b(1) >= a(1) && b(2) <= a(2)) || (a(1) >= b(1) && a(2) <= b(2));

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

    if any(A(:))
        G = graph(A);
        comp = conncomp(G);
    else
        comp = 1:N;
    end

    K = max(comp);
    merged_mask_clusters = cell(1,K);
    merged_bounds  = nan(K,2);
    merged_maxval  = nan(1,K);
    merged_maxchan = nan(1,K);
    merged_maxfreq = nan(1,K);

    for k = 1:K
        idk = find(comp == k);
        m = false(size(mask_clusters{1}));
        for t = idk(:)'
            m = m | mask_clusters{t};
        end
        merged_mask_clusters{k} = m;

        merged_bounds(k,1) = min(clust_bounds(idk,1));
        merged_bounds(k,2) = max(clust_bounds(idk,2));

        [~, ii] = max(abs(clust_maxval(idk)));
        pick = idk(ii);
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

% 9) Report and summarize
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
        fprintf('Cluster %d: %g to %g Hz (%g channels). Peak: %s at %g Hz (t = %.3f; %s = %.3f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), ...
            nElectrodes, chanlocs(merged_maxchan(i)).labels, peak_val, t, ES_name(i), ES_full(i));
    else
        fprintf('Cluster %d: %g to %g ms (%g channels). Peak: %s at %g ms (t = %.3f; %s = %.3f)\n', ...
            i, xaxis(merged_bounds(i,1)), xaxis(merged_bounds(i,2)), ...
            nElectrodes, chanlocs(merged_maxchan(i)).labels, peak_val, t, ES_name(i), ES_full(i));
    end
end

% 10) Summary table
nMerged = numel(merged_mask_clusters);
clust_summary = table;
clust_summary.Cluster          = (1:nMerged).';
clust_summary.Start            = xaxis(merged_bounds(:,1)).';
clust_summary.End              = xaxis(merged_bounds(:,2)).';
clust_summary.Peak             = xaxis(merged_maxfreq).';

% store full precision and rounded t
clust_summary.Tvalue_full      = merged_maxval(:);
clust_summary.Tvalue           = round(merged_maxval(:), 3);

% effect sizes
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
% metric: 'd' | 'g' | 'eta2' | 'r'
% When to use:
%   'd'   use for paired dz or independent d under equal variances
%   'g'  small-sample bias corrected d (paired or independent)
%   'eta2'      variance-explained for a single t test
%   'r'         correlation effect size for a single t test
if nargin < 5 || isempty(es_opts), es_opts = struct(); end
if ~isfield(es_opts,'df'), es_opts.df = []; end
if ~isfield(es_opts,'welch'), es_opts.welch = false; end
if ~isfield(es_opts,'s1'), es_opts.s1 = []; end
if ~isfield(es_opts,'s2'), es_opts.s2 = []; end
metric = lower(metric);

switch lower(grp)
    case 'dpt'   % paired
        n1 = n{1};
        df = n1 - 1;
        if ~isempty(es_opts.df), df = es_opts.df; end
        dz = t / sqrt(n1);               % Cohen dz from paired t
        J  = 1 - (3 / (4*df - 1));       % bias correction
        switch metric
            case 'd', es = dz;        name = "d";
            case 'g', es = dz * J;    name = "g";
            case 'eta2',     es = t^2 / (t^2 + df); name = "eta2";
            case 'r',        eta = t^2 / (t^2 + df); es = sign(t)*sqrt(eta); name = "r";
            otherwise, error('Unknown es_metric: %s', metric);
        end
        df_out = df;

    case 'idpt'  % independent
        n1 = n{1}; n2 = n{2};

        if es_opts.welch
            % Welch t path: need group SDs s1, s2 to back out the mean diff
            assert(~isempty(es_opts.s1) && ~isempty(es_opts.s2), ...
                'Welch conversion requires es_opts.s1 and es_opts.s2.');
            s1 = es_opts.s1; s2 = es_opts.s2;
            % Back out mean difference from t
            mdiff = t * sqrt(s1^2/n1 + s2^2/n2);
            % Pooled SD for standardizing d (Hedges suggested)
            sp = sqrt(((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2));
            d  = mdiff / sp;
            % Welch df for eta2 and r (or override)
            df_w = welch_df(s1, s2, n1, n2);
            if ~isempty(es_opts.df), df_w = es_opts.df; end
            J  = 1 - (3 / (4*df_w - 1));   % bias correction using Welch df
            switch metric
                case 'd', es = d;      name = "d";
                case 'g', es = d * J;  name = "g";
                case 'eta2',     es = t^2 / (t^2 + df_w); name = "eta2";
                case 'r',        eta = t^2 / (t^2 + df_w); es = sign(t)*sqrt(eta); name = "r";
                otherwise, error('Unknown es_metric: %s', metric);
            end
            df_out = df_w;

        else
            % Equal-variance pooled t
            df = n1 + n2 - 2;
            if ~isempty(es_opts.df), df = es_opts.df; end
            d  = t * sqrt(1/n1 + 1/n2);   % Cohen d from pooled t
            J  = 1 - (3 / (4*df - 1));
            switch metric
                case 'd', es = d;           name = "d";
                case 'g', es = d * J;      name = "g";
                case 'eta2',     es = t^2 / (t^2 + df); name = "eta2";
                case 'r',        eta = t^2 / (t^2 + df); es = sign(t)*sqrt(eta); name = "r";
                otherwise, error('Unknown es_metric: %s', metric);
            end
            df_out = df;
        end

    otherwise
        error("grp must be 'dpt' or 'idpt'.")
end
end

function dfw = welch_df(s1, s2, n1, n2)
% Welch–Satterthwaite degrees of freedom
v1 = s1^2 / n1; v2 = s2^2 / n2;
dfw = (v1 + v2)^2 / ( (v1^2 / (n1 - 1)) + (v2^2 / (n2 - 1)) );
end
