function hs = plot_clusters(summary_tbl, mask, tvals, map, f, chanlocs, line_label, varargin)
% PLOT_CLUSTERS  Scalp or cortex plot at cluster peak + frequency curve.
%
% Required
%   summary_tbl  table with columns Peak_Hz or Peak, and Channel
%   mask         cell {nClust}, each logical [nChan x nFreq]
%   tvals        [nChan x nFreq] t-values map (used for CI and cortex coloring)
%   map          GLM: [nChan x nFreq] beta map
%                PSD: [nChan x nFreq x nSub] per-subject values
%   f            [1 x nFreq] numeric frequency/scale vector
%   chanlocs     EEGLAB chanlocs struct
%   line_label   y-axis label for the curve
%
% Name-value
%   'LineNoiseHz'   scalar or vector freqs to mask in curves (default 60). [] disables
%   'LineNoiseBW'   bandwidth (Hz) around each line-noise freq (default 1)
%   'DataType'      'scalp' (default) or 'source'
%   'XTickLabels'   cellstr/strings for x-axis labels on the curve plot; if empty, numeric f is shown
%   'Domain'        'time' (ms) | 'frequency' (Hz) | 'nonlinear' (scale factors, default)
%
% Notes
%   - Axis computations always use numeric f
%   - XTick labels are applied only to the curve plot and never used for math

ip = inputParser;
addParameter(ip, 'LineNoiseHz', 60);
addParameter(ip, 'LineNoiseBW', 1);
addParameter(ip, 'DataType',    'scalp');
addParameter(ip, 'XTickLabels', []);
addParameter(ip, 'Domain',      'nonlinear');
parse(ip, varargin{:});
lnHz        = ip.Results.LineNoiseHz;
lnBW        = ip.Results.LineNoiseBW;
dataType    = ip.Results.DataType;
xTickLabels = ip.Results.XTickLabels;
domain      = lower(ip.Results.Domain);

% Unit strings derived from domain
switch domain
    case 'time'
        unit_label = 'ms';
        peak_fmt   = 'Peak: %.0f ms';
        xlabel_str = 'Time (ms)';
    case 'frequency'
        unit_label = 'Hz';
        peak_fmt   = 'Peak: %.2f Hz';
        xlabel_str = 'Frequency (Hz)';
    case 'nonlinear'
        unit_label = 'scale factors';
        peak_fmt   = 'Scale factor: %g';
        xlabel_str = 'Scale factors';
    otherwise
        error('plot_clusters: Domain must be ''time'', ''frequency'', or ''nonlinear''.');
end

% Basic sizes
[nChan, nFreq] = size(tvals);
assert(isnumeric(f) && numel(f) == nFreq, 'f must be numeric and match number of columns');

hs.topo  = cell(height(summary_tbl), 1);
hs.curve = cell(height(summary_tbl), 1);

% t-range for cortex coloring logic
if ~isempty(tvals)
    tmin = min(tvals(:), [], 'omitnan');
    tmax = max(tvals(:), [], 'omitnan');
else
    tmin = NaN; tmax = NaN;
end

for iClust = 1:height(summary_tbl)

    %% --- Peak frequency and channel ---
    if ismember('Peak', summary_tbl.Properties.VariableNames)
        peak_freq = summary_tbl.Peak(iClust);
    elseif ismember('Peak_Hz', summary_tbl.Properties.VariableNames)
        peak_freq = summary_tbl.Peak_Hz(iClust);
    else
        error('summary_tbl must have a Peak or Peak_Hz column.');
    end

    if ~isnumeric(peak_freq) || ~isscalar(peak_freq) || ~isfinite(peak_freq)
        warning('plot_clusters:BadPeak', 'Invalid Peak for cluster %d; using mid-bin.', iClust);
        fi = round(nFreq / 2);
    else
        [~, fi] = min(abs(f - peak_freq));
    end

    ch = summary_tbl.Channel{iClust};
    if isstring(ch), ch = char(ch); end

    %% --- Line-noise mask ---
    lnMask = false(1, nFreq);
    if ~isempty(lnHz)
        lnHz   = lnHz(:).';
        halfBW = max(lnBW, eps) / 2;
        for v = lnHz
            lnMask = lnMask | (abs(f - v) <= halfBW);
        end
    end

    %% --- Channel index for curve ---
    ci_curve = find(strcmpi({chanlocs.labels}, ch), 1);
    if isempty(ci_curve)
        warning('plot_clusters:CurveChannelMissing', ...
            'Curve skipped for cluster %d because channel "%s" was not found.', iClust, ch);
        hs.curve{iClust} = [];
        continue
    end

    %% --- Compute y, lo, hi ---
    sz = size(map);
    if numel(sz) == 2 && ~isempty(tvals) && all(sz == size(tvals))
        % GLM: beta with CI derived from t-values
        b_row  = map(ci_curve, :);
        t_row  = tvals(ci_curve, :);
        se_row = abs(b_row ./ t_row);
        se_row(~isfinite(se_row)) = NaN;
        y  = b_row;
        lo = b_row - 1.96 * se_row;
        hi = b_row + 1.96 * se_row;
    elseif numel(sz) == 3 && sz(1) == nChan && sz(2) == nFreq
        % PSD: mean across subjects with 95% CI from SE
        rows = squeeze(map(ci_curve, :, :)).';   % [nSub x nFreq]
        y    = mean(rows, 1, 'omitmissing');
        s    = std(rows,  0, 1, 'omitmissing');
        n    = sum(isfinite(rows), 1);
        se   = s ./ sqrt(max(n, 1));
        lo   = y - 1.96 * se;
        hi   = y + 1.96 * se;
    else
        if numel(sz) ~= 2 || ~all(sz == [nChan nFreq])
            error('map must be [nChan x nFreq] (GLM) or [nChan x nFreq x nSub] (PSD).');
        end
        % Fallback: no CI
        y  = map(ci_curve, :);
        lo = []; hi = [];
    end

    % Apply line-noise mask
    y(lnMask) = NaN;
    if ~isempty(lo), lo(lnMask) = NaN; end
    if ~isempty(hi), hi(lnMask) = NaN; end

    %% --- Significance mask for this channel ---
    sig_mask = false(1, nFreq);
    if iClust <= numel(mask) && ~isempty(mask{iClust})
        mc = mask{iClust};
        if size(mc, 1) >= ci_curve && size(mc, 2) == nFreq
            sig_mask = mc(ci_curve, :);
        end
    end
    sig_mask = sig_mask & ~lnMask;

    %% ===== CURVE PLOT =====
    hCurve = figure('Color', 'w', 'ToolBar', 'none', 'MenuBar', 'none');
    hs.curve{iClust} = hCurve;
    hold on
    set(gca, 'Color', 'w');  

    % Shaded CI
    if ~isempty(lo) && ~isempty(hi)
        good = isfinite(lo) & isfinite(hi);
        if any(good)
            dci    = [false, good, false];
            starts = find(diff(dci) ==  1);
            ends   = find(diff(dci) == -1) - 1;
            for s = 1:numel(starts)
                idx = starts(s):ends(s);
                fill([f(idx) fliplr(f(idx))], [lo(idx) fliplr(hi(idx))], ...
                    [0.4660 0.6740 0.1880], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            end
        end
    end

    % Baseline
    plot([f(1) f(end)], [0 0], '--', 'Color', 'k', 'LineWidth', 1.2);

    % Main line
    plot(f, y, 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2);

    % Significance bars along the bottom
    yl = ylim;
    H  = 0.015 * (yl(2) - yl(1));
    y0 = yl(1) - 1.6 * H;
    ylim([y0 yl(2)]);
    if numel(f) > 1
        d1    = f(2) - f(1);
        dN    = f(end) - f(end-1);
        edges = [f(1)-d1/2, 0.5*(f(1:end-1)+f(2:end)), f(end)+dN/2];
    else
        edges = [f(1)-0.5, f(1)+0.5];
    end
    dsig   = [false, sig_mask(:).', false];
    starts = find(diff(dsig) ==  1);
    ends   = find(diff(dsig) == -1);
    for k = 1:numel(starts)
        xs = edges(starts(k));
        xe = edges(ends(k));
        rectangle('Position', [xs, y0, max(eps, xe-xs), H], ...
            'FaceColor', 'k', 'EdgeColor', 'none');
    end

    % Optional x tick labels
    if ~isempty(xTickLabels)
        labs = string(xTickLabels(:).');
        if numel(labs) ~= nFreq
            stride = max(1, round(nFreq / 10));
            keep   = false(1, nFreq);
            keep(1:stride:end) = true;
            keep(end) = true;
            xticks(f(keep));
            xticklabels(labs(1:min(sum(keep), numel(labs))));
        else
            xticks(f);
            xticklabels(labs);
        end
        xtickangle(45);
        set(gca, 'TickLabelInterpreter', 'none');
    end

    xlabel(xlabel_str);
    ylabel(line_label);
    title(sprintf('Cluster %d - %s  (%s)', iClust, ch, unit_label));
    box on; grid off; axis tight

    % Apply formatting via explicit figure handle (avoids gcf ambiguity)
    set(findall(hCurve, 'type', 'axes'), 'FontSize', 14, 'FontWeight', 'bold', ...
        'XColor', 'k', 'YColor', 'k');
    set(get(gca, 'Title'), 'Color', 'k', 'FontSize', 16, 'FontWeight', 'bold');

    %% ===== TOPO / BRAIN PLOT =====
    if strcmpi(dataType, 'source')

        % 3-D cortex plot
        load cm17.mat
        load cortex; cortex = cortex_highres;  %#ok<NODEF>

        if tmin < 0 && tmax > 0
            allplots_cortex_BS(cortex, mean(tvals, 2, 'omitnan'), ...
                [-max(abs(tvals(:))) max(abs(tvals(:)))], cm17,  't-values', []);
        elseif tmax <= 0
            allplots_cortex_BS(cortex, mean(tvals, 2, 'omitnan'), ...
                [min(mean(tvals,2,'omitnan')) tmax], cm17b, 't-values', []);  %#ok<NODEF>
        elseif tmin >= 0
            allplots_cortex_BS(cortex, mean(tvals, 2, 'omitnan'), ...
                [tmin max(mean(tvals,2,'omitnan'))], cm17a, 't-values', []);  %#ok<NODEF>
        end
        hTopo = gcf;
        hs.topo{iClust} = hTopo;

    else

        % Scalp topography
        rdbu_cmap = interp1([1;128;256], ...
            [0.698 0.094 0.168; 1 1 1; 0.129 0.400 0.675], (1:256)', 'linear');
        dmap = flipud(max(0, min(1, rdbu_cmap)));

        topo_mask = false(nChan, 1);
        if iClust <= numel(mask) && ~isempty(mask{iClust})
            mc = mask{iClust};
            if size(mc, 2) == nFreq
                topo_mask = any(mc, 2);
            end
        end

        ci = find(strcmpi({chanlocs.labels}, ch), 1);
        if isempty(ci)
            warning('plot_clusters:MissingChan', ...
                'Skipping topo for cluster %d, channel %s not found.', iClust, ch);
            hs.topo{iClust} = [];
            continue
        end

        hTopo = figure('Color', 'w', 'ToolBar', 'none', 'MenuBar', 'none');
        hs.topo{iClust} = hTopo;

        topoplot(tvals(:, fi), chanlocs, 'colormap', dmap, 'pmask', topo_mask, ...
            'verbose', 'off', 'whitebk', 'on');
        c = colorbar;
        ylabel(c, 't-values', 'FontWeight', 'bold', 'FontSize', 12);

    end

    % Topo title and formatting via explicit figure handle
    figure(hTopo);
    title(sprintf(['Cluster %d (' peak_fmt ')'], iClust, f(fi)));
    set(findall(hTopo, 'type', 'axes'), 'FontSize', 14, 'FontWeight', 'bold', ...
        'XColor', 'k', 'YColor', 'k');
    set(get(gca, 'Title'), 'Color', 'k', 'FontSize', 16, 'FontWeight', 'bold');

end
end