function hs = plot_clusters(summary_tbl, mask_clusters, topo_map, curve_data, f, chanlocs, dmap, line_label, varargin)
% PLOT_CLUSTERS  Topography at peak + frequency curve (GLM or PSD) with your visual style.
%
% Required
%   summary_tbl    table with columns Peak (Hz) and Channel
%   mask_clusters  cell {nClust}, each logical [nChan x nFreq]
%   topo_map       [nChan x nFreq] values for topoplot at peak frequency
%   curve_data     GLM: [nChan x nFreq] beta map      OR
%                  PSD: [nChan x nFreq x nSub] per-subject values
%   f              [1 x nFreq] frequency vector
%   chanlocs       EEGLAB chanlocs struct
%   dmap           colormap for topoplot
%   line_label     y-axis label for the curve
%
% Name-value (minimal)
%   'TMap'         [nChan x nFreq] t map (only for GLM beta CI). Default []
%   'LineNoiseHz'  scalar or vector freqs to mask (default 60). [] disables
%   'LineNoiseBW'  bandwidth (Hz) around each line-noise freq (default 1)

ip = inputParser;
addParameter(ip,'TMap',[]);
addParameter(ip,'LineNoiseHz',60);
addParameter(ip,'LineNoiseBW',1);
parse(ip,varargin{:});
tmap = ip.Results.TMap;
lnHz = ip.Results.LineNoiseHz;
lnBW = ip.Results.LineNoiseBW;

% basic sizes
[nChan, nFreq] = size(topo_map);
assert(numel(f)==nFreq, 'f must match number of columns in topo_map');

hs.topo  = cell(height(summary_tbl),1);
hs.curve = cell(height(summary_tbl),1);

% line-noise mask
lnMask = false(1, nFreq);
if ~isempty(lnHz)
    lnHz = lnHz(:).';
    halfBW = max(lnBW, eps)/2;
    for v = lnHz
        lnMask = lnMask | (abs(f - v) <= halfBW);
    end
end

for iClust = 1:height(summary_tbl)
    % Peak freq and channel
    peak_freq = summary_tbl.Peak(iClust);
    [~, fi] = min(abs(f - peak_freq));
    ch = summary_tbl.Channel{iClust};
    if isstring(ch), ch = char(ch); end
    ci = find(strcmpi({chanlocs.labels}, ch), 1);
    if isempty(ci)
        warning('plot_clusters:MissingChan','Skipping cluster %d, channel %s not found.', iClust, ch);
        continue
    end

    % ----- Topography at peak freq -----
    topo_mask = false(nChan,1);
    if iClust <= numel(mask_clusters) && ~isempty(mask_clusters{iClust})
        mc = mask_clusters{iClust};
        if size(mc,2) == nFreq, topo_mask = mc(:, fi); end
    end

    hs.topo{iClust} = figure('Color','w');
    try
        topoplot(topo_map(:, fi), chanlocs, 'pmask', topo_mask, ...
                 'verbose','off','colormap', dmap, 'whitebk','on');
    catch
        topoplot(topo_map(:, fi), chanlocs, 'verbose','off','colormap', dmap, 'whitebk','on');
    end
    title(sprintf('Cluster %d (%.2f Hz)', iClust, f(fi)));
    colorbar;

    % ----- Frequency curve at peak channel -----
    % Decide GLM vs PSD from curve_data dims and presence of TMap
    sz = size(curve_data);
    if numel(sz)==2 && ~isempty(tmap) && all(size(curve_data)==size(tmap))
        % GLM: use beta line with CI from t
        b_row = curve_data(ci, :);
        t_row = tmap(ci, :);
        se_row = abs(b_row ./ t_row);
        se_row(~isfinite(se_row)) = NaN;
        lo = b_row - 1.96*se_row;
        hi = b_row + 1.96*se_row;
        y = b_row;
    elseif numel(sz)==3 && sz(1)==nChan && sz(2)==nFreq
        % PSD: mean across 3rd dim, CI from SE
        rows = squeeze(curve_data(ci, :, :)).';    % [nSub x nFreq]
        y  = nanmean(rows,1);
        s  = nanstd(rows,0,1);
        n  = sum(isfinite(rows),1);
        se = s ./ sqrt(max(n,1));
        lo = y - 1.96*se;
        hi = y + 1.96*se;
    else
        % Fallback: just plot the provided line with no CI
        if numel(sz)~=2 || ~all(size(curve_data)==[nChan nFreq])
            error('curve_data must be [nChan x nFreq] (GLM) or [nChan x nFreq x nSub] (PSD).');
        end
        y = curve_data(ci, :);
        lo = []; hi = [];
    end

    % apply line-noise mask
    if any(lnMask)
        y(lnMask) = NaN;
        if ~isempty(lo), lo(lnMask) = NaN; end
        if ~isempty(hi), hi(lnMask) = NaN; end
    end

    % significance mask for this channel
    sig_mask = false(1, nFreq);
    if iClust <= numel(mask_clusters) && ~isempty(mask_clusters{iClust})
        mc = mask_clusters{iClust};
        if size(mc,1) >= ci && size(mc,2) == nFreq
            sig_mask = mc(ci, :);
        end
    end
    sig_mask = sig_mask & ~lnMask;

    % plot (exactly your style + black dashed baseline)
    hs.curve{iClust} = figure('Color','w'); hold on

    % Shaded CI (blue), robust to NaNs
    if ~isempty(lo) && ~isempty(hi)
        good = isfinite(lo) & isfinite(hi);
        if any(good)
            d = [false, good, false];
            starts = find(diff(d)==1);
            ends   = find(diff(d)==-1) - 1;
            for s = 1:numel(starts)
                idx = starts(s):ends(s);
                fill([f(idx) fliplr(f(idx))], [lo(idx) fliplr(hi(idx))], ...
                     [0.4660, 0.6740, 0.1880], 'EdgeColor','none','FaceAlpha',0.5);
            end
        end
    end

    % Black dashed y=0
    plot([f(1) f(end)], [0 0], '--', 'Color', 'k', 'LineWidth', 1.2);

    % Main line (blue)
    plot(f, y, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 2);

    % Black significance bars at bottom, contiguous stretches
    ylims = ylim;
    ybar = ylims(1) - 0.05*(ylims(2)-ylims(1));
    d = [false, sig_mask, false];
    starts = find(diff(d)==1);
    ends   = find(diff(d)==-1) - 1;
    for s = 1:numel(starts)
        xs = f(starts(s)); xe = f(ends(s));
        line([xs xe], [ybar ybar], 'Color','k','LineWidth',4);
    end
    ylim([ybar ylims(2)]);

    xlabel('Frequency (Hz)');
    ylabel(line_label);
    title(sprintf('Cluster %d - %s', iClust, ch));
    set(gca,'FontSize',14,'LineWidth',1.2);
    box on; grid off; axis tight
end
end
