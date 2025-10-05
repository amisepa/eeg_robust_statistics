function hs = plot_clusters(summary_tbl, mask_clusters, tmap, map, f, chanlocs, line_label, varargin)
% PLOT_CLUSTERS  Scalp or cortex plot at cluster peak + frequency curve.
%
% Required
%   summary_tbl    table with columns Peak_Hz or Peak, and Channel
%   mask_clusters  cell {nClust}, each logical [nChan x nFreq]
%   tmap           [nChan x nFreq] t-values map (used for CI and cortex coloring)
%   map            GLM: [nChan x nFreq] beta map
%                  PSD: [nChan x nFreq x nSub] per-subject values
%   f              [1 x nFreq] frequency vector
%   chanlocs       EEGLAB chanlocs struct
%   line_label     y-axis label for the curve
%
% Name-value
%   'TopoMap'      [nChan x nFreq] to override the scalp topography source. Default []
%                  If empty, uses map (or mean(map,3) if PSD).
%   'LineNoiseHz'  scalar or vector freqs to mask in curves (default 60). [] disables
%   'LineNoiseBW'  bandwidth (Hz) around each line-noise freq (default 1)
%   'DataType'     'scalp' (default) or 'source'

% Silent loads for colormaps and cortex resources
load("colormap_bwr.mat"); dmap(1,:) = [.9 .9 .9]; % set NaNs to gray
load cm17.mat
load cortex; cortex = cortex_highres;

ip = inputParser;
addParameter(ip,'TopoMap',[]);
addParameter(ip,'LineNoiseHz',60);
addParameter(ip,'LineNoiseBW',1);
addParameter(ip,'DataType','scalp');
parse(ip,varargin{:});
topo_map = ip.Results.TopoMap;
lnHz     = ip.Results.LineNoiseHz;
lnBW     = ip.Results.LineNoiseBW;
dataType = ip.Results.DataType;

% Derive topo_map if not supplied
if isempty(topo_map)
    if ndims(map)==3
        topo_map = nanmean(map,3); % PSD case
    else
        topo_map = map;            % GLM case
    end
end

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

% Aliases for cortex snippet
f2    = f;
tvals = tmap;
mask  = mask_clusters;

% Precompute t-range for cortex logic
if ~isempty(tvals)
    tmin = min(tvals(:),[],'omitnan');
    tmax = max(tvals(:),[],'omitnan');
else
    tmin = NaN; tmax = NaN;
end

for iClust = 1:height(summary_tbl)
    % Peak freq and channel
    if ismember('Peak_Hz', summary_tbl.Properties.VariableNames)
        peak_freq = summary_tbl.Peak_Hz(iClust);
    elseif ismember('Peak', summary_tbl.Properties.VariableNames)
        peak_freq = summary_tbl.Peak(iClust);
    else
        error('summary_tbl must have Peak_Hz or Peak.');
    end
    [~, fi] = min(abs(f - peak_freq));
    ch = summary_tbl.Channel{iClust};
    if isstring(ch), ch = char(ch); end

    % ----- Topography or Cortex at peak freq -----
    topo_mask = false(nChan,1);
    if iClust <= numel(mask_clusters) && ~isempty(mask_clusters{iClust})
        mc = mask_clusters{iClust};
        if size(mc,2) == nFreq, topo_mask = mc(:, fi); end
    end

    if strcmpi(dataType,'source')
        % Your cortex plotting logic, colors and scaling preserved
        roiPSD = tvals(:, fi);
        if iClust <= numel(mask) && ~isempty(mask{iClust}) && size(mask{iClust},2) >= fi
            roiPSD(~mask{iClust}(:, fi)) = NaN;
        end

        % REMOVE this line to avoid the empty figure:
        % hs.topo{iClust} = figure('Color','w');

        if tmin < 0 && tmax > 0
            allplots_cortex_BS(cortex, mean(tvals,2,'omitnan'), [-max(abs(tvals(:))) max(abs(tvals(:)))], cm17, 't-values', []);
        elseif tmax <= 0
            allplots_cortex_BS(cortex, mean(tvals,2,'omitnan'), [min(mean(tvals,2,'omitnan')) tmax], cm17b, 't-values', []);
        elseif tmin >= 0
            allplots_cortex_BS(cortex, mean(tvals,2,'omitnan'), [tmin max(mean(tvals,2,'omitnan'))], cm17a, 't-values', []);
        end

        % Now capture the figure that allplots_cortex_BS created
        hs.topo{iClust} = gcf;
        title(sprintf('Cluster %d (%.2f Hz) - Cortex', iClust, f2(fi)));

    else
        % Scalp topography path
        ci = find(strcmpi({chanlocs.labels}, ch), 1);
        if isempty(ci)
            warning('plot_clusters:MissingChan','Skipping cluster %d, channel %s not found.', iClust, ch);
            continue
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
    end

    % ----- Frequency curve at peak channel -----
    sz = size(map);

    % Find channel index for curve
    ci_curve = [];
    if ~isempty(ch)
        ci_curve = find(strcmpi({chanlocs.labels}, ch), 1);
    end
    if isempty(ci_curve)
        warning('plot_clusters:CurveChannelMissing','Curve skipped for cluster %d because channel "%s" was not found.', iClust, ch);
        hs.curve{iClust} = [];
        continue
    end

    if numel(sz)==2 && ~isempty(tmap) && all(size(map)==size(tmap))
        % GLM: beta line with CI from t
        b_row = map(ci_curve, :);
        t_row = tmap(ci_curve, :);
        se_row = abs(b_row ./ t_row);
        se_row(~isfinite(se_row)) = NaN;
        lo = b_row - 1.96*se_row;
        hi = b_row + 1.96*se_row;
        y = b_row;
    elseif numel(sz)==3 && sz(1)==nChan && sz(2)==nFreq
        % PSD: mean across subjects, CI from SE
        rows = squeeze(map(ci_curve, :, :)).';    % [nSub x nFreq]
        y  = nanmean(rows,1);
        s  = nanstd(rows,0,1);
        n  = sum(isfinite(rows),1);
        se = s ./ sqrt(max(n,1));
        lo = y - 1.96*se;
        hi = y + 1.96*se;
    else
        % Fallback: provided line only
        if numel(sz)~=2 || ~all(size(map)==[nChan nFreq])
            error('map must be [nChan x nFreq] (GLM) or [nChan x nFreq x nSub] (PSD).');
        end
        y = map(ci_curve, :);
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
        if size(mc,1) >= ci_curve && size(mc,2) == nFreq
            sig_mask = mc(ci_curve, :);
        end
    end
    sig_mask = sig_mask & ~lnMask;

    % plot curve
    hs.curve{iClust} = figure('Color','w'); hold on

    % Shaded CI
    if ~isempty(lo) && ~isempty(hi)
        good = isfinite(lo) & isfinite(hi);
        if any(good)
            dci = [false, good, false];
            starts = find(diff(dci)==1);
            ends   = find(diff(dci)==-1) - 1;
            for s = 1:numel(starts)
                idx = starts(s):ends(s);
                fill([f(idx) fliplr(f(idx))], [lo(idx) fliplr(hi(idx))], ...
                     [0.4660, 0.6740, 0.1880], 'EdgeColor','none','FaceAlpha',0.5);
            end
        end
    end

    % Black dashed y=0
    plot([f(1) f(end)], [0 0], '--', 'Color', 'k', 'LineWidth', 1.2);

    % Main line
    plot(f, y, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 2);

    % Bottom significance bars
    ylims = ylim;
    ybar = ylims(1) - 0.05*(ylims(2)-ylims(1));
    dsig = [false, sig_mask, false];
    starts = find(diff(dsig)==1);
    ends   = find(diff(dsig)==-1) - 1;
    for s = 1:numel(starts)
        xs = f(starts(s)); xe = f(ends(s));
        line([xs xe], [ybar ybar], 'Color','k','LineWidth',10);
    end
    ylim([ybar ylims(2)]);

    xlabel('Frequency (Hz)');
    ylabel(line_label);
    title(sprintf('Cluster %d - %s', iClust, ch));
    set(gca,'FontSize',14,'LineWidth',1.2);
    box on; grid off; axis tight
end
end
