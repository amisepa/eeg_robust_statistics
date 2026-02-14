function hs = plot_clusters(summary_tbl, mask, tvals, map, f, chanlocs, line_label, varargin)
% PLOT_CLUSTERS  Scalp or cortex plot at cluster peak + frequency curve.
%
% Required
%   summary_tbl    table with columns Peak_Hz or Peak, and Channel
%   mask_clusters  cell {nClust}, each logical [nChan x nFreq]
%   tvals          [nChan x nFreq] t-values map (used for CI and cortex coloring)
%   map            GLM: [nChan x nFreq] beta map
%                  PSD: [nChan x nFreq x nSub] per-subject values
%   f              [1 x nFreq] numeric frequency/scale vector
%   chanlocs       EEGLAB chanlocs struct
%   line_label     y-axis label for the curve
%
% Name-value
%   'LineNoiseHz'   scalar or vector freqs to mask in curves (default 60). [] disables
%   'LineNoiseBW'   bandwidth (Hz) around each line-noise freq (default 1)
%   'DataType'      'scalp' (default) or 'source'
%   'XTickLabels'   cellstr/strings for x-axis labels on the curve plot; if empty, numeric f is shown
%
% Notes
%   - Axis computations always use numeric f
%   - XTick labels are applied only to the curve plot and never used for math

ip = inputParser;
addParameter(ip,'LineNoiseHz',60);
addParameter(ip,'LineNoiseBW',1);
addParameter(ip,'DataType','scalp');
addParameter(ip,'XTickLabels',[]);     % avoid shadowing xticks()
parse(ip,varargin{:});
lnHz        = ip.Results.LineNoiseHz;
lnBW        = ip.Results.LineNoiseBW;
dataType    = ip.Results.DataType;
xTickLabels = ip.Results.XTickLabels;

% Basic sizes
[nChan, nFreq] = size(tvals);
assert(isnumeric(f) && numel(f)==nFreq, 'f must be numeric and match number of columns');

hs.topo  = cell(height(summary_tbl),1);
hs.curve = cell(height(summary_tbl),1);

% t-range for cortex logic
if ~isempty(tvals)
    tmin = min(tvals(:),[],'omitnan');
    tmax = max(tvals(:),[],'omitnan');
else
    tmin = NaN; tmax = NaN;
end

for iClust = 1:height(summary_tbl)
    % Peak freq and channel
    if ismember('Peak', summary_tbl.Properties.VariableNames)
        peak_freq = summary_tbl.Peak(iClust);
    elseif ismember('Peak_Hz', summary_tbl.Properties.VariableNames)
        peak_freq = summary_tbl.Peak_Hz(iClust);
    else
        error('summary_tbl must have Peak or Peak_Hz.');
    end
    if ~isnumeric(peak_freq) || ~isscalar(peak_freq) || ~isfinite(peak_freq)
        warning('plot_clusters:BadPeak','Invalid Peak for cluster %d; using mid-bin.', iClust);
        fi = round(nFreq/2);
    else
        [~, fi] = min(abs(f - peak_freq));
    end
    ch = summary_tbl.Channel{iClust};
    if isstring(ch), ch = char(ch); end

        %%%%%%%%%%%%%%%%%% CURVE PLOT AT PEAK CHANNEL %%%%%%%%%%%%%%%%%%
    sz = size(map);

    % line-noise mask
    lnMask = false(1, nFreq);
    if ~isempty(lnHz)
        lnHz = lnHz(:).';
        halfBW = max(lnBW, eps)/2;
        for v = lnHz
            lnMask = lnMask | (abs(f - v) <= halfBW);
        end
    end

    % channel index for curve
    ci_curve = [];
    if ~isempty(ch)
        ci_curve = find(strcmpi({chanlocs.labels}, ch), 1);
    end
    if isempty(ci_curve)
        warning('plot_clusters:CurveChannelMissing','Curve skipped for cluster %d because channel "%s" was not found.', iClust, ch);
        hs.curve{iClust} = [];
        continue
    end

    if numel(sz)==2 && ~isempty(tvals) && all(size(map)==size(tvals))
        % GLM: beta with CI from t
        b_row = map(ci_curve, :);
        t_row = tvals(ci_curve, :);
        se_row = abs(b_row ./ t_row);
        se_row(~isfinite(se_row)) = NaN;
        lo = b_row - 1.96*se_row;
        hi = b_row + 1.96*se_row;
        y = b_row;
    elseif numel(sz)==3 && sz(1)==nChan && sz(2)==nFreq
        % PSD: mean across subjects, CI from SE
        rows = squeeze(map(ci_curve, :, :)).';    % [nSub x nFreq]
        y  = mean(rows,1,'omitmissing');
        s  = std(rows,0,1,'omitmissing');
        n  = sum(isfinite(rows),1);
        se = s ./ sqrt(max(n,1));
        lo = y - 1.96*se;
        hi = y + 1.96*se;
    else
        % fallback: provided line only
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

    % significance for this channel
    sig_mask = false(1, nFreq);
    if iClust <= numel(mask) && ~isempty(mask{iClust})
        mc = mask{iClust};
        if size(mc,1) >= ci_curve && size(mc,2) == nFreq
            sig_mask = mc(ci_curve, :);
        end
    end
    sig_mask = sig_mask & ~lnMask;

    % plot curve
    hs.curve{iClust} = figure('Color','w'); hold on

    % shaded CI
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

    % baseline
    plot([f(1) f(end)], [0 0], '--', 'Color', 'k', 'LineWidth', 1.2);

    % main line
    plot(f, y, 'color',[0.4660, 0.6740, 0.1880], 'LineWidth', 2);

    % bottom significance bars (robust to 1-bin clusters)
    yl = ylim; H = 0.015*(yl(2)-yl(1)); y0 = yl(1) - 1.6*H; ylim([y0 yl(2)]);
    if numel(f) > 1
        d1 = f(2)-f(1); dN = f(end)-f(end-1);
        edges = [f(1)-d1/2, 0.5*(f(1:end-1)+f(2:end)), f(end)+dN/2];
    else
        edges = [f(1)-0.5, f(1)+0.5];
    end
    dsig = [false, sig_mask(:).', false];
    starts = find(diff(dsig)== 1);
    ends   = find(diff(dsig)==-1);
    for k = 1:numel(starts)
        xs = edges(starts(k)); xe = edges(ends(k));
        rectangle('Position',[xs, y0, max(eps, xe-xs), H], ...
                  'FaceColor','k', 'EdgeColor','none');
    end

    % optional x tick labels
    if ~isempty(xTickLabels)
        % accept cellstr/string; otherwise convert
        labs = string(xTickLabels(:).');
        % if mismatch, sample to ~10 labels
        if numel(labs) ~= nFreq
            stride = max(1, round(nFreq/10));
            keep = false(1, nFreq); keep(1:stride:end) = true; keep(end) = true;
            xticks(f(keep));
            xticklabels(labs(1:min(sum(keep), numel(labs))));
        else
            xticks(f);
            xticklabels(labs);
        end
        xtickangle(45);
        set(gca,'TickLabelInterpreter','none');
    end

    xlabel('Time (ms)');     % or 'Frequency (Hz)' % 'Scale factor'
    ylabel(line_label);
    title(sprintf('Cluster %d - %s', iClust, ch));
    set(gca,'FontSize',14,'LineWidth',1.2);
    box on; grid off; axis tight
    set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');

    
    %%%%%%%%%%%% TOPO/BRAIN PLOT AT PEAK %%%%%%%%%%%%
    topo_mask = false(nChan,1);
    if iClust <= numel(mask) && ~isempty(mask{iClust})
        mc = mask{iClust};
        if size(mc,2) == nFreq
            % topo_mask = mc(:, fi);  % topo mask is extracted at only the peak time bin
            topo_mask = any(mc, 2);   % collapse across time
        end
    end
    
    % 3D BRAIN PLOT (SOURCE DATA)
    if strcmpi(dataType,'source')
        load cm17.mat
        load cortex; cortex = cortex_highres;

        roiPSD = tvals(:, fi);
        if iClust <= numel(mask) && ~isempty(mask{iClust}) && size(mask{iClust},2) >= fi
            roiPSD(~mask{iClust}(:, fi)) = NaN;
        end

        if tmin < 0 && tmax > 0
            allplots_cortex_BS(cortex, mean(tvals,2,'omitnan'), [-max(abs(tvals(:))) max(abs(tvals(:)))], cm17, 't-values', []);
        elseif tmax <= 0
            allplots_cortex_BS(cortex, mean(tvals,2,'omitnan'), [min(mean(tvals,2,'omitnan')) tmax], cm17b, 't-values', []);
        elseif tmin >= 0
            allplots_cortex_BS(cortex, mean(tvals,2,'omitnan'), [tmin max(mean(tvals,2,'omitnan'))], cm17a, 't-values', []);
        end
        hs.topo{iClust} = gcf;
    
    % SCALP TOPOGRAPHY
    else
        load("colormap_bwr.mat"); dmap(1,:) = [.9 .9 .9]; % NaNs to gray

        ci = find(strcmpi({chanlocs.labels}, ch), 1);
        if isempty(ci)
            warning('plot_clusters:MissingChan','Skipping cluster %d, channel %s not found.', iClust, ch);
            continue
        end

        hs.topo{iClust} = figure('Color','w');
        % try
            topoplot(tvals(:, fi), chanlocs, 'pmask', topo_mask, ...
                     'verbose','off','colormap', dmap, 'whitebk','on', ...
                     'gridscale', 300);

            % topoplot(tvals(:, fi), chanlocs, 'pmask', topo_mask, ...
            %     'verbose','off', 'colormap', dmap, 'whitebk','on', ...
            %     'gridscale', 300, 'interplimits', 'head','plotrad',0.5);

            % % After your topoplot call, add this:
            % h = findobj(gca, 'Type', 'surface');
            % if isempty(h), h = findobj(gca, 'Type', 'image'); end
            % if ~isempty(h)
            %     cd = get(h(1), 'CData');
            %     nanmask = isnan(cd);
            %     if any(nanmask(:))
            %         cd = regionfill(cd, nanmask);
            %         set(h(1), 'CData', cd);
            %     end
            % end


        % catch
            % topoplot(tvals(:, fi), chanlocs, 'verbose','off','colormap', dmap, 'whitebk','on');
        % end
        c = colorbar; ylabel(c,'t-values','FontWeight','bold','FontSize',11)
        set(gca,'CLim',[-max(abs(tvals(:))) max(abs(tvals(:)))])
        set(gcf, 'Color', 'white')
    end

    title(sprintf('Cluster %d (Time: %.0f ms)', iClust, f(fi)));
    % title(sprintf('Cluster %d (Scale factor: %g)', iClust, f(fi)));
    set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');

end
end
