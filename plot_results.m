function plot_results(datatype, chan_type, xaxis, stats, mask, chanlocs, plot_type, clust_tbl)
% PLOT_RESULTS Visualizes statistically significant EEG clusters
%
% Usage:
%   plot_results(datatype, chan_type, xaxis, stats, mask, chanlocs, plot_type, clust_tbl)
%
% Inputs:
%   - datatype     : 'time', 'frequency', or 'time-frequency'
%   - chan_type    : 'scalp' or 'source'
%   - xaxis        : time or frequency vector (1 x N)
%   - stats        : statistical map (channels x time/frequency), e.g., t-values
%   - mask         : binary mask (same size as stats); non-significant values set to NaN
%   - chanlocs     : EEGLAB-style channel structure
%   - plot_type    : 'all' to plot in one figure; otherwise, separate figures
%   - clust_tbl    : cluster summary table with fields:
%                     'Cluster', 'Start_Hz', 'End_Hz', 'Peak_Hz', 
%                     'Tvalue', 'Channel', 'Polarity'
%
% Assumes a custom diverging colormap file named 'colormap_bwr.mat' is available

if isempty(mask)
    disp("Nothing significant, nothing to plot")
    return
end

load colormap_bwr.mat  % provides variable `dmap`
dmap(1,:) = [0.9 0.9 0.9];  % gray for NaNs

figure('Color','w');

% Main statistical map plot
if strcmpi(plot_type, 'all')
    subplot(3,3,[1 2 4 5 7 8])
end

% Combine all cluster masks into one logical mask
global_mask = false(size(stats));
for i = 1:numel(mask)
    global_mask = global_mask | mask{i};
end

% Apply to stats
masked_stats = stats;
masked_stats(~global_mask) = NaN;

% Main plot
imagesc(xaxis, 1:size(stats,1), masked_stats);
colormap(gca, dmap);
c = colorbar; ylabel(c, 'T-values','FontWeight','bold','FontSize',11,'Rotation',-90)

% Y-axis labels
if ~isempty(chanlocs)
    Yticks = {chanlocs.labels};
    skip = 2;
    set(gca, 'YTick', 1:skip:length(Yticks), ...
             'YTickLabel', Yticks(1:skip:end), ...
             'FontWeight','normal');
    % ylabel(strcmpi(chan_type,'scalp') * 'EEG channels' + ...
    %        strcmpi(chan_type,'source') * 'Brain areas', ...
    %        'FontSize',12,'FontWeight','bold');
    if strcmpi(chan_type, 'scalp')
        ylab = 'EEG channels';
    elseif strcmpi(chan_type, 'source')
        ylab = 'Brain areas';
    else
        ylab = 'Channels';
    end
    ylabel(ylab, 'FontSize', 12, 'FontWeight', 'bold');

end

% X-axis label
if contains(datatype, 'time')
    xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
elseif strcmpi(datatype, 'frequency')
    xlabel('Frequency (Hz)', 'FontSize', 12, 'FontWeight', 'bold');
end

title('Corrected Statistical Map', 'FontSize', 12, 'FontWeight', 'bold');
clim([-max(abs(stats(:))) max(abs(stats(:)))]);
axis tight; box on; set(gca, 'LineWidth', 1);

% Topoplot for each cluster
if strcmpi(plot_type,'all')

    for iClust = 1:height(clust_tbl)
        peak_freq = clust_tbl.Peak_Hz(iClust);
        [~, freq_idx] = min(abs(xaxis - peak_freq));  % closest index
    
        peak_channel = clust_tbl.Channel{iClust};
        chan_idx = find(strcmpi({chanlocs.labels}, peak_channel));
        if isempty(chan_idx)
            warning('Channel %s not found in chanlocs', peak_channel);
            continue
        end
    
        topo_mask = mask{iClust}(:, freq_idx);  % extract the iClust-th logical mask, then select that frequency column
    
        % Sanity check: Make sure topo_mask is the same length as chanlocs
        if length(topo_mask) ~= length(chanlocs)
            warning('Topoplot skipped for cluster %d: mask size (%d) does not match chanlocs (%d)', ...
                iClust, length(topo_mask), length(chanlocs));
            continue
        end
    
        % Sanity check: Make sure data vector matches
        data_vec = stats(:, freq_idx);
        if length(data_vec) ~= length(chanlocs)
            warning('Topoplot skipped: stats length does not match chanlocs');
            continue
        end
    
        % Plot topography
        topoplot(data_vec, chanlocs, ...
            'colormap', dmap, 'pmask', topo_mask, ...
            'verbose','off', 'whitebk','on');
        % topoplot(stats(:, freq_idx), chanlocs, ...
        %     'colormap', dmap, 'pmask', topo_mask, ...
        %     'verbose','off', 'whitebk','on');
    
        title(sprintf('Cluster %d (%.1f Hz)', iClust, peak_freq), ...
            'FontSize', 13, 'FontWeight', 'bold');
        cb = colorbar;
        ylabel(cb, 'T-values', 'Rotation', 270, 'FontWeight', 'bold', 'FontSize', 11);
    end
end

% Course plot for strongest cluster (optional)
if strcmpi(plot_type, 'all') && ~isempty(clust_tbl)
    subplot(3,3,9)
    hold on
    [~, strongest_idx] = max(abs(clust_tbl.Tvalue));
    peak_freq = clust_tbl.Peak_Hz(strongest_idx);
    [~, freq_idx] = min(abs(xaxis - peak_freq));
    chan_idx = find(strcmpi({chanlocs.labels}, clust_tbl.Channel{strongest_idx}));
    plot(xaxis, stats(chan_idx, :), 'LineWidth', 2);
    ylabel('T-values','FontSize',11,'FontWeight','bold');
    xlabel('Frequency (Hz)','FontSize',11,'FontWeight','bold');
    title(sprintf('Course Plot: %s', clust_tbl.Channel{strongest_idx}), ...
        'FontSize', 11, 'FontWeight', 'bold');
    axis tight
end

% Final polish
set(gcf, 'Name', 'Results from mass-univariate analysis', ...
         'Color', 'w', 'NumberTitle', 'Off');
set(findall(gcf, 'type', 'axes'), 'FontSize', 10, 'FontWeight', 'bold');
