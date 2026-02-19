function plot_results(data_type, chan_type, xaxis, stats, mask, chanlocs, plot_type, clust_tbl)
% PLOT_RESULTS  Visualize significant EEG clusters
%
% plot_type: 'main', 'topo', 'course', or 'all'

if isempty(mask) || (iscell(mask) && ~any(cell2mat(mask), 'all')) || ...
   (~iscell(mask) && ~any(mask, 'all'))
    disp("Nothing significant, nothing to plot.")
    return
end

load colormap_bwr.mat; 
% load cm17.mat
% dmap = cm17a; 
dmap(1,:) = [0.9 0.9 0.9]; % set NaNs to grey


% Collapse mask if cell array
global_mask = false(size(stats));
if iscell(mask)
    for i = 1:numel(mask), global_mask = global_mask | mask{i}; end
else
    global_mask = mask;
end
masked_stats = stats; masked_stats(~global_mask) = NaN;

figure('Color','w');

% ---------------- MAIN PLOT ----------------
if strcmpi(plot_type,'main') || strcmpi(plot_type,'all')
    if strcmpi(plot_type,'all')
        subplot(2,2,[1 3]); % big panel
    end
    imagesc(xaxis, 1:size(stats,1), masked_stats);
    colormap(gca, dmap); c = colorbar;
    ylabel(c,'T-values','FontWeight','bold','FontSize',12)
    title('Corrected Statistical Map','FontSize',16,'FontWeight','bold')
    set(gca,'CLim',[-max(abs(stats(:))) max(abs(stats(:)))])
    if ~isempty(chanlocs)
        if length(chanlocs)>30
            set(gca,'YTick',1:2:numel(chanlocs), 'YTickLabel',{chanlocs(1:2:end).labels})
        else
            set(gca,'YTick',1:numel(chanlocs), 'YTickLabel',{chanlocs(1:2:end).labels})
        end
    end
    if contains(data_type,'time')
        xlabel('Time (ms)')
    elseif strcmpi(data_type,'frequency')
        xlabel('Frequency (Hz)')
    elseif strcmpi(data_type,'nonlinear')
        xlabel('Scale factor')
    end
    if strcmpi(chan_type, 'scalp')
        ylabel('Channels')
    elseif strcmpi(chan_type, 'source')
        ylabel("BRain Areas")
    end
end

% ---------------- TOPO PLOT ----------------
if strcmpi(plot_type,'topo') || strcmpi(plot_type,'all')
    if strcmpi(plot_type,'all')
        subplot(2,2,2);
    end
    % Pick strongest cluster
    [~, idx] = max(abs(clust_tbl.Tvalue));
    peak_val = clust_tbl.Peak(idx);
    chan_idx = find(strcmpi({chanlocs.labels}, clust_tbl.Channel{idx}));
    [~, freq_idx] = min(abs(xaxis - peak_val));
    topo_mask = mask{min(idx, numel(mask))}(:, freq_idx);
    topoplot(stats(:,freq_idx), chanlocs, 'colormap', dmap, 'pmask', topo_mask, ...
             'verbose','off','whitebk','on');
    title(sprintf('Cluster %d peak (%g)', idx, peak_val))
    colorbar;
end

% ---------------- COURSE PLOT ----------------
if strcmpi(plot_type,'course') || strcmpi(plot_type,'all')
    if strcmpi(plot_type,'all')
        subplot(2,2,4);
    end
    hold on
    [~, idx] = max(abs(clust_tbl.Tvalue));
    peak_val = clust_tbl.Peak(idx);
    chan_idx = find(strcmpi({chanlocs.labels}, clust_tbl.Channel{idx}));
    [~, freq_idx] = min(abs(xaxis - peak_val));
    plot(xaxis, stats(chan_idx,:), 'LineWidth',2);
    xlabel('X-axis'); ylabel('T-values')
    title(sprintf('Course Plot (%s)', clust_tbl.Channel{idx}))
    box on
end

set(findall(gcf, 'type', 'axes'), 'FontSize', 12, 'FontWeight', 'bold');


end
