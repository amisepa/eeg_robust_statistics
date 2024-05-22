%% Plots results using correction for multiple comparison masks
%
% Usage:
%    [peakClust, peakChan, peakLat, clustOrder] = plot_results(datatype, xaxis, stats, mask, pcorr, alpha, chanlocs, mcctype)
% 
% Example:
%   [peakClust, peakChan, peakLat, clustOrder] = plot_results(datatype, xaxis, stats, mask, pcorr, alpha, chanlocs, mcctype)
%
% Cedric Cannard, Sep 2022

function [peakClust, peakChan, peakLat, clustOrder] = plot_results(datatype, xaxis, stats, mask, pcorr, alpha, chanlocs, mcctype)
% Look up limo_display_image_tf for time-frequency data

peakClust = [];
peakChan = [];
peakLat = [];
clustOrder = [];

if isempty(mask)
    disp('Empty mask, computing one from corrected p-values.'); %return
    mask = pcorr < alpha;
end
if sum(mask,'all') == 0
    warning('No significant differences');
end

if sum(mask,'all') > 0

    % Get clusters properties (start/end, peak, channel, frame/freq)
    n_cluster     = max(mask(:));
    cluster_start = nan(1,n_cluster); % start of each cluster
    cluster_end   = nan(1,n_cluster); % end of each cluster
    cluster_maxv  = nan(1,n_cluster); % max value for each cluster
    cluster_maxe  = nan(1,n_cluster); % channel location of the max value of each cluster
    cluster_maxf  = nan(1,n_cluster); % frame location of the max value of each cluster
    for iClust = 1:n_cluster
        tmp = stats.*(mask==iClust);
        tmp(tmp==Inf) = NaN; tmp(tmp==-Inf) = NaN;
        sigframes = sum(tmp,1);
        cluster_start(iClust) = find(sigframes,1,'first');
        cluster_end(iClust) = find(sigframes,1,'last');
        [V,type] = max([abs(min(tmp(:))) max(tmp(:))]);
        if type == 2
            cluster_maxv(iClust) = V(1);
        else
            V = -V;
            cluster_maxv(iClust) = V(1);
        end
        [cluster_maxe(iClust), cluster_maxf(iClust)] = ind2sub(size(tmp),find(tmp==V(1)));
    end
    
    % Sort clusters from strongest to smallest t-value
    [~, maxEffect] = max(abs(cluster_maxv)); 
    [~,idx] = sort(cluster_maxv,'descend','ComparisonMethod','abs');
    % cluster_maxe = cluster_maxe(idx);
    % cluster_maxf = cluster_maxf(idx);
    % cluster_maxv = cluster_maxv(idx);
    % cluster_start = cluster_start(idx);
    % cluster_end = cluster_end(idx);
    % [~, maxEffect] = max(abs(cluster_maxv)); % should be 1 now, but to be sure
    clustOrder = idx;  % to store order of strength
    peakClust = [cluster_start(idx(1)) cluster_end(idx(1)) ];

    % Print in command window
    for iClust = 1:n_cluster
        if strcmpi(datatype, 'time')
            fprintf('Cluster %g: %g to %g ms. Peak effect: channel %s at %g ms (t = %g) \n', ...
                iClust, xaxis(cluster_start(idx(iClust))), xaxis(cluster_end(idx(iClust))), ...
                chanlocs(cluster_maxe(idx(iClust))).labels, xaxis(cluster_maxf(idx(iClust))), round(cluster_maxv(idx(iClust)),1) );
        elseif strcmpi(datatype, 'frequency')
            fprintf('Cluster %g: %g to %g Hz. Peak effect: channel %s at %g Hz (t = %g) \n', ...
                iClust, xaxis(cluster_start(idx(iClust))), xaxis(cluster_end(idx(iClust))), ...
                chanlocs(cluster_maxe(idx(iClust))).labels, xaxis(cluster_maxf(idx(iClust))), round(cluster_maxv(idx(iClust)),1) );
        end
    end

    %% MAIN PLOT

    figure('Color','w'); set(gcf,'InvertHardCopy','off');

    subplot(3,3,[1 2 4 5 7 8]);

    % Extract significant t-values for plotting
    effects = stats.*single(mask>0);
    effects(effects==0) = NaN;
    % v = max(effects(:));
    % [c, tf_stamp] = find(effects==v); % which channel and time/frequency frame
    % if length(c) > 1  % if we have multiple times the exact same max values
    %     c = c(1);
    %     tf_stamp = tf_stamp(1);  % then take the 1st (usually an artifact but allows to see it)
    % end

    % Image all channels and time/frequency data
    imagesc(xaxis,1:size(effects,1),effects);

    % Color palette
    load("colormap_bwr.mat");
    % dmap = colormap("viridis");
    % load("colormap_bgy.mat");
    % dmap = colormap("bone"); % "winter" "hot" 
    if sum(isnan(effects(:))) ~= 0
        dmap(1,:) = [.9 .9 .9]; % set NaNs to gray
    end
    colormap(gca, dmap); % colormap('parula')

    set(gca,'LineWidth',1)
    c = colorbar;
    ylabel(c, 'T-values','FontWeight','bold','FontSize',12,'Rotation',-90)
    % set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');

    % Y labels (EEG channels or brain areas)
    if contains(datatype,{'time','time-frequency'})
        xlabel('Time (ms)','FontSize',13,'FontWeight','bold')
    elseif strcmpi(datatype,'frequency')
        xlabel('Frequency (Hz)','FontSize',13,'FontWeight','bold')
    end

    % Y tick labels (electrode or area names)
    if strcmpi(datatype,'time-frequency')
        ylabel('Frequency (Hz)','FontSize',13,'FontWeight','bold')
    else
        % For scalp channels
        ylabel('EEG channels','FontSize',13,'FontWeight','bold');

        % For Brain areas
        % ylabel('Brain areas','FontSize',13,'FontWeight','bold');
        % for i = 1:length(chanlocs)
        %     chanlocs(i).labels = char(join(split(chanlocs(i).labels,'_')));
        % end
        
        % Add labels to plot
        Ylabels = {chanlocs.labels};
        % img_prop = get(gca);
        % newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        newticks = 1:2:length(Ylabels);
        newticks = unique(newticks);
        Ylabels  = Ylabels(newticks);
        set(gca,'YTick',newticks);
        set(gca,'YTickLabel', Ylabels,'FontWeight','bold');
    end
    correctoptions = {'Uncorrected' 'Max-corrected' 'Cluster-corrected' 'TFCE-corrected'};
    title(sprintf('%s (p<%g)', correctoptions{mcctype+1},alpha),'FontSize',13,'FontWeight','bold');

    % Clim
    maxval = max(abs(effects(:)));
    % if max(effects(:)) < 0
    %     clim([-maxval 0])
    % elseif min(effects(:)) > 0
    %     clim([0 maxval])
    % else
    clim([-maxval maxval])
    % end
   
    %% Course plot of t-values of strongest effect

    peakChan = cluster_maxe(maxEffect);

    subplot(3,3,9);
    plot(xaxis,stats(peakChan,:),'LineWidth',2);
    chanLabel = chanlocs(peakChan).labels;
    title(sprintf('Course plot: %s',chanLabel),'FontSize',11,'fontweight','bold')
    % plot(xaxis,stats(cluster_maxe,:),'LineWidth',2);  % plot peak effect of all clusters superimposed
    % chanLabel = {chanlocs(cluster_maxe).labels};
    % legend(chanLabel)
    grid on; axis tight;
    ylabel('t-values','FontSize',11,'fontweight','bold'); 
    xlabel('Frequency (Hz)','FontSize',11,'fontweight','bold')

    % Plot bars of significnace for peak electrode
    plotSigBar(mask(peakChan,:)~=0,xaxis);

    %% Scalp topography at peak latency/frequency (replace with 3D headplot?)

    peakLat = cluster_maxf(maxEffect);

    % figure

    subplot(3,3,6);

    % topoplot(stats(:, peakLat), chanlocs,'emarker',{'.','k',5,1}, ...             % normal electrodes
    %     'emarker2',{find(mask(:, peakLat)),'.',"red",7,1}, ...    % significant electrodes
    %     'verbose','off','colormap',dmap);                                                      % parameters
    topoplot(stats(:, peakLat), chanlocs,'emarker2',{find(mask(:, peakLat)),'.','k',10,1}, ...    % significant electrodes
        'verbose','off','colormap',dmap);                                                      % parameters
    if strcmpi(datatype, 'time')
        title(sprintf('Scalp topography: %g ms', xaxis(peakLat)), ...
            'FontSize',13,'fontweight','bold');
    elseif strcmpi(datatype, 'frequency')
        title(sprintf('Scalp topography: %g Hz', xaxis(peakLat)), ...
            'FontSize',13,'fontweight','bold');
    end
    set(gcf,'Name','Topography at peak frequency','color','w','Toolbar','none','Menu','none','NumberTitle','Off')

    %% Adjust labels and ticks font and weight for all plots

    set(gcf,'Name','Results','color','w','Toolbar','none','Menu','none','NumberTitle','Off')
    % set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

else
    disp('Nothing to plot (nothing is significant)')
end

