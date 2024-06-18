%% Plots results using correction for multiple comparison masks
%
% Usage:
%    [cluster_bounds, cluster_maxchan, cluster_maxfreq, cluster_order] = plot_results(datatype, xaxis, stats, mask, pcorr, alpha, chanlocs, mcctype)
% 
% Example:
%   [cluster_bounds, cluster_maxchan, cluster_maxfreq, cluster_order] = plot_results('time', times, tvals, mask, pcorr, 0.05, chanlocs, 0);
% 
% Cedric Cannard, Sep 2022

function [cluster_bounds, cluster_maxchan, cluster_maxfreq, cluster_maxval] = plot_results(datatype, xaxis, stats, mask, pcorr, alpha, chanlocs, mcctype)
% Look up limo_display_image_tf for time-frequency data

cluster_bounds = [];
cluster_maxchan = [];
cluster_maxfreq = [];
cluster_maxval = [];
% cluster_order = [];

if isempty(mask)
    disp('Empty mask, computing one from corrected p-values.'); %return
    mask = pcorr < alpha;
end

if sum(mask,'all') > 0

    % Get clusters properties (start/end, peak, channel, frame/freq)
    n_cluster     = max(mask(:));
    cluster_start = nan(1,n_cluster); % start of each cluster
    cluster_end   = nan(1,n_cluster); % end of each cluster
    cluster_maxval  = nan(1,n_cluster); % max value for each cluster
    cluster_maxchan  = nan(1,n_cluster); % channel location of the max value of each cluster
    cluster_maxfreq  = nan(1,n_cluster); % frame location of the max value of each cluster
    for iClust = 1:n_cluster
        tmp = stats.*(mask==iClust);
        tmp(tmp==Inf) = NaN; tmp(tmp==-Inf) = NaN;
        sigframes = sum(tmp,1);
        cluster_start(iClust) = find(sigframes,1,'first');
        cluster_end(iClust) = find(sigframes,1,'last');
        [V,type] = max([abs(min(tmp(:))) max(tmp(:))]);
        if type == 2
            cluster_maxval(iClust) = V(1);
        else
            V = -V;
            cluster_maxval(iClust) = V(1);
        end
        [cluster_maxchan(iClust), cluster_maxfreq(iClust)] = ind2sub(size(tmp),find(tmp==V(1)));
    end
    
    % Sort clusters from strongest to smallest t-value
    [~, maxEffect] = max(abs(cluster_maxval)); 
    % [~,idx] = sort(cluster_maxval,'descend','ComparisonMethod','abs');
    % cluster_maxe = cluster_maxe(idx);
    % cluster_maxf = cluster_maxf(idx);
    % cluster_maxv = cluster_maxv(idx);
    % cluster_start = cluster_start(idx);
    % cluster_end = cluster_end(idx);
    % [~, maxEffect] = max(abs(cluster_maxv)); % should be 1 now, but to be sure
    % cluster_order = idx;  % to store order of strength
    % peakClust = [cluster_start(idx(1)) cluster_end(idx(1)) ];
    cluster_bounds = [cluster_start; cluster_end ]';

    % Print in command window
    for iClust = 1:n_cluster
        if strcmpi(datatype, 'time')
            fprintf('Cluster %g: %g to %g ms. Peak effect: channel %s at %g ms (t = %g) \n', ...
                iClust, xaxis(cluster_start(iClust)), xaxis(cluster_end(iClust)), ...
                chanlocs(cluster_maxchan(iClust)).labels, xaxis(cluster_maxfreq(iClust)), round(cluster_maxval(iClust),1) );
        elseif strcmpi(datatype, 'frequency')
            fprintf('Cluster %g: %g to %g Hz. Peak effect: channel %s at %g Hz (t = %g) \n', ...
                iClust, xaxis(cluster_start(iClust)), xaxis(cluster_end(iClust)), ...
                chanlocs(cluster_maxchan(iClust)).labels, xaxis(cluster_maxfreq(iClust)), round(cluster_maxval(iClust),1) );
        end
    end

    %% MAIN PLOT

    figure('Color','w','InvertHardCopy','off');

    % subplot(3,3,[1 2 4 5 7 8]);

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
    ylabel(c, 'T-values','FontWeight','bold','FontSize',11,'Rotation',-90)
    % set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');

    % Y tick labels (electrode or area names)
    if strcmpi(datatype,'time-frequency')
        ylabel('Frequency (Hz)','FontSize',12,'FontWeight','bold')
    else
         % Add labels to plot
        Ylabels = {chanlocs.labels};
        % img_prop = get(gca);
        % newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        newticks = 1:2:length(Ylabels);
        newticks = unique(newticks);
        Ylabels  = Ylabels(newticks);
        set(gca,'YTick',newticks,'YTickLabel', Ylabels,'FontWeight','normal');

        % For scalp channels
        ylabel('EEG channels','FontSize',12,'FontWeight','bold');

        % For Brain areas
        % for i = 1:length(chanlocs)
        %     chanlocs(i).labels = char(join(split(chanlocs(i).labels,'_')));
        % end
        % ylabel('Brain areas','FontSize',13,'FontWeight','bold');
        
    end

    % X label
    if contains(datatype,{'time','time-frequency'})
        xlabel('Time (ms)','FontSize',12,'FontWeight','bold')
    elseif strcmpi(datatype,'frequency')
        xlabel('Frequency (Hz)','FontSize',12,'FontWeight','bold')
    end


    correctoptions = {'Uncorrected' 'Max-corrected' 'Cluster-corrected' 'TFCE-corrected'};
    title(sprintf('%s (p<%g)', correctoptions{mcctype+1},alpha),'FontSize',12,'FontWeight','bold');

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

    % subplot(3,4,6);
    % peakChan = cluster_maxchan(maxEffect);
    % plot(xaxis,stats(peakChan,:),'LineWidth',2);
    % chanLabel = chanlocs(peakChan).labels;
    % title(sprintf('Course plot: %s',chanLabel),'FontSize',11,'fontweight','bold')
    % % plot(xaxis,stats(cluster_maxe,:),'LineWidth',2);  % plot peak effect of all clusters superimposed
    % % chanLabel = {chanlocs(cluster_maxe).labels};
    % % legend(chanLabel)
    % grid on; axis tight;
    % ylabel('t-values','FontSize',11,'fontweight','bold'); 
    % xlabel('Frequency (Hz)','FontSize',11,'fontweight','bold')

    % Plot bars of significnace for peak electrode
    % plotSigBar(mask(peakChan,:)~=0,xaxis);

    %% Scalp topography at peak latency/frequency (replace with 3D headplot?)

    % figure

    % subplot(3,3,6);
    % peakLat = cluster_maxfreq(maxEffect);
    % topoplot(stats(:, peakLat), chanlocs,'emarker',{'.','k',5,1}, ...             % normal electrodes
    %     'emarker2',{find(mask(:, peakLat)),'.',"red",7,1}, ...    % significant electrodes
    %     'verbose','off','colormap',dmap);                                                      % parameters
    % topoplot(stats(:, peakLat), chanlocs,'emarker2',{find(mask(:, peakLat)),'.','k',10,1}, ...    % significant electrodes
    %     'verbose','off','colormap',dmap);                                                      % parameters
    % if strcmpi(datatype, 'time')
    %     title(sprintf('Scalp topography: %g ms', xaxis(peakLat)), ...
    %         'FontSize',13,'fontweight','bold');
    % elseif strcmpi(datatype, 'frequency')
    %     title(sprintf('Scalp topography: %g Hz', xaxis(peakLat)), ...
    %         'FontSize',13,'fontweight','bold');
    % end
    % set(gcf,'Name','Topography at peak frequency','color','w','Toolbar','none','Menu','none','NumberTitle','Off')

    %% Adjust labels and ticks font and weight for all plots

    set(gcf,'Name','Results','color','w','Toolbar','none','Menu','none','NumberTitle','Off')
    % set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

else
    disp('No significant differences, nothing to plot.')
end

