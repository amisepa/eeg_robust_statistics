%% Plots results using correction for multiple comparison masks
%
% Usage:
%    [cluster_bounds, cluster_maxchan, cluster_maxfreq, cluster_order] = plot_results(datatype, chan_type, xaxis, stats, mask, pcorr, alpha, chanlocs, mcctype)
%
% Example:
%   [cluster_bounds, cluster_maxchan, cluster_maxfreq, cluster_order] = plot_results('time', 'scalp', times, tvals, mask, pcorr, 0.05, chanlocs, 0);
%
% Adapted from limo_display_image.
% 
% Copyright (C) Cedric Cannard, Created on Sep 2022

function [cluster_bounds, cluster_maxchan, cluster_maxfreq, cluster_maxval, mask_clusters] = plot_results(datatype, chan_type, xaxis, stats, mask, pcorr, alpha, chanlocs, mcctype)

cluster_bounds = [];
cluster_maxchan = [];
cluster_maxfreq = [];
cluster_maxval = [];
mask_clusters = {};  % Initialize cell array to hold binary masks for each cluster

if isempty(mask)
    disp('Empty mask, computing one from corrected p-values.');
    mask = pcorr < alpha;
end

if sum(mask,'all') == 0
    disp('No significant differences, nothing to plot.')
    return
end

% Get clusters properties (start/end, peak, channel, frame/freq)
n_cluster     = max(mask(:));
cluster_start = nan(1,n_cluster);
cluster_end   = nan(1,n_cluster);
cluster_maxval  = nan(1,n_cluster);
cluster_maxchan = nan(1,n_cluster);
cluster_maxfreq = nan(1,n_cluster);

for iClust = 1:n_cluster
    % Binary mask for this cluster
    mask_i = (mask == iClust);
    mask_clusters{iClust} = mask_i;

    % Restrict stats to this cluster only
    tmp = stats .* mask_i;
    tmp(tmp == Inf | tmp == -Inf) = NaN;

    % Get bounds of the cluster
    sigframes = sum(tmp,1);
    cluster_start(iClust) = find(sigframes, 1, 'first');
    cluster_end(iClust) = find(sigframes, 1, 'last');

    % Get max value and location
    [V, type] = max([abs(min(tmp(:))) max(tmp(:))]);
    if type == 2
        cluster_maxval(iClust) = V(1);
    else
        V = -V;
        cluster_maxval(iClust) = V(1);
    end
    [cluster_maxchan(iClust), cluster_maxfreq(iClust)] = ind2sub(size(tmp), find(tmp == V(1), 1));
end

% Convert frequency or latency values at cluster peaks
cluster_peak_xvals = xaxis(cluster_maxfreq);

% Sort clusters by frequency or time
[~, idx] = sort(cluster_peak_xvals, 'ascend');
cluster_start     = cluster_start(idx);
cluster_end       = cluster_end(idx);
cluster_maxchan   = cluster_maxchan(idx);
cluster_maxfreq   = cluster_maxfreq(idx);
cluster_maxval    = cluster_maxval(idx);
mask_clusters     = mask_clusters(idx);  % sort masks accordingly

% Output cluster bounds
cluster_bounds = [cluster_start; cluster_end]';

% Print info
if ~isempty(chanlocs)
    for iClust = 1:n_cluster
        if strcmpi(datatype, 'time')
            fprintf('Cluster %g: %g to %g ms. Peak effect: %s at %g ms (t = %g) \n', ...
                iClust, xaxis(cluster_start(iClust)), xaxis(cluster_end(iClust)), ...
                chanlocs(cluster_maxchan(iClust)).labels, xaxis(cluster_maxfreq(iClust)), round(cluster_maxval(iClust),1));
        elseif strcmpi(datatype, 'frequency')
            fprintf('Cluster %g: %g to %g Hz. Peak effect: %s at %g Hz (t = %g) \n', ...
                iClust, xaxis(cluster_start(iClust)), xaxis(cluster_end(iClust)), ...
                chanlocs(cluster_maxchan(iClust)).labels, xaxis(cluster_maxfreq(iClust)), round(cluster_maxval(iClust),1));
        end
    end
end


%% MAIN PLOT

% figure('Color','w','InvertHardCopy','off');
figure('Color','w'); 
% subplot(3, 3, [1 2 4 5 7 8])

% Grey-out non-significant values
stats(mask==0) = NaN;

% Image all channels and time/frequency data
imagesc(xaxis, 1:size(stats,1), stats);

% Color map
load("colormap_bwr.mat");
if sum(isnan(stats(:))) ~= 0
    dmap(1,:) = [.9 .9 .9]; % set NaNs to gray
end
colormap(gca, dmap); % colormap('parula')
c = colorbar;
ylabel(c, 'T-values','FontWeight','bold','FontSize',11,'Rotation',-90)

% Y tick labels (electrode or area names)
if strcmpi(datatype,'time-frequency')
    ylabel('Frequency (Hz)','FontSize',12,'FontWeight','bold')
else

    % Edit area names when necessary
    if ~isempty(chanlocs)
        % for i = 1:length(chanlocs)
        %     chanlocs(i).labels = char(join(split(chanlocs(i).labels,'_')));  % to remove dashes in some atlas
        %     % chanlocs(i).labels = extractBefore(chanlocs(i).labels,' ');
        % end

        % % Add labels to plot
        % Yticks = {chanlocs.labels};
        % img_prop = get(gca);
        % newticks = round(linspace(1,length(Yticks),length(img_prop.YTick)*2));
        % % if length({chanlocs.labels}) < 20
        % %     newticks = 1:length(Yticks);
        % % else
        % %     newticks = 1:2:length(Yticks);
        % % end
        % % newticks = unique(newticks);
        % % Yticks  = Yticks(newticks);
        % % set(gca,'YTick',newticks,'YTickLabel', Yticks,'FontWeight','normal');
        % set(gca,'YTick',newticks,'YTickLabel', Yticks,'FontWeight','normal');

        % Add all labels, assuming each row of stats corresponds to a channel
        Yticks = {chanlocs.labels};
        skip = 2;  % 1 = no skip
        set(gca, 'YTick', 1:skip:length(Yticks), 'YTickLabel', Yticks(1:skip:end), 'FontWeight','normal');


        % Ylabel
        if strcmpi(chan_type, 'scalp')
            ylabel('EEG channels','FontSize',12,'FontWeight','bold');
        else
            ylabel('Brain areas','FontSize',12,'FontWeight','bold');
        end
    end
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
maxval = max(abs(stats(:)));
% if max(stats(:)) < 0
%     clim([-maxval 0])
% elseif min(stats(:)) > 0
%     clim([0 maxval])
% else
clim([-maxval maxval])
% end

axis tight; box on; set(gca,'LineWidth',1)

%% Scalp topography at peak latency/frequency (replace with 3D headplot?)

% if ~isempty(chanlocs) && strcmpi(chan_type, 'scalp')
%         subplot(3,3,6)
%         peakLat = cluster_maxfreq(maxEffect); % peak latency(for ERP) or frequency (power spectrum)
% 
%   BEST WITH MASK
    % figure('Color','w'); topoplot(tvals(:, peakFreqs(iClust)), chanlocs, ...
    %     'emarker2',{find(mask(:,peakFreqs(iClust))),'.','k',12,1}, ...    % significant electrodes
    %     'verbose','off','colormap',dmap, 'whitebk','on', ...
    %     'pmask', mask(:,peakFreqs(iClust))); 

%         topoplot(stats(:, peakLat), chanlocs, 'emarker', {'.','k',5,1}, ...             % normal electrodes
%             'emarker2', { find(mask(:, peakLat)==maxEffect),'.',"red",7,1 }, ...    % significant electrodes
%             'verbose','off','colormap',dmap);                                                      % parameters
%         % topoplot(stats(:, peakLat), chanlocs,'emarker2',{find(mask(:, peakLat)==maxEffect),'.','k',10,1}, ...    % significant electrodes
%         %     'verbose','off','colormap',dmap);                                                      % parameters
%         if strcmpi(datatype, 'time')
%             title(sprintf('Scalp topography: %g ms', xaxis(peakLat)), ...
%                 'FontSize',13,'fontweight','bold');
%         elseif strcmpi(datatype, 'frequency')
%             title(sprintf('Scalp topography: %g Hz', xaxis(peakLat)), ...
%                 'FontSize',13,'fontweight','bold');
%         end
%         % set(gcf,'Name','Topography at peak frequency','color','w','Toolbar','none','Menu','none','NumberTitle','Off')
%     % else
%     %     try
%     %         % Cortex surface plot at peak frequency
%     %         load cm17.mat; 
%     %         load cortex  % for cortex_highres
%     %         cortex = cortex_highres;
%     %         % cortex = tmp_upsight.roi.cortex;
%     % 
%     %         % figure('color','w','Name',sprintf('Cluster %g',iClust),'color','w','Toolbar','none','Menu','none','NumberTitle','Off');
%     %         roiPSD = stats(:, peakLat);  % t-values for peak fequency
%     %         roiPSD(mask(:, peakLat)~=maxEffect) = NaN;   % NaN non-significant values
%     %         allplots_cortex_BS(cortex_highres, roiPSD, [min(roiPSD) max(roiPSD)], cm17b, 't-values', 0.35, [], {});
%     %     catch
%     %         warning("Cortex plot failed. Make sure roiconnect plugin is installed and eeglab launched?")
%     %     end
%     % 
%     % end
% end

%% Course plot of t-values of strongest effect

% % if strcmpi(chan_type, 'scalp')
%     subplot(3,3,9)
% 
%     % figure('Color','w')
%     % hold on
%     % subplot(2,1,1)
% 
%     peakChan = cluster_maxchan(maxEffect);
%     plot(xaxis,stats(peakChan,:),'LineWidth',2);
%     if ~isempty(chanlocs)
%         chanLabel = chanlocs(peakChan).labels;
%         title(sprintf('Course plot: %s',chanLabel),'FontSize',11,'fontweight','bold')
%     end
%     % plot(xaxis,stats(cluster_maxe,:),'LineWidth',2);  % plot peak effect of all clusters superimposed
%     % chanLabel = {chanlocs(cluster_maxe).labels};
%     % legend(chanLabel)
%     grid off; axis tight;
%     ylabel('t-values','FontSize',11,'fontweight','bold');
%     xlabel('Frequency (Hz)','FontSize',11,'fontweight','bold')
% 
%     % Plot bars of significnace for peak electrode
%     plotSigBar(mask(peakChan,:)~=0,xaxis);
% % end


%% Adjust labels and ticks font and weight for all plots

% set(gcf,'Name','Results from mass-univariate analysis','color','w','Toolbar','none','Menu','none','NumberTitle','Off')
set(gcf,'Name','Results from mass-univariate analysis','color','w','NumberTitle','Off')
set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');



