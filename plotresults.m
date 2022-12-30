%% Plots results using correction for multiple comparison masks
%
% Cedric Cannard, Sep 2022

function plotresults(datatype, xaxis, stats, mask, pcorr, mcctype, chanlocs)

if isempty(mask)
    disp('Empty mask, computing one from corrected p-values.'); %return
    mask = pcorr < 0.05;
end
if sum(mask(:)) == 0
    warndlg('no values under threshold');
end

% if ndims(datatoplot) == 3
%     limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
% else
%     limo_display_image(LIMO,toplot,mask,mytitle,flag)
% end
%%%%%%% limo_display_image %%%%%%%%%%%
% Look at limo_display_image_tf for time-frequency data

if sum(sum(mask)) > 0
    figure('Color','w'); set(gcf,'InvertHardCopy','off');

    % Extract significant t-values for plotting
    effects = stats.*single(mask>0);
    effects(effects==0) = NaN;
%     v = max(effects(:));  
%     [c, tf_stamp] = find(effects==v); % which channel and time/frequency frame
%     if length(c) > 1  % if we have multiple times the exact same max values
%         c = c(1); 
%         tf_stamp = tf_stamp(1);  % then take the 1st (usually an artifact but allows to see it)
%     end
 
    % Image all channels and time/frequency data
    subplot(3,3,[1 2 4 5 7 8]);
    imagesc(xaxis,1:size(stats,1),effects);
    cc = color_images(effects);
    colormap(gca, cc);
%     colormap('parula')
    colorbar; % clim([0 0.1])
    % ylabel({chanlocs.labels})
    if strcmpi(datatype, 'time')
        set_imgaxes('Time (ms)', 'Channels', chanlocs, effects);
    elseif strcmpi(datatype, 'frequency')
        set_imgaxes('Frequency (Hz)', 'Channels', chanlocs, effects);
    end
    set(gca,'ytick',1:size(stats,1),'yticklabel',{chanlocs.labels});
    correctoptions = {'Uncorrected' 'Cluster-corrected' 'TFCE-corrected' 'Max-corrected'};
    title(['All electrodes and all bands (' correctoptions{mcctype+1} ')'],'Fontsize',12)
    
    % Return cluster info
    % if unthresholded, uncorrected, tfce or max = mask is made up of ones
    n_cluster     = max(mask(:));
    cluster_start = nan(1,n_cluster); % start of each cluster
    cluster_end   = nan(1,n_cluster); % end of each cluster
    cluster_maxv  = nan(1,n_cluster); % max value for each cluster
    cluster_maxe  = nan(1,n_cluster); % channel location of the max value of each cluster
    cluster_maxf  = nan(1,n_cluster); % frame location of the max value of each cluster
    warning off
    for c = 1:n_cluster
        tmp = stats.*(mask==c);
        tmp(tmp==Inf) = NaN;
        tmp(tmp==-Inf) = NaN;
        sigframes = sum(tmp,1);
        cluster_start(c) = find(sigframes,1,'first');
        cluster_end(c) = find(sigframes,1,'last');
        [V,type] = max([abs(min(tmp(:))) max(tmp(:))]);
        if type == 2
            cluster_maxv(c) = V(1);
        else
            V = -V;
            cluster_maxv(c) = V(1);
        end
        [cluster_maxe(c),cluster_maxf(c)] = ind2sub(size(tmp),find(tmp==V(1)));
    end
    for c = 1:n_cluster
        if strcmpi(datatype, 'time')
            fprintf('Cluster %g: from %g to %g ms. Peak effect at channel %s %g ms (t = %g) \n', ...
                c, xaxis(cluster_start(c)), xaxis(cluster_end(c)), ...
                chanlocs(cluster_maxe(c)).labels, xaxis(cluster_maxf(c)), round(cluster_maxv(c),1) );
        elseif strcmpi(datatype, 'frequency')
            fprintf('Cluster %g: from %g to %g Hz. Peak effect at channel %s %g Hz (t = %g) \n', ...
                c, xaxis(cluster_start(c)), xaxis(cluster_end(c)), ...
                chanlocs(cluster_maxe(c)).labels, xaxis(cluster_maxf(c)), round(cluster_maxv(c),1) );
        end
    end
    warning on
    
    % find max effect for topo and course plots
%     v = max(effects(:));                % find max from the 2D stats data
%     [c,tf_stamp] = find(effects == v);  % corresponding channel and time/frequency stamp

    % Take the 1st one if there are multiple same max values (usually an 
    %   artifact but allows to see it)
%     if length(c) > 1                
%         c = c(1); 
%         tf_stamp = tf_stamp(1);
%     end

    % Freq vector
%     xaxis = 1;
%     if size(freqvect,2) ~= size(stats,2)
%         xaxis = linspace(1,size(stats,2),size(stats,2));
%     end

    % Course plot at peak electrode from cluster 1
    if ~isnan(c)
        subplot(3,3,9);
        maxElec = cluster_maxe(c);
        plot(xaxis,stats(maxElec,:),'LineWidth',2);
        label = chanlocs(maxElec).labels;
%         plot(xaxis,stats(c,:),'LineWidth',2);
%         label = chanlocs(c).labels;
        if ~iscell(label)
            mytitle2 = sprintf('Channel %s', label);
        else
            mytitle2 = sprintf('Channel %s', char(label));
        end
        title(mytitle2,'FontSize',11)
        grid on; axis tight;
        ylabel('t-values')
    end

    % Topoplot at max time/frequency
    if ~isnan(c)
        subplot(3,3,6);
%         opt = {'maplimits','maxmin','verbose','off','colormap', limo_color_images(stats)};
        opt = {'maplimits','maxmin','verbose','off','colormap',color_images(effects)};
        topoplot(stats(:,cluster_maxf(c)),chanlocs,opt{:});
        % colormap('parula')
        % title(['topoplot @' num2str(round(xaxis(f))) 'Hz'],'FontSize',12);
        if size(stats,2) == 1
            title('Topoplot','FontSize',12)
        else
            if strcmpi(datatype, 'time')
                title(['Topography: ' num2str(xaxis(cluster_maxf(c))) ' ms'],'FontSize',10);
            elseif strcmpi(datatype, 'frequency')
                title(['Topography: ' num2str(xaxis(cluster_maxf(c))) ' Hz'],'FontSize',10);
            end
            % title(['Topography: ' freqlist(f)],'FontSize',10);
            set(gca,'XTickLabel', xaxis);
        end
    end

else
    disp('No significant differences.')
end
