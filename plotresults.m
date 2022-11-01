%% Plots results using correction for multiple comparison masks
%
% Cedric Cannard, Sep 2022

function plotresults(datatoplot, mask, pcorr, mcctype, chanlocs)


if isempty(mask)
    disp('no values computed'); return
elseif sum(mask(:)) == 0
    warndlg('  no values under threshold  ','no significant effect','modal');
else
    assignin('base','p_values', pcorr)
    assignin('base','mask', mask)
    assignin('base','stat_values', datatoplot)
end


% if ndims(datatoplot) == 3
%     limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
% else
%     limo_display_image(LIMO,toplot,mask,mytitle,flag)
% end


%%%%%%% limo_display_image %%%%%%%%%%%
% Look at limo_display_image_tf for time-frequency data

% Extract significant t-values for plotting
effects = datatoplot.*single(mask>0);
effects(effects==0) = NaN;
cc = color_images(effects); % get a color map commensurate to that
v = max(effects(:));       % from the 2D data to plot, find max
[e,f]=find(effects==v);    % which channel and time/frequency frame
if length(e)>1           % if we have multiple times the exact same max values
    e = e(1); f = f(1);  % then take the 1st (usually an artefact but allows to see it)
end

if sum(sum(mask)) > 0
    figure('Color','w'); %set(gcf,'InvertHardCopy','off');

    % find max effect for topo and course plots
    v = max(effects(:));            % find max from the 2D data to plot
    [e,f] = find(effects == v);     % corresponding channel and time/frequency
    if length(e) > 1                % take the 1st one if there are multiple same max values (usually an artefact but allows to see it)
        e = e(1); 
        f = f(1);
    end

    % Freq vector
    freqvect = 1;
    if size(freqvect,2) ~= size(datatoplot,2)
        freqvect = linspace(1,size(datatoplot,2),size(datatoplot,2));
    end

    % Course plot at max electrode
    if ~isnan(e)
        subplot(3,3,9);
        plot(freqvect,datatoplot(e,:),'LineWidth',2);
        label = chanlocs(e).labels;
        if ~iscell(label)
            mytitle2 = sprintf('Diff values: %s', label);
        else
            mytitle2 = sprintf('Diff values: %s', char(label));
        end
        title(mytitle2,'FontSize',11)
        %  title(mytitle,'FontSize',11)
        % set(gca,'xtick',freqvect,'xticklabel',freqlist);
        grid on; axis tight;
    end

    % Topoplot at max time/frequency
    if ~isnan(e)
        subplot(3,3,6);
        opt = {'maplimits','maxmin','verbose','off','colormap', limo_color_images(datatoplot)};
        topoplot(datatoplot(:,f),chanlocs,opt{:});
        % colormap('parula')
        % title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
        if size(datatoplot,2) == 1
            title('Topoplot','FontSize',12)
        else
            title(['Topography: ' num2str(freqvect(f)) ' Hz'],'FontSize',10);
            %     title(['Topography: ' freqlist(f)],'FontSize',10);
            set(gca,'XTickLabel', freqvect);
        end
    end

    % Images all
    subplot(3,3,[1 2 4 5 7 8]);
    imagesc(freqvect,1:size(datatoplot,1),effects);
    %     colormap('parula')
    colormap(gca, cc);
    set_imgaxes('Frequency', 'Channels',chanlocs, effects);
    correctoptions = {'Uncorrected' 'Cluster-corrected' 'TFCE-corrected' 'Max-corrected'};
    title(['All electrodes and all bands' correctoptions{mcctype+1}],'Fontsize',12)
    colorbar; %clim([0 0.1])
    % ylabel({chanlocs.labels})
    set(gca,'ytick',1:size(datatoplot,1),'yticklabel',{chanlocs.labels});

    % Return cluster info
    warning off
    % for each cluster, get start/end/max value
    % if unthresholded, uncorrected, tfce or max = mask is made up of ones
    n_cluster     = max(mask(:));
    cluster_start = NaN(1,n_cluster); % start of each cluster
    cluster_end   = NaN(1,n_cluster); % end of each cluster
    cluster_maxv  = NaN(1,n_cluster); % max value for each cluster
    cluster_maxe  = NaN(1,n_cluster); % channel location of the max value of each cluster
    cluster_maxf  = NaN(1,n_cluster); % frame location of the max value of each cluster
    for c=1:n_cluster
        tmp = datatoplot.*(mask==c);
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
    % report in command window
    % if contains(mytitle,'cluster')
    for c = 1:n_cluster
        fprintf('cluster %g starts at %gHz ends at %gHz, max %g @ %gHz channel %s \n', c, ...
            freqvect(cluster_start(c)),freqvect(cluster_end(c)), cluster_maxv(c), freqvect(cluster_maxf(c)), chanlocs(cluster_maxe(c)).labels);
    end
    %         end
    warning on

else
    disp('No significant differences.')
end
