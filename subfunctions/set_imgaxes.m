function set_imgaxes(analysis, type, chanlocs, effects)
    
%     analysis = 'Frequency';
%     type = 'Channels';

    img_prop = get(gca);
    set(gca,'LineWidth',2)

    % ----- X --------
    if strcmpi(analysis,'Time') || strcmpi(analysis,'Time-Frequency')
        xlabel('Time in ms','FontSize',10)
    elseif strcmpi(analysis,'Frequency')
        xlabel('Frequency in Hz','FontSize',10)
    end

    % ----- Y --------
    if strcmpi(analysis,'Time-Frequency')
        ylabel('Frequency in Hz','FontSize',10)
    else
        if strcmpi(type,'Components')
            if size(effects,1) == 1
                ylabel('Optimized component','FontSize',10);
            else
                ylabel('Components','FontSize',10);
            end
        else
            if size(effects,1) == 1
                ylabel('Optimized channel','FontSize',10);
            else
                ylabel('Channels','FontSize',10);
            end
        end

%         if isfield(LIMO.data, 'chanlocs')
        Ylabels = arrayfun(@(x)(x.labels), chanlocs, 'UniformOutput', false);
%         else
%             Ylabels = arrayfun(@(x)(x.labels), LIMO.data.expected_chanlocs, 'UniformOutput', false);
%         end

        newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        newticks = unique(newticks);
        Ylabels  = Ylabels(newticks);
        if size(effects,1) == 1
            set(gca,'YTick',1);
        else
            set(gca,'YTick',newticks);
            set(gca,'YTickLabel', Ylabels);
        end
    end

    % ----- Colormap --------
    try
        maxval = max(abs(max(effects(:))),abs(min(effects(:))));
        if max(effects(:)) < 0
            caxis([-maxval 0])
        elseif min(effects(:)) > 0
            caxis([0 maxval])
        else
            caxis([-maxval maxval])
        end
    catch caxiserror
        fprintf('axis issue: %s\n',caxiserror.message)
    end
