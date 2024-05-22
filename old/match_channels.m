function LIMO = match_channels(stattest,analysis_type,LIMO)

% note channel can be a singleton or a vector (channel/component optimized analysis)
if strcmpi(analysis_type,'1 channel/component only')
    if isempty(LIMO.design.electrode) && strcmp(LIMO.Type,'Channels')
        channel = inputdlg('which electrode to analyse [?]','Electrode option');
    elseif isempty(LIMO.design.electrode) && strcmp(LIMO.Type,'Components')
        channel = inputdlg('which component to analyse [?]','Component option'); % can be 1 nb or a vector of channels (channel optimized analysis)
    else
        if ischar(LIMO.design.electrode)
            LIMO.design.electrode = load(LIMO.design.electrode);
            LIMO.design.electrode = LIMO.design.electrode.(cell2mat(fieldnames(LIMO.design.electrode)));
        end
        channel = {num2str(LIMO.design.electrode)}; % reformat temporarilly as if from inputdlg
    end

    if isempty(cell2mat(channel))
        [file,filepath,index] = uigetfile('*.mat',['select a ' LIMO.Type ' file']);
        if isempty(file) || index == 0
            return
        else
            % check the vector has the same length as the number of files
            channel_vector = load(fullfile(filepath,file));
            tmpname        = fieldnames(channel_vector);
            channel_vector = getfield(channel_vector,tmpname{1}); %#ok<GFLD>
            clear tmpname
            if length(channel_vector) ~= length(Paths)
                if exist('errordlg2','file')
                    errordlg2(['the nb of ' LIMO.Type ' does not match the number of subjects'],'Error');
                else
                    errordlg(['the nb of ' LIMO.Type ' does not match the number of subjects'],'Error');
                end
                return
            end

            % add the name to LIMO if absent
            if ~isfield(LIMO.design,'name')
                if stattest == 1
                    LIMO.design.name = ['one sample t-test one ' LIMO.Type(1:end-1)];
                elseif stattest == 2
                    LIMO.design.name = ['two samples t-test one ' LIMO.Type(1:end-1)];
                elseif stattest == 3
                    LIMO.design.name = ['paired t-test one ' LIMO.Type(1:end-1)];
                elseif stattest == 4
                    LIMO.design.name = ['regression analysis one ' LIMO.Type(1:end-1)];
                elseif stattest == 5
                    LIMO.design.name = ['AN(C)OVA analysis one ' LIMO.Type(1:end-1)];
                end
            end

            % restric the channels
            if strcmp(LIMO.Type,'Channels')
                LIMO.data.chanlocs          = LIMO.data.expected_chanlocs;
                LIMO.data.expected_chanlocs = LIMO.data.expected_chanlocs(channel_vector);
                LIMO.design.electrode       = channel_vector;
            else
                LIMO.data.chanlocs    = [];
                LIMO.design.component = channel_vector;
            end
        end

    elseif numel(str2num(channel{1})) || ...
            max(size(cell2mat(channel))) == numel(LIMO.data.data) || ...
            max(size(cell2mat(channel))) == sum(cellfun(@numel,LIMO.data.data))

        if ~isfield(LIMO.design,'name')
            if stattest == 1
                LIMO.design.name = ['one sample t-test one ' LIMO.Type(1:end-1)];
            elseif stattest == 2
                LIMO.design.name = ['two samples t-test one ' LIMO.Type(1:end-1)];
            elseif stattest == 3
                LIMO.design.name = ['paired t-test one ' LIMO.Type(1:end-1)];
            elseif stattest == 4
                LIMO.design.name = ['regression analysis one ' LIMO.Type(1:end-1)];
            else
                LIMO.design.name = ['AN(C)OVA analysis one ' LIMO.Type(1:end-1)];
            end
        end

        if strcmp(LIMO.Type,'Channels')
            LIMO.design.electrode       = str2num(cell2mat(channel));
            LIMO.data.chanlocs          = LIMO.data.expected_chanlocs;
            LIMO.data.expected_chanlocs = LIMO.data.expected_chanlocs(LIMO.design.electrode);
        else
            LIMO.design.component = eval(cell2mat(channel));
            LIMO.data.chanlocs    = [];
        end
    else
        error(['the nb of ' LIMO.Type ' does not match the number of subjects'],[LIMO.Type(1:end-1) ' error']);
    end

    % ---------------
else % Full scalp
    % ---------------

    if ~isfield(LIMO.design,'name')
        if stattest == 1
            LIMO.design.name = ['One sample t-test all ' LIMO.Type];
        elseif stattest == 2
            LIMO.design.name = ['Two samples t-test all ' LIMO.Type];
        elseif stattest == 3
            LIMO.design.name = ['Paired t-test all ' LIMO.Type];
        elseif stattest == 4
            LIMO.design.name = ['Regression analysis all ' LIMO.Type];
        else
            LIMO.design.name = ['AN(C)OVA analysis all ' LIMO.Type];
        end
    end

    if isfield(LIMO.data,'expected_chanlocs')
        LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
    end

    if strcmpi(LIMO.Type,'Components')
        LIMO.design.component = [];
    else
        LIMO.design.electrode = [];
    end
end
end
