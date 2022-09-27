function [data,removed] = getdata(stattest,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO)

data = [];
removed = [];
disp('gathering data ...');
if stattest == 1 % one sample
    index = 1;
    if all(size(LIMO.data.data)==[1 1]) % cell of cell
        LIMO.data.data = LIMO.data.data{1};
    end

    for i=size(LIMO.data.data,2):-1:1 % for each subject
        tmp = load(char(LIMO.data.data{i}));

        % get indices to trim data
        if strcmp(LIMO.Analysis,'Time-Frequency')
            if contains(LIMO.data.data{i},'Betas')
                tmp = tmp.Betas;
            else
                tmp = tmp.con(:,:,:,1);
            end
            begins_at = fliplr((max(first_frame) - first_frame(i,:) + 1)); % returns time/freq/or freq-time
            ends_at(1) = size(tmp,2) - (last_frame(i,2) - min(last_frame(:,2)));
            ends_at(2) = size(tmp,3) - (last_frame(i,1) - min(last_frame(:,1)));
        else
            if contains(LIMO.data.data{i},'Betas')
                tmp = tmp.Betas;
            else
                tmp = tmp.con(:,:,1);
            end
            begins_at = max(first_frame) - first_frame(i) + 1;
            ends_at = size(tmp,2) - (last_frame(i) - min(last_frame));
        end

        % data dim [channel, freq/time, param, nb subjects]
        if isempty(LIMO.Type)
            if exist('questdlg2','file')
                LIMO.Type = questdlg2('Is the analysis on','Type is empty','Channels','Components','Channels');
            else
                LIMO.Type = questdlg('Is the analysis on','Type is empty','Channels','Components','Channels');
            end
        end

        if strcmp(analysis_type,'Full scalp analysis')
            if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                else
                    data(:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                end
            elseif strcmpi(LIMO.Type,'Components')
                try
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        data(:,:,:,:,index) = tmp(:,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                    else
                        data(:,:,:,index) = tmp(:,begins_at:ends_at,:);
                    end
                catch dim_error
                    if strcmp(dim_error,'Subscripted assignment dimension mismatch.')
                        disp('you are trying to match matrices of ICs of different size')
                        if isempty(expected_chanlocs)
                            disp('either cluster data are run 1st level batch, or input cluster file')
                        end
                    end
                end
            end
            index = index + 1;
            removed(i) = 0;

        elseif strcmp(analysis_type,'1 channel/component only') %&& size(subj_chanlocs(i).chanlocs,2) == size(tmp,1)

            % Use single channel
            if  length(LIMO.design.electrode) == 1
                if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(i).chanlocs) == size(tmp,1)
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta
                        else
                            data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                        else
                            data(1,:,index) = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                        end
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = tmp(LIMO.design.component,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        else
                            data(1,:,:,index) = tmp(LIMO.design.component,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = tmp(LIMO.design.component,begins_at:ends_at,:);
                        else
                            data(1,:,index) = tmp(LIMO.design.component,begins_at:ends_at,:);
                        end
                    end
                end
                index = index + 1;
                removed(i) = 0;

            else  % Use multiple single channels
                if strcmpi(LIMO.Type,'Channels')
                    out = limo_match_elec(subj_chanlocs(i).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = out(i,:,:,:);
                        else
                            data(1,:,:,index) = out(i,:,:,:);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                        else
                            data(1,:,index) = out(i,:,:); % matches the expected chanloc of the subject
                        end
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,:,index) = tmp(LIMO.design.component(i),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        else
                            data(1,:,:,index) = tmp(LIMO.design.component(i),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                        end
                    else
                        if contains(LIMO.data.data{i},'Betas')
                            data(1,:,:,index) = tmp(LIMO.design.component(i),begins_at:ends_at,:);
                        else
                            data(1,:,index) = tmp(LIMO.design.component(i),begins_at:ends_at,:);
                        end
                    end
                end
                index = index +1;
                removed(i) = 0;
            end
        else
            fprintf('subject %g discarded, channel description and data size don''t match \n',i);
            removed(i) = 1;
        end
        clear tmp
    end

elseif stattest == 2 % several samples
    subject_nb = 1;
    for igp = 1:length(LIMO.data.data)
        index = 1;
        for i=1:size(LIMO.data.data{igp},2) % for each subject per group
            tmp = load(cell2mat(LIMO.data.data{igp}(i)));

            % get indices to trim data
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                    tmp = tmp.Betas;
                else
                    tmp = tmp.con(:,:,:,1);
                end
                begins_at = fliplr((max(first_frame) - first_frame(subject_nb,:) + 1)); % returns time/freq/or freq-time
                ends_at(1) = size(tmp,2) - (last_frame(subject_nb,2) - min(last_frame(:,2)));
                ends_at(2) = size(tmp,3) - (last_frame(subject_nb,1) - min(last_frame(:,1)));
            else
                if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                    tmp = tmp.Betas;
                else
                    tmp = tmp.con(:,:,1);
                end
                begins_at = max(first_frame) - first_frame(subject_nb) + 1;
                ends_at = size(tmp,2) - (last_frame(subject_nb) - min(last_frame));
            end

            if strcmpi(analysis_type,'Full scalp analysis') %&& size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)
                if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(subject_nb).chanlocs) == size(tmp,1)
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        tmp_data(:,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                    else
                        tmp_data(:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp);
                    end
                elseif strcmpi(LIMO.Type,'Components')
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        tmp_data(:,:,:,:,index) = tmp(:,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:);
                    else
                        tmp_data(:,:,:,index) = tmp(:,begins_at:ends_at,:);
                    end

                end
                removed(igp,i) = 0;
                index = index + 1;

            elseif strcmpi(analysis_type,'1 channel/component only') %&& size(subj_chanlocs(subject_nb).chanlocs,2) == size(tmp,1)

                % Use single channel
                if length(LIMO.design.electrode) == 1
                    if strcmpi(LIMO.Type,'Channels') && length(subj_chanlocs(subject_nb).chanlocs) == size(tmp,1)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,:,1,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,1,index) = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % all param for beta, if con, adjust dim
                            end
                        end
                    elseif strcmpi(LIMO.Type,'Components')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = tmp(LIMO.design.electrode,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,:,1,index) = tmp(LIMO.design.electrode,begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % all param for beta, if con, adjust dim
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = tmp(LIMO.design.electrode,begins_at:ends_at,:); % all param for beta, if con, adjust dim
                            else
                                tmp_data(1,:,1,index) = tmp(LIMO.design.electrode,begins_at:ends_at,:); % all param for beta, if con, adjust dim
                            end
                        end
                    end
                    removed(igp,i) = 0;
                    index = index + 1;

                else % Use multiple single channels
                    if strcmpi(LIMO.Type,'Channels')
                        out = limo_match_elec(subj_chanlocs(subject_nb).chanlocs,LIMO.data.expected_chanlocs,begins_at,ends_at,tmp); % out is for all expected chanlocs, ie across subjects
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = out(subject_nb,:,:,:); % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,:,1,index) = out(subject_nb,:,:,:); % matches the expected chanloc of the subject
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = out(subject_nb,:,:);     % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,1,index) = out(subject_nb,:,:);     % matches the expected chanloc of the subject
                            end
                        end
                    elseif strcmpi(LIMO.Type,'Components')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,:,index) = tmp(LIMO.design.electrode(subject_nb),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,:,1,index) = tmp(LIMO.design.electrode(subject_nb),begins_at(1):ends_at(1),begins_at(2):ends_at(2),:); % matches the expected chanloc of the subject
                            end
                        else
                            if contains(cell2mat(LIMO.data.data{igp}(i)),'Betas')
                                tmp_data(1,:,:,index) = tmp(LIMO.design.electrode(subject_nb),begins_at:ends_at,:);     % matches the expected chanloc of the subject
                            else
                                tmp_data(1,:,1,index) = tmp(LIMO.design.electrode(subject_nb),begins_at:ends_at,:);     % matches the expected chanloc of the subject
                            end
                        end
                    end
                    removed(igp,i) = 0;
                    index = index +1;
                end
            else
                fprintf('subject %g of group %g discarded, channel description and data size don''t match \n',i, igp);
                removed(igp,i) = 1;
            end
            clear tmp
            subject_nb = subject_nb + 1;
        end
        data{igp} = tmp_data;
        clear tmp tmp_data
    end
end
end
