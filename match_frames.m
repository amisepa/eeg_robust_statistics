function [first_frame,last_frame,subj_chanlocs,channeighbstructmat,LIMO] = match_frames(Paths,LIMO)

% once we have all the files, we need to collect information to match the frames across subjects
% OUTPUT: first_frame and last _frame returns the beginning and end in terms of indices
%                                     for time-frequency these are vectors with time then frequency
%         subj_chanlocs the chanlocs per subjects
%         channeighbstructmat the neighbourg matrices
%
% the LIMO structure is also updated the reflect the smallest interval(s) across subjects,
% which is used for the second leve analysis

channeighbstructmat = [];
disp('match frames between subjects ...')
% check Paths format
if iscell(Paths{1})
    tmp = Paths; clear Paths
    index = 1;
    for gp=1:size(tmp,2)
        for s=1:size(tmp{gp},2)
            Paths{index} = tmp{gp}(s);
            index = index + 1;
        end
    end
end

% now loop loading the LIMO.mat for each subject to collect information
% ---------------------------------------------------------------------

ME = [];
for i=size(Paths,2):-1:1
    try
        if iscell(Paths{i})
            cd (cell2mat(Paths{i}));
        else
            cd (Paths{i});
        end
        limo = load('LIMO.mat');
        limo = limo.LIMO;
    catch
        if iscell(Paths{i})
            p=cell2mat(Paths{i});
        else
            p=Paths{i};
        end
        error('cannot load/find the LIMO.mat in %s',p)
    end

    if i==size(Paths,2)
        Analysis = limo.Analysis;
    else
        if ~strcmp(limo.Analysis,Analysis)
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
        end
    end

    sampling_rate(i) = limo.data.sampling_rate;
    if strcmpi(LIMO.Type,'Channels')
        subj_chanlocs(i).chanlocs = limo.data.chanlocs;
        if isfield(limo.data,'channeighbstructmat')
            channeighbstructmat = limo.data.channeighbstructmat;
        end
    else
        subj_chanlocs(i).chanlocs = [];
    end

    if strcmp(Analysis,'Time-Frequency')
        first_frame(i,1)            = limo.data.trim1;
        last_frame(i,1)             = limo.data.trim2;
        start(i,1)                  = limo.data.start;
        stop(i,1)                   = limo.data.end;

        first_frame(i,2)            = limo.data.trim_lowf;
        last_frame(i,2)             = limo.data.trim_highf;
        start(i,2)                  = limo.data.tf_freqs(1);
        stop(i,2)                   = limo.data.tf_freqs(end);

        tf_times{i}(1,:)            = limo.data.tf_times;
        tf_freqs{i}(1,:)            = limo.data.tf_freqs;
    else
        first_frame(i)              = limo.data.trim1;
        last_frame(i)               = limo.data.trim2;
        start(i)                    = limo.data.start;
        stop(i)                     = limo.data.end;

        if strcmp(Analysis,'Frequency')
            freqlist{i}(1,:)        = limo.data.freqlist;
        end
    end
end
clear limo

% quick check things are ok
if  strcmpi(LIMO.Type,'Channels') && ~isempty(ME) && isempty(LIMO.data.neighbouring_matrix)
    error('some subject(s) have a different channel structure \nplease load an expected chanloc when choosing a test');
end

if length(unique(sampling_rate)) ~= 1
    error('data have different sampling rates')
end

% match and return into LIMO
LIMO.Analysis           = Analysis;
LIMO.data.sampling_rate = sampling_rate(1);

% we need 1) to find the highest start in time and freq 2) the lowest end
% in time and freq and 3) match that on freqlist or tf_times/tf_freqs

[v,c] = max(first_frame);
if strcmp(Analysis,'Time-Frequency')
    LIMO.data.trim1 = v(1);
    LIMO.data.start = start(c(1),1);
    LIMO.data.trim_lowf = v(2);
    LIMO.data.lowf = start(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-LIMO.data.start));
        tf_times{i} = tf_times{i}(ind:end);
        [~,ind] = min(abs(tf_freqs{i}-LIMO.data.lowf));
        tf_freqs{i} = tf_freqs{i}(ind:end);
    end
else
    LIMO.data.trim1 = v;
    LIMO.data.start = start(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-LIMO.data.start));
            freqlist{i} = freqlist{i}(ind:end);
        end
    end
end

[v,c] = min(last_frame);
if strcmp(Analysis,'Time-Frequency')
    LIMO.data.trim2 = v(1);
    LIMO.data.end = stop(c(1),1);
    LIMO.data.trim_highf = v(2);
    LIMO.data.highf = stop(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-LIMO.data.end));
        tf_times{i} = tf_times{i}(1:ind);
        [~,ind] = min(abs(tf_freqs{i}-LIMO.data.highf));
        tf_freqs{i} = tf_freqs{i}(1:ind);
    end
else
    LIMO.data.trim2 = v;
    LIMO.data.end = stop(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-LIMO.data.end));
            freqlist{i} = freqlist{i}(1:ind);
        end
    end
end

% finally match everything
if strcmp(Analysis,'Frequency')
    % check all lists match ; if sampled the same across subject, all same sizes
    try
        freqlist = cell2mat(freqlist');
    catch list_issue
        assignin('base','freqlist',freqlist)
        error('the resolution of frequency lists doesn''t match between subjects \n%s',list_iisue.message)
    end
    LIMO.data.freqlist = mean(freqlist,1);
    LIMO.data.start    = LIMO.data.freqlist(1);
    LIMO.data.end      = LIMO.data.freqlist(end);

elseif strcmp(Analysis,'Time-Frequency')
    % check all lists match
    try
        tf_times = cell2mat(tf_times');
        tf_freqs = cell2mat(tf_freqs');
    catch list_issue
        error('the resolution of time/frequency lists doesn''t match between subjects \n%s',list_issue.message)
    end
    LIMO.data.tf_times = mean(tf_times,1);
    LIMO.data.tf_freqs = mean(tf_freqs,1);
    LIMO.data.start    = LIMO.data.tf_times(1);
    LIMO.data.lowf     = LIMO.data.tf_freqs(1);
    LIMO.data.end      = LIMO.data.tf_times(end);
    LIMO.data.highf    = LIMO.data.tf_freqs(end);
end
cd(LIMO.dir)
end