%% Hierarchichal linear modeling (HLM) with weighted least square (WLS)
% optimizaiton and threshold-free cluster enhancement (TFCE) prior to
% cluster-correction for multiple comparison
%
% Requirements:
%   - Does not support ICA mode.
%   - Channel loations must be already imported in each files.
%   - Bad channels must be already interpolated.
%
% Cedric Cannard, Sep 2022

clear; close all; clc;
dataDir = 'C:\Users\IONSLAB\Desktop\channeling_matlab\data\data_processed2';
outDir = 'G:\My Drive\HLM\limo';
addpath('G:\My Drive\HLM')
cd(dataDir)
tmpSub = dir;
idx = contains({tmpSub.name},'sub');
tmpSub = {tmpSub(idx).name}';
cd(outDir)
eeglab; close;
pop_editoptions('option_single',0); % ensure double precision

% add paths
root = fileparts(which('limo_eeg'));
addpath([root filesep 'limo_cluster_functions']);
addpath([root filesep 'external' filesep 'psom']);
addpath([root filesep 'external']);
addpath([root filesep 'help']);
ftPath = fileparts(which('ft_analysispipeline.m'));
addpath(fullfile(ftPath, 'external','signal'))
addpath(fullfile(ftPath, 'external','stats'))
addpath(fullfile(ftPath, 'external','images'))

% % Create study (PSD precomputed)
% commands = {};
% iFile = 1;
% for iSub = 1:4
%     disp('--------------------------------------------')
%     disp(['Processing subject ' num2str(iSub) ])
%     disp('--------------------------------------------')
%
%     EEG = pop_loadset('filename', sprintf('sub-%2.2d.set',iSub),'filepath',fullfile(dataDir,sprintf('sub-%2.2d',iSub)));
%
%     EEG.saved = 'no';
%     newpath = fullfile(outDir,EEG.subject); mkdir(newpath)
%     newname = sprintf('sub-%2.2d.set',iSub);
%     pop_saveset(EEG,'filepath',newpath,'filename',newname);
%     commands = [ commands(:)' 'index' iSub 'load' fullfile(newpath, newname) ];
%     [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','test','commands', commands, ...
%         'updatedat','on','savedat','off','rmclust','off');
%     [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
%
% end
% % [STUDY,EEG] = pop_savestudy(STUDY,EEG,'filename','test.study','filepath',outDir);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% STUDY = std_makedesign(STUDY,ALLEEG,1,'name','STUDY.design 1','delfiles','off', 'defaultdesign','off', ...
%     'variable1','type','values1',{'rest','trance'},'vartype1','categorical', ...
%     'subjselect',{'sub-01','sub-02','sub-03','sub-04'});
% [STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');
% [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on', ...
%     'rmicacomps','off', 'interp','off','recompute','on','spec','on', ...
%     'specparams',{'specmode','psd','logtrials','on','freqrange',[1 13]});
% [STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');
% cd(outDir)
% gong

% Load study
[STUDY, ALLEEG] = pop_loadstudy('filename','test.study','filepath',outDir);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

%% Channel neighbors

% find neighbors
neighbors = get_channelneighbors(ALLEEG(1).chanlocs,0);

% Generate LIMO chanlocs file
STUDY.etc.statistics.fieldtrip.channelneighbor = neighbors;
limoChanlocs.expected_chanlocs = ALLEEG(1).chanlocs;
limoChanlocs.channeighbstructmat = zeros(ALLEEG(1).nbchan);
for i = 1:length(neighbors)
    [tmp, posChan] = intersect({ALLEEG(1).chanlocs.labels}, neighbors(i).neighblabel);
    limoChanlocs.channeighbstructmat(i,posChan) = 1;
    limoChanlocs.channeighbstructmat(posChan,i) = 1;
end
% limoChanlocs = limoChanlocs;
mkdir([STUDY.filepath filesep 'derivatives']);
limoChanlocsFile = fullfile([STUDY.filepath filesep 'derivatives'], 'limo_gp_level_chanlocs.mat');
save('-mat', limoChanlocsFile, '-struct', 'limoChanlocs');
fprintf('Saving channel neighbors for correction for multiple comparisons in \n%s\n', limoChanlocsFile);

%% 1st level: Fixed effects

% HLM Parameters
studypath = STUDY.filepath;
cd(studypath);
opt.method = 'WLS';                 % 'OLS', 'IRLS', 'WLS' (default)
opt.measure = 'datspec';            % 'daterp', 'datspec', 'dattimef', 'icaerp', 'icaspec', 'icatimef'
opt.design = STUDY.currentdesign;   % default = 1
opt.erase = 'on';                   % default = 1
opt.splitreg = 'off';               % default = off
opt.interaction = 'off';             % default = off
opt.freqlim = [1 15];               % default = [], [1 30]
opt.timelim = [];                   % default = [], [-150 1000]
opt.chanloc = struct('no', {});     % default = struct('no', {})
% opt.neighboropt = {};               % default = {}
% opt.neighbormat = [];               % default = []
opt.zscore = 1;                     % 0 or 1 (default)
opt.ow_chanlocfile = 'no';          % 'yes' or 'no' (default)
% if strcmpi(opt.measure, 'datersp')
%     opt.measure = 'dattimef';
% end
opt.measureori = opt.measure;
Analysis = opt.measure;
design_index = opt.design;
model.defaults.datatype = Analysis(4:end);
model.defaults.type = 'Channels';

% LIMO 1st level
% pop_limo(STUDY,ALLEEG, ...
%     'method','WLS', ...         % GLM optimization: OLS, IRLS, WLS (default)
%     'measure','datspec', ...    % which EEG measure: daterp, datsec, datersp
%     'freqlim',[0 13] , ...      % freqs (e.g., [0 30]); or timelim for ERP
%     'erase','on', ...           % erase previous model
%     'splitreg','off', ...       % split regression for continuous indep. var.
%     'interaction','off');       % interaction model for categorical indep. var.

% model.set_files  = [];
% model.cat_files  = [];
% model.cont_files = [];

% Cleaning old files from the current design (Cleaning ALL)
if strcmp(opt.erase,'on')
    [~,filename] = fileparts(STUDY.filename);
    std_limoerase(STUDY.filepath, filename, STUDY.subject, num2str(STUDY.currentdesign));
    STUDY.limo = [];
end

% Check if the measures has been computed and if the channels are interpolated
% interpolated = zeros(1,length(STUDY.datasetinfo));
for iDat = 1:length(STUDY.datasetinfo)
    fileName = fullfile(STUDY.datasetinfo(iDat).filepath, [ STUDY.datasetinfo(iDat).subject '*.' opt.measure ]);
    % fileName should already match unless user moves / rename, hence using dir
    fileName = dir(fileName);
    if isempty(fileName)
        error('std_limo subject %s: Measures must be computed first',STUDY.datasetinfo(iDat).subject);
        %     else
        %         if strcmp(model.defaults.type,'Channels')
        %             tmpChans = load('-mat', fullfile(fileName(1).folder,fileName(1).name), 'labels');
        %             if length(tmpChans.labels) > ALLEEG(iDat).nbchan, interpolated(iDat) = 1; end
        %         end
    end
end

measureflags = struct('daterp','off',...
    'datspec','off',...
    'datersp','off',...
    'dattimef','off',...
    'datitc' ,'off',...
    'icaerp' ,'off',...
    'icaspec','off',...
    'icatimef','off',...
    'icaersp','off',...
    'icaitc','off');
measureflags.(lower(opt.measure)) = 'on';
STUDY.etc.measureflags = measureflags;
STUDY.etc.datspec = 'on';
mergedChanlocs = ALLEEG(1).chanlocs;

% generate temporary merged datasets needed by LIMO
fprintf('generating temporary files, pulling relevant trials ... \n')
allSubjects    = { STUDY.datasetinfo.subject };
allSessions    = { STUDY.datasetinfo.session };
uniqueSubjects = unique(allSubjects);
nb_subjects    = length(uniqueSubjects);
allSessions(cellfun(@isempty, allSessions)) = { 1 };
allSessions    = cellfun(@num2str, allSessions, 'uniformoutput', false);
uniqueSessions = unique(allSessions);

factors = pop_listfactors(STUDY.design(opt.design), 'gui', 'off', 'level', 'one');

% Create statistical models
for iSubj = 1:nb_subjects
    for iSess = 1:length(uniqueSessions)
        inds1 = strmatch( uniqueSubjects{iSubj}, allSubjects, 'exact');
        inds2 = strmatch( uniqueSessions{iSess}, allSessions, 'exact');
        inds  = intersect(inds1, inds2);

        if ~isempty(inds)
            if length(inds) ~= 1
                error([ 'Cannot calculate contrast because more than 1 dataset per session'  ...
                    'per subject. Merge datasets for each subject and try again.' ]);
            end

            % check subject name and session
            [~,subname] = fileparts(STUDY.datasetinfo(inds).filename);
            if isfield(ALLEEG,'filename')
                if ~strcmp(subname,ALLEEG(inds).filename(1:end-4))
                    error('STUDY and ALLEEG mismatch, can''t figure out which file to use')
                end
            else
                warning('No filename in ALLEEG, pulling data blindly from STUDY')
            end
            if strcmp(subname(1:4),'sub-')
                if contains(subname,'ses-')
                    filename = [subname '_design' num2str(design_index) '.set'];
                else
                    filename = [subname '_ses-' num2str(iSess) '_design' num2str(design_index) '.set'];
                end
            else
                if contains(subname,'ses-')
                    filename = ['sub-' subname '_design' num2str(design_index) '.set'];
                else
                    filename = ['sub-' subname '_ses-' num2str(iSess) '_design' num2str(design_index) '.set'];
                end
            end

            % Creating fields for limo
            fprintf('pulling trials for %s ... \n',filename)
            EEGTMP = std_lm_seteegfields(STUDY,ALLEEG(inds), inds,'datatype',model.defaults.type,'format', 'cell');
            ALLEEG = eeg_store(ALLEEG, EEGTMP, inds);

            % Return full path if 'filepath' is a relative path. The output format will
            % fit the one of 'filepath'. That means that if 'filepath' is a cell array,
            % then the output will a cell array too, and the same if is a string.
            filepath = ALLEEG(inds).filepath;
            %             nit = 1;
            %             if iscell(filepath)
            %                 nit = length(filepath);
            %             end
            %             for i = 1:nit
            %                 if iscell(filepath)
            %                     pathtmp = filepath{i};
            %                 else
            %                     pathtmp = filepath;
            %                 end

            %                 if strfind(pathtmp(end),filesep) %#ok<STRIFCND>
            %                     pathtmp = pathtmp(1:end-1);
            %                 end
            %                 if ~isempty(strfind(pathtmp(1:2),['.' filesep])) || ...
            %                         (isunix && pathtmp(1) ~= '/') || (ispc && pathtmp(2) ~= ':')
            %                     if iscell(filepath)
            %                         file_fullpath{i} = fullfile(studypath,pathtmp(1:end));
            %                     else
            %                         file_fullpath = fullfile(studypath,pathtmp(1:end));
            %                     end
            %                 else
            %                     if iscell(filepath)
            %                         file_fullpath{i} = pathtmp;
            %                     else
            %                         file_fullpath = pathtmp;
            %                     end
            %                 end
            %             end

            model.set_files{inds} = fullfile(filepath, filename);

            OUTEEG = [];
            if all([ALLEEG(inds).trials] == 1)
                OUTEEG.trials = 1;
            else
                OUTEEG.trials = sum([ALLEEG(inds).trials]);
            end

            %             filepath_tmp = rel2fullpath(STUDY.filepath,ALLEEG(inds).filepath);
            filepath_tmp        = filepath;
            OUTEEG.filepath     = filepath_tmp;
            OUTEEG.filename     = filename;
            OUTEEG.srate        = ALLEEG(inds).srate;
            OUTEEG.icaweights   = ALLEEG(inds).icaweights;
            OUTEEG.icasphere    = ALLEEG(inds).icasphere;
            OUTEEG.icawinv      = ALLEEG(inds).icawinv;
            OUTEEG.icachansind  = ALLEEG(inds).icachansind;
            OUTEEG.etc          = ALLEEG(inds).etc;
            OUTEEG.times        = ALLEEG(inds).times;
            %             if any(interpolated)
            %                 OUTEEG.chanlocs        = mergedChanlocs;
            %                 OUTEEG.etc.interpolatedchannels = setdiff(1:length(OUTEEG.chanlocs), std_chaninds(OUTEEG, { ALLEEG(inds).chanlocs.labels }));
            %             else
            OUTEEG.chanlocs = ALLEEG(inds).chanlocs;
            %             end
            OUTEEG.etc.merged{1}= ALLEEG(inds).filename;

            % Def fields
            %             OUTEEG.etc.datafiles.daterp   = [];
            %             OUTEEG.etc.datafiles.datspec  = [];
            %             OUTEEG.etc.datafiles.datersp  = [];
            %             OUTEEG.etc.datafiles.dattimef = [];
            %             OUTEEG.etc.datafiles.datitc   = [];
            %             OUTEEG.etc.datafiles.icaerp   = [];
            %             OUTEEG.etc.datafiles.icaspec  = [];
            %             OUTEEG.etc.datafiles.icaersp  = [];
            %             OUTEEG.etc.datafiles.icatimef = [];
            %             OUTEEG.etc.datafiles.icaitc   = [];

            % Filling fields
            %             single_trials_filename = ;
            %             if exist(single_trials_filename,'file')
            %                 if strcmpi(measureflags.daterp,'on')
            %                     OUTEEG.etc.datafiles.daterp = single_trials_filename;
            %                 elseif strcmpi(measureflags.datspec,'on')
            OUTEEG.etc.datafiles.datspec = EEGTMP.etc.datafiles.(opt.measureori);
            %                 elseif strcmpi(measureflags.datersp,'on')
            %                     OUTEEG.etc.datafiles.datersp = single_trials_filename;
            %                 elseif strcmpi(measureflags.datitc,'on')
            %                     OUTEEG.etc.datafiles.datitc = single_trials_filename;
            %                 elseif strcmpi(measureflags.dattimef,'on')
            %                     OUTEEG.etc.datafiles.dattimef = single_trials_filename;
            %                 elseif strcmpi(measureflags.icaerp,'on')
            %                     OUTEEG.etc.datafiles.icaerp = single_trials_filename;
            %                 elseif strcmpi(measureflags.icaspec,'on')
            %                     OUTEEG.etc.datafiles.icaspec = single_trials_filename;
            %                 elseif strcmpi(measureflags.icaersp,'on')
            %                     OUTEEG.etc.datafiles.icaersp = single_trials_filename;
            %                 elseif strcmpi(measureflags.icaitc,'on')
            %                     OUTEEG.etc.datafiles.icaitc = single_trials_filename;
            %                 elseif strcmpi(measureflags.icatimef,'on')
            %                     OUTEEG.etc.datafiles.icatimef = single_trials_filename;
            %                 end
            %             end

            % Save info
            EEG = OUTEEG;
            save('-mat', fullfile(filepath_tmp, OUTEEG.filename), 'EEG');
            clear OUTEEG filepath_tmp


            % Statistical model
            fprintf('making up statistical model for %s ... \n',filename)

            % create matices for continuous and categorical variables
            trialinfo = std_combtrialinfo(STUDY.datasetinfo, inds);
            [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, ...
                'splitreg','off','interaction',opt.interaction);
            if strcmpi(opt.splitreg,'on') && ~isempty(contMat)
                for c=size(contMat,2):-1:1
                    splitreg{c} = limo_split_continuous(catMat,contMat(:,c)); % std_limodesign does something else when splitting regressors
                end
                contMat = cell2mat(splitreg);
                opt.zscore = 0; % regressors are now zscored
            end

            % copy results
            model.cat_files{inds} = catMat;
            model.cont_files{inds} = contMat;
            if isfield(limodesign, 'categorical')
                STUDY.limo.categorical = limodesign.categorical;
            else
                STUDY.limo.categorical = {};
            end
            if isfield(limodesign, 'continuous')
                STUDY.limo.continuous = limodesign.continuous;
            else
                STUDY.limo.continuous = {};
            end
            STUDY.limo.subjects(inds).subject     = STUDY.datasetinfo(inds(1)).subject;
            STUDY.limo.subjects(inds).cat_file    = catMat;
            STUDY.limo.subjects(inds).cont_file   = contMat;
        end
    end
end

% Add contrasts for conditions that were merged during design selection
% i.e. multiple categorical variables (factors) and yet not matching the number
% of variables (contrasts are then a weighted sum of the crossed factors)
if ~isempty(factors) && isfield(factors, 'value') && ...
        sum(arrayfun(@(x) ~strcmpi(x.label,'group'), STUDY.design(opt.design).variable)) == 1 % only one non-continuous variable other than group
    if length(STUDY.design(opt.design).variable(1).value) ~= length(factors) % and this var has more values than the number of factors
        limocontrast = zeros(length(STUDY.design(opt.design).variable(1).value),length(factors)+1); % length(factors)+1 to add the constant
        for n=length(factors):-1:1
            factor_names{n} = factors(n).value;
        end

        index = find(arrayfun(@(x) ~strcmpi(x.label,'group'),STUDY.design(opt.design).variable)); % which one is not group
        for c=1:length(STUDY.design(opt.design).variable(index).value)
            limocontrast(c,1:length(factors)) = single(ismember(factor_names,STUDY.design(opt.design).variable(index).value{c}));
            limocontrast(c,1:length(factors)) = limocontrast(c,1:length(factors)) ./ sum(limocontrast(c,1:length(factors))); % scale by the number of variables
        end
    end
end

% transpose
% model.set_files  = model.set_files';
% model.cat_files  = model.cat_files';
% model.cont_files = model.cont_files';
if all(cellfun(@isempty, model.cat_files )), model.cat_files  = []; end
if all(cellfun(@isempty, model.cont_files)), model.cont_files = []; end

% set model.defaults - all conditions no bootstrap
% if strcmp(Analysis,'daterp') || strcmp(Analysis,'icaerp')
%     model.defaults.analysis = 'Time';
%     for s=nb_subjects:-1:1
%         vs(s) = ALLEEG(s).xmin*1000;
%         ve(s) = ALLEEG(s).xmax*1000;
%     end
%     model.defaults.start    = max(vs);
%     model.defaults.end      = min(ve);
%
%     if length(opt.timelim) == 2 && opt.timelim(1) < opt.timelim(end)
%         % start value
%         if opt.timelim(1) < model.defaults.start
%             fprintf('std_limo: Invalid time lower limit, using default value instead');
%         else
%             model.defaults.start = opt.timelim(1);
%         end
%         % end value
%         if opt.timelim(end) > model.defaults.end
%             fprintf('std_limo: Invalid time upper limit, using default value instead');
%         else
%             model.defaults.end = opt.timelim(end);
%         end
%     end
%     model.defaults.lowf  = [];
%     model.defaults.highf = [];
%
% elseif strcmp(Analysis,'datspec') || strcmp(Analysis,'icaspec')

model.defaults.analysis = 'Frequency';
if length(opt.freqlim) == 2
    model.defaults.lowf    = opt.freqlim(1);
    model.defaults.highf   = opt.freqlim(2);
else
    error('std_limo: Frequency limits need to be specified');
end
model.defaults.start = [];
model.defaults.end   = [];

% elseif strcmp(Analysis,'dattimef') || any(strcmp(Analysis,{'icatimef','icaersp'}))
%
%     model.defaults.analysis = 'Time-Frequency';
%     for s=nb_subjects:-1:1
%         vs(s) = ALLEEG(s).xmin*1000;
%         ve(s) = ALLEEG(s).xmax*1000;
%     end
%     model.defaults.start    = max(vs);
%     model.defaults.end      = min(ve);
%     model.defaults.lowf     = [];
%     model.defaults.highf    = [];
%
%     if length(opt.timelim) == 2
%         model.defaults.start    = opt.timelim(1);
%         model.defaults.end      = opt.timelim(2);
%     end
%     if length(opt.freqlim) == 2
%         model.defaults.lowf     = opt.freqlim(1);
%         model.defaults.highf    = opt.freqlim(2);
%     else
%         error('std_limo: Frequency limits need to be specified');
%     end
% end

model.defaults.fullfactorial    = 0;                 % all variables
model.defaults.zscore           = opt.zscore;        % done that already
model.defaults.bootstrap        = 0 ;                % only for single subject analyses - not included for studies
model.defaults.tfce             = 0;                 % only for single subject analyses - not included for studies
model.defaults.method           = opt.method;        % default is OLS - to be updated to 'WLS' once validated
model.defaults.Level            = 1;                 % 1st level analysis
model.defaults.type_of_analysis = 'Mass-univariate'; % option can be multivariate

% if ~exist('limocontrast','var')
[LIMO_files, procstatus] = limo_batch('model specification',model,[],STUDY);
% else
%     contrast.mat = limocontrast;
%     [LIMO_files, procstatus] = limo_batch('both',model,contrast,STUDY);
%     if exist(fullfile([STUDY.filepath filesep 'derivatives']),'dir')
%         save(fullfile([STUDY.filepath filesep 'derivatives'],[STUDY.design(opt.design).name '_contrast.mat']),'limocontrast');
%     else
%         save(fullfile(STUDY.filepath,[STUDY.design(opt.design).name '_contrast.mat']),'limocontrast');
%     end
% end

STUDY.limo.model = model;
STUDY.limo.datatype = Analysis;
STUDY.limo.chanloc = limoChanlocs;
% if exist('limocontrast','var')
%     STUDY.limo.contrast      = limocontrast;
% end

% Between session contrasts
index = 1;
for s = 1:nb_subjects
    sess_index = find(cellfun(@(x) strcmpi(x,uniqueSubjects{s}), allSubjects));
    % matches sess_index = find(contains(LIMO_files.mat,uniqueSubjects{s}))
    if length(sess_index) > 1
        fprintf('std_limo, computing additional between sessions contrasts for subject %s\n',uniqueSubjects{s})
        sess_name = allSessions(sess_index);
        pairs = nchoosek(1:length(sess_index),2); % do all session pairs
        parfor p=1:size(pairs,1)
            strpair = [cell2mat(sess_name(pairs(p,1))) cell2mat(sess_name(pairs(p,2)))];
            strpair(isspace(strpair)) = []; % remove spaces
            filesout{p} = limo_contrast_sessions(cell2mat(LIMO_files.mat(sess_index(pairs(p,1)))), ...
                cell2mat(LIMO_files.mat(sess_index(pairs(p,2)))),strpair);
        end

        for f=1:length(filesout)
            for ff=1:length(filesout{f})
                allcon{index} = filesout{f}{ff};
                index = index +1;
            end
        end
        clear filesout
    end
end

% use same glm_name as limo_batch
design_name = STUDY.design(STUDY.currentdesign).name;
design_name(isspace(design_name)) = [];
if contains(design_name,'STUDY.')
    design_name = design_name(7:end);
end
glm_name = [STUDY.filename(1:end-6) '_' design_name '_GLM_' model.defaults.type '_' model.defaults.analysis '_' model.defaults.method];

% Further split that list per regressor and group
% if exist('allcon','var')
%     maxcon = max(cellfun(@(x) str2double(x(strfind(x,'con_')+4:strfind(x,'sess_')-1)),allcon));
%     for con=1:maxcon
%         index = find(cellfun(@(x) ~isempty(x),cellfun(@(x) strfind(x,['con_' num2str(con)]),allcon','UniformOutput',false)));
%         cell2csv([LIMO_files.LIMO filesep 'Between_sessions_con_' num2str(con) '_' glm_name '.txt'], allcon(index)');
%         if length(STUDY.group) > 1
%             for g= 1:length(STUDY.group)
%                 % find subjects of group g
%                 subset = find(arrayfun(@(x)(strcmpi(x.group,STUDY.group{g})), STUDY.datasetinfo));
%                 for s=1:length(subset)
%                     sub{s} = STUDY.datasetinfo(subset(s)).subject;
%                 end
%                 sub = unique(sub);
%                 % find subjects of group g and contrast con
%                 subcon = [];
%                 for s = 1:length(sub)
%                     if strfind(sub{s},'sub-')
%                         subindex = find(cellfun(@(x) ~isempty(x),(cellfun(@(x) strfind(x,sub{s}),allcon','UniformOutput',false)))); % subject s group g
%                     else
%                         subindex = find(cellfun(@(x) ~isempty(x),(cellfun(@(x) strfind(x,['sub-' sub{s} ]),allcon','UniformOutput',false)))); % subject s group g
%                     end
%                     subcon = [subcon;intersect(index,subindex)];
%                 end
%                 % save
%                 if ~isempty(subcon)
%                     cell2csv([LIMO_files.LIMO filesep 'Between_sessions_con_' num2str(con) 'Gp' STUDY.group{g} '_' glm_name '.txt'], allcon(subcon)');
%                 end
%             end
%         end
%     end
% end

% Save STUDY
cd(STUDY.filepath);
STUDY = pop_savestudy(STUDY, [],'filepath', STUDY.filepath,'savemode','resave');
keep_files = 'no';
if all(procstatus)
    disp('All subjects have been successfully processed.')
else
    if sum(procstatus)==0 % not a WLS issue - limo_batch errors for that and tells the user
        errordlg2('all subjects failed to process, check limo batch report')
    else
        warndlg2('some subjects failed to process, check limo batch report','', 'non-modal')
    end
    % cleanup temp files - except for subjects without errors
    db = dbstack;
    if length(db) <= 2
        keep_files = questdlg('Do you want to keep temp files of unsuccessulfully processed subjects','option for manual debugging','yes','no','no');
    end

    % delete
    if isempty(keep_files) || strcmpi(keep_files,'no')
        for s = 1:nb_subjects
            delete(model.set_files{s});
        end
    else
        for s = find(procstatus)
            delete(model.set_files{s});
        end
    end
end

%% 2nd level: Random effects

% outDir = fullfile(dataDir, 'derivatives', 'LIMO_test', 'results');
% mkdir(outDir); cd(outDir)

%%%%%%%%%%%%%%%%%%%%% limo_random_select %%%%%%%%%%%%%%%%%%%%%%%%%
stattest = 'paired t-test';  %'one sample t-test', 'two-samples t-test', 'paired t-test', ...
%     'regression', 'N-Ways ANOVA', 'ANCOVA', 'Repeated Measures ANOVA'
expected_chanlocs = fullfile(studypath,'derivatives','limo_gp_level_chanlocs.mat');
% limo_random_select(stattest,expected_chanlocs,options)

LIMO.data = load(expected_chanlocs);
% if isfield(LIMO.data,'expected_chanlocs')
LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
% end
% if isfield(LIMO.data,'channeighbstructmat')
LIMO.data = renameStructField(LIMO.data, 'channeighbstructmat', 'neighbouring_matrix');
% end

% build LIMO.mat from inputs
% options = {'nboot' 1000 'tfce' 1 'type' };
LIMO.dir = 'G:\My Drive\HLM\limo\derivatives\LIMO_test\results';
mkdir(LIMO.dir); cd(LIMO.dir)
LIMO.Level = 2;
% LIMO.Analysis          = [];
% LIMO.Type              = [];
% LIMO.data.data         = [];
% LIMO.design.bootstrap  = 0;
% LIMO.design.tfce       = 0;
% LIMO.design.electrode  = [];
% LIMO.design.component  = [];
% LIMO.design.parameters = [];
LIMO.design.method = 'Trimmed mean';
% regressor_file = [];
analysis_type = 'Full scalp analysis'; % '1 channel/component only'
% zopt = [];
% skip_design_check = 'No';
LIMO.Type = 'Channels';
% idx = find(strcmpi(options,'nboot'));
LIMO.design.bootstrap = 1000;
% idx = find(strcmpi(options,'tfce'));
LIMO.design.tfce = 1;

warning on

% for in = 1:2:length(options)
%     if strcmpi(options{in},'LIMOfiles')
%         if ~iscell(options{in+1})
%             LIMO.data.data = {options{in+1}}; %#ok<CCAT1>
%         else
%             LIMO.data.data = options{in+1};
%         end
%     elseif strcmpi(options{in},'analysis type')
%     if strcmpi(options{in},'analysis type')
%         if any(strcmpi(options{in+1},{'Full scalp analysis','1 channel/component only'}))
%             analysis_type = options{in+1};
%         else
%             error('analysis type argument unrecognized')
%         end
%     elseif contains(options{in},'regressor','IgnoreCase',true)
%         regressor_file = options{in+1};
%     elseif strcmpi(options{in},'zscore')
%         zopt = options{in+1};
%     elseif strcmpi(options{in},'skip design check') || strcmpi(options{in},'skip_design_check')
%         skip_design_check = options{in+1};
%     elseif strcmpi(options{in},'channel')
% idx = strcmpi(options,'channel');
% LIMO.design.electrode = options{idx};

% idx = contains(options{:},'parameter','IgnoreCase',true)
%         if iscell(options{idx})
%             LIMO.design.parameters = options{in+1};
%         else
%             LIMO.design.parameters = {options{in+1}};
%         end
%     elseif strcmpi(options{in},'factor names')
%         LIMO.design.factor_names = options{in+1};
% idx = find(strcmpi(options,'type'));
% LIMO.Type = options{idx+1};
% idx = find(strcmpi(options,'nboot'));
% LIMO.design.bootstrap = options{idx+1};
% idx = find(strcmpi(options,'tfce'));
% LIMO.design.tfce = options{idx+1};


% Paired t-test
if strcmpi(stattest,'paired t-test')

    %     LIMO.design.X = [];
    %     [Names{1},Paths{1},LIMO.data.data{1}] = limo_get_files;
    %     Names = {};
    %     Paths = {};
    %     LIMO.data.data = {};
    %     LIMO.data = rmfield(LIMO.data, 'data');
    for i = 1:length(ALLEEG)
        Names(:,i) = {'Betas.mat'};
        tmppath = fullfile(studypath, 'derivatives', ALLEEG(i).subject);
        cd(tmppath); tmppath2 = dir; tmppath2 = tmppath2(contains({tmppath2.name},'GLM')).name;
        Paths(:,i) = {fullfile(tmppath, tmppath2)};
        %         LIMO.data.data(:,i) = char(fullfile(Paths(i), Names(i)));
        LIMO.data.data(:,i) = {fullfile(Paths(i), Names(i))};
    end
    Names = { Names };
    Paths = { Paths };
    LIMO.data.data = { LIMO.data.data };
    LIMO.data.data_dir = Paths{1};

    % Beta parameters to test
    parameters = [1 2]; % 1:3 if 3 conditions
    %     if size(parameters,2) == 1 % either con file, or command line beta with one parameter
    %
    %         % do 2nd pair of data
    %         if size(LIMO.data.data,2) == 1 % only 1st group of files
    %             [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],'select paired file');
    %             LIMO.data.data_dir{2} = Paths{2};
    %         else
    %             if list == 1 && ischar(LIMO.data.data{2})
    %                 [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],[],LIMO.data.data{2});
    %                 LIMO.data.data_dir{2} = Paths{2};
    %             else % Case when all paths are provided
    %                 [Names{2},Paths{2}] = breaklimofiles(LIMO.data.data(:,2));
    %                 LIMO.data = rmfield(LIMO.data,'data');
    %                 for sub = length(Paths{1}):-1:1; LIMO.data.data_dir{1}{sub} = Paths{1}{sub}; end
    %                 for sub = length(Paths{1}):-1:1; LIMO.data.data{1}{sub} = fullfile(Paths{1}{sub},Names{1}{sub}); end
    %                 for sub = length(Paths{2}):-1:1; LIMO.data.data_dir{2}{sub} = Paths{2}{sub}; end
    %                 for sub = length(Paths{2}):-1:1; LIMO.data.data{2}{sub} = fullfile(Paths{2}{sub},Names{2}{sub}); end
    %             end
    %         end
    %
    %         if ~isempty(LIMO.design.parameters)
    %                 % hack only availbale if beta files and command line argument // not allowed otherwise because it's a paired design
    %                 parameters(2) = check_files(Names{2},1,parameters(2));
    %             else
    %                 newparameters = check_files(Names{2},1);
    %                 if newparameters ~= 1
    %                     error('paired t-test second set must also be con files')
    %                 else
    %                     parameters = 1;
    %                 end
    %         end
    %
    %         con_parameters    = [str2double(unique(cellfun(@(x) x(5:end-4),Names{1}))) ...
    %             str2double(unique(cellfun(@(x) x(5:end-4),Names{2})))];
    %         if all(isnan(con_parameters))
    %             clear con_parameters % was betas from command line
    %         end
    %
    %         if size(Names{1},2) ~= size(Names{2},2)
    %             if exist('errordlg2','file')
    %                 errordlg2('the nb of files differs between pairs 1 and 2','Paired t-test error'); return
    %             else
    %                 errordlg('the nb of files differs between pairs 1 and 2','Paired t-test error'); return
    %             end
    %         end
    %     elseif size(parameters,2) ~= 2 % if it was beta file one needs a pair of parameters
    %         if exist('errordlg2','file')
    %             errordlg2('2 parameters must be selected for beta files','Paired t-test error'); return
    %         else
    %             errordlg('2 parameters must be selected for beta files','Paired t-test error'); return
    %         end
    %     else % check Betas match design
    for s = 1:length(Paths{1})
        sub_LIMO = load(fullfile(Paths{1}{s},'LIMO.mat'));
        if max(parameters) > size(sub_LIMO.LIMO.design.X,2)
            errordlg('invalid parameter(s)','Paired t-test error'); return
        end
    end
    cd(LIMO.dir);
    % end

    % Match frames and channels and update LIMO
    [first_frame,last_frame,subj_chanlocs,~,LIMO] = match_frames(Paths,LIMO);
    LIMO = match_channels(3,analysis_type,LIMO);

    % Get data
    if any(size(LIMO.data.data) == 2)
        data = getdata(2,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);
    else
        data = getdata(1,analysis_type,first_frame,last_frame,subj_chanlocs,LIMO);
    end

    % compute
    %     if strcmpi(LIMO.Analysis,'Time-Frequency')
    %         if strcmpi(analysis_type,'1 channel/component only')
    %             if size(parameters,2) == 2 % beta files
    %                 tmp                = squeeze(data(:,:,:,parameters(1),:));
    %                 tmp_data1          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 channel
    %                 tmp_data1(1,:,:,:) = tmp; clear tmp
    %                 tmp                = squeeze(data(:,:,:,parameters(2),:));
    %                 tmp_data2          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
    %                 tmp_data2(1,:,:,:) = tmp; clear tmp
    %             else % con files
    %                 tmp                = squeeze(data{1}(:,:,:,:,:));
    %                 tmp_data1          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3)); % add dim 1 = 1 channel
    %                 tmp_data1(1,:,:,:) = tmp; clear tmp
    %                 tmp                = squeeze(data{2}(:,:,:,:,:));
    %                 tmp_data2          = ones(1,size(tmp,1),size(tmp,2),size(tmp,3));
    %                 tmp_data2(1,:,:,:) = tmp; clear tmp
    %             end
    %         else % full scalp
    %             if size(parameters,2) == 2 % beta files
    %                 if iscell(data)
    %                     tmp_data1 = squeeze(data{1}(:,:,:,parameters(1),:));
    %                     tmp_data2 = squeeze(data{2}(:,:,:,parameters(2),:));
    %                 else
    %                     tmp_data1 = squeeze(data(:,:,:,parameters(1),:));
    %                     tmp_data2 = squeeze(data(:,:,:,parameters(2),:));
    %                 end
    %             else % con files
    %                 tmp_data1 = squeeze(data{1}(:,:,:,:));
    %                 tmp_data2 = squeeze(data{2}(:,:,:,:));
    %             end
    %         end
    %     else
    %         if strcmpi(analysis_type,'1 channel/component only')
    %             if size(parameters,2) == 2 % beta files
    %                 tmp              = squeeze(data(:,:,parameters(1),:));
    %                 tmp_data1        = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 channel
    %                 tmp_data1(1,:,:) = tmp; clear tmp
    %                 tmp              = squeeze(data(:,:,parameters(2),:));
    %                 tmp_data2        = ones(1,size(tmp,1),size(tmp,2));
    %                 tmp_data2(1,:,:) = tmp; clear tmp
    %             else % con files
    %                 tmp              = squeeze(data{1}(:,:,:,:));
    %                 tmp_data1        = ones(1,size(tmp,1),size(tmp,2)); % add dim 1 = 1 channel
    %                 tmp_data1(1,:,:) = tmp; clear tmp
    %                 tmp              = squeeze(data{2}(:,:,:,:));
    %                 tmp_data2        = ones(1,size(tmp,1),size(tmp,2));
    %                 tmp_data2(1,:,:) = tmp; clear tmp
    %             end
    %         else % full scalp
    %     if size(parameters,2) == 2 % beta files
    %         if iscell(data)
    %             data1 = squeeze(data{1}(:,:,parameters(1),:));
    %             data2 = squeeze(data{2}(:,:,parameters(2),:));
    %         else
    data1 = squeeze(data(:,:,parameters(1),:));
    data2 = squeeze(data(:,:,parameters(2),:));
    %         end
    %     else % con files
    %         tmp_data1 = squeeze(data{1}(:,:,:));
    %         tmp_data2 = squeeze(data{2}(:,:,:));
    %     end
    % end
    %     end

    %     if strcmp(LIMO.Analysis,'Time-Frequency')
    %         LIMO.data.size3D = [size(tmp_data1,1) size(tmp_data1,2)*size(tmp_data1,3) 5];
    %         LIMO.data.size4D = [size(tmp_data1,1) size(tmp_data1,2) size(tmp_data1,3) 5];
    %     end
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    Y1r = data1; save Y1r Y1r, clear Y1r
    Y2r = data2; save Y2r Y2r, clear Y2r
    %     if exist('con_parameters','var')
    %         parameters = con_parameters; % substitute param of data by con values for consistent naming
    %     end
    LIMO.design.parameters = parameters; % update in any cases

    %%%%%%%%%%%%%%%%%%% limo_random_robust %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function makes the result files for the random effects of various tests
    % as well as organizes and makes files for boostrap. It is interfaced with
    % limo_random_effect which itself interfaces with the user to select and pass
    % the data in the appropriate format. Limo_random_robust calls low level
    % functions to perform the actual computation and add the trimmed mean or mean
    % of parameters for the correponding test (helps for vizualizing effects)
    %
    % Input:
    %           limo_random_robust(3,y1,y2,parameter number,LIMO)
    %                    3 = paired t-test
    %                    y1 = data (dim channels, time or freq, subjects)
    %                       = data (dim channels, freq, time, subjects)
    %                       = the name of the Y1r file
    %                    y2 = data (dim channels, time or freq, subjects)
    %                       = data (dim channels, freq, time, subjects)
    %                       = the name of the Y2r file
    %                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
    %
    % Output:
    %           3 paired_samples_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
    %               H0_paired_samples_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)
    %
    % Called in limo_random_select with:
    %   tmpname = limo_random_robust(3,tmp_data1,tmp_data2,parameters,LIMO); % --------> calls limo_create_boot_table


    % Paired t-test // percentile bootstrap technique
    %     data1  = tmp_data1;
    %     data2  = tmp_data2;
    %     parameters = parameters;
    cd(LIMO.dir);
    %     if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
    %         data1 = limo_tf_4d_reshape(data1);
    %         data2 = limo_tf_4d_reshape(data2);
    %     end

    % check data structure
    %     if size(data1,1) ~= size(data2,1)
    %         error(['samples have a different number of ' LIMO.Type ', not a paired t-tests'])
    %     end
    %     for e = 1:size(data1,1)
    %         tmp = isnan(data1(e,1,:));
    %         tmp2 = isnan(data2(e,1,:));
    %         if length(tmp) ~= length(isnan(tmp2))
    %             errordlg([LIMO.Type ' ' num2str(e) ' has unpaired data - analysis aborded, not a paired t-test']);
    %             return
    %         elseif length(tmp) == sum(isnan(tmp))
    %             errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
    %             return
    %         elseif (length(tmp) - sum(isnan(tmp))) < 3
    %             errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects - analysis aborded']);
    %             return
    %         end
    %     end
    %     clear tmp tmp2

    %     if ~isfield(LIMO.design,'method')
    %         disp('Setting Trimmed mean as default method');
    %         LIMO.design.method = 'Trimmed Mean';
    %     end

    % Make a paired_samples file per parameter (channels, frames, [mean value, se, df, t, p])
    paired_samples = NaN(size(data1,1), size(data1,2),5);
    name = sprintf('paired_samples_ttest_parameter_%s',num2str(parameters')');

    %     array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
    for iChan = 1:size(data1,1)
        %         channel = array(e);
        fprintf('analyse parameter %s channel %g: %s',num2str(parameters')', iChan, LIMO.design.method);
        disp(' ');
        tmp = data1(iChan,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
        tmp = data2(iChan,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
        if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
            [   paired_samples(iChan,:,4), ...    % t-value
                paired_samples(iChan,:,1), ...    % diff
                paired_samples(iChan,:,2), ...    % se
                ~, ...                              % CI
                paired_samples(iChan,:,5), ...    % p-value
                ~, ...                              % tcrit
                paired_samples(iChan,:,3) ...     % df
                ] = limo_yuend_ttest(Y1,Y2);
        else % Mean
            [   paired_samples(iChan,:,1), ...    % mean
                paired_samples(iChan,:,3), ...    % dfe
                ~, ...                              % CI
                sd, ...
                n, ...
                paired_samples(iChan,:,4),...     % t-value
                paired_samples(iChan,:,5) ...     % p-value
                ] = limo_ttest(1,Y1,Y2,.05);
            paired_samples(iChan,:,2) = sd./sqrt(n);
        end
        clear Y1 Y2
    end

    %     if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
    %         paired_samples = limo_tf_4d_reshape(paired_samples);
    %     end
    save(name,'paired_samples','-v7.3')
    %     LIMOPath = LIMO.dir;

    %%%%%%%%%%%%%%%%% Bootstrap %%%%%%%%%%%%%%%%%%%%%%%%
    if LIMO.design.bootstrap > 0
        limo_check_ppool
        bootex = 1;

        % check if it was already computed
        boot_name = sprintf('H0_paired_samples_ttest_parameter_%s',num2str(parameters')');
        if exist(['H0', filesep, boot_name, '.mat'], 'file') && usejava('desktop')
            answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
            if ~strcmp(answer,'Yes')
                bootex = 0;
            end
        end

        if bootex == 1
            mkdir H0

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% limo_create_boot_table %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Builds a table of data to resample such as almost the same resampling is
            % applied across channels. At the 2nd level, not all subjects have the same
            % channels because missing data can be treated as NaN. The boot_table has
            % indexes which are common to all channels + some other indexes specific
            % to channels where there are some NaNs.
            %
            % For repeated measures (Time-Frequency, Repeated meaures ANOVA, etc)
            % call this function inputing 1 measure and apply the created table
            % to all measures
            %
            % FORMAT: boot_table = limo_create_boot_table(data,nboot)
            %
            % INPUT: data: a 3D matrix [channels x time frames x subjects];
            %        nboot: the number of bootstraps to do
            %
            % OUTPUT: boot_table is a cell array with one resampling matrix per
            %         channel

            disp('making boot table ...')
            % boot_table = limo_create_boot_table(data1,LIMO.design.bootstrap);
            data = data1;
            nboot = LIMO.design.bootstrap;

            % Minimum number of different trials/subjects.
            % lower variance means the stat values will be too high (see Pernet et al. 2014)
            Nmin = 3;
            if size(data,3) <= Nmin
                error(['Not enough subjects in dataset - need at least ' num2str(Nmin) ' subjects']);
            end

            % check data for NaNs
            if size(data,1) == 1, chdata = data(1,1,:); else, chdata = squeeze(data(:,1,:)); end
            if sum(sum(isnan(chdata),2)==size(data,ndims(data))) ~=0
                disp('some channels are empty (full of NaN) - still making the table but some cells will be empty')
            end
            if (sum((size(data,3) - sum(isnan(chdata),2))<=3)) ~=0
                disp('some cells have a very low count <=3 ; bootstrapping cannot work - still making the table but some cells will be empty')
            end

            % create boot_table
            B = 1;
            boot_index = zeros(size(data,3),nboot);
            while B~=nboot+1
                tmp = randi(size(data,3),size(data,3),1);
                if length(unique(tmp)) >= Nmin % at least Nmin different observations per boot
                    boot_index(:,B) = tmp;
                    B=B+1;
                end
            end
            clear chdata tmp

            % loop per channel, if no nan use boot_index else change it
            if size(data,1) > 1
                array = find(sum(squeeze(isnan(data(:,1,:))),2) < size(data,3)-3);
            else
                array = 1;
            end
            for iChan = size(array,1):-1:1
                channel       = array(iChan);
                tmp           = squeeze(data(channel,:,:)); % 2D
                bad_subjects  = find(isnan(tmp(1,:)));
                good_subjects = find(~isnan(tmp(1,:)));
                Y             = tmp(:,good_subjects); % remove NaNs

                if ~isempty(bad_subjects)
                    boot_index2 = zeros(size(Y,2),nboot);
                    for c=1:nboot
                        common  = ismember(boot_index(:,c),good_subjects');
                        current = boot_index(find(common),c); %#ok<FNDSB> % keep resampling of good subjets
                        % add or remove indices
                        add = size(Y,2) - size(current,1);
                        if add > 0
                            new_boot = [current ; good_subjects(randi(size(good_subjects),add,1))'];
                        else
                            new_boot = current(1:size(Y,2));
                        end
                        % change indices values
                        tmp_boot = new_boot;
                        for i=1:length(bad_subjects)
                            new_boot(find(tmp_boot > bad_subjects(i))) = tmp_boot(find(tmp_boot >  bad_subjects(i))) - i; %#ok<FNDSB> % change range
                        end
                        % new boot-index
                        boot_index2(:,c)  = new_boot;
                    end
                    boot_table{channel} = boot_index2;
                else
                    boot_table{channel} = boot_index;
                end
            end
            save(['H0', filesep, 'boot_table'], 'boot_table')

            % get results under H0
            disp('Centering data to estimate H0')
            if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                data1_centered = data1 - repmat(limo_trimmed_mean(data1),[1 1 size(data1,3)]);
                data2_centered = data2 - repmat(limo_trimmed_mean(data2),[1 1 size(data2,3)]);
            else % if strcmpi(LIMO.design.method,'Mean')
                data1_centered = data1 - repmat(mean(data1,3,'omitnan'),[1 1 size(data1,3)]);
                data2_centered = data2 - repmat(mean(data2,3,'omitnan'),[1 1 size(data2,3)]);
            end

            %%%%%%%%%% back to limo_random_robust %%%%%%%%%%%%%%%%%%
            H0_paired_samples = NaN(size(data1,1), size(data1,2), 2, LIMO.design.bootstrap); % stores T and p values for each boot
            for iChan = 1:size(array,1)
                channel = array(iChan);
                fprintf('bootstrapping channel %g/%g parameter %s \n',iChan,size(array,1),num2str(parameters')');
                tmp = data1_centered(channel,:,:);
                Y1 = tmp(1,:,find(~isnan(tmp(1,1,:))));
                clear tmp
                tmp = data2_centered(channel,:,:);
                Y2 = tmp(1,:,find(~isnan(tmp(1,1,:))));
                clear tmp
                %                 t = {}; p = {};
                if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
                    disp('Using Yuen paired t-test (i.e., trimmed means)')
                    parfor b = 1:LIMO.design.bootstrap
                        [t{b},~,~,~,p{b},~,~]=limo_yuend_ttest(Y1(1,:,boot_table{channel}(:,b)),Y2(1,:,boot_table{channel}(:,b)));
                    end
                else % if strcmpi(LIMO.design.method,'Mean')
                    disp('Using paired t-test (i.e., means)')
                    parfor b=1:LIMO.design.bootstrap
                        [~,~,~,~,~,t{b},p{b}]=limo_ttest(1,Y1(1,:,boot_table{channel}(:,b)),Y2(1,:,boot_table{channel}(:,b)));
                    end
                end
                for b=1:LIMO.design.bootstrap
                    H0_paired_samples(channel,:,1,b) = t{b};
                    H0_paired_samples(channel,:,2,b) = p{b};
                end
                clear t p Y1 Y2
            end
            %             if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            %                 H0_paired_samples = limo_tf_5d_reshape(H0_paired_samples);
            %             end
            save (['H0', filesep, boot_name],'H0_paired_samples','-v7.3');
        end
    end % closes if LIMO.design.bootstrap > 0

    %%%%%%%%%%%%%%%%%%%%% limo_tfce_handling %%%%%%%%%%%%%%%%%%%%%%
    % routine to create tfce files commensurate to boostrapped files
    %
    % FORMAT  [tfce_score,thresholded_maps] = limo_tfce_handling(filename,'checkfile','yes')
    %
    % INPUTS filename is the stat file that need to be tfced (if H0 exist it is done too)
    %        'checkfile' is 'yes' by default - if 'no' and tfce files already exist,
    %                     it overwrites without asking otherwise user is prompted
    %
    % OUTPUTS tfce_* files are saved on the drive in a tfce folder
    %         H0_tfce_* files are saved on the drive in the H0 folder
    %         tfce_score is the tfce_score for the data of the file in
    %         thresholded_maps is all the thresholded maps i.e. for each dh
    %                          value in sum(extent(h)^E*height^H*dh)
    %         if a boostrap file exist, thresholded_maps is a 1 * 2 cell array
    %         with {1} the maps for the observed data and {2} a cell of cells
    %         with the maps of each boostrap
    if LIMO.design.tfce ~= 0
        %         limo_tfce_handling(fullfile(LIMO.dir,name));
        file = fullfile(LIMO.dir,name);
        [filepath,filename,ext] = fileparts(file);
        if isempty(ext)
            file = dir(fullfile(filepath,[filename '*']));
            filename = file.name;
        end
        if exist(fullfile(filepath,[filename ext]),'file')
            if ~exist(fullfile(filepath,'LIMO.mat'),'file')
                error('no LIMO.mat found next to %s',filename)
            else
                filename = [filename ext];
                LIMO = load(fullfile(filepath,'LIMO.mat'));
                LIMO = LIMO.(cell2mat(fieldnames(LIMO)));
            end
        else
            error('can''t find %s',varargin{1})
        end
        checkfile = 'yes';

        % files to create
        tfce_file    = fullfile(LIMO.dir,['tfce' filesep 'tfce_' filename]);
        H0_tfce_file = fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' filename]);
        % given filename input, we expect H0 to be
        H0filename   = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);

        %     if strcmpi(checkfile,'yes')
        %         if exist(tfce_file,'file')
        %             answer = questdlg('tfce file already exist - overwrite?','data check','Yes','No','Yes');
        %             if strcmp(answer,'Yes')
        %                 LIMO.design.tfce = 1;
        %                 save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        %             else
        %                 LIMO.design.tfce = 0;
        %                 return
        %             end
        %         else
        %             LIMO.design.tfce = 1;
        %             if exist(LIMO.dir,'dir')
        %                 save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
        %                 if ~exist(fullfile(LIMO.dir,'tfce'),'dir')
        %                     mkdir(fullfile(LIMO.dir,'tfce'));
        %                 end
        %             else
        %                 save(fullfile(pwd,'LIMO.mat'),'LIMO','-v7.3')
        %                 if ~exist(fullfile(pwd,'tfce'),'dir')
        %                     mkdir(fullfile(pwd,'tfce'));
        %                 end
        %             end
        %         end
        %     end
        %     if isfield(LIMO.design,'bootstrap')
        %         nboot = LIMO.design.bootstrap;
        %     end

        % check if there is a neighbouring matrix
        % (since TFCE integrates over clusters)
        %         if ~isfield(LIMO.data,'neighbouring_matrix')
        %             warning('no neighbouring matrix found, this is required for TFCE')
        %             [~, LIMO.data.neighbouring_matrix] = limo_expected_chanlocs;
        %             if isempty(LIMO.data.neighbouring_matrix)
        %                 return
        %             else
        %                 save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        %             end
        %         end

        % create tfce folder
        if ~exist(fullfile(LIMO.dir,'tfce'),'dir')
            mkdir(fullfile(LIMO.dir,'tfce'))
        end
        fprintf('Thresholding %s using TFCE \n',filename);

        %         if contains(filename,'R2') || ...
        %                 contains(filename,'semi_partial') % these files last dimension is R2, F, p
        %         R2 = load(fullfile(LIMO.dir,'R2.mat'));
        %         R2 = R2.(cell2mat(fieldnames(R2)));
        %         if size(R2,1) == 1
        %             if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                 [tfce_score(1,:,:),thresholded_maps] = limo_tfce(2, squeeze(R2(:,:,:,2)),[]); % no neighbouring, time-freq cluster
        %             else
        %                 [tfce_score(1,:),thresholded_maps]   = limo_tfce(1, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
        %             end
        %         else
        %             if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                 [tfce_score,thresholded_maps] = limo_tfce(3, squeeze(R2(:,:,:,2)),LIMO.data.neighbouring_matrix);
        %             else
        %                 [tfce_score,thresholded_maps] = limo_tfce(2, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
        %             end
        %         end
        %         save(tfce_file,'tfce_score','-v7.3'); clear R2 ;
        %
        %         if exist(H0filename,'file')
        %             fprintf('Applying TFCE to null data ... \n')
        %             H0_R2 = load(H0filename);
        %             H0_R2 = H0_R2.(cell2mat(fieldnames(H0_R2)));
        %             tfce_H0_thmaps = cell(1,LIMO.design.bootstrap);
        %             if size(H0_R2,1) == 1
        %                 if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                     tfce_H0_score  = NaN(1,size(H0_R2,2),size(H0_R2,3),LIMO.design.bootstrap);
        %                     parfor b=1:nboot
        %                         [tfce_H0_score(1,:,:,b),tfce_H0_thmaps] = limo_tfce(2,squeeze(H0_R2(:,:,:,2,b)),[],0);
        %                     end
        %                 else
        %                     tfce_H0_score  = NaN(1,size(H0_R2,2),LIMO.design.bootstrap);
        %                     parfor b=1:nboot
        %                         [tfce_H0_score(1,:,b),tfce_H0_thmaps] = limo_tfce(1,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
        %                     end
        %                 end
        %             else
        %                 if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                     tfce_H0_score  = NaN(size(H0_R2,1),size(H0_R2,2),size(H0_R2,3),LIMO.design.bootstrap);
        %                     parfor b=1:nboot
        %                         [tfce_H0_score(:,:,:,b),tfce_H0_thmaps] = limo_tfce(3,squeeze(H0_R2(:,:,:,2,b)),LIMO.data.neighbouring_matrix,0);
        %                     end
        %                 else
        %                     tfce_H0_score  = NaN(size(H0_R2,1),size(H0_R2,2),LIMO.design.bootstrap);
        %                     parfor b=1:nboot
        %                         [tfce_H0_score(:,:,b),tfce_H0_thmaps] = limo_tfce(2,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
        %                     end
        %                 end
        %             end
        %             save(H0_tfce_file,'tfce_H0_score','-v7.3'); clear H0_R2 tfce_H0_score;
        %             tmp                 = thresholded_maps;     clear thresholded_maps;
        %             thresholded_maps{1} = tmp;                  clear tmp
        %             thresholded_maps{2} = tfce_H0_thmaps;       clear tfce_H0_thmaps;
        %         end
        %
        %         elseif contains(filename,'con') || contains(filename,'ess') || ...
        %                 contains(LIMO.design.name,'One sample','IgnoreCase',true) || ...
        %                 contains(LIMO.design.name,'Two samples','IgnoreCase',true) || ...
        %                 contains(LIMO.design.name,'Paired','IgnoreCase',true)  % these file last dimension is mean, se, df, t and p
        tval = load(filename);
        tval = tval.(cell2mat(fieldnames(tval)));
        %         if size(tval,1) == 1
        %             if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                 [tfce_score(1,:,:),thresholded_maps] = limo_tfce(2, squeeze(tval(:,:,:,end-1)),[]); % no neighbouring, time-freq cluster
        %             else
        %                 [tfce_score(1,:),thresholded_maps]   = limo_tfce(1, squeeze(tval(:,:,end-1)),LIMO.data.neighbouring_matrix);
        %             end
        %         else
        %             if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                 [tfce_score,thresholded_maps] = limo_tfce(3, squeeze(tval(:,:,:,end-1)),LIMO.data.neighbouring_matrix);
        %             else
        [tfce_score,thresholded_maps] = limo_tfce(2, squeeze(tval(:,:,4)),LIMO.data.neighbouring_matrix);
        %             end
        %         end
        save(tfce_file,'tfce_score','-v7.3'); clear tval ;

        % if exist(H0filename,'file')
        fprintf('Applying TFCE to null data ... \n')
        H0_tval = load(H0filename);
        H0_tval = H0_tval.(cell2mat(fieldnames(H0_tval)));
        tfce_H0_thmaps = cell(1,LIMO.design.bootstrap);
        %         if size(H0_tval,1) == 1
        %             if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                 tfce_H0_score  = NaN(1,size(H0_tval,2),size(H0_tval,3),LIMO.design.bootstrap);
        %                 parfor b=1:nboot
        %                     [tfce_H0_score(1,:,:,b),tfce_H0_thmaps{b}] = limo_tfce(2,squeeze(H0_tval(:,:,:,end-1,b)),[],0);
        %                 end
        %             else
        %                 tfce_H0_score  = NaN(1,size(H0_tval,2),LIMO.design.bootstrap);
        %                 parfor b=1:nboot
        %                     [tfce_H0_score(1,:,b),tfce_H0_thmaps{b}] = limo_tfce(1,squeeze(H0_tval(:,:,end-1,b)),LIMO.data.neighbouring_matrix,0);
        %                 end
        %             end
        %         else
        %             if strcmpi(LIMO.Analysis,'Time-Frequency')
        %                 tfce_H0_score  = NaN(size(H0_tval,1),size(H0_tval,2),size(H0_tval,3),LIMO.design.bootstrap);
        %                 parfor b=1:nboot
        %                     [tfce_H0_score(:,:,:,b),tfce_H0_thmaps{b}] = limo_tfce(3,squeeze(H0_tval(:,:,:,end-1,b)),LIMO.data.neighbouring_matrix,0);
        %                 end
        %             else
        tfce_H0_score  = NaN(size(H0_tval,1),size(H0_tval,2),LIMO.design.bootstrap);
        parfor b = 1:LIMO.design.bootstrap
            [tfce_H0_score(:,:,b),tfce_H0_thmaps{b}] = limo_tfce(2,squeeze(H0_tval(:,:,end-1,b)),LIMO.data.neighbouring_matrix,0);
        end
        %             end
        %         end
        save(H0_tfce_file,'tfce_H0_score','-v7.3'); clear H0_tval tfce_H0_score;
        tmp                 = thresholded_maps;     clear thresholded_maps;
        thresholded_maps{1} = tmp;                  clear tmp
        thresholded_maps{2} = tfce_H0_thmaps;       clear tfce_H0_thmaps;
    end

    %%%%%%%%%%%% end of limo_tfce_handling %%%%%%%%%%%%%%%

    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
    disp('Paired t-test done')

    %%%%%%%%%%%%%%%%%%%% end of limo_random_robust %%%%%%%%%%%%%%%%%

    %   tmpname = limo_random_robust(3,tmp_data1,tmp_data2,parameters,LIMO); % --------> calls limo_create_boot_table
    %     if nargout ~= 0
    %     LIMOPath = tmpname;
    %     end

    %%%%%%%%%%%%%%%%%%%% end of limo_random_select %%%%%%%%%%%%%%%%%

end

%% PLOT

%%%%%%%%%%%%%%%%%%%%%%%%% limo_display_results GROUP LEVEL ONLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% go to limo_display_results for subject level code
FileName = 'paired_samples_ttest_parameter_12.mat';
c = clock;
p = 0.05;
MCC = 1;  % 1 = uncorrected; 2 = cluster; 3 = tfce; 4 = max;
choice = 'use theoretical p values';
toplot = load(fullfile(LIMO.dir,FileName));
toplot = toplot.(cell2mat(fieldnames(toplot)));
Type = 1;   %  1 - 2D image with intensity as function of time/freq (x) and electrodes (y)
%              2 - scalp topography
%              3 - ERP data (original or modeled)
surfflag = 0; % to allow surfing the figure and click (1) or not (0)

%%%%%%%%%%%%%%%%%%%% limo_stat_values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,choice,[]);
matfile = load(FileName);
% M       = [];
% mask    = [];
% mytitle = [];

% if strcmpi(LIMO.Analysis,'Time-Frequency')
%     if strcmp(FileName,'R2.mat')
%         M         = squeeze(matfile.R2(:,:,:,2)); % F values
%         Pval      = squeeze(matfile.R2(:,:,:,3)); % P values
%         MCC_data  = 'H0_R2.mat';
%         titlename = 'R^2 Coef';
%     elseif strncmp(FileName,'Condition_effect',16)
%         effect_nb = eval(FileName(18:end-4));
%         M         = squeeze(matfile.Condition_effect(:,:,:,1));
%         Pval      = squeeze(matfile.Condition_effect(:,:,:,2));
%         MCC_data  = sprintf('H0_Condition_effect_%g.mat',effect_nb);
%         titlename = sprintf('Condition effect %g  F values',effect_nb);
%     elseif strncmp(FileName,'Covariate_effect',16)
%         effect_nb = eval(FileName(18:end-4));
%         M         = squeeze(matfile.Covariate_effect(:,:,:,1));
%         Pval      = squeeze(matfile.Covariate_effect(:,:,:,2));
%         MCC_data  = sprintf('H0_Covariate_effect_%g.mat',effect_nb);
%         titlename = sprintf('Covariate effect %g  F values',effect_nb);
%     elseif strncmp(FileName,'Interaction_effect',18)
%         effect_nb = eval(FileName(20:end-4));
%         M         = squeeze(matfile.Interaction_effect(:,:,:,1));
%         Pval      = squeeze(matfile.Interaction_effect(:,:,:,2));
%         MCC_data  = sprintf('H0_Interaction_effect_%g.mat',effect_nb);
%         titlename = sprintf('Interaction effect %g  F values',effect_nb);
%     elseif strncmp(FileName,'semi_partial_coef',17)
%         effect_nb = eval(FileName(19:end-4));
%         M         = squeeze(matfile.semi_partial_coef(:,:,:,2));
%         Pval      = squeeze(matfile.semi_partial_coef(:,:,:,3));
%         MCC_data  = sprintf('H0_semi_partial_coef_%g.mat',effect_nb);
%         titlename = sprintf('Semi Partial Coef %g',effect_nb);
%     elseif strncmp(FileName,'con_',4)
%         effect_nb = eval(FileName(5:end-4));
%         M         = squeeze(matfile.con(:,:,:,4));
%         Pval      = squeeze(matfile.con(:,:,:,5));
%         MCC_data  = sprintf('H0_con_%g.mat',effect_nb);
%         titlename = sprintf('Contrast %g T values',effect_nb);
%     elseif contains(FileName,'ttest') || contains(FileName,'LI_Map')
%         matfile   = matfile.(cell2mat(fieldnames(matfile)));
%         M         = matfile(:,:,:,4); % T values
%         Pval      = matfile(:,:,:,5);
%         MCC_data  = sprintf('H0_%s', FileName);
%         name      = FileName(1:strfind(FileName,'ttest')+4);
%         name(strfind(name,'_')) = ' ';
%         titlename = sprintf('%s t-test T values',name);
%     elseif strncmp(FileName,'ess_',4)
%         effect_nb = eval(FileName(5:end-4));
%         M         = squeeze(matfile.ess(:,:,:,end-1));
%         Pval      = squeeze(matfile.ess(:,:,:,end));
%         MCC_data  = sprintf('H0_ess_%g.mat',effect_nb);
%         titlename = sprintf('Contrast %g F values',effect_nb);
%     end
%
% else  % same with one less dimention
%     if strcmp(FileName,'R2.mat')
%         M         = squeeze(matfile.R2(:,:,2)); % F values
%         Pval      = squeeze(matfile.R2(:,:,3)); % P values
%         MCC_data  = 'H0_R2.mat';
%         titlename = 'R^2 Coef';
%     elseif strncmp(FileName,'Condition_effect',16)
%         effect_nb = eval(FileName(18:end-4));
%         M         = squeeze(matfile.Condition_effect(:,:,1));
%         Pval      = squeeze(matfile.Condition_effect(:,:,2));
%         MCC_data  = sprintf('H0_Condition_effect_%g.mat',effect_nb);
%         titlename = sprintf('Condition effect %g  F values',effect_nb);
%     elseif strncmp(FileName,'Covariate_effect',16)
%         effect_nb = eval(FileName(18:end-4));
%         M         = squeeze(matfile.Covariate_effect(:,:,1));
%         Pval      = squeeze(matfile.Covariate_effect(:,:,2));
%         MCC_data  = sprintf('H0_Covariate_effect_%g.mat',effect_nb);
%         titlename = sprintf('Covariate effect %g  F values',effect_nb);
%     elseif strncmp(FileName,'Interaction_effect',18)
%         effect_nb = eval(FileName(20:end-4));
%         M         = squeeze(matfile.Interaction_effect(:,:,1));
%         Pval      = squeeze(matfile.Interaction_effect(:,:,2));
%         MCC_data  = sprintf('H0_Interaction_effect_%g.mat',effect_nb);
%         titlename = sprintf('Interaction effect %g  F values',effect_nb);
%     elseif strncmp(FileName,'semi_partial_coef',17)
%         effect_nb = eval(FileName(19:end-4));
%         M         = squeeze(matfile.semi_partial_coef(:,:,2));
%         Pval      = squeeze(matfile.semi_partial_coef(:,:,3));
%         MCC_data  = sprintf('H0_semi_partial_coef_%g.mat',effect_nb);
%         titlename = sprintf('Semi Partial Coef %g',effect_nb);
%     elseif strncmp(FileName,'con_',4)
%         effect_nb = eval(FileName(5:end-4));
%         M         = squeeze(matfile.con(:,:,4));
%         Pval      = squeeze(matfile.con(:,:,5));
%         MCC_data  = sprintf('H0_con_%g.mat',effect_nb);
%         titlename = sprintf('Contrast %g T values',effect_nb);
%     elseif contains(FileName,'ttest') || contains(FileName,'LI_Map')
matfile   = matfile.(cell2mat(fieldnames(matfile)));
M         = matfile(:,:,4); % T values
Pval      = matfile(:,:,5);
MCC_data  = sprintf('H0_%s', FileName);
name      = FileName(1:strfind(FileName,'ttest')+4);
name(strfind(name,'_')) = ' ';
titlename = sprintf('%s T values',name);
%     elseif strncmp(FileName,'ess_',4)
%         effect_nb = eval(FileName(max(strfind(FileName,'_'))+1:end-4));
%         M         = squeeze(matfile.ess(:,:,end-1));
%         Pval      = squeeze(matfile.ess(:,:,end));
%         MCC_data  = sprintf('H0_ess_%g.mat',effect_nb);
%         titlename = sprintf('Contrast %g F values',effect_nb);
%     end
% end

% No correction for multiple testing
if ~isempty(M) && MCC == 1
    M       = Pval;
    mask    = Pval <= p;
    mytitle = sprintf('%s: uncorrected threshold',titlename);

    % Cluster-correction for multiple testing
elseif ~isempty(M) && MCC == 2
    %     if exist(['H0' filesep MCC_data],'file')
    %         try
    H0_data = load(['H0' filesep MCC_data]);
    H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
    %             if strcmpi(LIMO.Analysis,'Time-Frequency')
    %                 if contains(FileName,'R2') || contains(FileName,'semi_partial')
    %                     bootM = squeeze(H0_data(:,:,:,2,:)); % get all F values under H0
    %                     bootP = squeeze(H0_data(:,:,:,3,:)); % get all P values under H0
    %                 else
    %                     bootM = squeeze(H0_data(:,:,:,1,:));
    %                     bootP = squeeze(H0_data(:,:,:,2,:));
    %                 end
    %
    %                 if size(M,1) == 1
    %                     tmp = NaN(1,size(M,2),size(M,3),size(bootM,3));
    %                     tmp(1,:,:,:) = bootM; bootM = tmp;
    %                     tmp(1,:,:,:) = bootP; bootP = tmp;
    %                     clear tmp
    %                 end
    %             else
    %                 if contains(FileName,'R2') || contains(FileName,'semi_partial')
    %                     bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
    %                     bootP = squeeze(H0_data(:,:,3,:)); % get all P values under H0
    %                 else
    bootM = squeeze(H0_data(:,:,1,:));
    bootP = squeeze(H0_data(:,:,2,:));
    %                 end

    if size(M,1) == 1
        tmp = NaN(1,size(M,2),size(bootM,2));
        tmp(1,:,:,:) = bootM; bootM = tmp;
        tmp(1,:,:,:) = bootP; bootP = tmp;
        clear tmp
    end
    %             end

    %%%%%%%%%%%%%%%% limo_clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % finally get cluster mask and corrected p-values
%     if contains(FileName,'ttest') || contains(FileName,'LI_Map')
        % [mask,M] = limo_clustering(M.^2,Pval,bootM.^2,bootP,LIMO,MCC,p); % mask and cluster p values
%     end

    function mask = correct_cluster(M, Pval, bootM, bootP, neighbormatrix, MCC, p)

    %%%%%%%%%%%%%%% end of limo_clustering %%%%%%%%%%%%%%%%%%%%%%%%
    Nclust   = unique(mask);
    Nclust   = length(Nclust)-1; % mask = mask>0;
    if Nclust <= 1
        Mclust = 'cluster';
    else
        Mclust = 'clusters';
    end
    mytitle = sprintf('%s cluster correction (%g %s)', titlename, Nclust, Mclust);
    %     end

    % Max-correction
elseif ~isempty(M) && MCC == 4 % Stat max
    %     if exist(['H0' filesep MCC_data],'file')
    %         try
    H0_data = load(['H0' filesep MCC_data]);
    H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
    %             if strcmpi(LIMO.Analysis,'Time-Frequency')
    %                 if contains(FileName,'R2') || contains(FileName,'semi_partial')
    %                     bootM = squeeze(H0_data(:,:,:,2,:)); % get all F values under H0
    %                 else
    %                     bootM = squeeze(H0_data(:,:,:,1,:));
    %                 end
    %             else
    %                 if contains(FileName,'R2') || contains(FileName,'semi_partial')
    %                     bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
    %                 else
    bootM = squeeze(H0_data(:,:,1,:));
    %                 end
    %             end
    clear H0_data;
    [mask,M] = limo_max_correction(abs(M),abs(bootM),p);
    mytitle  = sprintf('%s: correction by max',titlename);
    %         catch ME
    %             errordlg(sprintf('error log: %s \n',ME.message),'max correction failure')
    %             return
    %         end
    %     else
    %         errordlg(['H0' filesep MCC_data ' not found'],'max correction failure')
    %     end

    % TFCE-correction
elseif ~isempty(M) && MCC == 3 % Stat max
    %     if exist(fullfile(LIMO.dir,['tfce' filesep 'tfce_' FileName]),'file')
    %         try
    score    = load(fullfile(LIMO.dir,['tfce' filesep 'tfce_' FileName]));
    score    = score.(cell2mat(fieldnames(score)));
    H0_score = load(fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' FileName]));
    H0_score = H0_score.(cell2mat(fieldnames(H0_score)));
    [mask,M] = limo_max_correction(score,H0_score,p);
    mytitle  = sprintf('%s: correction using TFCE',titlename);
    %         catch ME
    %             errordlg(sprintf('error log: %s \n',ME.message),'tfce correction failure')
    %             return
    %         end
    %     else
    %         errordlg('no tfce tfce file was found','missing data')
    %     end
end

% % Repeated measures ANOVA
% if contains(FileName,'Rep_ANOVA')
%
%     % all files have dim electrode x [freq/time] frames x F/p
%     if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
%         M    = matfile.(cell2mat(fieldnames(matfile)))(:,:,:,1);
%         PVAL = matfile.(cell2mat(fieldnames(matfile)))(:,:,:,2);
%     else
%         M    = matfile.(cell2mat(fieldnames(matfile)))(:,:,1);
%         PVAL = matfile.(cell2mat(fieldnames(matfile)))(:,:,2);
%     end
%     MCC_data = fullfile(LIMO.dir,['H0' filesep 'H0_' FileName]);
%
%     % no correction for multiple testing
%     % -----------------------------------
%     if MCC == 1
%         mask = PVAL <= p;
%         M    = PVAL;
%         if contains(FileName,'Rep_ANOVA_Interaction')
%             mytitle = sprintf('Interaction F-values uncorrected threshold');
%         elseif contains(FileName,'Rep_ANOVA_Gp_effect')
%             mytitle = sprintf('Gp effect F-values uncorrected threshold');
%         elseif contains(FileName,'Rep_ANOVA_Main')
%             mytitle = sprintf('Main Effect F-values uncorrected threshold');
%         end
%
%         % cluster correction for multiple testing
%         % ---------------------------------------
%     elseif MCC == 2
%
%         if exist(MCC_data,'file')
%             try
%                 H0_data = load(MCC_data);
%                 H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
%                 if strcmpi(LIMO.Analysis,'Time-Frequency')
%                     bootT = squeeze(H0_data(:,:,:,1,:));
%                     bootP = squeeze(H0_data(:,:,:,2,:));
%                     if size(M,1) == 1
%                         tmp = NaN(1,size(M,2),size(M,3),size(bootT,4));
%                         tmp(1,:,:,:) = bootT; bootT = tmp;
%                         tmp(1,:,:,:) = bootP; bootP = tmp;
%                         clear tmp
%                     end
%                 else
%                     bootT = squeeze(H0_data(:,:,1,:));
%                     bootP = squeeze(H0_data(:,:,2,:));
%                     if size(M,1) == 1
%                         tmp = NaN(1,size(M,2),size(bootT,4));
%                         tmp(1,:,:) = bootT; bootT = tmp;
%                         tmp(1,:,:) = bootP; bootP = tmp;
%                         clear tmp
%                     end
%                 end
%
%                 if size(M,1) == 1
%                     [mask,M] = limo_clustering(M,PVAL,bootT,bootP,LIMO,3,p); % temporal clustering
%                 else
%                     [mask,M] = limo_clustering(M,PVAL,bootT,bootP,LIMO,2,p); % spatial-temporal clustering
%                 end
%                 Nclust = unique(mask); Nclust = length(Nclust)-1; % mask = mask>0;
%                 if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
%                 if contains(FileName,'Rep_ANOVA_Interaction')
%                     mytitle = sprintf('Interaction F-values cluster correction (%g %s)', Nclust, Mclust);
%                 elseif contains(FileName,'Rep_ANOVA_Gp_effect')
%                     mytitle = sprintf('Gp effect F-values cluster correction (%g %s)', Nclust, Mclust);
%                 elseif contains(FileName,'Rep_ANOVA_Main')
%                     mytitle = sprintf('Main effect F-values cluster correction (%g %s)', Nclust, Mclust);
%                 end
%
%             catch ME
%                 errordlg(sprintf('error log: %s \n',ME.message),'cluster correction failure')
%                 return
%             end
%         else
%             errordlg(['H0' filesep MCC_data ' not found'],'cluster correction failure')
%             return
%         end
%
%
%         % T max correction for multiple testing
%         % -------------------------------------
%     elseif MCC == 4 % Stat max
%         if exist(MCC_data,'file')
%             try
%                 H0_data = load(MCC_data);
%                 H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
%                 if strcmpi(LIMO.Analysis,'Time-Frequency')
%                     bootT = squeeze(H0_data(:,:,:,1,:));
%                     if size(M,1) == 1
%                         tmp = NaN(1,size(M,2),size(M,3),size(bootT,4));
%                         tmp(1,:,:,:) = bootT; bootT = tmp;
%                         clear tmp
%                     end
%                 else
%                     bootT = squeeze(H0_data(:,:,1,:));
%                     if size(M,1) == 1
%                         tmp = NaN(1,size(M,2),size(bootT,3));
%                         tmp(1,:,:) = bootT; bootT = tmp;
%                         clear tmp
%                     elseif size(M,2) == 1 %% for Weight bias testing
%                         tmp = NaN(size(M,1),1,size(bootT,2));
%                         tmp(:,1,:) = bootT; bootT = tmp;
%                         clear tmp
%                     end
%                 end
%
%                 [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
%                 if strncmp(FileName,'Rep_ANOVA_Interaction',21)
%                     mytitle = sprintf('Interaction correction by T max');
%                 elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
%                     mytitle = sprintf('Gp effect correction by T max');
%                 elseif strncmp(FileName,'Rep_ANOVA',9)
%                     mytitle = sprintf('Main Effect correction by T max');
%                 end
%             catch ME
%                 errordlg(sprintf('error log: %s \n',ME.message),'max correction failure')
%                 return
%             end
%
%         else
%             errordlg(['H0' filesep MCC_data ' not found'],'max correction failure')
%             return
%         end
%
%         % Correction using TFCE
%         % -------------------------------------
%     elseif MCC == 3 % Stat tfce
%         tfce_data    = sprintf('tfce%stfce_%s',filesep, FileName);
%         H0_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
%         if exist(tfce_data,'file') && exist(H0_tfce_data,'file')
%             try
%                 tfce_data    = load(tfce_data);
%                 tfce_data    = tfce_data.(cell2mat(fieldnames(tfce_data)));
%                 H0_tfce_data = load(H0_tfce_data);
%                 H0_tfce_data = H0_tfce_data.(cell2mat(fieldnames(H0_tfce_data)));
%                 [mask,M]     = limo_max_correction(tfce_data, H0_tfce_data,p);
%                 if strncmp(FileName,'Rep_ANOVA_Interaction',21)
%                     mytitle = sprintf('Interaction correction using TFCE');
%                 elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
%                     mytitle = sprintf('Gp effect correction using TFCE');
%                 elseif strncmp(FileName,'Rep_ANOVA',9)
%                     mytitle = sprintf('Main Effect correction using TFCE');
%                 end
%             catch ME
%                 errordlg(sprintf('error log: %s \n',ME.message),'tfce correction failure')
%                 return
%             end
%         else
%             errordlg('no tfce file or bootstrap file was found to compute the max distribution','missing data')
%             return
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%% end of limo_stat_values back to limo_display_results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assignin('base','p_values',squeeze(M))
assignin('base','mask',squeeze(mask))

% if strcmpi(LIMO.Analysis,'Time-Frequency') || strcmpi(LIMO.Analysis,'ITC')
%     if contains(FileName,'R2') || ...
%             contains(FileName,'semi_partial')
%         toplot = squeeze(toplot(:,:,:,1));
%     elseif contains(FileName,'ttest','IgnoreCase',true) || ...
%             contains(FileName,'LI_Map','IgnoreCase',true)
%         toplot = squeeze(toplot(:,:,:,4));
%     elseif strncmp(FileName,'con_',4)
%         toplot = squeeze(toplot(:,:,:,4));
%     elseif strncmp(FileName,'ess_',4)
%         if ~exist('ess','var')
%             effect_nb = eval(FileName(5:end-4)); %#ok<NASGU>
%         end
%         toplot = squeeze(toplot(:,:,:,end-1));
%     elseif contains(FileName,'Condition') || ...
%             contains(FileName,'Covariate') || ...
%             contains(FileName,'Rep_ANOVA')
%         toplot = squeeze(toplot(:,:,:,1));
%     else
%         disp('file no supported'); return
%     end
% else
%     if contains(FileName,'R2') || ...
%             contains(FileName,'semi_partial')
%         toplot = squeeze(toplot(:,:,1));
%     elseif contains(FileName,'ttest','IgnoreCase',true) || ...
%             contains(FileName,'LI_Map','IgnoreCase',true)
toplot = squeeze(toplot(:,:,4));
%     elseif strncmp(FileName,'con_',4)
%         toplot = squeeze(toplot(:,:,4));
%     elseif strncmp(FileName,'ess_',4)
%         toplot = squeeze(toplot(:,:,4));
%     elseif contains(FileName,'Condition') || ...
%             contains(FileName,'Covariate') || ...
%             contains(FileName,'Rep_ANOVA')
%         toplot = squeeze(toplot(:,:,1));
%     else
%         disp('file no supported'); return
%     end
% end
assignin('base','stat_values',toplot)
% data_cached = 0;

%      Image and topoplot
if Type == 1 || Type == 2

    %     % cache the results for next time
    %     % ------------------------------
    %     if data_cached == 0
    %         LIMO.cache.fig.name       = FileName;
    %         LIMO.cache.fig.MCC        = MCC;
    %         LIMO.cache.fig.stats      = toplot;
    %         LIMO.cache.fig.threshold  = p;
    %         LIMO.cache.fig.pval       = squeeze(M);
    %         LIMO.cache.fig.mask       = squeeze(mask);
    %         LIMO.cache.fig.title      = mytitle;
    %         if exist(LIMO.dir,'dir')
    %             save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
    %         else
    %             disp('Cached data in LIMO.mat cannot be updated - LIMO dir doesn''t exist (likely moved files)')
    %         end
    %     end

    %%%% Image all results %%%%
    if Type == 1 && ~strcmpi(LIMO.Analysis,'Time-Frequency') && ~strcmpi(LIMO.Analysis,'ITC')
        %         limo_display_image(LIMO,toplot,mask,mytitle,surfflag)
        %%%%% limo_display_image %%%%%%%%%%%%5

        % what do we plot?  the data (toplot) masked (tpically of significance)
        scale           = toplot.*single(mask>0);
        scale(scale==0) = NaN;
        cc              = limo_color_images(scale); % get a color map commensurate to that

        v = max(scale(:));       % from the 2D data to plot, find max
        [e,f] = find(scale == v);    % which channel and time/frequency frame
        if length(e) > 1           % if we have multiple times the exact same max values
            e = e(1); f = f(1);  % then take the 1st (usually an artefact but allows to see it)
        end

        % for each cluster, get start/end/max value
        % if unthresholded, uncorrected, tfce or max = mask is made up of ones
        n_cluster     = max(mask(:));
        cluster_start = NaN(1,n_cluster); % start of each cluster
        cluster_end   = NaN(1,n_cluster); % end of each cluster
        cluster_maxv  = NaN(1,n_cluster); % max value for each cluster
        cluster_maxe  = NaN(1,n_cluster); % channel location of the max value of each cluster
        cluster_maxf  = NaN(1,n_cluster); % frame location of the max value of each cluster

        for c=1:n_cluster
            tmp                               = toplot.*(mask==c);
            tmp(tmp==Inf)                     = NaN;
            tmp(tmp==-Inf)                    = NaN;
            sigframes                         = sum(tmp,1);
            cluster_start(c)                  = find(sigframes,1,'first');
            cluster_end(c)                    = find(sigframes,1,'last');
            [V,type]                          = max([abs(min(tmp(:))) max(tmp(:))]);
            if type == 2
                cluster_maxv(c)               = V(1);
            else
                V = -V;
                cluster_maxv(c)               = V(1);
            end
            [cluster_maxe(c),cluster_maxf(c)] = ind2sub(size(tmp),find(tmp==V(1)));
        end

        %% get frame information
        %         if strcmpi(LIMO.Analysis,'Time')
        %             if isfield(LIMO.data,'timevect')
        %                 timevect = LIMO.data.timevect;
        %                 if size(timevect,2) == 1; timevect = timevect'; end
        %             else
        %                 timevect = [];
        %             end
        %
        %             if size(timevect,2) ~= size(toplot,2)
        %                 timevect           = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
        %                 LIMO.data.timevect =  timevect;
        %                 if exist(LIMO.dir,'dir')
        %                     save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
        %                 end
        %             end
        %
        %             ratio =  abs(timevect(end)-timevect(1)) / length(timevect); % this the diff in 'size' between consecutive frames
        %             if LIMO.data.start < 0
        %                 frame_zeros = find(timevect == 0);
        %                 if isempty(frame_zeros)
        %                     frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
        %                 end
        %             else
        %                 frame_zeros = -round(min(timevect)/ ratio);
        %             end
        %
    elseif strcmpi(LIMO.Analysis,'Frequency')
        %         if isfield(LIMO.data,'freqlist')
        freqvect = LIMO.data.freqlist;
        %             freqvect = 1:15;
        %             if size(freqvect,2) == 1
        %                 freqvect = freqvect';
        %             end
        %         else
        %             freqvect = [];
        %         end

        if size(freqvect,2) ~= size(toplot,2)
            freqvect           = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
            LIMO.data.freqlist = freqvect;
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        end

        ratio =  abs(freqvect(end)-freqvect(1)) / length(freqvect);
        %         if LIMO.data.start < 0
        %             frame_zeros = find(freqvect == 0);
        %             if isempty(frame_zeros)
        %                 frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
        %             end
        %         else
        %             frame_zeros = -round(min(freqvect)/ ratio);
        %         end
        frame_zeros = 2;
        %     elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        %         if isfield(LIMO.data,'tf_times')
        %             timevect = LIMO.data.tf_times;
        %             if size(timevect,2) == 1; timevect = timevect'; end
        %         else
        %             timevect = [];
        %         end
        %
        %         if size(timevect,2) ~= size(toplot,2)
        %             timevect           = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
        %             LIMO.data.tf_times = timevect;
        %             save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        %         end
        %
        %         ratio =  abs(timevect(end)-timevect(1)) / length(timevect); % this the diff in 'size' between consecutive frames
        %         if LIMO.data.start < 0
        %             frame_zeros = find(timevect == 0);
        %             if isempty(frame_zeros)
        %                 frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
        %             end
        %         else
        %             frame_zeros = -round(min(timevect)/ ratio);
        %         end
        %
        %         if isfield(LIMO.data,'tf_freqs')
        %             freqvect = LIMO.data.tf_freqs;
        %             if size(freqvect,2) == 1; freqvect = freqvect'; end
        %         else
        %             freqvect = [];
        %         end
        %
        %         if size(freqvect,2) ~= size(toplot,1)
        %             freqvect           = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,1));
        %             LIMO.data.tf_freqs =  freqvect;
        %             save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        %         end
        %
        %     else
        %         error('LIMO.Analysis unspecfied')
        %     end

        %     if isempty(mytitle)
        %         if isfield(LIMO.design,'name')
        %             mytitle = LIMO.design.name;
        %         else
        mytitle = ' ';
        %         end
        %     end

        % make the main figure
        figure; set(gcf,'Color','w','InvertHardCopy','off');

        % course plot at best electrode
        ax(3) = subplot(3,3,9);
        %     if ~isfield(LIMO.data, 'chanlocs') || isfield(LIMO.data,'expected_chanlocs')
        %         LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
        %     end

        %     if size(toplot,1) == 1
        %         if strcmpi(LIMO.Analysis,'Time')
        %             plot(timevect,toplot,'LineWidth',3);
        %         elseif strcmpi(LIMO.Analysis,'Frequency')
        %             plot(freqvect,toplot,'LineWidth',3);
        %         end
        %         grid on; ylabel('stat value'); axis tight
        %
        %         if isfield(LIMO,'Type')
        %             if strcmpi(LIMO.Type,'Components')
        %                 mytitle2 = 'Average component';
        %             elseif strcmpi(LIMO.Type,'Channels')
        %                 mytitle2 = 'Average channel';
        %             end
        %         else
        %             mytitle2 = 'Average channel';
        %         end
        %     else
        if strcmpi(LIMO.Analysis,'Time')
            plot(timevect,toplot(e,:),'LineWidth',3);
        elseif strcmpi(LIMO.Analysis,'Frequency')
            plot(freqvect,toplot(e,:),'LineWidth',3);
        elseif strcmpi(LIMO.Analysis,'Time-Frequency')
            plot(timevect,toplot(e,:),'LineWidth',3);
            mytitle2 = sprintf('stat values @ %g Hz', freqvect(e));
        end
        grid on; axis tight

        %         if ~strcmpi(LIMO.Analysis,'Time-Frequency')
        %             if isfield(LIMO,'Type')
        %                 if strcmpi(LIMO.Type,'Components')
        %                     mytitle2 = sprintf('stat values @ \n component %g', e);
        %                 elseif strcmpi(LIMO.Type,'Channels')
        label = LIMO.data.chanlocs(e).labels;
        mytitle2 = sprintf('stat values @ \n channel %s (%g)', label,e);
        %                 end
        %             else
        %                 try
        %                     label = LIMO.data.chanlocs(e).labels;
        %                     mytitle2 = sprintf('stat values @ \n channel %s (%g)', label.labels,e);
        %                 catch
        %                     mytitle2 = sprintf('stat values @ y=%g', e);
        %                 end
        %             end
        %         end
    end
    title(mytitle2,'FontSize',12)

    % topoplot at max time
    if size(toplot,1) ~= 1 && ~strcmpi(LIMO.Analysis,'Time-Frequency')

        ax(2) = subplot(3,3,6);
        opt   = {'maplimits','maxmin','verbose','off','colormap', limo_color_images(toplot)};

        %         if isfield(LIMO,'Type')
        %             if strcmpi(LIMO.Type,'Components')
        %                 opt = {'maplimits','absmax','electrodes','off','verbose','off','colormap', limo_color_images(toplot)};
        %                 topoplot(toplot(:,f),LIMO.data.chanlocs,opt{:});
        %             else
        %                 topoplot(toplot(:,f),LIMO.data.chanlocs,opt{:});
        %             end

        %             if size(toplot,2) == 1
        %                 title('Topoplot','FontSize',12)
        %             else
        %                 if strcmpi(LIMO.Analysis,'Time')
        %                     title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
        %                     set(gca,'XTickLabel', timevect);
        %                 elseif strcmpi(LIMO.Analysis,'Frequency')
        title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
        set(gca,'XTickLabel', LIMO.data.freqlist);
        %                 end
        %             end

    elseif ~isempty(LIMO.data.chanlocs)
        topoplot(toplot(:,f),LIMO.data.chanlocs,opt{:});
        if size(toplot,2) == 1
            title('Topoplot','FontSize',12)
        else
            %                 if strcmpi(LIMO.Analysis,'Time')
            %                     title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
            %                     set(gca,'XTickLabel', timevect);
            %                 elseif strcmpi(LIMO.Analysis,'Frequency')
            title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
            set(gca,'XTickLabel', LIMO.data.freqlist);
            %                 end
        end
    end
end

% images toplot
ax(1) = subplot(3,3,[1 2 4 5 7 8]);
%     if strcmpi(LIMO.Analysis,'Time')
%         imagesc(timevect,1:size(toplot,1),scale);
%         colormap(gca, cc);
%     elseif strcmpi(LIMO.Analysis,'Frequency')
imagesc(freqvect,1:size(toplot,1),scale);
colormap(gca, cc);
%     elseif strcmpi(LIMO.Analysis,'Time-Frequency')
%         imagesc(timevect,freqvect,scale);
%         colormap(gca, cc);
%     end
set_imgaxes(LIMO,scale);
title(mytitle,'Fontsize',12)

% return cluster info
warning off
if contains(mytitle,'cluster')
    for c=1:n_cluster
        %             if strcmpi(LIMO.Analysis,'Time')
        %                 fprintf('cluster %g starts at %gms ends at %gms, max %g @ %gms channel %s \n', c, ...
        %                     timevect(cluster_start(c)),timevect(cluster_end(c)), cluster_maxv(c), timevect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
        %             elseif strcmpi(LIMO.Analysis,'Frequency')
        fprintf('cluster %g starts at %gHz ends at %gHz, max %g @ %gHz channel %s \n', c, ...
            freqvect(cluster_start(c)),freqvect(cluster_end(c)), cluster_maxv(c), freqvect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
        %             elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        %                 f1 = scale(:,cluster_start(c)); f2 = scale(:,cluster_end(c)); f3 = scale(:,cluster_maxf(c));
        %                 fprintf('cluster %g at %gms %gHz, ends at %gms %gHz, max %g @ %gms %gHz \n', c, ...
        %                     timevect(cluster_start(c)),freqvect(find(f1==max(f1))), ...
        %                     timevect(cluster_end(c)), freqvect(find(f2==max(f2))), ...
        %                     cluster_maxv(c), timevect(cluster_maxf(c)), freqvect(find(f3==max(f3))));
        %             end
    end
else % no clusters (we have make one)
    %         if strcmpi(LIMO.Analysis,'Time')
    %             fprintf('1st significant frame at %gms, last signifiant frame at %gms, max %g @ %gms channel %s \n', ...
    %                 timevect(cluster_start(c)),timevect(cluster_end(c)), cluster_maxv(c), timevect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
    %         elseif strcmpi(LIMO.Analysis,'Frequency')
    fprintf('1st significant frame at %gHz, last signifiant frame at %gHz, max %g @ %gHz channel %s \n', ...
        freqvect(cluster_start(c)),freqvect(cluster_end(c)), cluster_maxv(c), freqvect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
    %         elseif strcmpi(LIMO.Analysis,'Time-Frequency')
    %             f1 = scale(:,cluster_start(c)); f2 = scale(:,cluster_end(c)); f3 = scale(:,cluster_maxf(c));
    %             fprintf('1st significant frame at %gms %gHz, last signifiant frame at %gms %gHz, max %g @ %gms %gHz \n', ...
    %                 timevect(cluster_start(c)),freqvect(find(f1==max(f1))), ...
    %                 timevect(cluster_end(c)), freqvect(find(f2==max(f2))), ...
    %                 cluster_maxv(c), timevect(cluster_maxf(c)), freqvect(find(f3==max(f3))));
    %         end
end
warning on

%     % update with mouse clicks
%     if dynamic == 1
%         if size(toplot,1) > 1
%             update = 0;
%             while update ==0
%                 try
%                     [x,y,button] = ginput(1);
%                 catch
%                     break
%                 end
%
%                 if button > 1
%                     update = 1; % right click to come out of the dynamic figure
%                 end
%
%                 clickedAx = gca;
%                 if clickedAx ~=ax(1)
%                     disp('right click to exit')
%                 else
%                     % topoplot at new time or freq
%                     % because x axis is already with the correct value
%                     frame = frame_zeros + round(x / ratio);
%                     if frame<=0; frame = 1; end
%                     if frame>=size(toplot,2); frame=size(toplot,2); end
%
%                     if strcmpi(LIMO.Analysis,'Time-Frequency')
%                         % course plot at selected frequency
%                         y          = round(y);
%                         [~,yframe] = min(abs(freqvect - y));
%                         if size(toplot,1)> 1 && yframe>size(toplot,1)
%                             yframe = size(toplot,1);
%                             y      = freqvect(yframe);
%                         elseif size(toplot,1)> 1 && y<1
%                             yframe = 1;
%                             y      = freqvect(yframe);
%                         end
%                     else
%                         % course plot at selected  channel
%                         y = round(y);
%                         if size(toplot,1)> 1 && y>size(toplot,1)
%                             y = size(toplot,1);
%                         elseif size(toplot,1)> 1 && y<1
%                             y = 1;
%                         end
%                     end
%
%                     if strcmpi(LIMO.Analysis,'Time')
%                         if ~contains(LIMO.design.name, ['one ' LIMO.Type(1:end-1)]) && ~isempty(LIMO.data.chanlocs)
%                             subplot(3,3,6,'replace');
%                             if size(toplot,2) == 1
%                                 topoplot(toplot(:,1),LIMO.data.chanlocs,opt{:});
%                                 title('topoplot','FontSize',12)
%                             else
%                                 topoplot(toplot(:,frame),LIMO.data.chanlocs,opt{:});
%                                 title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',12)
%                             end
%                         end
%
%                     elseif strcmpi(LIMO.Analysis,'Frequency')
%                         if ~contains(LIMO.design.name, ['one ' LIMO.Type(1:end-1)]) && ~isempty(LIMO.data.chanlocs)
%                             subplot(3,3,6,'replace');
%                             if size(toplot,2) == 1
%                                 topoplot(toplot(:,1),LIMO.data.chanlocs,opt{:});
%                                 title('topoplot','FontSize',12)
%                             else
%                                 topoplot(toplot(:,frame),LIMO.data.chanlocs,opt{:});
%                                 title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
%                             end
%                         end
%                     end
%
%                     subplot(3,3,9,'replace');
%                     if size(toplot,2) == 1
%                         bar(toplot(y,1)); grid on; axis([0 2 0 max(toplot(:))+0.2]); ylabel('stat value')
%                         if isfield(LIMO,'Type')
%                             if strcmpi(LIMO.Type,'Components')
%                                 mytitle2 = sprintf('component %g', y);
%                             elseif strcmpi(LIMO.Type,'Channels')
%                                 mytitle2 = sprintf('channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
%                             end
%                         else
%                             mytitle2 = sprintf('channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
%                         end
%                     else
%                         if strcmpi(LIMO.Analysis,'Time')
%                             plot(timevect,toplot(y,:),'LineWidth',3);
%                         elseif strcmpi(LIMO.Analysis,'Frequency')
%                             plot(freqvect,toplot(y,:),'LineWidth',3);
%                         elseif strcmpi(LIMO.Analysis,'Time-Frequency')
%                             plot(timevect,toplot(yframe,:),'LineWidth',3);
%                             mytitle2 = sprintf('stat values @ %g Hz', y);
%                         end
%                         grid on; axis tight
%
%                         if ~strcmpi(LIMO.Analysis,'Time-Frequency')
%                             if isfield(LIMO,'Type')
%                                 if strcmpi(LIMO.Type,'Components')
%                                     mytitle2 = sprintf('stat values @ \n component %g', y);
%                                 elseif strcmpi(LIMO.Type,'Channels')
%                                     mytitle2 = sprintf('stat values @ \n channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
%                                 end
%                             else
%                                 try
%                                     mytitle2 = sprintf('stat values @ \n channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
%                                 catch
%                                     mytitle2 = sprintf('stat values @ \n y=%g)', y);
%                                 end
%                             end
%                         end
%                     end
%                     title(mytitle2,'FontSize',12);
%
%                     subplot(3,3,[1 2 4 5 7 8]); colormap(gca, cc);
%
%                     try
%                         p_values = evalin('base','p_values');
%                         if strcmpi(LIMO.Analysis,'Time-Frequency')
%                             if ~isnan(p_values(yframe,frame))
%                                 fprintf('Stat value: %g, p_value %g \n',toplot(yframe,frame),p_values(yframe,frame));
%                             else
%                                 fprintf('Stat value: %g, p_value is a NaN? \n',toplot(yframe,frame));
%                             end
%                         else
%                             if ~isnan(p_values(round(y),frame))
%                                 fprintf('Stat value: %g, p_value %g \n',toplot(round(y),frame),p_values(round(y),frame));
%                             else
%                                 fprintf('Stat value: %g, p_value is a NaN? \n',toplot(round(y),frame));
%                             end
%                         end
%                     catch pvalerror
%                         fprintf('couldn''t figure the stats values?? %s \n',pvalerror.message)
%                     end
%                 end
%             end
%         end
%     end

%     elseif Type == 1 && strcmpi(LIMO.Analysis,'Time-Frequency') || ...
%             Type == 1 && strcmpi(LIMO.Analysis,'ITC')
%         if ndims(toplot)==3
%             limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
%         else
%             limo_display_image(LIMO,squeeze(toplot),squeeze(mask),mytitle,flag)
%         end

% elseif Type == 2
%--------------------------
% topoplot
%--------------------------

%         if strcmpi(LIMO.Analysis,'Time-Frequency')
%             errordlg('topoplot not supported for time-frequency analyses')
%         else
%             if isfield(LIMO.design,'channel')  % not full scalp
%                 if ~isempty(LIMO.design.electrode)
%                     msgbox('Only one channel found','No topoplot')
%                     return
%                 end
%             end
%         end

%         if sum(mask(:)) == 0
%             warndlg('no values under threshold','no significant effect');
%         else
EEG.data     = toplot;
EEG.setname  = mytitle;
EEG.chanlocs = LIMO.data.chanlocs;

if size(toplot,2) == 1
    opt = {'maplimits','maxmin','verbose','off'};
    if isfield(LIMO,'Type')
        if strcmpi(LIMO.Type,'Components')
            opt = {'maplimits','absmax','electrodes','off','verbose','off'};
        end
    end
    figure; set(gcf,'Color','w','InvertHardCopy','off');
    topoplot(toplot(:,1),EEG.chanlocs,opt{:});
    title('Topoplot','FontSize',12)
else
    if strcmpi(LIMO.Analysis,'Time')
        EEG.xmin  = LIMO.data.start/1000; % in sec
        EEG.xmax  = LIMO.data.end/1000;   % in sec
        EEG.times = (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
        call_topolot(EEG,FileName,LIMO.Analysis)
        % pop_topoplot(EEG);
    elseif strcmpi(LIMO.Analysis,'Frequency')
        EEG.xmin  = LIMO.data.freqlist(1);
        EEG.xmax  = LIMO.data.freqlist(end);
        freqlist  = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
        %                     if isempty(freqlist)
        %                         return
        %                     else
        EEG.freq = str2double(cell2mat(freqlist));
        if min(EEG.freq)<EEG.xmin || max(EEG.freq)>EEG.xmax
            errordlg('selected frequency out of bound'); return
        end
        %                     end
        %                     call_topolot(EEG,FileName,LIMO.Analysis)
        %                 function call_topolot(EEG,FileName,Domain)
        Domain = LIMO.Analysis;
        EEG.pnts     = size(EEG.data,2);
        EEG.nbchan   = size(EEG.data,1);
        EEG.trials   = 1;

        %                 if strcmpi(FileName,'R2.mat') || strcmpi(FileName,'R2')
        %                     newname = 'R^2';
        %                 else
        newname = [FileName(1:min(strfind(FileName,'_'))-1),FileName(max(strfind(FileName,'_'))+1:end-4)];
        %                 end

        if strcmpi(Domain,'Time')
            if ~isfield(EEG,'setname')
                if contains(FileName,'con','IgnoreCase',true)
                    EEG.setname = sprintf('%s - T values',newname);
                else
                    EEG.setname = sprintf('%s - F values',newname);
                end
            end
            pop_topoplot(EEG);
            % set(gca,'Colormap',limo_color_images(EEG.data),'CLim',[min(EEG.data(:)),max(EEG.data(:))])
        else % freq
            N = size(EEG.freq,2);
            figure;
            for f=1:N
                if N<=6
                    subplot(1,N,f)
                else
                    subplot(ceil(N/6),6,f);
                end
                [~,ind] = min(abs(EEG.freq-EEG.freq(f)));
                opt = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', limo_color_images(EEG.data(:,ind))};
                topoplot(EEG.data(:,ind),EEG.chanlocs,opt{:});
                if isfield(EEG,'setname')
                    title(sprintf('Frequency %g Hz from \n%s',round(EEG.freq(ind)),EEG.setname));
                else
                    title(sprintf('Frequency %g Hz from %s - F values',round(EEG.freq(ind)),newname));
                end
            end
        end
        %             end
        assignin('base','Plotted_data',EEG.data)
    end
end
% end

% elseif Type == 3
%
%     %--------------------------
%     % Course plot
%     %--------------------------
%
%     if contains(FileName,'one_sample','IgnoreCase',true) || contains(FileName,'two_samples','IgnoreCase',true) || ...
%             contains(FileName,'paired_samples','IgnoreCase',true) || contains(FileName,'con_','IgnoreCase',true) || ...
%             contains(FileName,'ess_','IgnoreCase',true)
%         % ------------------------------------------------------------------------------------------------------------
%         % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
%         % H0 file dim = (electrodes,frames,[t, p],nboot)
%
%         data = load(fullfile(PathName,FileName));
%         data = data.(cell2mat(fieldnames(data)));
%         if strcmpi(LIMO.Analysis,'Time-Frequency')
%             [~,channel,freq,time] = limo_display_reducedim(data(:,:,:,[4 5]),LIMO,g.channels,g.restrict,g.dimvalue);
%             data                  = squeeze(data(channel,freq,time,:,:)); % 2D
%             sig                   = squeeze(single(mask(channel,freq,time))); %1D
%         else
%             [~,channel,freq,time] = limo_display_reducedim(data(:,:,[4 5]),LIMO,g.channels);
%             data                  = squeeze(data(channel,time,:));
%             sig                   = single(mask(channel,:));
%         end
%         sig(sig==0)=NaN;
%
%         % compute
%         trimci      = NaN(size(data,1),3);
%         trimci(:,2) = data(:,1); % mean values
%         if contains(FileName,'ess','IgnoreCase',true)
%             start_at = max(strfind(FileName,'_'))+1;
%             C = LIMO.contrast{eval(FileName(start_at:end-4))}.C;
%             df = rank(C); % rank of the relevant contrast
%             trimci(:,1) = squeeze(trimci(:,2))-(finv(1-p./2*size(C,1),df,data(:,3)).*data(:,2));
%             trimci(:,3) = squeeze(trimci(:,2))+(finv(1-p./2*size(C,1),df,data(:,3)).*data(:,2));
%         else
%             trimci(:,1) = squeeze(trimci(:,2))-(tinv(1-p./2,data(:,3)).*data(:,2));
%             trimci(:,3) = squeeze(trimci(:,2))+(tinv(1-p./2,data(:,3)).*data(:,2));
%         end
%
%         % plot
%         if strcmpi(LIMO.Analysis,'Time')
%             if isfield(LIMO.data,'timevect')
%                 xvect = LIMO.data.timevect;
%             else
%                 xvect = [];
%             end
%
%             if size(xvect,2) ~= size(toplot,2)
%                 xvect              = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
%                 LIMO.data.timevect = xvect;
%                 save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
%             end
%
%         elseif strcmpi(LIMO.Analysis,'Frequency')
%             if isfield(LIMO.data,'freqlist')
%                 xvect=LIMO.data.freqlist;
%             else
%                 xvect = [];
%             end
%
%             if size(xvect,2) ~= size(toplot,2)
%                 xvect              = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
%                 LIMO.data.freqlist = xvect;
%                 save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
%             end
%
%         elseif strcmpi(LIMO.Analysis,'Time-Frequency')
%             if length(time) > 1 && isfield(LIMO.data,'tf_times')
%                 xvect = LIMO.data.tf_times;
%             elseif length(freq) > 1 && isfield(LIMO.data,'tf_freqs')
%                 xvect = LIMO.data.tf_freqs;
%             else
%                 xvect = [];
%             end
%
%             if length(time) > 1&& size(xvect,2) ~= size(data,1)
%                 xvect              = linspace(LIMO.data.start,LIMO.data.end,size(data,1));
%                 LIMO.data.tf_times =  xvect;
%                 save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
%             elseif length(freq) > 1 && size(xvect,2) ~= size(data,1)
%                 xvect              = linspace(LIMO.data.lowf,LIMO.data.highf,size(data,1));
%                 LIMO.data.tf_freqs =  xvect;
%                 save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
%             end
%         end
%
%         figure;
%         set(gcf,'Color','w')
%         plot(xvect,squeeze(trimci(:,2)),'LineWidth',3);
%         fillhandle = patch([xvect,fliplr(xvect)], [trimci(:,1)' fliplr(trimci(:,3)')], [1 0 0]);
%         set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);% set edge color
%         grid on; box on; axis tight
%         h = axis;  hold on;
%         plot(xvect,(sig./10+1).*h(3),'k.','MarkerSize',20)
%         if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
%             xlabel('Time in ms','FontSize',14)
%             ylabel('Amplitude (A.U.)','FontSize',14)
%         elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
%             xlabel('Frequency in Hz','FontSize',14)
%             ylabel('Spectral Power (A.U.)','FontSize',14)
%         end
%         if isempty(LIMO.design.electrode)
%             title(sprintf('%s \n%s %s %s (%g)',mytitle,'Mean values',LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel),'FontSize',16); drawnow;
%         else
%             title(sprintf('%s \n%s virtual %s',mytitle,'Mean values',LIMO.Type(1:end-1)),'FontSize',16); drawnow;
%         end
%         assignin('base','Plotted_data',trimci);
%
%
%     elseif contains(LIMO.design.name,'regression','IgnoreCase',true) && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
%             contains(LIMO.design.name,'ANOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
%             contains(LIMO.design.name,'ANCOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
%         % --------------------------------------------------------------------------------
%
%         % which variable(s) to plot
%         % ----------------------
%         if size(LIMO.design.X,2) >= 2
%             if contains(FileName,'Condition_effect_') || ...
%                     contains(FileName,'Covariate_effect_')
%                 regressor = eval(FileName(18:end-4));
%                 if contains(FileName,'Covariate_effect_')
%                     regressor = regressor+sum(LIMO.design.nb_conditions);
%                 end
%             else
%                 input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2)-1);
%                 regressor = inputdlg(input_title,'Plotting option');
%             end
%
%             if isempty(regressor)
%                 warning on
%                 warning('couldn''t figure out the regressor number/column, plot aborded')
%                 return
%             end
%
%             try
%                 if iscell(regressor)
%                     regressor = sort(eval(cell2mat(regressor)));
%                 end
%                 if max(regressor) > size(LIMO.design.X,2)
%                     errordlg('invalid regressor number');
%                 end
%             catch reginput_error
%                 fprintf('error with regressor numbers/columns line 1373:\n %s',reginput_error.message)
%                 return
%             end
%         else
%             regressor = 1;
%         end
%
%         categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
%         if max(regressor) == size(LIMO.design.X,2)
%             tmp = regressor(1:end-1);
%         else
%             tmp = regressor;
%         end
%
%         if sum(tmp<=categorical) >=1 && sum(tmp>categorical) >=1
%             errordlg('you can''t plot categorical and continuous regressors together'); return
%         end
%
%         % load the effect
%         % --------------
%         data = load(FileName);
%         data = data.(cell2mat(fieldnames(data)));
%         if numel(size(data)) == 3 && size(data,2) == 1
%             errordlg2('single time point detected, plot aborded'); return
%         elseif numel(size(data)) == 4 && size(data,3) == 1
%             errordlg2('single time point detected, plot aborded'); return
%         end
%
%         % which course plot to make
%         % -------------------------
%         if isempty(g.plot3type)
%             extra = questdlg('Plotting data','Options','Original data','Modelled data','Adjusted data','Modelled data');
%         else
%             extra = g.plot3type;
%             % allow typos
%             if contains(extra,'orig','Ignorecase',true)
%                 extra = 'Original';
%             elseif contains(extra,'Model','Ignorecase',true)
%                 extra = 'Modelled';
%             elseif contains(extra,'Adj','Ignorecase',true)
%                 extra = 'Adjusted';
%             else
%                 if exist(errodlg2,'file')
%                     errordlg2(sprintf('input option ''%s'' invalid',extra)); return
%                 else
%                     errordlg(sprintf('input option ''%s'' invalid',extra)); return
%                 end
%             end
%         end
%
%         if isempty(extra)
%             return
%         elseif strcmpi(extra,'Original data')
%             if regressor == size(LIMO.design.X,2)
%                 errordlg('you can''t plot adjusted mean for original data'); return
%             end
%         end
%
%         if strcmpi(LIMO.Analysis,'Time-Frequency')
%             [~,channel,freq,time] = limo_display_reducedim(data,LIMO,g.channels,g.restrict,g.dimvalue);
%             if length(freq) == 1; g.restrict = 'time';
%             else; g.restrict = 'frequency'; end
%             sig         = squeeze(single(mask(channel,freq,time))); %1D
%         else
%             [~,channel] = limo_display_reducedim(data,LIMO,g.channels);
%             sig         = single(mask(channel,:));
%         end
%         sig(sig==0)=NaN;
%         clear data
%
%         % down to business
%         % ----------------------
%         probs = [p/2; 1-p/2];
%         z     = norminv(probs);
%         Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
%         Yr    = Yr.Yr;
%
%         if contains(extra,'Original','Ignorecase',true)
%             if regressor <= length(LIMO.design.nb_conditions) && ...
%                     LIMO.design.nb_conditions ~= 0 % for categorical variables
%                 if length(LIMO.design.nb_conditions) == 1
%                     start = 1;
%                 else
%                     start = sum(LIMO.design.nb_conditions(1:regressor-1));
%                 end
%
%                 for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
%                     index{i} = find(LIMO.design.X(:,i));
%                     if strcmpi(LIMO.Analysis,'Time-Frequency')
%                         data = squeeze(Yr(channel,freq,:,index{i}));
%                     else
%                         data = squeeze(Yr(channel,:,index{i}));
%                     end
%                     average(i,:) = nanmean(data,2);
%                     se           = (nanstd(data,0,2) ./ sqrt(numel(index{i})));
%                     ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Yr,2));
%                 end
%                 clear mytitle
%                 if isempty(LIMO.design.electrode)
%                     mytitle = sprintf('Original subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
%                 else
%                     mytitle = sprintf('Original subjects'' parameters at optimized channel');
%                 end
%             else % continuous variable
%                 for i=max(regressor):-1:min(regressor)
%                     index{i}           = find(LIMO.design.X(:,i));
%                     [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
%                     reg_values(i,:)    = LIMO.design.X(sorting_values,i);
%                     if strcmpi(LIMO.Analysis,'Time-Frequency')
%                         if strcmpi(g.restrict,'time')
%                             continuous(i,:,:) = Yr(channel,freq,:,sorting_values);
%                         else
%                             continuous(i,:,:) = Yr(channel,:,time,sorting_values);
%                         end
%                     else
%                         continuous(i,:,:) = Yr(channel,:,sorting_values);
%                     end
%                     clear mytitle
%                     if isempty(LIMO.design.electrode)
%                         mytitle{i} = sprintf('Original subjects'' parameters \n sorted by regressor %g channel %s (%g)', i, LIMO.data.chanlocs(channel).labels, channel);
%                     else
%                         mytitle{i} = sprintf('Original subjects'' parameters, \n sorted by regressor %g at optimized channel', i);
%                     end
%                 end
%                 remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
%                 continuous(remove,:,:) = [];
%                 reg_values(remove,:)   = [];
%             end
%
%         elseif contains(extra,{'Modelled','Modeled'},'Ignorecase',true)
%             if exist('Betas.mat','file') % OLS & IRLS GLM
%                 if strcmpi(LIMO.Analysis,'Time-Frequency')
%                     Betas = load('Betas.mat');
%                     if strcmpi(g.restrict,'time')
%                         Betas = squeeze(Betas.Betas(channel,freq,:,:));
%                     else
%                         Betas = squeeze(Betas.Betas(channel,:,time,:));
%                     end
%                     R     = eye(size(Yr,4)) - (LIMO.design.X*pinv(LIMO.design.X));
%                 else
%                     Betas = load('Betas.mat');
%                     Betas = squeeze(Betas.Betas(channel,:,:));
%                     R     = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
%                 end
%                 Yh = (LIMO.design.X*Betas')'; % modelled data
%             else % strcmpi(LIMO.design.method,'Generalized Welch's method')
%                 [~,~,Yh,~,dfe] = limo_robust_1way_anova(squeeze(Yr(channel,:,:)),LIMO.design.X);
%                 Res            = squeeze(Yr(channel,:,:))-Yh;
%             end
%
%             if regressor <= length(LIMO.design.nb_conditions) && ...
%                     LIMO.design.nb_conditions ~= 0 % for categorical variables
%                 if length(LIMO.design.nb_conditions) == 1
%                     start = 1;
%                 else
%                     start = sum(LIMO.design.nb_conditions(1:regressor-1));
%                 end
%
%                 for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
%                     index{i}     = find(LIMO.design.X(:,i));
%                     data         = squeeze(Yh(:,index{i}));
%                     average(i,:) = nanmean(data,2);
%                     index{i}     = index{i}(find(~isnan(squeeze(Yr(channel,1,index{i}))))); %#ok<FNDSB>
%                     if exist('R','var')
%                         var      = diag(((R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')) / LIMO.model.model_df(2));
%                     else
%                         var      = diag(Res*Res')./dfe;
%                     end
%                     CI           = sqrt(var/size(index{i},1))*z';
%                     ci(i,:,:)    = (repmat(nanmean(data,2),1,2)+CI)';
%                 end
%                 clear mytitle
%                 if isempty(LIMO.design.electrode)
%                     mytitle = sprintf('Modelled subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
%                 else
%                     mytitle = sprintf('Modelled subjects'' parameters at optimized channel');
%                 end
%                 clear Yr
%             else % continuous variable
%                 for i=max(regressor):-1:min(regressor)
%                     index{i}           = find(LIMO.design.X(:,i));
%                     [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
%                     reg_values(i,:)    = LIMO.design.X(sorting_values,i);
%                     continuous(i,:,:)  = Yh(:,sorting_values);
%                     clear mytitle
%                     if isempty(LIMO.design.electrode)
%                         mytitle = sprintf('Modelled subjects'' parameters \n sorted by regressor %g channel %s (%g)', ...
%                             i, LIMO.data.chanlocs(channel).labels, channel);
%                     else
%                         mytitle = sprintf('Modelled subjects'' parameters \n sorted by regressor %g at optimized channel', i);
%                     end
%                 end
%                 remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
%                 continuous(remove,:,:) = [];
%                 reg_values(remove,:)   = [];
%             end
%
%         elseif contains(extra,'Adjusted','Ignorecase',true)
%             if length(LIMO.design.nb_conditions) == 1 && LIMO.design.nb_continuous == 0
%                 warning on;
%                 if exist('warndlg2','file')
%                     warndlg2('Only one condition detected, no adjusted data possible');return
%                 else
%                     warndlg('Only one condition detected, no adjusted data possible');return
%                 end
%             end
%
%             allvar = 1:size(LIMO.design.X,2)-1;
%             allvar(regressor)=[]; % all but constant and columns of interest
%             if strcmpi(LIMO.Analysis,'Time-Frequency')
%                 Betas     = load('Betas.mat');
%                 if strcmpi(g.restrict,'time')
%                     Yr    = squeeze(Yr(channel,freq,:,:));
%                     Betas = squeeze(Betas.Betas(channel,freq,:,:));
%                 else
%                     Yr    = squeeze(Yr(channel,:,time,:));
%                     Betas = squeeze(Betas.Betas(channel,:,time,:));
%                 end
%                 confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
%             else
%                 Yr        = squeeze(Yr(channel,:,:));
%                 Betas     = load('Betas.mat');
%                 Betas     = squeeze(Betas.Betas(channel,:,:));
%                 confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
%             end
%             Ya = Yr - confounds;
%             clear Yr Betas confounds;
%
%             if regressor <= length(LIMO.design.nb_conditions) && ...
%                     LIMO.design.nb_conditions ~= 0 % for categorical variables
%                 if length(LIMO.design.nb_conditions) == 1
%                     start = 1;
%                 else
%                     start = sum(LIMO.design.nb_conditions(1:regressor-1));
%                 end
%
%                 for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
%                     index{i}     = find(LIMO.design.X(:,i));
%                     data         = squeeze(Ya(:,index{i})); % use adjusted data
%                     average(i,:) = nanmean(data,2);
%                     se           = nanstd(data,0,2) ./ sqrt(numel(index{i}));
%                     ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Ya,1));
%                 end
%                 clear mytitle
%                 if isempty(LIMO.design.electrode)
%                     mytitle = sprintf('Adjusted subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
%                 else
%                     mytitle = sprintf('Adjusted subjects'' parameters at  at optimized channel');
%                 end
%                 clear Yr
%
%             else % continuous variable ; regressor value already + sum(LIMO.design.nb_conditions)
%                 for i=max(regressor):-1:min(regressor)
%                     index{i}           = find(LIMO.design.X(:,i));
%                     [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
%                     reg_values(i,:)    = LIMO.design.X(sorting_values,i);
%                     continuous(i,:,:)  = Ya(:,sorting_values);
%                     clear mytitle
%                     if isempty(LIMO.design.electrode)
%                         mytitle = sprintf('Adjusted subjects'' parameters \n sorted by regressor %g channel %s (%g)', i, LIMO.data.chanlocs(channel).labels, channel);
%                     else
%                         mytitle = sprintf('Adjusted subjects'' parameters \n sorted by regressor %g at optimized channel', i);
%                     end
%                 end
%                 remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
%                 continuous(remove,:,:) = [];
%                 reg_values(remove,:)   = [];
%             end
%         else
%             error('unspecified data type to plot ''Original'',''Modelled'' or ''Adjusted'' ')
%         end
%
%         % make the figure(s)
%         % ------------------
%         figure;set(gcf,'Color','w')
%         if regressor <= length(LIMO.design.nb_conditions) && ...
%                 LIMO.design.nb_conditions ~= 0 % for categorical variables
%             for i=1:size(average,1)
%                 if strcmpi(LIMO.Analysis,'Time') || strcmpi(g.restrict,'time')
%                     timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
%                     plot(timevect,squeeze(average(i,:)),'LineWidth',1.5); hold on
%                     xlabel('Time in ms','FontSize',14)
%                     ylabel('Amplitude (A.U.)','FontSize',14)
%                 elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(g.restrict,'frequency')
%                     freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
%                     plot(freqvect,squeeze(average(i,:)),'LineWidth',1.5); hold on
%                     xlabel('Frequency in Hz','FontSize',14)
%                     ylabel('Spectral Power (A.U.)','FontSize',14)
%                 else
%                     error('couldn''t figure out what dimension to plot')
%                 end
%
%                 if i==1
%                     colorOrder = get(gca, 'ColorOrder');
%                     colorOrder = repmat(colorOrder,ceil(size(average,1)/size(colorOrder,1)),1);
%                 end
%                 x = squeeze(ci(i,1,:)); y = squeeze(ci(i,2,:));
%                 fillhandle = patch([timevect fliplr(timevect)], [x',fliplr(y')], colorOrder(i,:));
%                 set(fillhandle,'EdgeColor',colorOrder(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
%             end
%
%             h = axis;
%             if strcmpi(LIMO.Analysis,'Time')
%                 plot(timevect,(sig./10+1).*h(3),'r*','LineWidth',2)
%             else
%                 plot(freqvect,(sig./10+1).*h(3),'r*','LineWidth',2)
%             end
%
%             axis tight; grid on; box on
%             title(mytitle,'FontSize',16); drawnow;
%             assignin('base','Plotted_data', average)
%             v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
%             set(gca,'FontSize',14);
%             if strcmpi(LIMO.Analysis,'Time')
%                 ylabel('Amplitude (A.U.)','FontSize',16)
%                 xlabel('Time in ms','FontSize',16)
%             else
%                 ylabel('Spectral Power (A.U.)','FontSize',16)
%                 xlabel('Frequency in Hz','FontSize',16)
%             end
%
%         else % 3D plots
%             for i=1:size(continuous,1)
%                 if i > 1; figure;set(gcf,'Color','w'); end
%                 index = find(~isnan(squeeze(continuous(i,1,:))));
%                 if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency')
%                     if strcmpi(LIMO.Analysis,'Time')
%                         timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in msec
%                     else
%                         timevect = LIMO.data.tf_times;
%                     end
%                     surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
%                     ylabel('Time in ms','FontSize',16)
%                     zlabel('Amplitude (A.U.)','FontSize',16)
%                 else
%                     surf(index,LIMO.data.freqlist,squeeze(continuous(i,:,index)));shading interp
%                     ylabel('Frequency in Hz','FontSize',16)
%                     zlabel('Spectral Power (A.U.)','FontSize',16)
%                 end
%                 % --
%                 axis tight; title(mytitle,'FontSize',14); drawnow;
%                 xlabel('Sorted variable','FontSize',14)
%                 try
%                     set(gca,'XTick',index, 'XTickLabels', reg_values(index));
%                 catch label_err
%                     warning on; warning('could not set X-labels:\n%s',label_err)
%                 end
%             end
%         end
%
%
%     elseif contains(LIMO.design.name,'Repeated','IgnoreCase',true)   % All stuffs for repeated measures ANOVA
%         if contains(FileName,'LIMO')
%             error('Select summary stat file, nothing to infer from LIMO file')
%         end
%
%         % which summary stat
%         if ~isempty(g.sumstats) && any(strcmpi(g.sumstats,{'Mean','Trimmed'}))
%             extra = g.sumstats;
%         else
%             if contains(LIMO.design.name,'robust','Ignorecase',true)
%                 extra = 'Trimmed Mean';
%             else
%                 extra = 'Mean';
%             end
%             % let's not give GUI option and follow the design
%             % extra = questdlg('Summarize data using:','Data plot option','Mean','Trimmed Mean','Mean');
%             % if isempty(extra)
%             %     return
%             % end
%         end
%
%         if contains(FileName,'Rep_ANOVA_Main')
%
%             % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
%             if contains(FileName,'Main_effect','IgnoreCase',true)
%                 index1 = strfind(FileName,'Main_effect')+12;
%             elseif contains(FileName,'Interaction','IgnoreCase',true)
%                 index1 = strfind(FileName,'Interaction')+12;
%             end
%             index2                = max(strfind(FileName,'_'))-1;
%             effect_nb             = eval(FileName(index1:index2));
%             C                     = LIMO.design.C{effect_nb};
%             Data                  = load(fullfile(LIMO.dir,'Yr.mat'));
%             Data                  = Data.(cell2mat(fieldnames(Data)));
%             if strcmpi(LIMO.Analysis,'Time-Frequency')
%                 [~,channel,freq,time] = limo_display_reducedim(Data,LIMO,g.channels,g.restrict,g.dimvalue);
%                 Data              = squeeze(Data(channel,freq,time,:,:));
%                 sig               = squeeze(single(mask(channel,freq,time)));
%             else
%                 [~,channel,freq,time] = limo_display_reducedim(Data,LIMO,g.channels);
%                 Data              = squeeze(Data(channel,time,:,:)); % note freq/time variables have the same values
%                 sig               = single(mask(channel,:));
%             end
%             sig(sig==0)=NaN;
%
%             % compute differences between pairs using C and Cov
%             n = size(Data,2);
%             if strcmpi(extra,'Mean')
%                 for time_or_freq = size(Data,1):-1:1
%                     avg(time_or_freq,:) = nanmean(C*squeeze(Data(time_or_freq,:,:))',2);
%                     S(time_or_freq,:,:) = nancov(squeeze(Data(time_or_freq,:,:)));
%                 end
%                 if isempty(LIMO.design.electrode)
%                     mytitle = sprintf('Original %s \n %s %s (%g)',mytitle,LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel);
%                 else
%                     mytitle = sprintf('Original %s \n virtual %s',mytitle,LIMO.Type(1:end-1));
%                 end
%             else
%                 g=floor((20/100)*n); %% compute for 20% trimmed mean
%                 for time_or_freq = size(Data,1):-1:1
%                     [v,indices]          = sort(squeeze(Data(time_or_freq,:,:))); % sorted data
%                     TD(time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
%                     avg(time_or_freq,:)  = nanmean(C*squeeze(TD(time_or_freq,:,:))',2);
%                     v(1:g+1,:)           = repmat(v(g+1,:),g+1,1);
%                     v(n-g:end,:)         = repmat(v(n-g,:),g+1,1); % winsorized data
%                     [~,reorder]          = sort(indices);
%                     for j = size(Data,3):-1:1
%                         SD(:,j)          = v(reorder(:,j),j);
%                     end % restore the order of original data
%                     S(time_or_freq,:,:)  = cov(SD); % winsorized covariance
%                 end
%                 clear mytitle
%                 if isempty(LIMO.design.electrode)
%                     mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
%                 else
%                     mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
%                 end
%             end
%
%             % CI
%             dfe = size(Data,2)-size(Data,3)+1;
%             % c = avg + 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C'))); % uses Bonferoni inequality
%             % b = avg - 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C')));
%             bound = (abs(tinv(p./(2*size(C,1)),dfe)).*diag((sqrt(C*squeeze(S(time_or_freq,:,:))*C'))));
%             c = avg + repmat(bound', [length(avg),1]);
%             b = avg - repmat(bound', [length(avg),1]);
%
%             % do the figure
%             colours = limo_color_images(size(avg,2));
%             if strcmpi(LIMO.Analysis,'Time')
%                 if isfield(LIMO.data,'timevect')
%                     xvect = LIMO.data.timevect;
%                 else
%                     LIMO.data.timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
%                     save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.timevect;
%                 end
%             elseif strcmpi(LIMO.Analysis,'Frequency')
%                 if isfield(LIMO.data,'timevect')
%                     xvect = LIMO.data.freqlist;
%                 else
%                     LIMO.data.freqlist = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,2));
%                     save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.freqlist;
%                 end
%             elseif strcmpi(LIMO.Analysis,'Time-Frequency')
%                 if length(time) > 1
%                     if isfield(LIMO.data,'tf_times')
%                         xvect = LIMO.data.tf_times;
%                     else
%                         LIMO.data.tf_times = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
%                         save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.tf_times;
%                     end
%                 elseif length(freq) > 1
%                     if isfield(LIMO.data,'tf_freqs')
%                         xvect = LIMO.data.tf_freqs;
%                     else
%                         LIMO.data.tf_freqs = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,2));
%                         save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.tf_freqs;
%                     end
%                 end
%             end
%
%             figure;set(gcf,'Color','w')
%             for cond = 1:size(c,2)
%                 plot(xvect,avg(:,cond)','LineWidth',3,'color',colours(cond,:));
%                 fillhandle = patch([xvect fliplr(xvect)], [c(:,cond)',fliplr(b(:,cond)')], colours(cond,:));
%                 set(fillhandle,'EdgeColor',colours(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
%                 hold on
%             end
%             grid on; box on; axis tight; hold on;
%             h = axis;  plot(xvect,(sig./10+1).*h(3),'k.','MarkerSize',20)
%             if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
%                 xlabel('Time in ms','FontSize',14)
%                 ylabel('Amplitude (A.U.)','FontSize',14)
%             elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
%                 xlabel('Frequency in Hz','FontSize',14)
%                 ylabel('Spectral Power (A.U.)','FontSize',14)
%             end
%             set(gca,'FontSize',14,'layer','top');
%             title(mytitle,'FontSize',16); drawnow;
%             trimci = [c ; avg ; b];
%             assignin('base','Plotted_data',trimci);
%
%             % ----------------------
%         elseif contains(FileName,'Rep_ANOVA_Gp')  %% plot pairs of gp differences
%
%             % -------------------
%             % which ERP to make
%             % ------------------
%             extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Original');
%             if isempty(extra)
%                 return;
%             end
%             % -----------------------
%             Rep_ANOVA_Gp_effect = load(FileName);
%             Rep_ANOVA_Gp_effect = Rep_ANOVA_Gp_effect.(cell2mat(fieldnames(Rep_ANOVA_Gp_effect)));
%             if strcmpi(LIMO.Analysis,'Time-Frequency')
%                 [~,channel,freq,time] = limo_display_reducedim(Rep_ANOVA_Gp_effect,LIMO,g.channels,g.restrict,g.dimvalue);
%                 Rep_ANOVA_Gp_effect   = squeeze(Rep_ANOVA_Gp_effect(channel,freq,time,:,:));
%                 sig                   = squeeze(single(mask(channel,freq,time)));
%             else
%                 [~,channel,freq,time] = limo_display_reducedim(Rep_ANOVA_Gp_effect,LIMO,g.channels); %#ok<ASGLU>
%                 Rep_ANOVA_Gp_effect   = squeeze(Rep_ANOVA_Gp_effect(channel,time,:,:)); % note freq/time variables have the same values
%                 sig                   = single(mask(channel,:));
%             end
%
%             % check channel to plot
%             if size(Rep_ANOVA_Gp_effect,1) > 1
%                 if isempty(channel)
%                     [v,e] = max(Rep_ANOVA_Gp_effect(:,:,1));
%                     [~,c] = max(v); channel = e(c);
%                 else
%                     if ischar(channel)
%                         channel = str2double(channel);
%                     end
%
%                     if length(channel) > 1
%                         error('1 channel only can be plotted')
%                     elseif channel > size(Rep_ANOVA_Gp_effect,1)
%                         error('channel number invalid')
%                     end
%                 end
%             else
%                 if length(LIMO.design.electrode) == 1
%                     channel = LIMO.design.electrode;
%                 else
%                     channel = 1;  % accomodates the fact that all matrices have the channel dim (even = 1)
%                 end
%             end
%             clear Rep_ANOVA_Gp_effect
%
%             % compute the pair-wise differences and plot
%             Yr = load(fullfile(LIMO.dir,'Yr.mat'));
%             Yr = Yr.Yr;
%             if strcmpi(extra,'Original')
%                 Data = mean(squeeze(Yr(channel,:,:,:)),3);
%                 combinations = nchoosek(1:LIMO.design.nb_conditions,2);
%                 for d = size(combinations,1):-1:1
%                     Effect(:,d) = nanmean(Data(:,find(LIMO.data.Cat == combinations(d,1))),2) - nanmean(Data(:,find(LIMO.data.Cat == combinations(d,2))),2); %#ok<FNDSB>
%                 end
%             end
%
%             % design matrix
%             X = zeros(size(Yr,3),LIMO.design.nb_conditions+1);
%             X(:,end) = 1;
%             for i=1:LIMO.design.nb_conditions
%                 X(find(LIMO.data.Cat == i),i) = 1; %#ok<FNDSB>
%             end
%
%             % data again
%             Y = nanmean(squeeze(Yr(channel,:,:,:)),3);
%             X = X(find(~isnan(Y(1,:))),:); %#ok<FNDSB>
%             Y = Y(:,find(~isnan(Y(1,:))))'; %#ok<FNDSB>
%             if strcmpi(extra,'Modelled')
%                 beta = pinv(X)*Y; Yhat = X*beta;
%                 combinations = nchoosek(1:LIMO.design.nb_conditions,2);
%                 for d = 1:size(combinations,1)
%                     Effect(:,d) = nanmean(Yhat(find(X(:,combinations(d,1))),:),1)' - nanmean(Yhat(find(X(:,combinations(d,2))),:),1)'; %#ok<FNDSB>
%                 end
%             end
%
%             Res    = (Y'*(eye(size(Y,1)) - (X*pinv(X)))*Y);
%             df     = size(Y,1)-rank(X); t = tcdf(1-p,df);
%             sigma2 = sum((Res.^2./df),2);
%             v      = t.*sqrt(sigma2 ./ norm(X(:,1:end-1)).^2);
%             b      = Effect - v(channel,:);
%             c      = Effect + v(channel,:);
%             if channel == 1
%                 if length(LIMO.design.electrode) == 1
%                     if strcmpi(extra,'Original')
%                         mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
%                     else
%                         mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
%                     end
%                 else
%                     if strcmpi(extra,'Original')
%                         mytitle = sprintf('Mean parameter difference between groups \n optimized channel');
%                     else
%                         mytitle = sprintf('Modelled parameter difference between groups \n optimized channel');
%                     end
%                 end
%             else
%                 if strcmpi(extra,'Original')
%                     mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
%                 else
%                     mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
%                 end
%             end
%
%             figure;set(gcf,'Color','w'); hold on
%             RGB = limo_color_images(size(Effect,2));
%             if strcmpi(LIMO.Analysis,'Time')
%                 xvect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
%             else
%                 xvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
%             end
%
%             for d = size(combinations,1):-1:1
%                 plot(xvect,Effect(:,d),'LineWidth',3,'Color',RGB(d,:));
%                 fillhandle = patch([xvect fliplr(xvect)], [c(:,d)',fliplr(b(:,d)')], RGB(d,:));
%                 set(fillhandle,'EdgeColor',RGB(d,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
%                 Gp_difference{d} = [c(:,d)' ; Effect(:,d)' ; b(:,d)'];
%             end
%             grid on; box on; axis tight; hold on;
%             h = axis; sig(sig==0)=NaN;
%             plot(xvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
%             set(gca,'FontSize',14,'layer','top');
%             title(mytitle,'FontSize',16); drawnow;
%             assignin('base','Plotted_data',Gp_difference);
%
%             % -------------------------
%         elseif contains(FileName,'Rep_ANOVA_Interaction') % Gp * Repeated measures - plot differences btween condition per gp
%             % ------------------------
%
%             Rep_ANOVA_Interaction_with_gp = load(FileName);
%             Rep_ANOVA_Interaction_with_gp = Rep_ANOVA_Interaction_with_gp.(cell2mat(fieldnames(Rep_ANOVA_Interaction_with_gp)));
%             if strcmpi(LIMO.Analysis,'Time-Frequency')
%                 [~,channel,freq,time]         = limo_display_reducedim(Rep_ANOVA_Interaction_with_gp,LIMO,g.channels,g.restrict,g.dimvalue);
%                 Rep_ANOVA_Interaction_with_gp = squeeze(Rep_ANOVA_Interaction_with_gp(channel,freq,time,:,:));
%                 sig                           = squeeze(single(mask(channel,freq,time)));
%             else
%                 [~,channel,freq,time]         = limo_display_reducedim(Rep_ANOVA_Interaction_with_gp,LIMO,g.channels); %#ok<ASGLU>
%                 Rep_ANOVA_Interaction_with_gp = squeeze(Rep_ANOVA_Interaction_with_gp(channel,time,:,:)); % note freq/time variables have the same values
%                 sig                           = single(mask(channel,:));
%             end
%             sig(sig==0)=NaN;
%
%             if size(Rep_ANOVA_Interaction_with_gp,1) > 1
%                 if isempty(channel)
%                     [v,e] = max(Rep_ANOVA_Interaction_with_gp(:,:,1));
%                     [~,c] = max(v); channel = e(c);
%                 else
%                     if ischar(channel)
%                         channel = str2double(channel);
%                     end
%
%                     if length(channel) > 1
%                         error('1 channel only can be plotted')
%                     elseif channel > size(Rep_ANOVA_Interaction_with_gp,1)
%                         error('channel number invalid')
%                     end
%                 end
%             else
%                 channel = 1;
%             end
%
%             % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
%             effect_nb = str2double(FileName(strfind(FileName,'Factor_')+7:strfind(FileName,'Factor_')+6+strfind(FileName(strfind(FileName,'Factor_')+6:end),'_')));
%             C         = LIMO.design.C{effect_nb};
%             Yr        = load(fullfile(LIMO.dir,'Yr.mat'));
%             Yr        = Yr.Yr;
%             Data      = squeeze(Yr(channel,:,:,:));
%
%             % compute differences between pairs using C and Cov
%             if strcmpi(extra,'Mean')
%                 for gp = 1:LIMO.design.nb_conditions
%                     index = find(LIMO.data.Cat==gp);
%                     for time=size(Data,1):-1:1
%                         avg(gp,time,:) = nanmean(C*squeeze(Data(time,index,:))',2);
%                         S(gp,time,:,:) = cov(squeeze(Data(time,index,:)));
%                     end
%                 end
%                 if ~isempty(LIMO.design.electrode)
%                     mytitle = sprintf('Original %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
%                 else
%                     mytitle = sprintf('Original %s \n optimized channel',mytitle);
%                 end
%             else
%                 for gp = 1:LIMO.design.nb_conditions
%                     index = find(LIMO.data.Cat==gp);
%                     n     = length(index);
%                     g     = floor((20/100)*n); %% compute for 20% trimmed mean
%                     for time_or_freq=size(Data,1):-1:1
%                         [v,indices]             = sort(squeeze(Data(time_or_freq,index,:))); % sorted data
%                         TD(gp,time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
%                         avg(gp,time_or_freq)    = nanmean(C*squeeze(TD(gp,time_or_freq,:,:))',2);
%                         v(1:g+1,:)              = repmat(v(g+1,:),g+1,1);
%                         v(n-g:end,:)            = repmat(v(n-g,:),g+1,1); % winsorized data
%                         [~,reorder]             = sort(indices);
%                         for j = 1:size(Data,3)
%                             SD(:,j)             = v(reorder(:,j),j); % restore the order of original data
%                         end
%                         S(gp,time_or_freq,:,:)  = cov(SD); % winsorized covariance
%                     end
%                     clear SD
%                 end
%                 if ~isempty(LIMO.design.electrode)
%                     mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
%                 else
%                     mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
%                 end
%             end
%
%             figure; set(gcf,'Color','w'); hold on
%             if strcmpi(LIMO.Analysis,'Time')
%                 xvect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
%             else
%                 xvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
%             end
%
%             colorindex = 1;
%             trimci     = cell(size(avg,1)*size(avg,3),1);
%             RGB        = limo_color_images(size(avg,1)*size(avg,3));
%             for gp = 1:LIMO.design.nb_conditions
%                 % get the variance per comparison
%                 for frame = size(avg,2):-1:1
%                     varc(:,frame) = diag(sqrt(C*squeeze(S(gp,frame,:,:))*C'));
%                 end
%                 % plot each comparison
%                 for c = 1:size(avg,3)
%                     plot(xvect,squeeze(avg(gp,:,c)),'Color',RGB(colorindex,:),'LineWidth',3);
%                     % there is an error in the following 2 formulas I cannot fix -- @disbeat
%                     index      = find(LIMO.data.Cat==gp);
%                     dfe        = size(Data,2)-length(index)+1;
%                     up         = avg(gp,:,c) + tinv(p./(2*size(C,1)),dfe).* varc(c,:);
%                     down       = avg(gp,:,c) - tinv(p./(2*size(C,1)),dfe).* varc(c,:);
%                     fillhandle = patch([xvect fliplr(xvect)], [up,fliplr(down)], RGB(colorindex,:));
%                     set(fillhandle,'EdgeColor',RGB(colorindex,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
%                     trimci{colorindex} = [up ; squeeze(avg(gp,:,c)); down];
%                     colorindex = colorindex + 1;
%                 end
%             end
%
%             grid on; box on; axis tight
%             h = axis;  hold on;
%             plot(xvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
%             ylabel('Amplitude (A.U.)','FontSize',14)
%             if strcmpi(LIMO.Analysis,'Time')
%                 xlabel('Time in ms','FontSize',14)
%             else
%                 ylabel('Frequency in Hz','FontSize',14)
%             end
%             set(gca,'FontSize',14,'layer','top');
%             title(mytitle,'FontSize',16); drawnow;
%             assignin('base','Plotted_data',trimci);
%         end
%     else
%         errordlg('this file is not supported for this kind of plot','Nothing plotted')
%     end
% end % closes type

% elseif strcmpi(LIMO.Level,'LI')
%
%     [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,[],[]);
%
%     if Type == 1
%         %--------------------------
%         % imagesc of the results
%         %--------------------------
%         if sum(mask(:)) == 0
%             warndlg('no values under threshold parameter','no significant effect');
%         else
%             scale = M.*mask;
%             if min(scale(:))<0
%                 scale(scale==0)=min(scale(:))+(min(scale(:))/10);
%             else
%                 scale(scale==0)=NaN;
%             end
%
%             figure; set(gcf,'Color','w');
%             timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(M,2));
%             imagesc(timevect,1:size(M,1),scale);
%             title(mytitle,'FontSize',18);
%             color_images_(scale,LIMO);
%             assignin('base','Plotted_data',scale)
%             assignin('base','Mask_of_sig',mask)
%         end
%
%     elseif Type == 2
%         %--------------------------
%         % topoplot
%         %--------------------------
%         if sum(mask(:)) == 0
%             warndlg('no values under threshold','no significant effect');
%         else
%             EEG.data = M.*mask;
%             EEG.setname = 'Lateralization Map';
%             EEG.pnts = size(EEG.data,2);
%             EEG.xmin = LIMO.data.start/1000; % in sec
%             EEG.xmax = LIMO.data.end/1000;   % in sec
%             EEG.times =  (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
%             EEG.trials = 1;
%             EEG.chanlocs = LIMO.data.chanlocs;
%             EEG.nbchan = size(EEG.data,1);
%             pop_topoplot(EEG);
%             assignin('base','Plotted_data',EEG.data)
%         end
%
%     elseif Type == 3
%         disp('no ERP Plots to Lateralization'); % we could but nobody asked for it
%     end
% end

%%%%%%%%%%%% end of limo_display_results %%%%%%%%%%%%%%%%%
