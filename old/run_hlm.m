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
mkdir(outDir); cd(outDir)
eeglab; close;
pop_editoptions('option_single',0); % ensure double precision
pop_editoptions('option_parallel',1); % turn parrallel pool off
            
% add paths
% root = fileparts(which('limo_eeg'));
% addpath([root filesep 'limo_cluster_functions']);
% addpath([root filesep 'external' filesep 'psom']);
% addpath([root filesep 'external']);
% addpath([root filesep 'help']);
% ftPath = fileparts(which('ft_analysispipeline.m'));
% addpath(fullfile(ftPath, 'external','signal'))
% addpath(fullfile(ftPath, 'external','stats'))
% addpath(fullfile(ftPath, 'external','images'))

% % Create study and compute PSD
% commands = {};
% for iSub = 1:5
%     disp('--------------------------------------------')
%     disp(['Processing subject ' num2str(iSub) ])
%     disp('--------------------------------------------')
% 
%     EEG = pop_loadset('filename', sprintf('sub-%2.2d.set',iSub),'filepath',fullfile(dataDir,sprintf('sub-%2.2d',iSub)));
%     EEG.saved = 'no';
%     newpath = fullfile(outDir,EEG.subject); mkdir(newpath)
%     newname(iSub,:) = { sprintf('sub-%2.2d.set',iSub) };
%     pop_saveset(EEG,'filepath',newpath,'filename', newname{iSub});
%     commands = [ commands(:)' 'index' iSub 'load' fullfile(newpath, newname{iSub}) ];
%     [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','test','commands', commands, ...
%         'updatedat','on','savedat','off','rmclust','off');
%     [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% end
% [STUDY, EEG] = pop_savestudy(STUDY,EEG,'filename','test.study','filepath',outDir);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% STUDY = std_makedesign(STUDY,ALLEEG,1,'name','STUDY.design 1','delfiles','off', 'defaultdesign','off', ...
%     'variable1','type','values1',{'rest','trance'},'vartype1','categorical', ...
%     'subjselect',newname);
% [STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');
% [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on', ...
%     'rmicacomps','off', 'interp','off','recompute','on','spec','on', ...
%     'specparams',{'specmode','psd','logtrials','on','freqrange',[1 15]});
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% [STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');
% cd(outDir)
% gong

% Load study
[STUDY, ALLEEG] = pop_loadstudy('filename','test.study','filepath',outDir);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% STUDY = pop_statparams(STUDY,'condstats','on','method','perm','mcorrect','fdr','alpha',0.05);
% STUDY = pop_specparams(STUDY, 'plotconditions','together','freqrange',[1 15] ,'averagechan','on');
% [STUDY, specdata, specfreqs, pgroup, pcond] = std_specplot(STUDY,ALLEEG, ...
%     'channels',{ALLEEG(1).urchanlocs.labels},'design',1,'noplot','on');
% cond1 = cell2mat(specdata(1)); % condition 1
% cond2 = cell2mat(specdata(2)); % condition 2
% plotHDI(specfreqs',cond1,cond2,'Mean','dependent',.05,cell2mat(pcond)','rest','trance','Mean PSD + 95% HDI');

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

% LIMO 1st level
% pop_limo(STUDY,ALLEEG, ...
%     'method','WLS', ...         % GLM optimization: OLS, IRLS, WLS (default)
%     'measure','datspec', ...    % which EEG measure: daterp, datsec, datersp
%     'freqlim',[0 13] , ...      % freqs (e.g., [0 30]); or timelim for ERP
%     'erase','on', ...           % erase previous model
%     'splitreg','off', ...       % split regression for continuous indep. var.
%     'interaction','off');       % interaction model for categorical indep. var.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% std_limo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% model.set_files  = [];
% model.cat_files  = [];
% model.cont_files = [];

% Cleaning old files from the current design (Cleaning ALL)
% if strcmp(opt.erase,'on')
%     [~,filename] = fileparts(STUDY.filename);
%     std_limoerase(STUDY.filepath, filename, STUDY.subject, num2str(STUDY.currentdesign));
%     STUDY.limo = [];
% end

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
            %             single_trials_filename = EEGTMP.etc.datafiles.(opt.measureori);
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

% % Add contrasts for conditions that were merged during design selection
% % i.e. multiple categorical variables (factors) and yet not matching the number
% % of variables (contrasts are then a weighted sum of the crossed factors)
% if ~isempty(factors) && isfield(factors, 'value') && ...
%         sum(arrayfun(@(x) ~strcmpi(x.label,'group'), STUDY.design(opt.design).variable)) == 1 % only one non-continuous variable other than group
%     if length(STUDY.design(opt.design).variable(1).value) ~= length(factors) % and this var has more values than the number of factors
%         limocontrast = zeros(length(STUDY.design(opt.design).variable(1).value),length(factors)+1); % length(factors)+1 to add the constant
%         for n=length(factors):-1:1
%             factor_names{n} = factors(n).value;
%         end
% 
%         index = find(arrayfun(@(x) ~strcmpi(x.label,'group'),STUDY.design(opt.design).variable)); % which one is not group
%         for c=1:length(STUDY.design(opt.design).variable(index).value)
%             limocontrast(c,1:length(factors)) = single(ismember(factor_names,STUDY.design(opt.design).variable(index).value{c}));
%             limocontrast(c,1:length(factors)) = limocontrast(c,1:length(factors)) ./ sum(limocontrast(c,1:length(factors))); % scale by the number of variables
%         end
%     end
% end

% Transpose
model.set_files  = model.set_files';
model.cat_files  = model.cat_files';
model.cont_files = model.cont_files';
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% limo_batch %%%%%%%%%%%%%%%%%%%%%%%%%%% 

% [LIMO_files, procstatus] = limo_batch('model specification',model,[],STUDY);
% [LIMO_files, procstatus] = limo_batch(option,model,contrast)

% run several 1st level analyses select directories and files - possibly enter contrasts of
% interests and let it run. The batch relies on PSOM (see Ref). 
% see opt.mode for parallel computing on grid using qsub or msub

folder = fileparts(which('run_hlm.m'));
addpath(genpath(folder))

opt.mode                = 'session';    % run in the current session -- see psom for other options // in batch we use parfor
opt.max_queued          = Inf;          % with a maximum of possible sessions
opt.time_between_checks = 3;            % and x sec between job submission
opt.flag_pause          = false;        % don't bother asking to start jobs
opt.flag_debug          = false;        % report a bit more of issues
option = 'model specification';
% model = model;
% batch_contrast = varargin{3};

% psom_gb_vars                            % initialize PSOM variables

mkdir(['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6)]);
mkdir(['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6) filesep 'limo_batch_report']);
LIMO_files.LIMO = [pwd filesep ['derivatives' filesep 'LIMO_' STUDY.filename(1:end-6)]];

% Quick check
if ~isempty(model.cat_files)
    if size(model.cat_files,1) ~= size(model.set_files,1)
        error('the number of set and cat files disagree')
    end
end
if ~isempty(model.cont_files)
    if size(model.cont_files,1) ~= size(model.set_files,1)
        error('the number of set and cat files disagree')
    end
end

% Build pipelines
for subject = 1:size(model.set_files,1)

    % build LIMO.mat files from import
    command = 'limo_batch_import_data(files_in,opt.cat,opt.cont,opt.defaults)';
    pipeline(subject).import.command = command; 
    pipeline(subject).import.files_in = model.set_files{subject};
    pipeline(subject).import.opt.defaults = model.defaults;

%     if isfield(model.defaults,'type')
    pipeline(subject).import.opt.defaults.type = model.defaults.type;
%     else
%         pipeline(subject).import.opt.defaults.type = 'Channels';
%     end

%     if isfield(model.defaults,'method')
    pipeline(subject).import.opt.defaults.method = model.defaults.method;
%     else
%         pipeline(subject).import.opt.defaults.method = 'WLS';
%     end

%     if isfield(model.defaults,'type_of_analysis')
    pipeline(subject).import.opt.defaults.type_of_analysis = model.defaults.type_of_analysis;
%     else
%         pipeline(subject).import.opt.defaults.type_of_analysis = 'Mass-univariate';
%     end

%     if exist('STUDY','var')
%         if ~contains(STUDY.datasetinfo(subject).filename,{'sub-'}) && ...
%                 ~contains(STUDY.datasetinfo(subject).filename,{'_task-'}) % not bids
%             root = [fileparts(LIMO_files.LIMO) filesep 'sub-' STUDY.datasetinfo(subject).subject];
%         else
        subname = STUDY.datasetinfo(subject).subject;
        extra   = STUDY.datasetinfo(subject).filepath(strfind(STUDY.datasetinfo(subject).filepath,subname)+length(subname):end);
        root    = [fileparts(LIMO_files.LIMO) filesep subname extra]; % still in derivatives via LIMO_files.LIMO
%         end

        % if session and data are not in a derivatives/sess, make subdir
%         if ~isempty(STUDY.datasetinfo(subject).session)
            nsess = sum(strcmp(STUDY.datasetinfo(subject).subject,{STUDY.datasetinfo.subject}));
            if ~contains(root,'ses-') && nsess >= 1
%                 if ischar(STUDY.datasetinfo(subject).session)
%                     reuse = dir(fullfile(root,['ses-*' STUDY.datasetinfo(subject).session]));
%                     if ~isempty(reuse)
%                         index = find(arrayfun(@(x) STUDY.datasetinfo(subject).session == eval(x.name(5:end)), reuse));
%                         root = fullfile(reuse(index).folder,reuse(index).name);
%                     else
%                         root  = fullfile(root,['ses-' STUDY.datasetinfo(subject).session]);
%                     end
%                 else
                    reuse = dir(fullfile(root,['ses-*' num2str(STUDY.datasetinfo(subject).session)]));
                    if ~isempty(reuse)
                        index = find(arrayfun(@(x) STUDY.datasetinfo(subject).session == eval(x.name(5:end)), reuse));
                        root = fullfile(reuse(index).folder,reuse(index).name);
                    else
                        root = fullfile(root,['ses-' num2str(STUDY.datasetinfo(subject).session)]);
                    end
%                 end
            end
%         end

        % [root filesep eeg] - case of bids without ses-
%         if exist(fullfile(root,'eeg'),'dir')
%             root = fullfile(root,'eeg');
%         end

%         if exist(root,'dir') ~= 7
        mkdir(root);
%         end
        design_name = STUDY.design(STUDY.currentdesign).name;
        design_name(isspace(design_name)) = [];
        if strfind(design_name,'STUDY.') %#ok<STRIFCND>
            design_name = design_name(7:end);
        end
        glm_name = [STUDY.filename(1:end-6) '_' design_name '_GLM_' model.defaults.type '_' model.defaults.analysis '_' model.defaults.method];
        batch_contrast.LIMO_files{subject} = [root filesep glm_name filesep 'LIMO.mat'];
        % pipeline(subject).import.opt.defaults.studyinfo = STUDY.design_info;
%     else
%         [root,~,~] = fileparts(model.set_files{subject});
%         for l=min(length(LIMO_files.LIMO),length(root)):-1:1
%             common(l) = root(l) == LIMO_files.LIMO(l);
%         end
%         root = fullfile(LIMO_files.LIMO,root(min(find(diff(common))):end)); %#ok<MXFND>
%         glm_name = ['GLM_' model.defaults.method '_' model.defaults.analysis '_' model.defaults.type];
%     end
    pipeline(subject).import.files_out = [root filesep glm_name filesep 'LIMO.mat'];

%     if strcmp(option,'both') && ~isfield(batch_contrast,'LIMO_files')
%         batch_contrast.LIMO_files{subject} = [root filesep glm_name filesep 'LIMO.mat'];
%         batch_contrast.LIMO_files = batch_contrast.LIMO_files';
%     end
% 
%     if ~isempty(model.cat_files)
    pipeline(subject).import.opt.cat = model.cat_files{subject};
%     else
%         pipeline(subject).import.opt.cat = [];
%     end
    if ~isempty(model.cont_files)
        pipeline(subject).import.opt.cont = model.cont_files{subject};
    else
        pipeline(subject).import.opt.cont = [];
    end
    pipeline(subject).import.opt.defaults.name = fileparts(pipeline(subject).import.files_out);
    LIMO_files.mat{subject}  = [root filesep glm_name filesep 'LIMO.mat'];
    LIMO_files.Beta{subject} = [root filesep glm_name filesep 'Betas.mat'];

    % make design and evaluate
    command = 'limo_batch_design_matrix(files_in)';
    pipeline(subject).design.command = command;
    pipeline(subject).design.files_in = pipeline(subject).import.files_out;
    pipeline(subject).design.files_out = [root filesep glm_name filesep 'Yr.mat'];

    % run GLM
    command = 'limo_eeg(4,files_in)';
    pipeline(subject).glm.command = command;
    pipeline(subject).glm.files_in = pipeline(subject).import.files_out;
    pipeline(subject).glm.files_out = [root filesep glm_name filesep 'Betas.mat'];
end

% if strcmp(option,'contrast only') || strcmp(option,'both')
%     
%     if ~exist('model','var')
%         model.defaults.bootstrap = 0;
%         model.defaults.tfce      = 0;
%     end
%     
%     for subject = 1:length(batch_contrast.LIMO_files)
%         command = 'limo_batch_contrast(files_in,opt.C)';
%         pipeline(subject).n_contrast.command = command;
%         pipeline(subject).n_contrast.files_in = batch_contrast.LIMO_files{subject};
%         if iscell(batch_contrast.mat)
%             pipeline(subject).n_contrast.opt.C = cell2mat(batch_contrast.mat);
%         else
%             pipeline(subject).n_contrast.opt.C = batch_contrast.mat;
%         end
%         
%         if exist(batch_contrast.LIMO_files{subject},'file')
%             sub_LIMO = load(batch_contrast.LIMO_files{subject});
%             if ~isfield(sub_LIMO.LIMO,'contrast')
%                 start = 0;
%             else
%                 start = length(sub_LIMO.LIMO.contrast);
%             end
%         else
%             start = 0;
%         end
%         
%         for c=1:size(batch_contrast.mat,1)
%             name{c} = [fileparts(batch_contrast.LIMO_files{subject}) filesep 'con_' num2str(c+start) '.mat'];
%         end
%         pipeline(subject).n_contrast.files_out = name; % name{1};
%         LIMO_files.con{subject} = name;
%     end
% end


% if strcmp(option,'model specification') || strcmp(option,'both')
N = size(model.set_files,1);
LIMO_files.mat = LIMO_files.mat';
LIMO_files.Beta = LIMO_files.Beta';
remove_limo = zeros(1,N);
% else
%     N = length(batch_contrast.LIMO_files);
% end
procstatus = zeros(1,N);

% if isfield(LIMO_files,'con')
%     LIMO_files.con = LIMO_files.con';
%     remove_con     = zeros(1,N);
% else
remove_con = 0;
% end

% allocate names
for subject = 1:N
    limopt{subject} = opt;
    limopt{subject}.path_logs = [LIMO_files.LIMO filesep 'limo_batch_report' filesep glm_name filesep 'subject' num2str(subject)];
end

% limo_settings_script
limo_settings.psom = 1;

% if model.defaults.bootstrap ~= 0 || ~limo_settings.psom % debugging mode, serial analysis    
%     for subject = 1:N 
%         disp('--------------------------------')
%         fprintf('processing model %g/%g \n',subject,N)
%         disp('--------------------------------')
%         
%         psom_pipeline_debug(pipeline(subject));
%         if strcmp(option,'contrast only')
%             name = fileparts(batch_contrast.LIMO_files{subject}); %#ok<PFBNS,PFTUSW>
%         else
%             [~,name]=fileparts(model.set_files{subject}); %#ok<PFBNS>
%         end
%         sub = min(strfind(name,'sub-'));
%         ses = min(strfind(name,'ses-'));
%         und = strfind(name,'_');
%         if ~isempty(sub) && ~isempty(ses) && ~isempty(und)
%             try
%                 sub_und = und(und>sub); ses_und = und(und>ses);
%                 report{subject} = ['subject ' name(sub+4:sub+min(abs(sub_und-sub))-1) ' session ' name(ses+4:ses+min(abs(ses_und-ses))-1) ' processed'];
%             catch
%                 report{subject} = ['subject ' num2str(subject) ' processed'];
%             end
%         else
%             report{subject} = ['subject ' num2str(subject) ' processed'];
%         end
%         procstatus(subject) = 1;
%     end
%     
% else % parallel call to the pipeline
%     limo_check_ppool

%% PSOM 
for subject = 1:N 
    disp('--------------------------------')
    fprintf('processing model %g/%g \n',subject,N)
    disp('--------------------------------')
        
%         try
            %%%%%%%%%%%%%%%%%% psom_run_pipeline %%%%%%%%%%%%%%%%%%%%%%%%%
%             psom_run_pipeline(pipeline(subject),limopt{subject})
        psompipeline = pipeline(subject);
        psomopt = limopt{subject};
        psom_gb_vars
        name_pipeline = 'PIPE';
        gb_name_structure = 'psomopt';
        gb_list_fields    = {'flag_short_job_names' , 'nb_resub'       , 'type_restart' , 'flag_pause' , 'init_matlab'       , 'flag_update' , 'flag_debug' , 'path_search'       , 'restart' , 'shell_options'       , 'path_logs' , 'command_matlab' , 'flag_verbose' , 'mode'       , 'mode_pipeline_manager' , 'max_queued'       , 'qsub_options'       , 'time_between_checks' , 'nb_checks_per_point' , 'time_cool_down' };
        gb_list_defaults  = {true                   , gb_psom_nb_resub , 'substring'    , true         , gb_psom_init_matlab , true          , false        , gb_psom_path_search , {}        , gb_psom_shell_options , NaN         , ''               , true           , gb_psom_mode , gb_psom_mode_pm         , gb_psom_max_queued , gb_psom_qsub_options , []                    , []                    , []               };
        psom_set_defaults

        if ~strcmp(psomopt.path_logs(end),filesep)
            psomopt.path_logs = [psomopt.path_logs filesep];
            path_logs = psomopt.path_logs;
        end

        if isempty(path_search)
            path_search = path;
            psomopt.path_search = path_search;
        end

        if isempty(psomopt.command_matlab)
            if strcmp(gb_psom_language,'matlab')
                psomopt.command_matlab = gb_psom_command_matlab;
            else
                psomopt.command_matlab = gb_psom_command_octave;
            end
        end

        if strcmp(psomopt.mode,'session')
            psomopt.max_queued = 1;
            max_queued = 1;
        end

%             if max_queued == 0
%                 switch psomopt.mode
%                     case {'batch','background'}
%                         if isempty(gb_psom_max_queued)
%                             psomopt.max_queued = 1;
%                             max_queued = 1;
%                         else
%                             psomopt.max_queued = gb_psom_max_queued;
%                             max_queued = gb_psom_max_queued;
%                         end
%                     case {'session','qsub','msub','condor'}
%                         if isempty(gb_psom_max_queued)
%                             psomopt.max_queued = Inf;
%                             max_queued = Inf;
%                         else
%                             psomopt.max_queued = gb_psom_max_queued;
%                             max_queued = gb_psom_max_queued;
%                         end
%                 end % switch action
%             end % default of max_queued

            if ~ismember(psomopt.mode,{'session','background','batch','qsub','msub','condor'})
                error('%s is an unknown mode of pipeline execution. Sorry dude, I must quit ...',psomopt.mode);
            end

%             switch psomopt.mode
%                 case 'session'
                if isempty(time_between_checks)
                    time_between_checks = 0;
                end
                if isempty(nb_checks_per_point)
                    nb_checks_per_point = Inf;
                end
%                 otherwise
%                     if isempty(time_between_checks)
%                         time_between_checks = 0;
%                     end
%                     if isempty(nb_checks_per_point)
%                         nb_checks_per_point = 60;
%                     end
%             end

            file_pipe_running = cat(2,path_logs,filesep,name_pipeline,'.lock');
            file_logs = cat(2,path_logs,filesep,name_pipeline,'_history.txt');

            % Initialize the logs folder
            opt_init.path_logs      = psomopt.path_logs;
            opt_init.path_search    = psomopt.path_search;
            opt_init.command_matlab = psomopt.command_matlab;
            opt_init.flag_verbose   = psomopt.flag_verbose;
            opt_init.restart        = psomopt.restart;
            opt_init.flag_update    = psomopt.flag_update;    
            opt_init.flag_pause     = psomopt.flag_pause;
            opt_init.type_restart   = psomopt.type_restart;
            
            if flag_debug
                opt_init
            end
            
            [tmp,flag_start] = psom_pipeline_init(psompipeline,opt_init);   
            if ~flag_start, return; end
            
            % Run the pipeline manager
            file_pipeline = cat(2,path_logs,filesep,name_pipeline,'.mat');
            
            opt_proc.mode                  = psomopt.mode;
            opt_proc.mode_pipeline_manager = psomopt.mode_pipeline_manager;
            opt_proc.max_queued            = psomopt.max_queued;
            opt_proc.qsub_options          = psomopt.qsub_options;
            opt_proc.shell_options         = shell_options;
            opt_proc.command_matlab        = psomopt.command_matlab;
            opt_proc.time_between_checks   = psomopt.time_between_checks;
            opt_proc.nb_checks_per_point   = psomopt.nb_checks_per_point;
            opt_proc.flag_short_job_names  = psomopt.flag_short_job_names;
            opt_proc.flag_debug            = psomopt.flag_debug;
            opt_proc.flag_verbose          = psomopt.flag_verbose;
            opt_proc.init_matlab           = psomopt.init_matlab;
            opt_proc.nb_resub              = psomopt.nb_resub;
            
            if flag_debug, opt_proc; end
            
            % Read the number of characters that are currently in the history
%             if flag_verbose&&~strcmp(psomopt.mode_pipeline_manager,'session')
%                 hf = fopen(file_logs,'r');
%                 if hf~=-1
%                     str_logs = fread(hf,Inf,'uint8=>char')';
%                     nb_chars = ftell(hf);
%                     fclose(hf);
%                 else
%                     nb_chars = 0;
%                 end
%             end
            
%%%%%%%%%%%%%%%%%%%%%%% psom_pipeline_process %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             psom_pipeline_process(file_pipeline,opt_proc);
        psompipeline2 = file_pipeline;
        opt2 = opt_proc;

        gb_name_structure = 'psomopt2';
        gb_list_fields    = { 'flag_short_job_names' , 'nb_resub'       , 'flag_verbose' , 'init_matlab'       , 'flag_debug' , 'shell_options'       , 'command_matlab' , 'mode'    , 'mode_pipeline_manager' , 'max_queued' , 'qsub_options'       , 'time_between_checks' , 'nb_checks_per_point' , 'time_cool_down' };
        gb_list_defaults  = { true                   , gb_psom_nb_resub , true           , gb_psom_init_matlab , true         , gb_psom_shell_options , ''               , 'session' , ''                      , 0            , gb_psom_qsub_options , []                    , []                    , []               };
        psom_set_defaults

        flag_verbose = flag_verbose || flag_debug;

        if isempty(opt2.nb_resub)
            switch opt2.mode
                case {'session','batch','background'}
                    opt2.nb_resub = 0;
                    nb_resub = 0;
                otherwise
                    opt2.nb_resub = 1;
                    nb_resub = 1;
            end % switch action
        end % default of max_queued

        if isempty(time_between_checks)
            opt2.time_between_checks = 0;
            time_between_checks = 0;
        end
        if isempty(nb_checks_per_point)
            opt2.nb_checks_per_point = Inf;
            nb_checks_per_point = Inf;
        end
        if isempty(time_cool_down)
            opt2.time_cool_down = 0;
            time_cool_down = 0;
        end

        % generic messages
        hat_qsub_o = sprintf('\n\n*****************\nOUTPUT QSUB\n*****************\n');
        hat_qsub_e = sprintf('\n\n*****************\nERROR QSUB\n*****************\n');

        % generating file names
        [path_logs,name_pipeline,ext_pl] = fileparts(psompipeline2);
        file_pipe_running   = [ path_logs filesep name_pipeline '.lock'               ];
        file_pipe_log       = [ path_logs filesep name_pipeline '_history.txt'        ];
        file_news_feed      = [ path_logs filesep name_pipeline '_news_feed.csv'        ];
        file_manager_opt    = [ path_logs filesep name_pipeline '_manager_opt.mat'    ];
        file_logs           = [ path_logs filesep name_pipeline '_logs.mat'           ];
        file_logs_backup    = [ path_logs filesep name_pipeline '_logs_backup.mat'    ];
        file_status         = [ path_logs filesep name_pipeline '_status.mat'         ];
        file_status_backup  = [ path_logs filesep name_pipeline '_status_backup.mat'  ];
        file_jobs           = [ path_logs filesep name_pipeline '_jobs.mat'           ];
        file_profile        = [ path_logs filesep name_pipeline '_profile.mat'        ];
        file_profile_backup = [ path_logs filesep name_pipeline '_profile_status.mat' ];

        logs    = load( file_logs    );
        status  = load( file_status  );
        profile = load( file_profile );

        path_tmp = [path_logs filesep 'tmp'];
        if exist(path_tmp,'dir')
            delete([path_tmp '*']);
        else
            mkdir(path_tmp);
        end

        % Create a running tag on the pipeline (if not done during the initialization phase)
        str_now = datestr(clock);
        save(file_pipe_running,'str_now');

        % open the log file
        hfpl = fopen(file_pipe_log,'a');

        % Open the news feed file
        hfnf = fopen(file_news_feed,'w');

        % Print general info about the pipeline
%         msg_line1 = sprintf('The pipeline %s is now being processed.',name_pipeline);
%         msg_line2 = sprintf('Started on %s',datestr(clock));
%         msg_line3 = sprintf('user: %s, host: %s, system: %s',gb_psom_user,gb_psom_localhost,gb_psom_OS);
%         size_msg = max([size(msg_line1,2),size(msg_line2,2),size(msg_line3,2)]);
%         msg = sprintf('%s\n%s\n%s',msg_line1,msg_line2,msg_line3);
%         stars = repmat('*',[1 size_msg]);    
%         sub_add_line_log(hfpl,sprintf('\n%s\n%s\n%s\n',stars,msg,stars),flag_verbose);
    
        % Load the pipeline
        load(psompipeline2,'list_jobs','graph_deps','files_in');                
        
        % Update dependencies
        mask_finished = false([length(list_jobs) 1]);
        for num_j = 1:length(list_jobs)
            mask_finished(num_j) = strcmp(status.(list_jobs{num_j}),'finished');
        end
        graph_deps(mask_finished,:) = 0;
        mask_deps = max(graph_deps,[],1)>0;
        mask_deps = mask_deps(:);
        
        % Track number of submissions
        nb_sub = zeros([length(list_jobs) 1]);
    
        % Initialize the to-do list
        mask_todo = false([length(list_jobs) 1]);
        for num_j = 1:length(list_jobs)
            mask_todo(num_j) = strcmp(status.(list_jobs{num_j}),'none');
            if ~mask_todo(num_j)
                sub_add_line_log(hfnf,sprintf('%s , finished\n',list_jobs{num_j}),false);
            end
        end    
        mask_done = ~mask_todo;
    
        mask_failed = false([length(list_jobs) 1]);
        for num_j = 1:length(list_jobs)
            mask_failed(num_j) = strcmp(status.(list_jobs{num_j}),'failed');
            if mask_failed(num_j)
                sub_add_line_log(hfnf,sprintf('%s , failed\n',list_jobs{num_j}),false);
            end
        end    
        list_num_failed = find(mask_failed);
        list_num_failed = list_num_failed(:)';
        for num_j = list_num_failed
            mask_child = false([1 length(mask_todo)]);
            mask_child(num_j) = true;
            mask_child = sub_find_children(mask_child,graph_deps);
            mask_todo(mask_child) = false; % Remove the children of the failed job from the to-do list
        end    
        mask_running = false(size(mask_done));
    
        % Initialize miscallenaous variables
        nb_queued   = 0;                   % Number of queued jobs
        nb_todo     = sum(mask_todo);      % Number of to-do jobs
        nb_finished = sum(mask_finished);  % Number of finished jobs
        nb_failed   = sum(mask_failed);    % Number of failed jobs
        nb_checks   = 0;                   % Number of checks to print a points
        nb_points   = 0;                   % Number of printed points
    
        lmax = 0;
        for num_j = 1:length(list_jobs)
            lmax = max(lmax,length(list_jobs{num_j}));
        end   
    
    % GLM HAPPENS HERE
    while (any(mask_todo) || any(mask_running)) && psom_exist(file_pipe_running)

        % Update logs & status
        save(file_logs           ,'-struct','logs');
        copyfile(file_logs,file_logs_backup,'f');        
        save(file_status         ,'-struct','status');
        copyfile(file_status,file_status_backup,'f');
        save(file_profile        ,'-struct','profile');
        copyfile(file_profile,file_profile_backup,'f');
        flag_nothing_happened = true;
        
        % Update the status of running jobs
        list_num_running = find(mask_running);
        list_num_running = list_num_running(:)';
        list_jobs_running = list_jobs(list_num_running);
        new_status_running_jobs = psom_job_status(path_logs,list_jobs_running,opt2.mode);        
        pause(time_cool_down); % pause for a while to let the system finish to write eqsub and oqsub files (useful in 'qsub' mode).
        
        % Loop over running jobs to check the new status
        num_l = 0;
        for num_j = list_num_running
            num_l = num_l+1;
            name_job = list_jobs{num_j};
            flag_changed = ~strcmp(status.(name_job),new_status_running_jobs{num_l});
            
            if flag_changed
                
                if flag_nothing_happened % if nothing happened before...
                    % Reset the 'dot counter'
                    flag_nothing_happened = false;
                    nb_checks = 0;
                    if nb_points>0
                        sub_add_line_log(hfpl,sprintf('\n'),flag_verbose);
                    end
                    nb_points = 0;
                end
                
                % update status of the job in the status file                
                status.(name_job) = new_status_running_jobs{num_l};
                
                if strcmp(status.(name_job),'exit') % the script crashed ('exit' tag)
                    sub_add_line_log(hfpl,sprintf('%s - The script of job %s terminated without generating any tag, I guess we will count that one as failed.\n',datestr(clock),name_job),flag_verbose);
                    status.(name_job) = 'failed';
                    nb_failed = nb_failed + 1;
                end
                
                if strcmp(status.(name_job),'failed')||strcmp(status.(name_job),'finished')
                    % for finished or failed jobs, transfer the individual
                    % test log files to the matlab global logs structure
                    nb_queued = nb_queued - 1;
                    text_log    = sub_read_txt([path_logs filesep name_job '.log']);
                    text_qsub_o = sub_read_txt([path_logs filesep name_job '.oqsub']);
                    text_qsub_e = sub_read_txt([path_logs filesep name_job '.eqsub']);                    
                    if isempty(text_qsub_o)&&isempty(text_qsub_e)
                        logs.(name_job) = text_log;                        
                    else
                        logs.(name_job) = [text_log hat_qsub_o text_qsub_o hat_qsub_e text_qsub_e];
                    end
                    % Update profile for the jobs
                    file_profile_job = [path_logs filesep name_job '.profile.mat'];
                    if psom_exist(file_profile_job)
                        profile.(name_job) = load(file_profile_job);
                    end
                    profile.(name_job).nb_submit = nb_sub(num_j);
                    sub_clean_job(path_logs,name_job); % clean up all tags & log                    
                end
                
                switch status.(name_job)
                    
                    case 'failed' % the job has failed, too bad !

                        if nb_sub(num_j) > nb_resub % Enough attempts to submit the jobs have been made, it failed !
                            nb_failed = nb_failed + 1;   
                            msg = sprintf('%s %s%s failed   ',datestr(clock),name_job,repmat(' ',[1 lmax-length(name_job)]));
                            sub_add_line_log(hfpl,sprintf('%s (%i run / %i fail / %i done / %i left)\n',msg,nb_queued,nb_failed,nb_finished,nb_todo),flag_verbose);
                            sub_add_line_log(hfnf,sprintf('%s , failed\n',name_job),false);
                            mask_child = false([1 length(mask_todo)]);
                            mask_child(num_j) = true;
                            mask_child = sub_find_children(mask_child,graph_deps);
                            mask_todo(mask_child) = false; % Remove the children of the failed job from the to-do list
                        else % Try to submit the job one more time (at least)
                            mask_todo(num_j) = true;
                            status.(name_job) = 'none';
                            new_status_running_jobs{num_l} = 'none';
                            nb_todo = nb_todo+1;
                            msg = sprintf('%s %s%s reset    ',datestr(clock),name_job,repmat(' ',[1 lmax-length(name_job)]));
                            sub_add_line_log(hfpl,sprintf('%s (%i run / %i fail / %i done / %i left)\n',msg,nb_queued,nb_failed,nb_finished,nb_todo),flag_verbose);
                        end

                    case 'finished'

                        nb_finished = nb_finished + 1;                        
                        msg = sprintf('%s %s%s completed',datestr(clock),name_job,repmat(' ',[1 lmax-length(name_job)]));
                        sub_add_line_log(hfpl,sprintf('%s (%i run / %i fail / %i done / %i left)\n',msg,nb_queued,nb_failed,nb_finished,nb_todo),flag_verbose);
                        sub_add_line_log(hfnf,sprintf('%s , finished\n',name_job),false);
                        graph_deps(num_j,:) = 0; % update dependencies

                end
                
            end % if flag changed
        end % loop over running jobs
        
        if ~flag_nothing_happened % if something happened ...
            
            % update the to-do list
            mask_done(mask_running) = ismember(new_status_running_jobs,{'finished','failed','exit'});
            mask_todo(mask_running) = mask_todo(mask_running)&~mask_done(mask_running);
            
            % Update the dependency mask
            mask_deps = max(graph_deps,[],1)>0;
            mask_deps = mask_deps(:);
            
            % Finally update the list of currently running jobs
            mask_running(mask_running) = mask_running(mask_running)&~mask_done(mask_running);
            
        end
        
        % Time to (try to) submit jobs !!
        list_num_to_run = find(mask_todo&~mask_deps);
        num_jr = 1;
        
        while (nb_queued < max_queued) && (num_jr <= length(list_num_to_run))
            
            if flag_nothing_happened % if nothing happened before...
                % Reset the 'dot counter'
                flag_nothing_happened = false;
                nb_checks = 0;
                if nb_points>0                    
                    sub_add_line_log(hfpl,sprintf('\n'),flag_verbose);
                end
                nb_points = 0;
            end
            
            % Pick up a job to run
            num_job = list_num_to_run(num_jr);
            num_jr = num_jr + 1;
            name_job = list_jobs{num_job};
            
            mask_todo(num_job) = false;
            mask_running(num_job) = true;
            nb_queued = nb_queued + 1;
            nb_todo = nb_todo - 1;
            nb_sub(num_job) = nb_sub(num_job)+1;
            status.(name_job) = 'submitted';
            msg = sprintf('%s %s%s submitted',datestr(clock),name_job,repmat(' ',[1 lmax-length(name_job)]));            
            sub_add_line_log(hfpl,sprintf('%s (%i run / %i fail / %i done / %i left)\n',msg,nb_queued,nb_failed,nb_finished,nb_todo),flag_verbose);
            sub_add_line_log(hfnf,sprintf('%s , submitted\n',name_job),false);
            
            % Execute the job in a "shelled" environment
            file_job        = [path_logs filesep name_job '.mat'];
            opt_logs.txt    = [path_logs filesep name_job '.log'];
            opt_logs.eqsub  = [path_logs filesep name_job '.eqsub'];
            opt_logs.oqsub  = [path_logs filesep name_job '.oqsub'];
            opt_logs.failed = [path_logs filesep name_job '.failed'];
            opt_logs.exit   = [path_logs filesep name_job '.exit'];
            opt_script.path_search    = psompipeline2;
            opt_script.name_job       = name_job;
            opt_script.mode           = opt2.mode;
            opt_script.init_matlab    = opt2.init_matlab;
            opt_script.flag_debug     = opt2.flag_debug;        
            opt_script.shell_options  = opt2.shell_options;
            opt_script.command_matlab = opt2.command_matlab;
            opt_script.qsub_options   = opt2.qsub_options;
            opt_script.flag_short_job_names = opt2.flag_short_job_names;
            opt_script.file_handle    = hfpl;
            cmd = sprintf('psom_run_job(''%s'')',file_job);
                
            if ispc % this is windows
                script = [path_tmp filesep name_job '.bat'];
            else
                script = [path_tmp filesep name_job '.sh'];
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%% psom_run_script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             [flag_failed,errmsg] = psom_run_script(cmd,script,opt_script,opt_logs);
%             [flag_failed,msg] = psom_run_script(cmd,script,opt,logs)
            cmd3 = cmd;
            script3 = script;
            opt3 = opt_script;
            logs3 = opt_logs;

            list_fields    = { 'flag_short_job_names' , 'path_search'       , 'file_handle' , 'name_job'    , 'init_matlab'       , 'flag_debug' , 'shell_options'       , 'command_matlab' , 'mode' , 'qsub_options'       };
            list_defaults  = { true                   , gb_psom_path_search , []            , 'psom_script' , gb_psom_init_matlab , false        , gb_psom_shell_options , ''               , NaN    , gb_psom_qsub_options };
            opt3 = psom_struct_defaults(opt3,list_fields,list_defaults);
            
%             if opt3.flag_debug
%                 msg = sprintf('\n    The execution mode is %s\n',opt3.mode);
%                 fprintf(msg);
%                 if ~isempty(opt3.file_handle)
%                     fprintf(opt3.file_handle,'%s',msg);
%                 end
%             end
            
%             if isempty(opt3.path_search)
%                 opt3.path_search = path;
%             end
            
%             if ~isempty(opt3.init_matlab)&&~ismember(opt3.init_matlab(end),{',',';'})
%                 opt3.init_matlab = [opt3.init_matlab ','];
%             end

%             if isempty(opt3.command_matlab)
%                 if strcmp(gb_psom_language,'matlab')
%                     opt3.command_matlab = gb_psom_command_matlab;
%                 else
%                     opt3.command_matlab = gb_psom_command_octave;
%                 end
%             end

%             % Logs
%             if nargin < 4
%                 opt_logs = [];
%             else
%                 list_fields   = { 'txt' , 'eqsub' , 'oqsub' , 'failed' , 'exit' };
%                 if ismember(psomopt3.mode,{'qsub','msub','condor'})
%                     list_defaults = { NaN   , NaN     , NaN     , NaN    , ''     };
%                 else
%                     list_defaults = { NaN   , ''      , ''      , ''     , ''     };
%                 end
%                 opt_logs = psom_struct_defaults(opt_logs,list_fields,list_defaults);
%             end

% % Check that the execution mode exists
% if ~ismember(opt3.mode,{'session','background','batch','qsub','msub','condor'})
%     error('%s is an unknown mode of command execution. Sorry dude, I must quit ...',opt3.mode);
% end

% Generate the script
% Generate some OS-appropriate options to start Matlab/Octave
% switch gb_psom_language
%     case 'matlab'
%         if ispc
            opt_matlab = '-automation -nodesktop -r';
%         else
%             opt_matlab = '-nosplash -nodesktop -r';
%         end        
%     case 'octave'
%         opt_matlab = '--silent --eval';       
% end

% Set-up the search path for the job
% if ~strcmp(opt3.mode,'session')&&~isempty(cmd)
%     if (length(opt3.path_search)>4)&&(strcmp(opt3.path_search(end-3:end),'.mat'))
%         file_path = opt3.path_search;
%     else
%         [path_f,name_f] = fileparts(script);
%         file_path = fullfile(path_f,[name_f '_path.mat']);
%         path_work = opt3.path_search;
%         save(file_path,'path_work');
%     end 
%     opt3.init_matlab = [sprintf('load(''%s'',''path_work''), if ~ismember(path_work,{''gb_niak_omitted'',''gb_psom_omitted''}), path(path_work), end,',file_path) opt3.init_matlab];
% else
    file_path = '';
% end
        
% Add an appropriate call to Matlab/Octave
% if ~isempty(cmd)            
    instr_job = sprintf('"%s" %s "%s %s,exit"',opt3.command_matlab,opt_matlab,opt3.init_matlab,cmd);
%     if ~isempty(opt_logs)
%         if opt3.flag_debug
%             instr_job = sprintf('%s >"%s" 2>&1\n',instr_job,opt_logs.txt);
%         else
            instr_job = sprintf('%s >"%s"\n',instr_job,opt_logs.txt);
%         end
%     else
%         instr_job = sprintf('%s\n',instr_job);
%     end
% else
%     instr_job = '';
% end
        
% Add shell options
% if ~isempty(opt3.shell_options)
%     instr_job = sprintf('%s\n%s',opt3.shell_options,instr_job);
% end    

% Add a .exit tag file
% if ~isempty(opt_logs)&&~isempty(opt_logs.exit)
%     if ispc % this is windows
        instr_job = sprintf('%stype nul > "%s"\nexit\n',instr_job,opt_logs.exit);
%     else
%         instr_job = sprintf('%stouch "%s"',instr_job,opt_logs.exit);
%     end
% else
%     if ispc
%         instr_job = sprintf('%sexit\n',instr_job);
%     end
% end

% Write the script
% if ~strcmp(opt3.mode,'session')            
%     if opt3.flag_debug
%         msg = sprintf('    This is the content of the script used to run the command :\n"\n%s\n"\n',instr_job);
%         if ~strcmp(opt3.path_search,'gb_niak_omitted')&&~isempty(cmd)
%             msg = sprintf('%s    The following matlab search path is used (may be truncated):\n%s\n',msg,opt3.path_search(1:min(100,length(opt3.path_search))));
%             msg = sprintf('%s    The search path will be loaded from the following file:\n%s\n',msg,file_path);
%         end
%         fprintf('%s',msg);
%         if ~isempty(opt3.file_handle)
%             fprintf(opt3.file_handle,'%s',msg);
%         end
%     end
%     
%     hf = fopen(script,'w');
%     fprintf(hf,'%s',instr_job);
%     fclose(hf);
% else
%     if opt3.flag_debug
%         msg = sprintf('    The following command is going to be executed :\n%s\n\n',cmd);
%         fprintf('%s',msg);
%         if ~isempty(opt3.file_handle)
%             fprintf(opt3.file_handle,'%s',msg);
%         end
%     end
% end

% Execute the script 
% switch opt3.mode
% 
%     case 'session'

        try
            if ~isempty(opt_logs)
                diary(opt_logs.txt);
                eval([ cmd ';' ])  %%%%%%%%%% ====> HAPPENS HERE %%%%%%%
                %%%%%%%%%%%%%%%%%% psom_run_job %%%%%%%%%%%%%%%%%%%%%%%%
%                 psom_run_job('G:\My Drive\HLM\limo\derivatives\LIMO_test\limo_batch_report\test_design1_GLM_Channels_Frequency_WLS\subject1\\import.mat')
                
%          1) IMPORT line 128 with sub_eval --> limo_batch_import_data(set_files,opt.cat,opt.cont,opt.defaults)
%           find what is called for 2) design and 3) GLM (same line in psom_run_job)


                %%%%%%%%%%%%%%%%% end of psom_run_job %%%%%%%%%%%%%%%%%%
                diary off;
            else
                sub_eval([ cmd ';' ]);
            end
            flag_failed = false;
            msg = '';
        catch
            flag_failed = true;
            errmsg = lasterror;
            msg = errmsg.message;
            if isfield(errmsg,'stack')
                for num_e = 1:length(errmsg.stack)
                    msg = sprintf('%s\nFile %s at line %i\n',msg,errmsg.stack(num_e).file,errmsg.stack(num_e).line);
                end           
            end
        end
        if ~isempty(opt_logs.exit)
            save(opt_logs.exit,'flag_failed')
        end

%     case 'background'
% 
%        if ispc
%             cmd_script = ['"' script '"']; % /min instead of /b ?
%        else
%             cmd_script = ['. "' script '"'];
%        end
% 
%        if opt3.flag_debug
%            if strcmp(gb_psom_language,'octave')
%                cmd_script = [cmd_script ' 2>&1']; % In octave, the error stream is lost. Redirect it to standard output
%            end
%            msg = sprintf('    The script is executed using the command :\n%s\n\n',cmd_script);
%            fprintf('%s',msg);
%            if ~isempty(opt3.file_handle)
%                fprintf(opt3.file_handle,'%s',msg);
%            end
%            [flag_failed,msg] = system(cmd_script);
%        else
%            if strcmp(gb_psom_language,'octave')
%                system(cmd_script,false,'async');
%                flag_failed = 0;
%            else 
%                flag_failed = system([cmd_script ' &']);
%            end
%            msg = '';
%        end
% 
%     case 'batch'
% 
%         if ispc
%             instr_batch = sprintf('start /min "%s" "%s"',opt3.name_job,script); 
%         else
%             instr_batch = ['at -f "' script '" now'];
%         end
%         if strcmp(gb_psom_language,'octave')
%             instr_batch = [instr_batch ' 2>&1']; % In octave, the error stream is lost. Redirect it to standard output
%         end
%         if opt3.flag_debug 
%             msg = sprintf('    The script is executed using the command :\n%s\n\n',instr_batch);
%             fprintf('%s',msg);
%             if ~isempty(opt3.file_handle)
%                 fprintf(opt3.file_handle,'%s',msg);
%             end
%         end
%         [flag_failed,msg] = system(instr_batch);    
%         
%     case {'qsub','msub','condor'}
%         script_submit = [gb_psom_path_psom 'psom_submit.sh'];
%         switch opt3.mode
%             case {'qsub','msub'}
%                 sub = opt3.mode;
%             case 'condor'
%                 sub = [gb_psom_path_psom 'psom_condor.sh'];
%         end
%         if ~isempty(opt_logs)
%             qsub_logs = [' -e \"' opt_logs.eqsub '\" -o \"' opt_logs.oqsub '\"'];
%         else 
%             qsub_logs = '';
%         end
%         if opt3.flag_short_job_names
%             name_job = opt3.name_job(1:min(length(opt3.name_job),8));
%         else
%             name_job = opt3.name_job;
%         end
%         instr_qsub = sprintf('%s%s -N %s %s %s',sub,qsub_logs,name_job,opt3.qsub_options,['\"' script '\"']);            
%         if ~isempty(opt_logs)
%             instr_qsub = [script_submit ' "' instr_qsub '" ' opt_logs.failed ' ' opt_logs.exit ' ' opt_logs.oqsub ];
%         end
%         if opt3.flag_debug
%             if strcmp(gb_psom_language,'octave')
%                 instr_qsub = [instr_qsub ' 2>&1']; % In octave, the error stream is lost. Redirect it to standard output
%             end
%             msg = sprintf('    The script is executed using the command :\n%s\n\n',instr_qsub);
%             fprintf('%s',msg);
%             if ~isempty(opt3.file_handle)
%                 fprintf(opt3.file_handle,'%s',msg);
%             end
%             [flag_failed,msg] = system(instr_qsub);
%         else 
%             if strcmp(gb_psom_language,'octave')
%                 system([instr_qsub ' > /dev/null'],false,'async');
%                 flag_failed = 0;
%             else
%                 flag_failed = system([instr_qsub ' > /dev/null &']);
%             end
%             msg = '';
%         end
% end





            %%%%%%%%%%%%%%%%%%%%%%% end of psom_run_script %%%%%%%%%%%%%%%%%%%%%%%

            if flag_failed~=0
                msg = fprintf('\n    The execution of the job %s failed.\n The feedback was : %s\n',name_job,errmsg);
                sub_add_line_log(hfpl,msg,true);
                error('Something went bad with the execution of the job.')
            elseif flag_debug
                msg = fprintf('\n    The feedback from the execution of job %s was : %s\n',name_job,errmsg);
                sub_add_line_log(hfpl,msg,true);
            end            
        end % submit jobs
        
        if (any(mask_todo) || any(mask_running)) && psom_exist(file_pipe_running)
            pause(time_between_checks); % To avoid wasting resources, wait a bit before re-trying to submit jobs
        end
                
        if nb_checks >= nb_checks_per_point
            nb_checks = 0;
            if flag_verbose
                fprintf('.');
            end
            sub_add_line_log(hfpl,sprintf('.'),flag_verbose);
            nb_points = nb_points+1;
        else
            nb_checks = nb_checks+1;
        end
        
    end % While there are jobs to do
    
% catch
%     
%     errmsg = lasterror;        
%     sub_add_line_log(hfpl,sprintf('\n\n******************\nSomething went bad ... the pipeline has FAILED !\nThe last error message occured was :\n%s\n',errmsg.message),flag_verbose);
%     if isfield(errmsg,'stack')
%         for num_e = 1:length(errmsg.stack)
%             sub_add_line_log(hfpl,sprintf('File %s at line %i\n',errmsg.stack(num_e).file,errmsg.stack(num_e).line),flag_verbose);
%         end
%     end
%     if exist('file_pipe_running','var')
%         if exist(file_pipe_running,'file')
%             delete(file_pipe_running); % remove the 'running' tag
%         end
%     end
%     
%     %% Close the log file
%     if strcmp(gb_psom_language,'matlab')
%         fclose(hfpl);
%         fclose(hfnf);
%     end
%     return
% end

Update the final status
save(file_logs           ,'-struct','logs');
copyfile(file_logs,file_logs_backup,'f');
save(file_status         ,'-struct','status');
copyfile(file_status,file_status_backup,'f');
save(file_profile        ,'-struct','profile');
copyfile(file_profile,file_profile_backup,'f');

Print general info about the pipeline
msg_line1 = sprintf('The processing of the pipeline is terminated.');
msg_line2 = sprintf('See report below for job completion status.');
msg_line3 = sprintf('%s',datestr(now));
size_msg = max([size(msg_line1,2),size(msg_line2,2)]);
msg = sprintf('%s\n%s\n%s',msg_line1,msg_line2,msg_line3);
stars = repmat('*',[1 size_msg]);
sub_add_line_log(hfpl,sprintf('\n%s\n%s\n%s\n',stars,msg,stars),flag_verbose);

Report if the lock file was manually removed
if exist('file_pipe_running','var')
    if ~exist(file_pipe_running,'file')
        sub_add_line_log(hfpl,sprintf('The pipeline manager was interrupted because the .lock file was manually deleted.\n'),flag_verbose);
    end    
end

% Print a list of failed jobs
mask_failed = false([length(list_jobs) 1]);
for num_j = 1:length(list_jobs)
    mask_failed(num_j) = strcmp(status.(list_jobs{num_j}),'failed');
end
mask_todo = false([length(list_jobs) 1]);
for num_j = 1:length(list_jobs)
    mask_todo(num_j) = strcmp(status.(list_jobs{num_j}),'none');
end
list_num_failed = find(mask_failed);
list_num_failed = list_num_failed(:)';
list_num_none = find(mask_todo);
list_num_none = list_num_none(:)';
flag_any_fail = ~isempty(list_num_failed);

if flag_any_fail
    if length(list_num_failed) == 1
        sub_add_line_log(hfpl,sprintf('The execution of the following job has failed :\n\n    '),flag_verbose);
    else
        sub_add_line_log(hfpl,sprintf('The execution of the following jobs have failed :\n\n    '),flag_verbose);
    end
    for num_j = list_num_failed
        name_job = list_jobs{num_j};        
        sub_add_line_log(hfpl,sprintf('%s ; ',name_job),flag_verbose);
    end    
    sub_add_line_log(hfpl,sprintf('\n\n'),flag_verbose);
    sub_add_line_log(hfpl,sprintf('More infos can be found in the individual log files. Use the following command to display these logs :\n\n    psom_pipeline_visu(''%s'',''log'',JOB_NAME)\n\n',path_logs),flag_verbose);
end

Print a list of jobs that could not be processed
if ~isempty(list_num_none)
    if length(list_num_none) == 1
        sub_add_line_log(hfpl,sprintf('The following job has not been processed due to a dependence on a failed job or the interruption of the pipeline manager :\n\n    '),flag_verbose);
    else
        sub_add_line_log(hfpl,sprintf('The following jobs have not been processed due to a dependence on a failed job or the interruption of the pipeline manager :\n\n    '),flag_verbose);
    end
    for num_j = list_num_none
        name_job = list_jobs{num_j};
        sub_add_line_log(hfpl,sprintf('%s ; ',name_job),flag_verbose);
    end    
    sub_add_line_log(hfpl,sprintf('\n\n'),flag_verbose);
end

% Give a final one-line summary of the processing
if flag_any_fail    
    sub_add_line_log(hfpl,sprintf('All jobs have been processed, but some jobs have failed.\nYou may want to restart the pipeline latter if you managed to fix the problems.\n'),flag_verbose);
else
    if isempty(list_num_none)
        sub_add_line_log(hfpl,sprintf('All jobs have been successfully completed.\n'),flag_verbose);
    end
end

if ~strcmp(opt2.mode_pipeline_manager,'session')&& strcmp(gb_psom_language,'octave')   
    sub_add_line_log(hfpl,sprintf('Press CTRL-C to go back to Octave.\n'),flag_verbose);
end

% Close the log file
if strcmp(gb_psom_language,'matlab')
    fclose(hfpl);
    fclose(hfnf);
end

if exist('file_pipe_running','var')
    if exist(file_pipe_running,'file')
        delete(file_pipe_running); % remove the 'running' tag
    end
end

if strcmp(opt2.mode,'session')&&strcmp(opt2.mode_pipeline_manager,'session')&&flag_any_fail
    error('All jobs have been processed, but some jobs have failed. You may want to restart the pipeline latter if you managed to fix the problems.')
end





%%%%%%%%%%%%%%%%%%% end of psom_pipeline_process %%%%%%%%%%%%%%%%

            %% If not in session mode, monitor the output of the pipeline
            if flag_verbose&&~strcmp(opt2.mode_pipeline_manager,'session')
                psom_pipeline_visu(path_logs,'monitor',nb_chars);
            end            
%%%%%%%%%%%%%%%%%%% end of psom_run_pipeline %%%%%%%%%%%%%%%%


            % example of debugging
            % ---------------------
            % psom reported with function failed, eg limo_batch_import
            % pipeline(subject).import tells you the command line to test
            % put the point brack where needed and call e.g.
            % limo_batch_import_data(pipeline(subject).import.files_in,pipeline(subject).import.psomopt2.cat,pipeline(subject).import.psomopt2.cont,pipeline(subject).import.psomopt2.defaults)
            % limo_batch_design_matrix(pipeline(subject).design.files_in)
            % limo_eeg(4,fileparts(pipeline(subject).glm.files_in))
            % limo_batch_contrast(pipeline(subject).n_contrast.files_in,pipeline(subject).n_contrast.psomopt2.C)
            
            if strcmp(option,'contrast only')
                name = fileparts(batch_contrast.LIMO_files{subject}); %#ok<PFBNS,PFTUSW>
            else
                [~,name]=fileparts(model.set_files{subject}); %#ok<PFBNS>
            end
            sub = min(strfind(name,'sub-'));
            ses = min(strfind(name,'ses-'));
            und = strfind(name,'_');
            
            if ~isempty(sub) && ~isempty(ses) && ~isempty(und)
                try
                    sub_und = und(und>sub); ses_und = und(und>ses);
                    if strcmp(option,'contrast only')
                        report{subject} = ['subject ' name(sub:sub+min(abs(sub_und-sub))-1) ' processed'];
                    else
                        report{subject} = ['subject ' name(sub+4:sub+min(abs(sub_und-sub))-1) ' session ' name(ses+4:ses+min(abs(ses_und-ses))-1) ' processed'];
                    end
                catch
                    report{subject} = ['subject ' num2str(subject) ' processed'];
                end
            else
                report{subject} = ['subject ' num2str(subject) ' processed'];
            end
            procstatus(subject) = 1;
%         catch ME
%             report{subject} = sprintf('subject %g failed: %s',subject,ME.message');
%             if strcmp(option,'model specification')
%                 remove_limo(subject) = 1;
%             elseif strcmp(option,'both')
%                 remove_limo(subject) = 1;
%                 remove_con(subject) = 1;
%             elseif strcmp(option,'contrast only')
%                 remove_con(subject) = 1;
%             end
%         end
%     end
% end

% Save txt files
% save as txt file the list of .set, Betas, LIMO and con
% these lists can then be used in second level analyses
cd(LIMO_files.LIMO)

if strcmp(option,'model specification') || strcmp(option,'both')
    if ~all(remove_limo)
        cell2csv([LIMO_files.LIMO filesep 'LIMO_files_' glm_name '.txt'], LIMO_files.mat(find(~remove_limo),:))
        cell2csv([LIMO_files.LIMO filesep 'Beta_files_' glm_name '.txt'], LIMO_files.Beta(find(~remove_limo),:))
    end
end

if strcmp(option,'contrast only') || strcmp(option,'both')
    for c=1:size(batch_contrast.mat,1)
        index = 1; clear name
        for subject = 1:N
            if strcmp(option,'contrast only')
                LIMO = load([fileparts(pipeline(subject).n_contrast.files_in) filesep 'LIMO.mat']); LIMO = LIMO.LIMO;
                if isfield(LIMO,'contrast')
                    con_num = max(find(cellfun(@(x) isequal(x.C,limo_contrast_checking(LIMO.dir,LIMO.design.X,batch_contrast.mat(c,:))),LIMO.contrast))); % if several identical contrasts, take max
                else
                    con_num = c;
                end
                name{index} = [fileparts(pipeline(subject).n_contrast.files_in) filesep 'con_' num2str(con_num) '.mat'];
            else
                name{index} = [fileparts(pipeline(subject).glm.files_out) filesep 'con_' num2str(c) '.mat'];
                con_num = c;
            end
            index = index + 1;
        end
        name = name';
        
        if ~all(remove_con)
            cell2csv([LIMO_files.LIMO filesep 'con_' num2str(con_num) '_files_' glm_name '.txt'], name(find(~remove_con),:));
        end
    end
end

% save the report from psom
cell2csv([LIMO_files.LIMO filesep 'limo_batch_report' filesep 'batch_report_' glm_name '.txt'], report')

cd(current);
failed = zeros(1,N);
for subject=1:N
    if strfind(report{subject},'failed')
        failed(subject) = 1;
    end
end

if sum(failed) == 0
    disp('LIMO batch processing finished succesfully')
else
    if sum(failed) == N % all subjects
        warning('LIMO batch done but all subjects failed')
    else
        warning('LIMO batch done, some errors where detected\ncheck limo batch report subjects %s',num2str(find(failed)))
    end
end

% if EEGLAB STUDY check for groups and sessions 
% and further export txt files
if exist('STUDY','var')
    try
        if isfield(model, 'set_files')
            cell2csv([LIMO_files.LIMO filesep 'EEGLAB_set_' glm_name '.txt'],model.set_files)
        end
        
        if ~isempty(STUDY.datasetinfo(subject).session)
            sesvalues = unique(arrayfun(@(x) x.session, STUDY.datasetinfo));
        else
            sesvalues = 1;
        end
        
        % split txt files if more than 1 group or session
        if length(STUDY.group) > 1 || length(sesvalues)>1
            for s=1:length(sesvalues)
                for g= 1:length(STUDY.group)
                    if length(STUDY.group) > 1
                        subset = arrayfun(@(x)(strcmpi(x.group,STUDY.group{g})), STUDY.datasetinfo);
                    end
                    
                    if length(sesvalues) > 1
                        sesset = arrayfun(@(x) x.session==s, STUDY.datasetinfo);
                    end
                    
                    if isfield(LIMO_files,'mat') && isfield(LIMO_files,'Beta')
                        if length(STUDY.group) > 1 && length(sesvalues)==1 % only groups
                            if any(subset)
                                cell2csv(fullfile(LIMO_files.LIMO, ['LIMO_files_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.mat(subset));
                                cell2csv(fullfile(LIMO_files.LIMO, ['Beta_files_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.Beta(subset));
                            end
                        elseif length(STUDY.group) == 1 && length(sesvalues) > 1 % only sessions
                            if any(sesset)
                                cell2csv(fullfile(LIMO_files.LIMO, ['LIMO_files_ses-' num2str(s) '_' glm_name '.txt']), LIMO_files.mat(sesset));
                                cell2csv(fullfile(LIMO_files.LIMO, ['Beta_files_ses-' num2str(s) '_' glm_name '.txt']), LIMO_files.Beta(sesset));
                            end
                        else % groups and sessions
                            if any(subset.*sesset)
                                cell2csv(fullfile(LIMO_files.LIMO, ['LIMO_files_ses-' num2str(s) '_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.mat(logical(subset.*sesset)));
                                cell2csv(fullfile(LIMO_files.LIMO, ['Beta_files_ses-' num2str(s) '_Gp-' STUDY.group{g} '_' glm_name '.txt']), LIMO_files.Beta(logical(subset.*sesset)));
                            end
                        end
                    end
                    
                    if isfield(LIMO_files,'con')
                        if length(STUDY.group) > 1 && length(sesvalues)==1 % only groups
                            tmpcell = LIMO_files.con(subset);
                            if ~isempty(tmpcell{1})
                                for c=1:length(tmpcell{1})
                                    [~,con_name,~] = fileparts(LIMO_files.con{1}{c});
                                    cell2csv(fullfile(LIMO_files.LIMO, [con_name '_files_Gp-' STUDY.group{g} '_' glm_name '.txt']),cellfun(@(x) x(c), tmpcell));
                                end
                            end
                        elseif length(STUDY.group) == 1 && length(sesvalues) > 1 % only sessions
                            tmpcell = LIMO_files.con(sesset);
                            if ~isempty(tmpcell{1})
                                for c=1:length(tmpcell{1})
                                    [~,con_name,~] = fileparts(LIMO_files.con{1}{c});
                                    cell2csv(fullfile(LIMO_files.LIMO, [con_name '_files_ses-' num2str(s) '_' glm_name '.txt']),cellfun(@(x) x(c), tmpcell));
                                end
                            end
                        else
                            tmpcell = LIMO_files.con(logical(subset.*sesset));
                            if ~isempty(tmpcell)
                                for c=1:length(tmpcell{1})
                                    [~,con_name,~] = fileparts(LIMO_files.con{1}{c});
                                    cell2csv(fullfile(LIMO_files.LIMO, [con_name '_files_ses-' num2str(s) '_Gp-' STUDY.group{g} '_' glm_name '.txt']),cellfun(@(x) x(c), tmpcell));
                                end
                            end
                        end
                    end
                end
            end
        end
    catch writtingerr
        if sum(failed) == 0
            warning(writtingerr.identifier,'all LIMO files created but failing to write some metadata txt files ''%s''\n ',writtingerr.message);
        else
            warning(writtingerr.identifier,'also failing to write some metadata txt files ''%s''\n ',writtingerr.message);
        end
    end
end

disp('LIMO batch works thanks to PSOM by Bellec et al. (2012)')
disp('The Pipeline System for Octave and Matlab. Front. Neuroinform. 6:7')


%%%%%%%%%%%%%%%%%%%%%%%%%%% end of limo_batch %%%%%%%%%%%%%%%%%%%%%%5
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
            if size(data,3)-1 <= Nmin
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

    mask = correct_cluster(M, Pval, bootM, bootP, neighbormatrix, MCC, p);

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
