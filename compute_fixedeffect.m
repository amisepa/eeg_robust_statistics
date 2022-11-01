%% Analyze data using hierarchichal linear modeling (HLM)
% 
% Cedric Cannard, Oct 2022

clear; close all; clc;
dataDir = 'C:\Users\IONSLAB\Desktop\channeling_matlab\data\data_processed2';
studypath = 'G:\My Drive\HLM\limo';
addpath('G:\My Drive\HLM')
cd(dataDir)
tmpSub = dir;
idx = contains({tmpSub.name},'sub');
tmpSub = {tmpSub(idx).name}';
mkdir(studypath); cd(studypath)
eeglab; close;
pop_editoptions('option_single',0); % double precision on (1) or off (0)
pop_editoptions('option_parallel',0); %  parallel computing on (1) or off (0)

% % Create study
% commands = {};
% for iSub = 1:5
%     disp('--------------------------------------------')
%     disp(['Subject ' num2str(iSub) ])
%     disp('--------------------------------------------')
% 
%     EEG = pop_loadset('filename', sprintf('sub-%2.2d.set',iSub),'filepath',fullfile(dataDir,sprintf('sub-%2.2d',iSub)));
%     EEG.saved = 'no';
%     EEG.session = 1;
%     allpaths(iSub,:) = { fullfile(studypath,EEG.subject) }; 
%     mkdir(allpaths{iSub})
%     allsub(iSub,:) = { sprintf('sub-%2.2d_ses-%d',iSub,EEG.session) };
%     pop_saveset(EEG,'filepath',allpaths{iSub},'filename', [allsub{iSub} '.set']);
%     commands = [ commands(:)' 'index' iSub 'load' fullfile(allpaths{iSub}, [allsub{iSub} '.set']) ];
%     [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','test','commands', commands, ...
%         'updatedat','on','savedat','off','rmclust','off');
%     [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% end
% [STUDY, EEG] = pop_savestudy(STUDY,EEG,'filename','test.study','filepath',studypath);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% % Make study design
% conditions = unique({ALLEEG(1).event.type});
% STUDY = std_makedesign(STUDY,ALLEEG,1,'name','STUDY.design 1','delfiles','off', 'defaultdesign','off', ...
%     'variable1','type','values1',conditions,'vartype1','categorical', ...
%     'subjselect',allsub);
% [STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% % Compute PSD
% [STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on', ...
%     'rmicacomps','off', 'interp','off','recompute','on','spec','on', ...
%     'specparams',{'specmode','psd','logtrials','on','freqrange',[1 15]});
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
% [STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');
% % cd(studypath)

% Load study
[STUDY, ALLEEG] = pop_loadstudy('filename','test.study','filepath',studypath);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% % Basic stats (FDR-correction, permutations, 95% HDI)
% STUDY = pop_statparams(STUDY, 'condstats','on','method','perm','mcorrect','fdr','alpha',0.05);
% STUDY = pop_specparams(STUDY, 'plotconditions','together','freqrange',[1 10] ,'averagechan','on');
% for iChan = 1:ALLEEG(1).nbchan
%     disp(num2str(iChan))
%     [STUDY, specdata, specfreqs, pgroup(iChan,:,:), pcond(iChan,:,:)] = std_specplot(STUDY,ALLEEG, ...
%         'channels',{ALLEEG(1).chanlocs(iChan).labels},'design',1,'noplot','on');   
%     cond1(iChan,:,:) = specdata{1}; % condition 1
%     cond2(iChan,:,:) = specdata{2}; % condition 2
% end
% chantoplot = 30;
% plotHDI(specfreqs', squeeze(cond1(chantoplot,:,:)), squeeze(cond2(chantoplot,:,:)), ...
%     'Trimmed mean', 'dependent', .05, pcond{chantoplot}', 'rest', 'trance', 'Trimmed mean & 95% HDI');

%% 1st level: Fixed effects

folder = fileparts(which('run_hlm.m'));
addpath(genpath(folder))

freqlim = [1 15];                                   % default = [], [1 30]
model.measure = 'datspec';                          % 'daterp', 'datspec', 'dattimef', 'icaerp', 'icaspec', 'icatimef'
model.method = 'WLS';                               % 'OLS', 'IRLS', 'WLS' (default)
model.datatype = 'Channels';
model.defaults.method = 'WLS';              
model.defaults.analysistype = 'Mass-univariate';    % 'Mass-univariate' or 'multivariate'
model.defaults.level = 1;                           % 1st level analysis
model.defaults.fullfactorial = 0;                   % all variables
model.defaults.zscore = 1;                          % 0 or 1 (default)
model.splitreg = 'off';
model.interaction = 'off';
model.design = STUDY.currentdesign;
model.defaults.type = 'Channels';
% model.defaults.datatype = model.measure(4:end);

% Parameters
% design = 1;                     % STUDY design to analyze
% method = 'WLS';                 % 'OLS', 'IRLS', 'WLS' (default)
% measure = 'datspec';            % 'daterp', 'datspec', 'dattimef', 'icaerp', 'icaspec', 'icatimef'
% datatype = 'Channels';          % 'Channels' or 'Components'
% erase = 'on';                   % default = 1
% level = 1;                      % subject (1) or group level (2)
% splitreg = 'off';               % default = off
% interaction = 'off';             % default = off
% timelim = [];                   % default = [], [-150 1000]
% opt.chanloc = struct('no', {});     % default = struct('no', {})
% zscoreDat = 1;                     % 0 or 1 (default)
% opt.ow_chanlocfile = 'no';          % 'yes' or 'no' (default)
% opt.measureori = opt.measure;
% Analysis = measure;
% design_index = design;

if model.defaults.level == 1
    model.defaults.bootstrap        = 0 ;                   % only for single subject analyses
    model.defaults.tfce             = 0;                    % only for single subject analyses
end
if strcmpi(model.measure, 'datspec')
    model.defaults.analysis = 'Frequency';
    model.defaults.lowf    = freqlim(1);
    model.defaults.highf   = freqlim(2);
    model.defaults.start = [];
    model.defaults.end   = [];
end

% Get STUDY info for LIMO
% measureflags = struct('daterp','off',...
%     'datspec','off', ...
%     'datersp','off',...
%     'dattimef','off',...
%     'datitc' ,'off',...
%     'icaerp' ,'off',...
%     'icaspec','off',...
%     'icatimef','off',...
%     'icaersp','off',...
%     'icaitc','off');
% measureflags.(lower(measure)) = 'on';
% STUDY.etc.measureflags = measureflags;
% STUDY.etc.datspec = 'on';
% mergedChanlocs = ALLEEG(1).chanlocs;

% Generate temporary merged datasets needed by LIMO
allSub = { STUDY.datasetinfo.subject };
allSess = { STUDY.datasetinfo.session };
uniqueSub = unique(allSub);
nSub = length(uniqueSub);
allSess(cellfun(@isempty, allSess)) = { 1 };
allSess = cellfun(@num2str, allSess, 'uniformoutput', false);
uniqueSess = unique(allSess);
factors = pop_listfactors(STUDY.design(model.design), 'gui', 'off', 'level', 'one');

% Channel neighbors
% STUDY.etc.statistics.fieldtrip.channelneighbor = neighbors; % for plot
chanlocs = ALLEEG(1).chanlocs;
nChan = ALLEEG(1).nbchan;
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs,0);

% Build statistical models
for iSub = 1:nSub
    for iSess = 1:length(uniqueSess)
        inds1 = strcmp( uniqueSub{iSub}, allSub);
        inds2 = strcmp( uniqueSess{iSess}, allSess);
        idx  = intersect(inds1, inds2);

        if ~isempty(idx)

            % Check all data per session are in one file
            if length(idx) ~= 1
                error([ 'Cannot calculate contrast because more than 1 dataset per session'  ...
                    'per subject. Merge datasets for each subject and try again.' ]);
            end

            % Check subject name and session
            filepath = ALLEEG(iSub).filepath;
            filename = [uniqueSub{iSub} '_ses-' num2str(iSess) '_design' num2str(model.design) '.set'];

            % Creating fields for limo
%             EEGTMP = std_lm_seteegfields(STUDY, ALLEEG(idx), idx, 'datatype', datatype, 'format', 'cell');
%             ALLEEG = eeg_store(ALLEEG, EEGTMP, idx);
            model.set_files{idx} = fullfile(filepath, filename);
% 
%             OUTEEG = [];
%             if all([ALLEEG(idx).trials] == 1)
%                 OUTEEG.trials = 1;
%             else
%                 OUTEEG.trials = sum([ALLEEG(idx).trials]);
%             end
%             OUTEEG.filepath     = filepath;
%             OUTEEG.filename     = filename;
%             OUTEEG.srate        = ALLEEG(idx).srate;
%             OUTEEG.icaweights   = ALLEEG(idx).icaweights;
%             OUTEEG.icasphere    = ALLEEG(idx).icasphere;
%             OUTEEG.icawinv      = ALLEEG(idx).icawinv;
%             OUTEEG.icachansind  = ALLEEG(idx).icachansind;
%             OUTEEG.etc          = ALLEEG(idx).etc;
%             OUTEEG.times        = ALLEEG(idx).times;
%             OUTEEG.chanlocs = ALLEEG(idx).chanlocs;
%             OUTEEG.etc.merged{1}= ALLEEG(idx).filename;
%             single_trials_filename = char(fullfile(filepath, [filename(1:6) '.' model.measure ]));
%             if exist(single_trials_filename,'file')
%                 if strcmpi(measure(4:end),'erp')
%                     OUTEEG.etc.datafiles.daterp = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'spec')
%                     OUTEEG.etc.datafiles.datspec = EEGTMP.etc.datafiles.(measure);
%                 elseif strcmpi(measure(4:end),'ersp')
%                     OUTEEG.etc.datafiles.datersp = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'datitc')
%                     OUTEEG.etc.datafiles.datitc = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'dattimef')
%                     OUTEEG.etc.datafiles.dattimef = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'icaerp')
%                     OUTEEG.etc.datafiles.icaerp = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'icaspec')
%                     OUTEEG.etc.datafiles.icaspec = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'icaersp')
%                     OUTEEG.etc.datafiles.icaersp = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'icaitc')
%                     OUTEEG.etc.datafiles.icaitc = single_trials_filename;
%                 elseif strcmpi(measure(4:end),'icatimef')
%                     OUTEEG.etc.datafiles.icatimef = single_trials_filename;
%                 end
%             end
%             % Save info
%             EEG = OUTEEG;
%             save('-mat', fullfile(filepath, OUTEEG.filename), 'EEG');
%             clear OUTEEG filepath_tmp

            % Statistical model
            fprintf('Building statistical model for %s ... \n',filename)

            % Create matices for continuous and categorical variables
            trialinfo = std_combtrialinfo(STUDY.datasetinfo, idx);
            [catMat,contMat,limodesign] = std_limodesign(factors, trialinfo, ...
                'splitreg',model.splitreg,'interaction', model.interaction);
            if strcmpi(model.splitreg,'on') && ~isempty(model.contMat) % std_limodesign does something else when splitting regressors
                for c = size(contMat,2):-1:1
                    splitreg{c} = limo_split_continuous(catMat,contMat(:,c)); 
                end
                contMat = cell2mat(splitreg);
                zscoreDat = 0; % regressors are now zscored
            end
            model.cat_files{idx} = catMat;
            model.cont_files{idx} = contMat;
%             if isfield(limodesign, 'categorical')
%                 limo.categorical = limodesign.categorical;
%             else
%                 limo.categorical = {};
%             end
%             if isfield(limodesign, 'continuous')
%                 limo.continuous = limodesign.continuous;
%             else
%                 limo.continuous = {};
%             end
%             limo.subjects(idx).subject     = STUDY.datasetinfo(idx(1)).subject;
%             limo.subjects(idx).cat_file    = catMat;
%             limo.subjects(idx).cont_file   = contMat;
        end
    end
end

% Add contrasts for conditions that were merged during design selection
% i.e. multiple categorical variables (factors) and yet not matching the number
% of variables (contrasts are then a weighted sum of the crossed factors)
if ~isempty(factors) && isfield(factors, 'value') && sum(arrayfun(@(x) ~strcmpi(x.label,'group'), STUDY.design(model.design).variable)) == 1 % only one non-continuous variable other than group
    if length(STUDY.design(model.design).variable(1).value) ~= length(factors) % and this var has more values than the number of factors
        limocontrast = zeros(length(STUDY.design(opt.design).variable(1).value),length(factors)+1); % length(factors)+1 to add the constant
        for n=length(factors):-1:1
            factor_names{n} = factors(n).value;
        end

        index = find(arrayfun(@(x) ~strcmpi(x.label,'group'),STUDY.design(opt.design).variable)); % which one is not group
        for c = 1:length(STUDY.design(opt.design).variable(index).value)
            limocontrast(c,1:length(factors)) = single(ismember(factor_names,STUDY.design(opt.design).variable(index).value{c}));
            limocontrast(c,1:length(factors)) = limocontrast(c,1:length(factors)) ./ sum(limocontrast(c,1:length(factors))); % scale by the number of variables
        end
    end
end

model.set_files  = model.set_files';
model.cat_files  = model.cat_files';
model.cont_files = model.cont_files';
if all(cellfun(@isempty, model.cat_files )), model.cat_files  = []; end
if all(cellfun(@isempty, model.cont_files)), model.cont_files = []; end

[LIMO_files, procstatus] = limo_batch('model specification',model,[],STUDY);
% limo_batch 
% psom_run_pipeline
% psom_pipeline_process
% psom_run_script
% psom_run_job
% import: line 128 - limo_batch_import_data
% design: line 128 - 
% gm: line128 - 

STUDY.limo.model = model;
STUDY.limo.datatype = model.defaults;
STUDY.limo.chanloc = limoChanlocs;

% Between session contrasts
index = 1;
for s = 1:nb_subjects
    sess_index = find(cellfun(@(x) strcmpi(x,uniqueSubjects{s}), allSubjects));
    if length(sess_index) > 1
        fprintf('std_limo, computing additional between sessions contrasts for subject %s\n',uniqueSubjects{s})
        sess_name = allSessions(sess_index);
        pairs = nchoosek(1:length(sess_index),2); % do all session pairs
        parfor p = 1:size(pairs,1)
            strpair = [cell2mat(sess_name(pairs(p,1))) cell2mat(sess_name(pairs(p,2)))];
            strpair(isspace(strpair)) = []; % remove spaces
            filesout{p} = limo_contrast_sessions(cell2mat(LIMO_files.mat(sess_index(pairs(p,1)))), ...
                cell2mat(LIMO_files.mat(sess_index(pairs(p,2)))),strpair);
        end

        for f = 1:length(filesout)
            for ff = 1:length(filesout{f})
                allcon{index} = filesout{f}{ff};
                index = index +1;
            end
        end
        clear filesout
    end
end

% GLM name
design_name = STUDY.design(STUDY.currentdesign).name;
design_name(isspace(design_name)) = [];
if contains(design_name,'STUDY.')
    design_name = design_name(7:end);
end
glm_name = [STUDY.filename(1:end-6) '_' design_name '_GLM_' model.defaults.type '_' model.defaults.analysis '_' model.defaults.method];


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

