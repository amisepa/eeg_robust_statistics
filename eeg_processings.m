%% Create study from cleaned data and compute PSD
% 
%  Cedric Cannard, Sep 2022

clear; close all; clc;
folder = 'G:\My Drive\HLM';
dataDir = 'C:\Users\IONSLAB\Documents\MATLAB\precog_eeg\data_processed';
eeglab; close;

cd(dataDir)
tmp = dir;
tmpSub = {tmp(3:end).name};

cd(folder)
mkdir(fullfile(folder, 'limo'))

commands = {};
count = 1;
for iSub = 1:6
    filepath = fullfile(dataDir, tmpSub{iSub}, 'ses-01', 'eeg');
    filename = [tmpSub{iSub} '_ses-01_task-precognition_eeg.set'];

    EEG = pop_loadset('filename',filename,'filepath',filepath);
    EEG.subject = sprintf('sub-%2.2d',iSub);
    EEG.run = [];
    EEG.saved = 'no';
    
    newpath = fullfile(folder, 'limo', sprintf('sub-%2.2d',iSub));
    newname = sprintf('sub-%2.2d.set',iSub);
    mkdir(newpath)
    pop_saveset(EEG,'filepath',newpath,'filename',newname);

    commands = [ commands(:)' 'index' count 'load' fullfile(newpath, newname) ];
    [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG, ...
        'name','test', ...
        'commands', commands, ...
        'updatedat','on','savedat','off','rmclust','off');
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); 
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
    count = count+1;

end

STUDY = pop_savestudy(STUDY,EEG,'filename','test.study','filepath',fullfile(folder,'limo'));
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
eeglab redraw

% Study design 1
STUDY = std_makedesign(STUDY, ALLEEG, 1, ...
    'name','design 1','delfiles','off', 'defaultdesign','off', ...
    'variable1','type','values1',{'2','4','8'},'vartype1','categorical', ...
    'subjselect',{'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06'}); 
[STUDY, EEG] = pop_savestudy(STUDY,EEG,'savemode','resave');

% Compute PSD
[STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on','recompute','on', ...
    'spec','on','specparams',{'specmode','psd','logtrials','off','freqrange',[1 15]});


%% LIMO 1st level

pop_limo(STUDY,ALLEEG, ... 
    'method','WLS', ...         % GLM optimization: OLS, IRLS, WLS (default)
    'measure','datspec', ...    % which EEG measure: daterp, datsec, datersp
    'freqlim',[0 30] , ...      % freqs (e.g., [0 30]); or timelim for ERP
    'erase','on', ...           % erase previous model
    'splitreg','off', ...       % split regression for continuous indep. var.
    'interaction','on');        % interaction model for categorical indep. var.

gong

% LIMO 2nd level and plot results --> eegh doesn't work
limo_eeg

%% LIMO PLANNED CONTRASTS (STATS at the subject level)

limoFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_eeg\analysis_eeg\derivatives\LIMO_precognition_eeg';
cd(limoFolder);

%1) Create a contrast file
% contrast = [1 -1 0 0; 0 -1 1 0; 1 -2 1 0; 1 0 -1 0; 0 1 -1 0; -1 0 1 0; -1 1 0 0] ;   %Pleasant vs Neutral; Unpleasant vs Neutral; Emotional vs Neutral ; Pleasant vs Unpleasant; %Neutral vs Unpleasant; Unpleasant vs Pleasant; Neutral vs Pleasant
% contrast = [1 -1 0 0; 0 -1 1 0; 1 -2 1 0] ;   %Pleasant vs Neutral; Unpleasant vs Neutral; Emotional vs Neutral ; Pleasant vs Unpleasant; %Neutral vs Unpleasant; Unpleasant vs Pleasant; Neutral vs Pleasant
% save('C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\bids_data\LIMO_precognition\contrast.mat', 'contrast'); 

% 2)call limo_batch in the command window, 
limo_batch

% 4) "contrast only"
% 5) select the list of LIMO files (txt)
% 6) select the contrast file 
% 7) use the generated con files in the ANOVA/t-test


















%% Create study from .set files already processed by Arno (ERROR DURING PRECOMPUTE)

dataDir = 'C:\Users\IONSLAB\Desktop\channeling_matlab\data\bids';
cd(dataDir)

count = 1;
for iSub = 1:13
    for iSes = 1:2
        filepath = fullfile(dataDir, sprintf('sub-%3.3d',iSub), sprintf('ses-%2.2d',iSes), 'eeg');
        cd(filepath)
        tmp = dir;
        tmp = {tmp(3:8).name}';
        for iFile = 3:length(tmp)
            [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','trance_eeg','commands',...
                {{ 'index', count, 'load', char(fullfile(filepath,tmp(iFile))) }},...
                'updatedat','on','savedat','off','rmclust','off');
            [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
            count = count+1;
        end
    end
end
[STUDY,EEG] = pop_savestudy(STUDY,EEG,'filename','trance_eeg.study','filepath',fullfile(mainDir,'limo'));
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% Study design 1
STUDY = std_makedesign(STUDY,ALLEEG,1,'name','Design 1','delfiles','off','defaultdesign','off',...
    'variable1','condition','values1',{'rest','tran'},'vartype1','categorical',...
    'variable2','session','values2',{1,2},'vartype2','continuous',...
    'subjselect',{'sub-001','sub-002','sub-003','sub-004','sub-005','sub-006',...
    'sub-007','sub-008','sub-009','sub-010','sub-011','sub-012','sub-013'});

% Study design 2 (without session)
STUDY = std_makedesign(STUDY,ALLEEG,2,'name','Design 2','delfiles','off','defaultdesign','off',...
    'variable1','condition','values1',{'rest','tran'},'vartype1','categorical',...
    'subjselect',{'sub-001','sub-002','sub-003','sub-004','sub-005','sub-006',...
    'sub-007','sub-008','sub-009','sub-010','sub-011','sub-012','sub-013'});

[STUDY, EEG] = pop_savestudy(STUDY,EEG,'savemode','resave');

% Precompute FFT
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','rmicacomps','on',...
    'interp','on','recompute','on','spec','on','specparams',{'specmode','fft','logtrials','off'});



%% Reprocess everything for LIMO (cannot have differnet files for sessions and trials)
commands = {};
count = 1;
for iSub = 1:13
    for iSes = 1:2

        disp('--------------------------------------------')
        disp(['Processing subject ' num2str(iSub) ' session ' num2str(iSes)])
        disp('--------------------------------------------')

        % load data
        filepath = fullfile(dataDir, sprintf('sub-%3.3d',iSub), sprintf('ses-%2.2d',iSes), 'eeg');
        cd(filepath)
        tmp = dir; tmp = {tmp(3:8).name}';
        filename = tmp(contains(tmp, '.set'));
        EEG = pop_loadset('filename',filename,'filepath',filepath);
        
        % preprocessing
        EEG = pop_chanedit(EEG,'rplurchanloc',1,'lookup',fullfile(chanLocs,'standard_BEM','elec','standard_1005.elc'));
        EEG = pop_resample(EEG,256);
%         EEG = pop_eegfiltnew(EEG,'locutoff',1);
        EEG = pop_eegfiltnew(EEG,'locutoff',0.75,'filtorder',1690);
        EEG = pop_reref(EEG,[]);
        % EEGcsd = csd_transform(EEG);
        % figure; plot(EEGav.data(1,1:256)); hold on; plot(EEGcsd.data(1,1:256)); legend('average', 'csd')
        
        % keep only rest and trance data (exclude first 10 s of each run) and separate them
        evIdx = contains({EEG.event.type}, 'rest');
        rest = nan(11,2); trance = nan(11,2);
        for iEv = 1:2:length(evIdx)
            if evIdx(iEv) == 1 && evIdx(iEv+1) == 1
                rest(iEv,:) = [ EEG.event(iEv).latency + EEG.srate*10
                    EEG.event(iEv+1).latency - EEG.srate*1 ];
                %                 rest(iEv,:) = [ EEG.event(iEv).latency-1
                %                     EEG.event(iEv+1).latency+1 ];
            elseif evIdx(iEv) == 0 && evIdx(iEv+1) == 0
                trance(iEv,:) = [ EEG.event(iEv).latency + EEG.srate*10
                    EEG.event(iEv+1).latency - EEG.srate*1 ];
            else
                error('check events')
            end
        end
        rest(isnan(rest(:,1)),:) = []; trance(isnan(trance(:,1)),:) = [];
        if iSes == 1
            restdata1 = pop_select(EEG, 'point', rest);       % rest data
            trancedata1 = pop_select(EEG, 'point', trance);     % trance data
        else
            restdata2 = pop_select(EEG, 'point', rest);       % rest data
            trancedata2 = pop_select(EEG, 'point', trance);     % trance data
        end
    end

    % balance channels across sessions (take file with least channels as ref)
%     if restdata1.nbchan ~= restdata2.nbchan
%         if restdata1.nbchan > restdata2.nbchan
%             idx = ~ismember({restdata1.chanlocs.labels}, {restdata2.chanlocs.labels});
%             restdata1 = pop_select(restdata1,'nochannel',{restdata1.chanlocs(idx).labels});
%             trancedata1 = pop_select(trancedata1,'nochannel',{trancedata1.chanlocs(idx).labels});
%         else
%             idx = ~ismember({restdata2.chanlocs.labels}, {restdata1.chanlocs.labels});
%             restdata2 = pop_select(restdata2,'nochannel',{restdata2.chanlocs(idx).labels});
%             trancedata2 = pop_select(trancedata2,'nochannel',{trancedata2.chanlocs(idx).labels});
%         end
%     end
%     restdata1 = eeg_checkset(restdata1);
%     restdata2 = eeg_checkset(restdata2);
%     trancedata1 = eeg_checkset(trancedata1);
%     trancedata2 = eeg_checkset(trancedata2);

    % each condition merged across the two sessions
    EEGR = pop_mergeset(restdata1, restdata2);
    EEGT = pop_mergeset(trancedata1, trancedata1, 1);
    EEGR = eeg_checkset(EEGR);
    EEGT = eeg_checkset(EEGT);

    %reject bad channels
    EEGR = clean_flatlines(EEGR,5,20); 
    EEGT = clean_flatlines(EEGT,5,20); 
    [EEGR,remchansR] = clean_channels(EEGR,.75,4,5,.4,50,.25);
    [EEGT,remchansT] = clean_channels(EEGT,.75,4,5,.4,50,.25);
%     EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,'ChannelCriterion',0.75, ...
%         'LineNoiseCriterion',4,'Highpass','off','BurstCriterion','off', ...
%         'WindowCriterion','off','BurstRejection','off','Distance','Euclidian');
%         EEG = pop_select(EEG,'nochannel',{'POz'});
    if EEGR.nbchan ~= EEGT.nbchan
        warning('different number of channels between the two conditions.')
    end

    % ASR
    EEGR = clean_asr(EEGR,20,max(0.5,1.5*EEGR.nbchan/EEGR.srate),[],2/3,0.075,[-3.5 5.5],1,[],false); 
    EEGT = clean_asr(EEGT,20,max(0.5,1.5*EEGT.nbchan/EEGT.srate),[],2/3,0.075,[-3.5 5.5],1,[],false); 

    % Line noise
    %     [psd, f] = get_psd(EEGR.data(1,:),EEG.srate*2,'hamming',50,[],EEG.srate,[1 80],'psd');
    %     figure; plot(f,psd)
    %     EEGR = pop_eegfiltnew(EEGR,'hicutoff',50,'filtorder',170); % lowpass after ASR for better detection of bursts
    EEGR = pop_cleanline(EEGR,'bandwidth',2,'chanlist',1,'computepower',1,...
        'linefreqs',60,'newversion',1,'normSpectrum',0,'p',0.01,'pad',2,...
        'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2, ...
        'tau',100,'verb',1,'winsize',4,'winstep',1);
%     [psd, f] = get_psd(EEGR.data(1,:),EEG.srate*2,'hamming',50,[],EEG.srate,[1 80],'psd');
%     hold on; plot(f,psd); legend('before', 'after')
    EEGT = pop_cleanline(EEGT,'bandwidth',2,'chanlist',1,'computepower',1,...
        'linefreqs',60,'newversion',1,'normSpectrum',0,'p',0.01,'pad',2,...
        'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2, ...
        'tau',100,'verb',1,'winsize',4,'winstep',1);

    % ICA
    EEGR = pop_runica(EEGR,'icatype','runica','extended',1);
    EEGR = pop_iclabel(EEGR,'default');
    EEGR = pop_icflag(EEGR, [NaN NaN;0.95 1;0.95 1;0.95 1;0.95 1;0.95 1;NaN NaN]);
    
    EEGT = pop_runica(EEGT,'icatype','runica','extended',1);
    EEGT = pop_iclabel(EEGT,'default');
    EEGT = pop_icflag(EEGT, [NaN NaN;0.95 1;0.95 1;0.95 1;0.95 1;0.95 1;NaN NaN]);

    % check
    EEGR = eeg_checkset(EEGR);
    EEGT = eeg_checkset(EEGT);

    % display each segment length for quick visual check
    disp(['REST Run 1: ' num2str((EEGR.event(2).latency - EEGR.event(1).latency) / 512 /60) ' min']);
    disp(['REST Run 2: ' num2str((EEGR.event(3).latency - EEGR.event(2).latency) / 512 /60) ' min']);
    disp(['REST Run 3: ' num2str((EEGR.event(4).latency - EEGR.event(3).latency) / 512 /60) ' min']);
    disp(['REST Run 4: ' num2str((EEGR.event(5).latency - EEGR.event(4).latency) / 512 /60) ' min']);
    disp(['REST Run 5: ' num2str((EEGR.event(6).latency - EEGR.event(5).latency) / 512 /60) ' min']);
    disp(['REST Run 6: ' num2str((EEGR.pnts - EEGR.event(6).latency) / 512 /60) ' min']);
    disp(['TRANCE Run 1: ' num2str((EEGT.event(2).latency - EEGT.event(1).latency) / 512 /60) ' min']);
    disp(['TRANCE Run 2: ' num2str((EEGT.event(3).latency - EEGT.event(2).latency) / 512 /60) ' min']);
    disp(['TRANCE Run 3: ' num2str((EEGT.event(4).latency - EEGT.event(3).latency) / 512 /60) ' min']);
    disp(['TRANCE Run 4: ' num2str((EEGT.event(5).latency - EEGT.event(4).latency) / 512 /60) ' min']);
    disp(['TRANCE Run 5: ' num2str((EEGT.event(6).latency - EEGT.event(5).latency) / 512 /60) ' min']);
    disp(['TRANCE Run 6: ' num2str((EEGT.pnts - EEGT.event(6).latency) / 512 /60) ' min']);

    % save REST and load into STUDY
    EEGR.subject = sprintf('sub-%2.2d',iSub);
    EEGR.condition = 'rest';
    EEGR.run = [];
    EEGR.saved = 'no';
    newpath = fullfile(outDir, 'data_processed', EEGR.subject);
    newname = sprintf('sub-%2.2d_rest.set',iSub);
    mkdir(newpath)
    pop_saveset(EEGR,'filepath',newpath,'filename',newname);
    commands = [ commands(:)' 'index' count 'load' fullfile(newpath, newname) ];
    count = count+1;
    [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','trance_eeg','commands', commands, ...
        'updatedat','on','savedat','off','rmclust','off');
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

    % save TRANCE and load into STUDY
    EEGT.subject = sprintf('sub-%2.2d',iSub);
    EEGT.condition = 'trance';
    EEGT.run = [];
    EEGT.saved = 'no';
    newname = sprintf('sub-%2.2d_trance.set',iSub);
    pop_saveset(EEGT,'filepath',newpath,'filename',newname);
    commands = [ commands(:)' 'index' count 'load' fullfile(newpath, newname) ];
    count = count+1;
    [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','trance_processed','commands', commands, ...
        'updatedat','on','savedat','off','rmclust','off');
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

    disp('--------------------------------------------')
    disp(['Subject ' num2str(iSub) ' done.'])
    disp('--------------------------------------------')
end

[STUDY,EEG] = pop_savestudy(STUDY,EEG,'filename','trance_eeg.study','filepath',outDir);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% study design
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','Design 1','delfiles','off','defaultdesign','off', ...
    'variable1','condition','values1',{'rest','trance'},'vartype1','categorical', ...
    'subjselect',{'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06', ...
    'sub-07','sub-08','sub-09','sub-10','sub-11','sub-12','sub-13'});
[STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');

% Precompute FFT
[STUDY, ALLEEG] = std_precomp(STUDY,ALLEEG,{},'savetrials','on','rmicacomps','on', ...
    'interp','on','recompute','on','spec','on', ...
    'specparams',{'specmode','fft','logtrials','off','freqrange',[1 50]});
[STUDY, EEG] = pop_savestudy(STUDY, EEG, 'savemode','resave');

% Plot average RMS all freqs
STUDY = pop_statparams(STUDY,'condstats','on','method','perm','mcorrect','fdr','alpha',0.05);
STUDY = pop_specparams(STUDY,'freqrange',[1 50],'averagechan','rms');
STUDY = std_specplot(STUDY,ALLEEG,'channels',{EEG(1).urchanlocs.labels},'design',1);

% Topography for each band
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[1 3],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{EEG(1).urchanlocs.labels},'design',1);
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[3 7],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{EEG(1).urchanlocs.labels},'design',1);
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[8 13],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{EEG(1).urchanlocs.labels},'design',1);
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[14 30],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{EEG(1).urchanlocs.labels},'design',1);
STUDY = pop_specparams(STUDY,'plotconditions','apart','averagechan','off','topofreq',[30 50],'freqrange',[]);
std_specplot(STUDY,ALLEEG,'channels',{EEG(1).urchanlocs.labels},'design',1);

gong
