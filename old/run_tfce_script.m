%% attempt to run TFCE
% 
% Cedric Cannard, Sep 2022

clear; close all; clc;
dataDir = 'G:\My Drive\data_processed2';
outDir = 'G:\My Drive\HLM\limo';
addpath('G:\My Drive\HLM')
cd(dataDir)
tmpSub = dir;
idx = contains({tmpSub.name},'sub');
tmpSub = {tmpSub(idx).name}';
cd(outDir)
eeglab; close;
pop_editoptions('option_single',0); % ensure double precision

% Load study
% [STUDY, ALLEEG] = pop_loadstudy('filename','test.study','filepath',outDir);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
EEG = pop_loadset('filename','sub-01.set','filepath','G:\\My Drive\\data_processed2\\sub-01\\');

%% get channel neighbors
% load a file to get channel locations
% eeglab; close;
% pop_editoptions('option_single',0); % ensure double precision
% filepath = 'G:\Shared drives\Grants\Post Award Grants\(736) Bial Full-trance 2017\Research\Data\EEG\trance_bids_2021\sub-001\ses-01\eeg';
% filename = 'sub-001_ses-01_task-trance_eeg.set';
% EEG = pop_loadset('filename',filename, 'filepath',filepath);
% chanlocPath = fileparts(which('dipfitdefs.m'));
% EEG = pop_chanedit(EEG,'lookup',fullfile(chanlocPath,'standard_BEM','elec','standard_1005.elc'));
% chanlocs = EEG.chanlocs;
load(fullfile(dataDir,'chanlocs.mat'))
neighbors = get_channelneighbors(chanlocs,1);
folder = 'C:\Users\IONSLAB\Desktop\channeling_matlab\plot_results\power_electrode';

%%
% data1: elec x freq x subject (e.g., 64 x 45 x 13)
% data2: elec x freq x subject (e.g., 64 x 45 x 13)

stattest = 'paired t-test';
% expected_chanlocs = fullfile(studypath, 'derivatives', 'limo_gp_level_chanlocs.mat');
% options = {'nboot' 1000 'tfce' 1 'type' 'Channels'};
% LIMO.data = load(expected_chanlocs);
% LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
% LIMO.dir = 'G:\My Drive\HLM\limo\derivatives\LIMO_test\results';
method = 'Trimmed mean';
analysis_type = 'Full scalp analysis';
% parameters = [1 2];

for iChan = 1:size(data1,1)
    tmp = isnan(data1(iChan,1,:));
    tmp2 = isnan(data2(iChan,1,:));
    if length(tmp) ~= length(isnan(tmp2))
        errordlg([LIMO.Type ' ' num2str(iChan) ' has unpaired data - analysis aborded, not a paired t-test']);
        return
    elseif length(tmp) == sum(isnan(tmp))
        errordlg([LIMO.Type ' ' num2str(iChan) ' is empty - analysis aborded']);
        return
    elseif (length(tmp) - sum(isnan(tmp))) < 3
        errordlg([LIMO.Type ' ' num2str(iChan) ' has less than 3 subjects - analysis aborded']);
        return
    end
end

% make a paired_samples file per parameter (channels, frames, [mean value, se, df, t, p])
paired_samples = NaN(size(data1,1), size(data1,2),5);
name = sprintf('paired_samples_ttest_parameter_%s',num2str(parameter')');
% array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
% array = 1:EEG.nbchan;
channels = 1:64;

for iChan = 1:length(channels)
    % remove NaNs
    nans = isnan(data1(iChan,1,:));
    y1 = data1(iChan,:,~nans);
    nans = isnan(data2(iChan,1,:));
    y = data2(iChan,:,~nans);
    
    % t-test
    if contains(method,'Trimmed Mean')
        disp(['Applying Yuen t-test on real data: channel ' num2str(iChan)]);
        [   paired_samples(iChan,:,4), ...    % t-value
            paired_samples(iChan,:,1), ...    % diff
            paired_samples(iChan,:,2), ...    % se
            ~, ...                              % CI
            paired_samples(iChan,:,5), ...    % p-value
            ~, ...                              % tcrit
            paired_samples(iChan,:,3) ...     % df
            ] = limo_yuend_ttest(y1,y);
    else 
        disp(['Applying normal t-test on real data: channel ' num2str(iChan)]);
        [   paired_samples(iChan,:,1), ...    % mean
            paired_samples(iChan,:,3), ...    % dfe
            ~, ...                              % CI
            sd, ...
            n, ...
            paired_samples(iChan,:,4),...     % t-value
            paired_samples(iChan,:,5) ...     % p-value
            ] = limo_ttest(1,y1,y,p);
        paired_samples(iChan,:,2) = sd./sqrt(n);
    end
end

%% Bootstrap
nboot = 1000;

% prep data
H0_paired_samples = NaN(size(data1,1), size(data1,2),2,nboot); % stores T and p values for each boot
data1_centered = data1 - repmat(mean(data1,3,'omitnan'),[1 1 size(data1,3)]);
data2_centered = data2 - repmat(mean(data2,3,'omitnan'),[1 1 size(data2,3)]);

% Build a table of data to resample such as almost the same resampling is applied across channels.
disp('Making boot table ...')
data = data1;
Nmin = 3; % minimum 3 subjects
B = 1;
boot_index = zeros(size(data,3),nboot);
while B ~= nboot + 1
    tmp = randi(size(data,3),size(data,3),1);
    if length(unique(tmp)) >= Nmin % at least Nmin different observations per boot
        boot_index(:,B) = tmp;
        B=B+1;
    end
end
clear chdata tmp

% loop per channel, if no nan use boot_index else change it
array = find(sum(squeeze(isnan(data(:,1,:))),2) < size(data,3)-3);
for iChan = size(array,1):-1:1
    channel       = array(iChan);
    tmp           = squeeze(data(channel,:,:)); % 2D
    nSub = find(~isnan(tmp(1,:)));
    Y             = tmp(:,nSub); % remove NaNs

    boot_index2 = zeros(size(Y,2),nboot);
    for c = 1:nboot
        common  = ismember(boot_index(:,c),nSub');
        current = boot_index(find(common),c); %#ok<FNDSB> % keep resampling of good subjets
        % add or remove indices
        add = size(Y,2) - size(current,1);
        if add > 0
            new_boot = [current ; nSub(randi(size(nSub),add,1))'];
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
end
% save(['H0', filesep, 'boot_table'], 'boot_table')

%% Get results under H0

for iChan = 1:size(array,1)
    channel = array(iChan);
    fprintf('bootstrapping channel %g/%g parameter %s \n',iChan,size(array,1),num2str(parameter')');
    tmp = data1_centered(channel,:,:); 
    y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); 
    clear tmp
    tmp = data2_centered(channel,:,:); 
    y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); 
    clear tmp
    if contains(LIMO.design.method,'Trimmed Mean','IgnoreCase',true)
        disp('Using Yuen paired t-test (i.e., trimmed means)')
        parfor b=1:LIMO.design.bootstrap
            [t{b},~,~,~,p{b},~,~] = limo_yuend_ttest(y1(1,:,boot_table{channel}(:,b)),y2(1,:,boot_table{channel}(:,b)));
        end
    else % if strcmpi(LIMO.design.method,'Mean')
        disp('Using paired t-test (i.e., means)')
        parfor b=1:LIMO.design.bootstrap
            [~,~,~,~,~,t{b},p{b}] = limo_ttest(1,y1(1,:,boot_table{channel}(:,b)),y2(1,:,boot_table{channel}(:,b)));
        end
    end
    for b=1:LIMO.design.bootstrap
        H0_paired_samples(channel,:,1,b) = t{b};
        H0_paired_samples(channel,:,2,b) = p{b};
    end
    clear t p y1 y2
end
% save (['H0', filesep, boot_name],'H0_paired_samples','-v7.3');

%% TFCE

% neighboring matrix since TFCE integrates over clusters
% [~, LIMO.data.neighbouring_matrix] = limo_expected_chanlocs;

% calculate TFCE threshold
fprintf('Calculating threshold-free cluster enhancement (TFCE) for %s \n',filename);
tval = paired_samples;
[tfce_score,thresholded_maps] = limo_tfce(2, squeeze(tval(:,:,4)),LIMO.data.neighbouring_matrix);
H0_tval = H0_paired_samples;
tfce_H0_thmaps = cell(1,nboot);
tfce_H0_score  = NaN(size(H0_tval,1),size(H0_tval,2),nboot);
parfor b=1:nboot
    [tfce_H0_score(:,:,b),tfce_H0_thmaps{b}] = limo_tfce(2,squeeze(H0_tval(:,:,end-1,b)),LIMO.data.neighbouring_matrix,0);
end
% save(H0_tfce_file,'tfce_H0_score','-v7.3'); clear H0_tval tfce_H0_score;
tmp = thresholded_maps; clear thresholded_maps;
thresholded_maps{1} = tmp; clear tmp
thresholded_maps{2} = tfce_H0_thmaps; clear tfce_H0_thmaps;
% save(tfce_file,'tfce_score','-v7.3'); clear Fval;

%% PLOT

c = clock;
p = 0.05;
MCC = 1;  % 1 = uncorrected; 2 = cluster; 3 = tfce; 4 = max;
choice = 'use theoretical p values';
FileName = 'paired_samples_ttest_parameter_12.mat';
toplot = load(fullfile(LIMO.dir,FileName));
toplot = toplot.(cell2mat(fieldnames(toplot)));
Type = 1;   %  1 - 2D image with intensity as function of time/freq (x) and electrodes (y)
%              2 - scalp topography
%              3 - ERP data (original or modeled)
surfflag = 0; % to allow surfing the figure and click (1) or not (0)
