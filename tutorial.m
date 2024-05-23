%% Load sample ERP data containing:
% - 19 subjects with 65-channel EEG and 238 time points
% - three conditions: pleasant images, neutral images, unpleasant images
%   (from the IAPS dataset)

clear; close all; clc
folder_dir = fileparts(which('run_stats_bootstrap'));
disp("Loading ERP data and electrode locations...")
load(fullfile(folder_dir,'ERP_data.mat'))
load(fullfile(folder_dir,'chanlocs.mat'))

%% Compute permutaton statistics

% nPerm = 1000;              % number of permutations to do
% alphalvl = 0.05;            % alpha level
% method =  'trimmed mean';   % estimator: 'mean', 'trimmed mean' (more robust)
% dpt = 'dpt';                % whether the two data variables are dependent ('dpt') or not ('idpt')
% 
% [tvals, pvals] = run_stats_perm(data1, data2, alphalvl, nPerm, dpt);

%% Compute bootstrap statistics implementing mass-univariate analysis to 
% compare brain response to Unpleasant images relative to Neutral ones.

nBoots = 1000;              % number of bootstrap iterations to do
alphalvl = 0.05;            % alpha level
method =  'trimmed mean';   % estimator: 'mean', 'trimmed mean' (more robust)
dpt = 'dpt';                % whether the two data variables are dependent ('dpt') or not ('idpt')

tstart = tic;
[tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(UNPLEASANT, NEUTRAL, nBoots, method, dpt);
fprintf('Time to run: %g min\n', round(toc(tstart)/60,2))

%% Apply robust corrections for type 1 error (family-wise error) and plot 
% results

% Get electrode neighbors for spatiotemporal cluster-based corrections
[~, neighbormatrix] = get_channelneighbors(chanlocs);

% Uncorrected
corr_method = 0;  % 0 = uncorrected; 1 = max; 2 = cluster; 3 = tfce
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, corr_method, alphalvl, neighbormatrix);
[peakClust, peakChan, peakLat] = plot_results('time', times, tvals, mask, pcorr, alphalvl, chanlocs, corr_method);
set(gcf,'Name','Uncorrected results','Toolbar','none','Menu','none','NumberTitle','Off'); 

% Maximum likelihood correction
corr_method = 1;  % 0 = uncorrected; 1 = max; 2 = cluster; 3 = tfce
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, corr_method, alphalvl, neighbormatrix);
[peakClust, peakChan, peakLat] = plot_results('time', times, tvals, mask, pcorr, alphalvl, chanlocs, corr_method);
set(gcf,'Name','Maximum likelihood correction','Toolbar','none','Menu','none','NumberTitle','Off'); 

% Spatiotemporal cluster correction
corr_method = 2;  % 0 = uncorrected; 1 = max; 2 = cluster; 3 = tfce
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, corr_method, alphalvl, neighbormatrix);
[peakClust, peakChan, peakLat] = plot_results('time', times, tvals, mask, pcorr, alphalvl, chanlocs, corr_method);
set(gcf,'Name','Spatiotemporal cluster correction','Toolbar','none','Menu','none','NumberTitle','Off'); 

% Threshold-free cluster enhancement (TFCE) correction --> best
corr_method = 3;  % 0 = uncorrected; 1 = max; 2 = cluster; 3 = tfce
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, corr_method, alphalvl, neighbormatrix);
[peakClust, peakChan, peakLat] = plot_results('time', times, tvals, mask, pcorr, alphalvl, chanlocs, corr_method);
set(gcf,'Name','Threshold-free cluster enhancement (TFCE) correction','Toolbar','none','Menu','none','NumberTitle','Off'); 

%% Visualize ERP at the peak electrode in more details using 20% trimmed mean 
% and 95% confidence intervals (CIs)

plotDiff(times, squeeze(UNPLEASANT(peakChan,:,:)), squeeze(NEUTRAL(peakChan,:,:)),'Trimmed mean','CI',mask(peakChan,:),'unpleasant','neutral')
subplot(2,1,1); title(sprintf('Channel %s (Mean + 95%% CIs)',chanlocs(peakChan).labels))

%% Visualize ERP at the peak electrode in more details using 20% trimmed mean 
% and 95% Bayesian high-density intervals (HDI)

plotHDI(times, squeeze(UNPLEASANT(peakChan,:,:)), squeeze(NEUTRAL(peakChan,:,:)),'Trimmed mean',alphalvl,mask(peakChan,:),'unpleasant','neutral')
subplot(2,1,1); title(sprintf('Channel %s (Mean + 95%% HDIs)',chanlocs(peakChan).labels))


%% Perform GLM

clear; close all; clc;
dataDir = 'C:\Users\Tracy\Desktop\data';
outDir = fullfile(dataDir,'glm_poststim'); mkdir(outDir)
eeglab; close;

tlims = [-100 950];     % time window (in ms)
optimization = 'WLS';   % 'OLS' 'IRLS' or 'WLS'
weightmethod = 'Hubert';   % 'PCP' 'Hubert' 'Tukey'

% Get filenames
cd(dataDir)
filenames = dir; filenames = {filenames.name}';
filenames(~contains(filenames,'sub')) = [];

F = [];
R2 = [];
P  = [];
BETAS = [];
PLEASANT = [];
NEUTRAL = [];
UNPLEASANT = [];
nSub = length(filenames);
progressbar('Running GLM on all subjects')
for iSub = 1:nSub

    disp('-----------------------------------------')
    fprintf('               Subject %g\n', iSub)
    disp('-----------------------------------------')

    filepath = fullfile(dataDir,filenames{iSub});
    tmp = dir(filepath); tmp = {tmp.name}';
    idx = ~contains(tmp,'.set');
    for i = 1:sum(idx)
        if idx(i)
            delete(fullfile(filepath,tmp{i}))
        end
    end
    filename = sprintf('%s.set',filenames{iSub});
    EEG = pop_loadset(fullfile(filepath,filename));
    times = EEG.times;

    % Run GLM
    [betas,rsquared,fstat,pvals,times] = run_glm(EEG.data,times,tlims,EEG.event,optimization,weightmethod);
    
    % Save the results
    subject = EEG.subject;
    chanlocs = EEG.chanlocs;
    save(fullfile(filepath, sprintf('%s_GLM_%s.mat', subject, optimization)), 'times','chanlocs','betas','rsquared','pvals','fstat');

    % Group level
    BETAS(:,:,:,iSub) = betas;
    F(:,:,iSub) = fstat;
    R2(:,:,iSub) = rsquared;
    P(:,:,iSub) = pvals;
    
    % ERP raw for comparison
    idx = EEG.times>=tlims(1) & EEG.times<=tlims(2);
    EEG.data = EEG.data(:,idx,:);
    if length({EEG.event.type})~=size(EEG.data,3)
        warning("Event structure has different number of trials than EEG data. Trying to correct.")
        lats = [EEG.event.latency];
        idx = diff(lats)<mean(diff(lats));
        warning("Removing %g events with wrong latency", sum(idx))
        EEG.event(idx) = [];
    end
    idx = strcmp({EEG.event.type},'2');
    PLEASANT(:,:,iSub) =  trimmean(EEG.data(:,:,idx),20,3);     % for comparison
    idx = strcmp({EEG.event.type},'4');
    NEUTRAL(:,:,iSub) =  trimmean(EEG.data(:,:,idx),20,3);      % for comparison
    idx = strcmp({EEG.event.type},'8');
    UNPLEASANT(:,:,iSub) =  trimmean(EEG.data(:,:,idx),20,3);   % for comparison

    progressbar(iSub/nSub)

end

% Save GLM outputs
chanlocs = EEG.chanlocs;
save(fullfile(outDir,sprintf('GLM_%s.mat', optimization)),'times','chanlocs','BETAS','F','R2','P')
save(fullfile(outDir,'RAW_ERP.mat'),'times','chanlocs','PLEASANT','NEUTRAL','UNPLEASANT')

% Compare Raw ERP with GLM Betas for each condition
% idx = times>=tlims(1) & times<=tlims(2);
% times = times(idx);
plotHDI(times, squeeze(trimmean(PLEASANT,20,1)), squeeze(trimmean(BETAS(:,:,1,:),20,1)),'Trimmed mean',0.05,[],'raw','beta')
subplot(2,1,1); title(sprintf('Pleasant (%s)',optimization))
plotHDI(times, squeeze(trimmean(NEUTRAL,20,1)), squeeze(trimmean(BETAS(:,:,2,:),1)),'Trimmed mean',0.05,[],'raw','beta')
subplot(2,1,1); title(sprintf('Neutral (%s)',optimization))
plotHDI(times, squeeze(trimmean(UNPLEASANT,20,1)), squeeze(trimmean(BETAS(:,:,3,:),1)),'Trimmed mean',0.05,[],'raw','beta')
subplot(2,1,1); title(sprintf('Unpleasant (%s)',optimization))

%% Run mass-univariate analysis and nonparametric statistics with spatiotemporal
% cluster correction, to compare the effects at the group level, as if we
% ran a study. 
