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


