%% Run t-tests comparing two EEG datasets (e.g., conditions or groups)
% on real data and H0 data (i.e., bootstrap).
% 
% INPUTS:
%	data1 	- 3D data, computing bootstrap and stats from 3rd dimension (e.g., channels, samples/frequency, subjects)
%	data2 	- 3D data to compare (another condition or group)
%	nBoot 	- number of iterations for the bootstrap (default = 1000)
% 	method 	- 'mean' to use paired t-test and 'trimmed mean' to use Yuen t-test default; 20% trim)
%				Yuen t-test better accounts for outliers and non-normal distributions.
%   dpt     - variables are dependent ('dpt', paired t-test) or not ('idpt', two-sample t-test)
%
% OUTPUTS
% 	tvals 	    - t-values for real data
% 	pvals 	    - p-values for real data
%	tvals_H0    - t-values for H0 data
%	pvals_H0    - p-values for H0 data
%
% EXAMPLE
%   [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(data1,data2,nBoot,method,dpt)
% 
% Cedric Cannard, Sep 2022
% Modified Jan 2025: Fixed independent samples with unequal N

function [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(data1,data2,nBoot,method,dpt)

% add path to subfunctions
tmp = fileparts(which('run_stats_bootstrap'));
addpath(fullfile(tmp,'functions'))

% Check input parameters
if nargin < 3 || isempty(nBoot)
    nBoot = 1000; % Default number of bootstraps
end
if nargin < 4 || isempty(method)
    method = 'trimmed mean'; % Default statistical method
end
if nargin < 5 || isempty(dpt)
    error("Variables not defined. Please define if data area paired or not to select adequate statistical test. ")
end
if ndims(data1) == 2 % if one channel and squeezed, add an extra channel dimension to be compatible
    data1 = reshape(data1, 1, size(data1, 1), size(data1, 2));
    data2 = reshape(data2, 1, size(data2, 1), size(data2, 2));
end

% Data sizes
nChan = size(data1,1);  % number of channels
nTimes = size(data1,2); % number of time/frequency points
nSub1 = size(data1,3);  % number of subjects/trials in group 1
nSub2 = size(data2,3);  % number of subjects/trials in group 2

% For paired design, require equal N
if strcmpi(dpt, 'dpt') && nSub1 ~= nSub2
    error('Paired design requires equal number of subjects in data1 and data2');
end

% Run stats on real data (All electrodes)
tvals = nan(size(data1,1),size(data1,2));
pvals = nan(size(data1,1),size(data1,2));
disp('Performing statistical test on observed data (all channels)...');

parfor iChan = 1:nChan    
    % Initialize temporaries for parfor
    tval = nan(1, nTimes);
    pval = nan(1, nTimes);
    
    x1 = data1(iChan,:,:);
    x2 = data2(iChan,:,:);

    % Remove subjects that contain any NaN
    if strcmpi(dpt, 'dpt')
        % Paired: remove NaN subjects from both
        nanSubj = squeeze(any(any(isnan(x1), 1), 2)) | squeeze(any(any(isnan(x2), 1), 2));
        if any(nanSubj)
            warning('%g NaN subject(s) detected and removed from both variables!', sum(nanSubj));
            x1(:,:,nanSubj) = [];
            x2(:,:,nanSubj) = [];
        end
    else
        % Independent: remove NaNs separately
        nanSubj1 = squeeze(any(any(isnan(x1), 1), 2));
        nanSubj2 = squeeze(any(any(isnan(x2), 1), 2));
        if any(nanSubj1)
            warning('%g NaN subject(s) detected and removed from data1!', sum(nanSubj1));
            x1(:,:,nanSubj1) = [];
        end
        if any(nanSubj2)
            warning('%g NaN subject(s) detected and removed from data2!', sum(nanSubj2));
            x2(:,:,nanSubj2) = [];
        end
    end
    
    if ~isempty(x1) && ~isempty(x2) && size(x1,3)>0 && size(x2,3)>0
        if strcmpi(method,'trimmed mean')
            if strcmpi(dpt, 'dpt')
                [tval,~,~,~,pval,~,~] = limo_yuend_ttest(x1,x2,20,0.05);
            elseif strcmpi(dpt, 'idpt')
                [tval,~,~,~,pval,~,~] = limo_yuen_ttest(x1,x2,20,0.05);
            else
                error("'dpt' input must be 'dpt' (paired data) or 'idpt' (unpaired data)")
            end
        elseif strcmpi(method,'mean')
            if strcmpi(dpt, 'dpt')
                [~,~,~,~,~,tval,pval] = limo_ttest(1,x1,x2,.05);
            elseif strcmpi(dpt, 'idpt')
                [~,~,~,~,~,tval,pval] = limo_ttest(2,x1,x2,.05);
            else
                error("'dpt' input must be 'dpt' (paired data) or 'idpt' (unpaired data)")
            end
        else
            error('The method input must be ''mean'' or ''trimmed mean'' ')
        end
    else
        warning('No data left for channel %d after removing NaNs', iChan);
    end
    
    tvals(iChan,:) = tval;
    pvals(iChan,:) = pval;
end
clear tval; clear pval

% minimum number of degrees of freedom
if strcmpi(dpt, 'dpt')
    if nSub1-1 < 4
        error('Not enough subjects, minimum is 4 for degrees of freedom (i.e. n = 3)')
    end
else
    if min(nSub1, nSub2) < 3
        error('Not enough subjects in at least one group, minimum is 3')
    end
end

% Generate boot tables (H0) - SEPARATE for independent groups
disp('Generating boot table(s) (H0)...')
if strcmpi(dpt, 'dpt')
    % Paired: single boot table
    boot_index1 = zeros(nSub1, nBoot);
    for b = 1:nBoot
        boot_index1(:,b) = randi(nSub1, nSub1, 1);
    end
    % Replicate for all channels
    for iChan = nChan:-1:1
        boot_table1{iChan} = boot_index1;
    end
    boot_table2 = boot_table1;  % Same indices for paired data
else
    % Independent: separate boot tables for each group
    boot_index1 = zeros(nSub1, nBoot);
    boot_index2 = zeros(nSub2, nBoot);
    for b = 1:nBoot
        boot_index1(:,b) = randi(nSub1, nSub1, 1);
        boot_index2(:,b) = randi(nSub2, nSub2, 1);
    end
    % Replicate for all channels
    for iChan = nChan:-1:1
        boot_table1{iChan} = boot_index1;
        boot_table2{iChan} = boot_index2;
    end
end

% Center data to estimate H0
if strcmpi(method,'trimmed mean')
    data1_centered = data1 - repmat(limo_trimmed_mean(data1),[1 1 nSub1]);
    data2_centered = data2 - repmat(limo_trimmed_mean(data2),[1 1 nSub2]);
elseif strcmpi(method,'mean')
    data1_centered = data1 - repmat(mean(data1,3,'omitnan'),[1 1 nSub1]);
    data2_centered = data2 - repmat(mean(data2,3,'omitnan'),[1 1 nSub2]);
else
    error('The method input must be ''mean'' or ''trimmed mean'' ')
end

% Estimate H0 for each channel using ttests on null data
% Preallocate outputs
tvals_H0 = nan(nChan, nTimes, nBoot);
pvals_H0 = nan(nChan, nTimes, nBoot);
disp('Running statistical tests on H0 data...');
parfor iChan = 1:nChan
    fprintf('Computing bootstrap statistics under H0 for channel %g/%g\n', iChan, nChan)

    % Extract data for this channel
    x1 = data1_centered(iChan,:,:);
    x2 = data2_centered(iChan,:,:);

    % Remove NaN subjects
    if strcmpi(dpt, 'dpt')
        % Paired: remove from both
        nanSubj = squeeze(isnan(x1(:,1,:)) | isnan(x2(:,1,:)))';
        if any(nanSubj)
            x1(:,:,nanSubj) = [];
            x2(:,:,nanSubj) = [];
        end
    else
        % Independent: remove separately
        nanSubj1 = squeeze(isnan(x1(:,1,:)))';
        nanSubj2 = squeeze(isnan(x2(:,1,:)))';
        if any(nanSubj1), x1(:,:,nanSubj1) = []; end
        if any(nanSubj2), x2(:,:,nanSubj2) = []; end
    end

    % Preallocate temporary arrays for this channel
    tval_tmp = nan(nTimes, nBoot);
    pval_tmp = nan(nTimes, nBoot);
    
    for b = 1:nBoot
        % Initialize temporaries for inner loop
        tval_b = nan(1, nTimes);
        pval_b = nan(1, nTimes);
        
        % Get bootstrap indices (different for each group in independent case)
        idx1 = boot_table1{iChan}(:,b);
        idx2 = boot_table2{iChan}(:,b);

        if strcmpi(method,'trimmed mean')
            if strcmpi(dpt, 'dpt')
                [tval_b, ~, ~, ~, pval_b, ~, ~] = limo_yuend_ttest(x1(:,:,idx1), x2(:,:,idx2));
            elseif strcmpi(dpt, 'idpt')
                [tval_b, ~, ~, ~, pval_b, ~, ~] = limo_yuen_ttest(x1(:,:,idx1), x2(:,:,idx2));
            else
                error("'dpt' input must be 'dpt' or 'idpt'")
            end

        elseif strcmpi(method,'mean')
            if strcmpi(dpt, 'dpt')
                [~,~,~,~,~,tval_b,pval_b] = limo_ttest(1, x1(:,:,idx1), x2(:,:,idx2), 0.05);
            elseif strcmpi(dpt, 'idpt')
                [~,~,~,~,~,tval_b,pval_b] = limo_ttest(2, x1(:,:,idx1), x2(:,:,idx2), 0.05);
            else
                error("'dpt' input must be 'dpt' or 'idpt'")
            end
        else
            error("The method input must be 'mean' or 'trimmed mean'")
        end

        tval_tmp(:,b) = tval_b;
        pval_tmp(:,b) = pval_b;
    end

    % Store results for this channel
    tvals_H0(iChan,:,:) = tval_tmp;
    pvals_H0(iChan,:,:) = pval_tmp;
end

disp('Bootstrap statistics completed.')
end