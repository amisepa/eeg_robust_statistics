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
% Cedric Cannard, Sep 2022

function [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(data1,data2,nBoot,method,dpt)

% Check input parameters
if nargin < 3 || isempty(nBoot)
    nBoot = 1000; % Default number of bootstraps
end
if nargin < 4 || isempty(method)
    method = 'trimmed mean'; % Default statistical method
end
if nargin < 5 || isempty(dpt)
    error("Variables not defined. Please define if data area paired or not to select adequate statistical test. ")
    % dpt = 'dpt'; % Default to independent data
end
if ndims(data1) == 2 % Add an extra dimension if needed (one channel squeezed)
    data1 = reshape(data1, 1, size(data1, 1), size(data1, 2));
    data2 = reshape(data2, 1, size(data2, 1), size(data2, 2));
end

% Data sizes
nChan = size(data1,1);  % number of channels
nTimes = size(data1,2); % number of time/frequency points
nSub = size(data1,3);   % number of subjects

% Parpool with max number of workers
addons = ver;
if any(contains({addons.Name}, 'Parallel'))
    ps = parallel.Settings;
    fprintf('Parallel computing set to ON. \n')
    ps.Pool.AutoCreate = true;
    p = gcp('nocreate');
    % delete(gcp('nocreate')) % shut down opened parpool
    if isempty(p) % if not already on, launch it
        disp('Initiating parrallel computing (all cores and threads -2)...')
        c = parcluster; % cluster profile
        % N = feature('numcores');          % only physical cores
        N = getenv('NUMBER_OF_PROCESSORS'); % all processor (cores + threads)
        if ischar(N), N = str2double(N); end
        c.NumWorkers = N-2;  % update cluster profile to include all workers
        c.parpool();
    end
end

% Run stats on real data (All electrodes)
tvals = nan(size(data1,1),size(data1,2));
pvals = nan(size(data1,1),size(data1,2));
disp('Performing statistical test on observed data (all channels)...');
progressbar('EEG channels')
for iChan = 1:nChan

    x1 = data1(iChan,:,:);
    x2 = data2(iChan,:,:);

    nanSubj = squeeze(isnan(x1(:,1,:)) | isnan(x2(:,1,:)))';
    if any(nanSubj)
        warning('%g NaN subject(s) detected and removed from both variables!',sum(nanSubj))
        x1(:,:,nanSubj) = [];
        x2(:,:,nanSubj) = [];        
    end

    if strcmpi(method,'trimmed mean')
        if strcmpi(dpt, 'dpt')
            [tval,~,~,~,pval,~,~] = limo_yuend_ttest(x1,x2,20,0.05);
%             [tval,~,~,~,~,pval] = yuend(x1,x2,20,0.05);   % for 2D vecotrs
        elseif strcmpi(dpt, 'idpt')
            [tval,~,~,~,pval,~,~] = limo_yuen_ttest(x1,x2,20,0.05);
%             [tval,~,~,~,~,~,pval] = yuen(x1,x2,20,0.05);  % for 2D vectors
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
    tvals(iChan,:) = tval;
    pvals(iChan,:) = pval;

    progressbar(iChan/nChan)
end
clear tval; clear pval

% Generate boot table (H0)
b = 1;
boot_index = zeros(nSub-sum(nanSubj),nBoot);
disp('Generating boot table (H0)...')
while b ~= nBoot + 1
    tmp = randi(nSub-sum(nanSubj), nSub-sum(nanSubj), 1);
    if length(unique(tmp)) >= 4   % minimum number of subjects/trials
        boot_index(:,b) = tmp;
        b = b + 1;
    else
        error('Not enough subjects, minimum is 4 for degrees of freedom (i.e. n = 3)')
    end
end
clear tmp
for iChan = size(data1,1):-1:1
    boot_table{iChan} = boot_index;
end

% Center data to estimate H0
if strcmpi(method,'trimmed Mean')
    data1_centered = data1 - repmat(limo_trimmed_mean(data1),[1 1 nSub]);
    data2_centered = data2 - repmat(limo_trimmed_mean(data2),[1 1 nSub]);
elseif strcmpi(method,'mean')
    data1_centered = data1 - repmat(mean(data1,3,'omitnan'),[1 1 nSub]);
    data2_centered = data2 - repmat(mean(data2,3,'omitnan'),[1 1 nSub]);
else
    error('The method input must be ''mean'' or ''trimmed mean'' ')
end

% Estimate H0 for each channel using ttests on null data
tvals_H0 = nan(nChan, nTimes, nBoot);
pvals_H0 = nan(nChan, nTimes, nBoot);
disp('Running statistical tests on H0 data...');
progressbar('Estimating H0 on all channels')
for iChan = 1:nChan
    fprintf('Estimating H0 for channel %g/%g\n', iChan, nChan)
    % nanSubj = squeeze(isnan(data1_centered(iChan,1,:)));
    % x1 = data1_centered(iChan,:,~nanSubj);
    % nanSubj = squeeze(isnan(data2_centered(iChan,1,:)));
    % x2 = data2_centered(iChan,:,~nanSubj);

    x1 = data1_centered(iChan,:,:);
    x2 = data2_centered(iChan,:,:);

    nanSubj = squeeze(isnan(x1(:,1,:)) | isnan(x2(:,1,:)))';
    if any(nanSubj)
        % warning('%g NaN subject(s) detected and removed from both variables!',sum(nanSubj))
        x1(:,:,nanSubj) = [];
        x2(:,:,nanSubj) = [];        
    end

    parfor b = 1:nBoot
        if strcmpi(method,'trimmed Mean')
            if strcmpi(dpt, 'dpt')
                [tval{b}, ~, ~, ~, pval{b}, ~, ~] = limo_yuend_ttest( x1(:,:,boot_table{iChan}(:,b)), x2(:,:,boot_table{iChan}(:,b)) ); % Yuen t-test for depedent variables
            elseif strcmpi(dpt, 'idpt')
                [tval{b}, ~, ~, ~, pval{b}, ~, ~] = limo_yuen_ttest( x1(:,:,boot_table{iChan}(:,b)), x2(:,:,boot_table{iChan}(:,b)) ); % Yuen t-test for indepedent variables
            else
                error("'dpt' input must be 'dpt' (paired data) or 'idpt' (unpaired data)")
            end

        elseif strcmpi(method,'mean')
            if strcmpi(dpt, 'dpt')
			    [~,~,~,~,~,tval{b},pval{b}] = limo_ttest(1,x1(:,:,boot_table{iChan}(:,b)), x2(:,:,boot_table{iChan}(:,b)), 0.05); % paired t-test for depedent variables
            elseif strcmpi(dpt, 'idpt')
			    [~,~,~,~,~,tval{b},pval{b}] = limo_ttest(2,x1(:,:,boot_table{iChan}(:,b)), x2(:,:,boot_table{iChan}(:,b)), 0.05); % paired t-test for independent variables
            else
                error("'dpt' input must be 'dpt' (paired data) or 'idpt' (unpaired data)")
            end
        else
            error("The method input must be 'mean' or 'trimmed mean'")
        end
    end

    parfor b = 1:nBoot
        tvals_H0(iChan,:,b) = tval{b};
        pvals_H0(iChan,:,b) = pval{b};
    end

    progressbar(iChan/nChan)

end

disp('Boostrap statistics completed.')
