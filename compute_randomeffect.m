%% Run t-tests comparing two EEG datasets (e.g., conditions or groups)
% on real data and H0 data (i.e., bootstrap).
% 
% INPUTS:
%	data1 	- 3D data, computing bootstrap and stats from 3rd dimension (e.g., channels, samples/frequency, subjects)
%	data2 	- 3D data to compare (another condition or group)
%	nboot 	- number of iterations for the bootstrap (default = 1000)
% 	method 	- 'mean' to use paired t-test and 'trimmed mean' to use Yuen t-test default; 20% trim)
%				Yuen t-test better accounts for outliers and non-normal distributions.
%
% OUTPUTS
% 	results 	- t and p values for real data
%	results_H0 	- t and p values for H0 data
%
% Cedric Cannard, Sep 2022

function [results, results_H0] = compute_randomeffect(data1,data2,nboot,method)

if ~exist(nboot,'var') || isempty(nboot)
	nboot = 1000;
end

if ~exist(method,'var') || isempty(method)
	method = 'trimmed mean';
end


% Run stats on real data
results = nan(size(data1,1),size(data1,2),2); % 2 for tval and pval
disp('Running statistical tests on real data...');
progressbar('EEG channels')
for iChan = 1:size(data1,1)
    idx = isnan(data1(iChan,1,:));
    y1 = data1(iChan,:,~idx);
    idx = isnan(data2(iChan,1,:));
    y2 = data2(iChan,:,~idx);

    if strcmpi(method,'trimmed Mean')
		% ADD DETECTION IF DATA ARE DEPENDENT OR INDEPENDENT (DIFFERENT SIZE)
        [tval, ~, ~, ~, pval, ~, ~] = yuend(y1,y2,20); 		% 20% trimmed means
    elseif strcmpi(method,'mean')
		% ADD DETECTION IF DATA ARE DEPENDENT OR INDEPENDENT (DIFFERENT SIZE)
        [~, ~, ~, ~, ~, tval, pval] = limo_ttest(1,y1,y2,.05);
    else
        errordlg('The method input must be ''mean'' or ''trimmed mean'' ')
    end
    results(iChan,:,1) = tval;
    results(iChan,:,2) = pval;

    progressbar(iChan / size(data1,1))
end

% Generate boot table (surrogate)
b = 1;
boot_index = zeros(size(data1,3),nboot);
disp('Generating boot table of H0 data...')
while b ~= nboot + 1
    tmp = randi(size(data1,3),size(data1,3),1);
    if length(unique(tmp))-1 >= 3   % minimum number of subjects 
        boot_index(:,b) = tmp;
        b = b + 1;
    end
end
clear tmp
for iChan = size(data1,1):-1:1
    boot_table{iChan} = boot_index;
end

% Center data to estimate H0
disp('Centering H0 data.')
if strcmpi(method,'trimmed Mean')
    data1_centered = data1 - repmat(limo_trimmed_mean(data1),[1 1 size(data1,3)]);
    data2_centered = data2 - repmat(limo_trimmed_mean(data2),[1 1 size(data2,3)]);
elseif strcmpi(method,'mean')
    data1_centered = data1 - repmat(mean(data1,3,'omitnan'),[1 1 size(data1,3)]);
    data2_centered = data2 - repmat(mean(data2,3,'omitnan'),[1 1 size(data2,3)]);
else
    errordlg('The method input must be ''mean'' or ''trimmed mean'' ')
end

% Estimate H0 for each channel using ttests on null data
results_H0 = NaN(size(data1,1), size(data1,2), 2, nboot);
disp('Running statistical tests on H0 data...');
progressbar('EEG channels','Boots')
for iChan = 1:size(data1,1)
    disp(['Estimating H0 for channel ' num2str(iChan) '/' num2str(size(data1,1))])
    idx = isnan(data1_centered(iChan,1,:));
    y1 = data1_centered(iChan,:,~idx);
    idx = isnan(data2_centered(iChan,1,:));
    y2 = data2_centered(iChan,:,~idx);

    for b = 1:nboot
        if strcmpi(method,'trimmed Mean')
            % disp(['Estimating H0 using Yuen t-test for channel ' num2str(iChan)]);
            [t, ~, ~, ~, p, ~, ~] = limo_yuend_ttest(y1(:,:,boot_table{iChan}(:,b)),y2(:,:,boot_table{iChan}(:,b)));
        elseif strcmpi(method,'mean')
            % disp(['Estimating H0 using paired t-test on channel ' num2str(iChan)  '/' num2str(size(data1,1))]);
			[~,~,~,~,~,t,p] = limo_ttest(1,y1(:,:,boot_table{iChan}(:,b)), y2(:,:,boot_table{iChan}(:,b)), 0.05);
        else
            errordlg('The method input must be ''mean'' or ''trimmed mean'' ')
        end

        results_H0(iChan,:,1,b) = t;
        results_H0(iChan,:,2,b) = p;

        frac2 = b/nboot;
        frac1 = (iChan-1 + frac2) / size(data1,1);
        progressbar(frac1, frac2)
    end
end

disp('Random effects conputed.')