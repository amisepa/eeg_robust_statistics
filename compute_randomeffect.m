%% Run t-tests comparing two EEG datasets (e.g., conditions or groups)
% on real data and H0 data (i.e., bootstrap).
% 
% Cedric Cannard, Sep 2022

function [results, results_H0] = compute_randomeffect(data1,data2,nboot,method)

% Run stats on real data
results = nan(size(data1,1),size(data1,2),2); % 2 for tval and pval
progressbar('EEG channels')
for iChan = 1:size(data1,1)
    idx = isnan(data1(iChan,1,:));
    y1 = data1(iChan,:,~idx);
    idx = isnan(data2(iChan,1,:));
    y2 = data2(iChan,:,~idx);

    if strcmpi(method,'Trimmed Mean')
        disp(['Applying Yuen t-test for channel ' num2str(iChan)]);
        [tval, ~, ~, ~, pval, ~, ~] = limo_yuend_ttest(y1,y2);
    elseif strcmpi(method,'Mean')
        disp(['Applying paired t-test on channel ' num2str(iChan)  '/' num2str(size(data1,1))]);
        [~, ~, ~, ~, ~, tval, pval] = limo_ttest(1,y1,y2,.05);
    end
    results(iChan,:,1) = tval;
    results(iChan,:,2) = pval;

    progressbar(iChan / size(data1,1))
end

% Generate boot table (surrogate)
b = 1;
boot_index = zeros(size(data1,3),nboot);
disp('Generating boot table...')
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
if strcmpi(method,'Trimmed Mean')
    data1_centered = data1 - repmat(limo_trimmed_mean(data1),[1 1 size(data1,3)]);
    data2_centered = data2 - repmat(limo_trimmed_mean(data2),[1 1 size(data2,3)]);
elseif strcmpi(method,'Mean')
    data1_centered = data1 - repmat(mean(data1,3,'omitnan'),[1 1 size(data1,3)]);
    data2_centered = data2 - repmat(mean(data2,3,'omitnan'),[1 1 size(data2,3)]);
end

% Estimate H0 for each channel using ttests on null data
results_H0 = NaN(size(data1,1), size(data1,2), 2, nboot);
progressbar('EEG channels','Boots')
for iChan = 1:size(data1,1)
    disp(['Estimating H0 for channel ' num2str(iChan) '/' num2str(size(data1,1))])
    idx = isnan(data1_centered(iChan,1,:));
    y1 = data1_centered(iChan,:,~idx);
    idx = isnan(data2_centered(iChan,1,:));
    y2 = data2_centered(iChan,:,~idx);

    for b = 1:nboot
        if strcmpi(method,'Trimmed Mean')
            disp(['Estimating H0 using Yuen t-test for channel ' num2str(iChan)]);
            [t, ~, ~, ~, p, ~, ~] = limo_yuend_ttest(y1(:,:,boot_table{iChan}(:,b)),y2(:,:,boot_table{iChan}(:,b)));
        elseif strcmpi(method,'Mean')
            disp(['Estimating H0 using paired t-test on channel ' num2str(iChan)  '/' num2str(size(data1,1))]);
        [~,~,~,~,~,t,p] = limo_ttest(1,y1(:,:,boot_table{iChan}(:,b)), y2(:,:,boot_table{iChan}(:,b)), 0.05);
        end

        results_H0(iChan,:,1,b) = t;
        results_H0(iChan,:,2,b) = p;

        frac2 = b/nboot;
        frac1 = (iChan-1 + frac2) / size(data1,1);
        progressbar(frac1, frac2)
    end
end
disp('Random effects done.')

