%% RUN GLM on ERP DATA - Vectorized approach to significantly increase
% computation cost. Vectorized code performs matrix operations directly,
% which are generally more efficient and less error-prone than iterating
% through individual elements. Matrix multiplication and other linear
% algebra operations in MATLAB are highly optimized.
%
% The core computations, such as the design matrix multiplication (XTWX, XTWY),
% fitting the GLM (pinv(XTWX) * XTWY), and computing residuals (data_reshaped - Yhat),
% are the same in both approaches.
%
% The update of weights in IRLS (and the use of weights in WLS) is performed
% in a manner that considers all data points simultaneously, just as in the
% looped approach.
%
% Some considerations:
% Ensure the numerical precision and stability are maintained, especially
% with operations like pseudoinverse (pinv). In some cases, the use of
% regularization or other stabilization techniques may be necessary.
% Ensure that edge cases (e.g., handling of zero or negative weights, rank
% deficiency) are managed consistently in both approaches.

clear; close all; clc;
dataDir = 'C:\Users\Tracy\Desktop\data';
outDir = fullfile(dataDir,'glm_prestim'); mkdir(outDir)
eeglab; close;

optimization = 'WLS'; % 'OLS' 'IRLS' or 'WLS'
tlims = [-1000 10];     % time window (in ms)

% Get filenames
cd(dataDir)
filenames = dir; filenames = {filenames.name}';
filenames(~contains(filenames,'sub')) = [];

% Convergence criteria for IRLS optimization
if strcmp(optimization,'IRLS')
    max_iter = 100;
    tol = 1e-6;
    k = 1.345; % Tuning parameter for Huber's function
    % Combines squared error for small residuals and absolute error for large residuals, offering a balance between efficiency and robustness.
    % Tukey's Biweight Function: Reduces the influence of outliers more aggressively by setting large residuals to zero weight.
end

% BETAS = nan(65,750,3,20);
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
    EEG = csd_transform(EEG,'C:\Users\Tracy\Documents\MATLAB\PAA_EEG\paa_eeg_code\chanlocs_BEM.ced');
    data = EEG.data;

    % trim time window of interest
    idx = EEG.times>=tlims(1) & EEG.times<=tlims(2);
    data = data(:,idx,:);
    times = EEG.times(idx);

    % ensure data and events have same number of trials
    event_types = {EEG.event.type};
    event_types = cellfun(@num2str, event_types, 'UniformOutput', false);  % ensure events are strings


    if length(event_types)~=size(data,3)
        warning("Event structure has different number of trials than EEG data. Trying to correct.")
        lats = [EEG.event.latency];
        idx = diff(lats)<mean(diff(lats));
        warning("Removing %g events with wrong latency", sum(idx))
        event_types(idx) = [];

        % Check again
        if length(event_types)~=size(data,3)
            error("Still not the same amount of trials betwen event_types and EEG data")
        end
        % data(:,:,idx) = [];
    end

    % Remove unecessary conditions
    idx = strcmp(event_types,'1');
    warning("removing %g trials from unecessary condition", sum(idx))
    event_types(idx) = [];
    data(:,:,idx) = [];

    % Define the categorical design matrix
    % conditions = cell(size(event_types));
    % conditions(strcmp(event_types,'2')) = {'pleasant'};
    % conditions(strcmp(event_types,'4')) = {'neutral'};
    % conditions(strcmp(event_types,'8')) = {'unpleasant'};
    % conditions = categorical(conditions);  % convert to categorical

    % % Define one-hot encoded design matrix
    % % Assuming '2' for pleasant stimuli, '4' for neutral stimuli, and '8' for unpleasant stimuli
    % % event_types = cellfun(@str2double, event_types); % Convert to numeric if needed
    design_matrix = zeros(length(event_types), 3);
    design_matrix(strcmp(event_types,'2'), 1) = 1;  % Column 1: Pleasant stimuli
    design_matrix(strcmp(event_types,'4'), 2) = 1;  % Column 2: Neutral stimuli
    design_matrix(strcmp(event_types,'8'), 3) = 1;  % Column 3: Unpleasant stimuli

    % Count the number of trials per condition
    nTrials_cond1 = sum(design_matrix(:, 1));
    nTrials_cond2 = sum(design_matrix(:, 2));
    nTrials_cond3 = sum(design_matrix(:, 3));
    fprintf('Number of trials for condition 1: %d\n', nTrials_cond1);
    fprintf('Number of trials for condition 2: %d\n', nTrials_cond2);
    fprintf('Number of trials for condition 3: %d\n', nTrials_cond3);

    % Visualize design matrix
    % figure; imagesc(design_matrix);
    % title('Design Matrix'); xlabel('Conditions'); ylabel('Trials');

    % Calculate the correlation matrix of the design matrix
    corr_matrix = corr(design_matrix);
    disp('Correlation matrix of the design matrix:');
    disp(corr_matrix);

    % Perform singular value decomposition (SVD)
    [U, S, V] = svd(design_matrix);
    singular_values = diag(S);
    disp('Singular values of the design matrix:');
    disp(singular_values);

    % Determine the rank using a tolerance
    % tolerance = 1E-7;
    % effective_rank = sum(singular_values > tolerance);
    effective_rank = sum(eig(cov(double(design_matrix')))>1E-7);
    fprintf('Effective rank of the design matrix: %d\n',effective_rank);
    ncol = size(design_matrix, 2);
    if effective_rank == ncol
        disp('The design matrix is full rank.');
    else
        warning('The design matrix is rank deficient.');
    end

    % Convert the 3D EEG data into a 2D matrix for vectorized processing.
    [nChan, nFrames, nTrials] = size(data);
    data_reshaped = reshape(data, nChan*nFrames, nTrials)';
    
    % Run GLM
    switch optimization

            case 'OLS'

            disp("Running GLM using OLS optimization...")

            % Fit the GLM using OLS
            XTX = design_matrix' * design_matrix;
            XTY = design_matrix' * data_reshaped;

            % If the design matrix is rank deficient, handle it appropriately
            % Tikhonov regularization or "ridge regression" (instead of pinv)
            % to handle potential numerical instability. We add a very small
            % value (here 1e-5) to the diagonal elements of the design matrix
            % before computing the pseudoinverse. regularization parameter 1e-5
            % may need to be adjusted.
            effective_rank = sum(eig(cov(double(XTX')))>1E-7);
            if effective_rank < size(XTX, 2)
                warning('Design matrix is rank deficient. Applying Tikhonov regularization.');
                XTX = XTX + 1e-5 * eye(size(XTX));
            end

            % Calculate betas
            betas = pinv(XTX) * XTY;

            % Compute fitted values and residuals
            Yhat = design_matrix * betas;
            residuals = data_reshaped - Yhat;

        case 'IRLS'

            disp("Running GLM using IRLS optimization...")

            % IRLS initial weights
            w = ones(nTrials, 1);

            % Initialize output variables
            betas = nan(nChan * nFrames, size(design_matrix, 2));
            residuals = nan(nTrials, nChan * nFrames);
            % rsquared = nan(nChan * nFrames, 1);
            % fstat = nan(nChan * nFrames, 1);
            % pvals = nan(nChan * nFrames, 1);

            for iter = 1:max_iter

                % Ensure w are a column vector and non-negative and non-zero (replace
                % with a small positive value)
                % w = w(:);
                w(w <= 0) = 1e-5;

                % Fit the GLM with current weights using matrix multiplication
                W = diag(w);
                XTWX = design_matrix' * W * design_matrix;
                XTWY = design_matrix' * W * data_reshaped;

                % If the design matrix is rank deficient, handle it appropriately
                % Tikhonov regularization or "ridge regression" (instead of pinv)
                % to handle potential numerical instability. We add a very small
                % value (here 1e-5) to the diagonal elements of the design matrix
                % before computing the pseudoinverse. regularization parameter 1e-5
                % may need to be adjusted.
                effective_rank = sum(eig(cov(double(XTWX')))>1E-7);
                if effective_rank < size(XTWX, 2)
                    warning('Design matrix is rank deficient. Applying Tikhonov regularization.');
                    XTWX = XTWX + 1e-5 * eye(size(XTWX));
                end

                % betas
                betas = pinv(XTWX) * XTWY;

                % Compute fitted values and residuals
                Yhat = design_matrix * betas;
                residuals = data_reshaped - Yhat;

                % Update weights using Huber's function
                abs_residuals = abs(residuals);
                w_new = ones(size(w));
                for i = 1:nTrials
                    if max(abs_residuals(i, :)) > k
                        w_new(i) = k / max(abs_residuals(i, :));
                    end
                end

                % Check for convergence
                if max(abs(w_new - w)) < tol, break; end

                % Update w for the next iteration
                w = w_new;
            end

        case 'WLS'
            disp("Running GLM using WLS optimization...")

            % Calculate w using a principal components projection for all channels
            disp("Calculating w using principal components projection...")
            [betas_limo, w, ~] = limo_WLS(design_matrix, data_reshaped);

            % Ensure w are a column vector and non-negative and non-zero (replace
            % with a small positive value)
            % w = w(:);
            w(w <= 0) = 1e-5;

            % create diagonal matrix of weights to fit the GLM using matrix multiplication
            W = diag(w);

            % Fit the GLM with current w
            XTWX = design_matrix' * W * design_matrix;
            XTWY = design_matrix' * W * data_reshaped;

            % If the design matrix is rank deficient, handle it appropriately
            % Tikhonov regularization or "ridge regression" (instead of pinv)
            % to handle potential numerical instability. We add a very small
            % value (here 1e-5) to the diagonal elements of the design matrix
            % before computing the pseudoinverse. regularization parameter 1e-5
            % may need to be adjusted.
            effective_rank = sum(eig(cov(double(XTWX')))>1E-7);
            if effective_rank < size(XTWX, 2)
                warning('Design matrix is rank deficient. Applying Tikhonov regularization.');
                XTWX = XTWX + 1e-5 * eye(size(XTWX));
            end

            % Get the betas
            betas = pinv(XTWX) * XTWY;

            % Compute fitted values and residuals
            Yhat = design_matrix * betas;
            residuals = data_reshaped - Yhat;

    end

    % Calculate R2 statistics
    Ymean = mean(data_reshaped);
    SS_tot = sum((data_reshaped - Ymean).^2, 1);
    SS_res = sum(residuals.^2, 1);
    rsquared = 1 - SS_res ./ SS_tot;

    % Adjust degrees of freedom for F-statistics
    df1 = size(design_matrix, 2) - 1;
    df2 = nTrials - size(design_matrix, 2);

    % Calculate F statistics
    fstat = (rsquared / df1) ./ ((1 - rsquared) / df2);

    % Reshape betas, R2, and F statistics back to the original dimensions
    betas = reshape(betas, [size(design_matrix, 2), nChan, nFrames]);
    % betas_limo = reshape(betas_limo, [size(design_matrix, 2), nChan, nFrames]);
    rsquared = reshape(rsquared, [nChan, nFrames]);
    fstat = reshape(fstat, [nChan, nFrames]);

    % Compute p-values (for simplicity, assuming F-distribution)
    pvals = 1 - fcdf(fstat, df1, df2);
    pvals = reshape(pvals, [nChan, nFrames]);   % reshape p-values back to the original dimensions
    % mask = pvals<0.05;
    % F(~mask) = NaN;
    % figure; imagesc(EEG.times,1:size(mask,1),F); colorbar

    % reshape betas to have channels x frames x parameters
    betas = permute(betas, [2, 3, 1]);
    % betas_limo = permute(betas_limo, [2, 3, 1]);

    % Save the results
    % times = EEG.times;
    subject = EEG.subject;
    save(fullfile(filepath, sprintf('%s_GLM_%s.mat', subject, optimization)), 'betas', 'residuals', 'rsquared', 'pvals', 'fstat', 'times');

    % Group level
    BETAS(:,:,:,iSub) = betas;
    F(:,:,iSub) = fstat;
    R2(:,:,iSub) = rsquared;
    P(:,:,iSub) = pvals;

    % ERP raw for comparison
    idx = strcmp(event_types,'2');
    PLEASANT(:,:,iSub) =  trimmean(data(:,:,idx),20,3); % for comparison
    idx = strcmp(event_types,'4');
    NEUTRAL(:,:,iSub) =  trimmean(data(:,:,idx),20,3); % for comparison
    idx = strcmp(event_types,'8');
    UNPLEASANT(:,:,iSub) =  trimmean(data(:,:,idx),20,3); % for comparison

    progressbar(iSub/nSub)

end

% Save GLM outputs
save(fullfile(outDir,sprintf('GLM_%s.mat', optimization)),'times','chanlocs','BETAS','F','R2','P')
save(fullfile(outDir,'RAW_ERP.mat'),'times','chanlocs','PLEASANT','NEUTRAL','UNPLEASANT')

% % Visualize to compare raw ERP with Beta (one subject)
% figure; hold on
% plot(times,squeeze(betas(65,:,1)))
% plot(times,squeeze(betas_limo(65,:,1)))
% plot(times,squeeze(PLEASANT(65,:,iSub)))
% legend('my_WLS','limo_WLS','raw'); axis tight

% figure('color','w')
% subplot(2,1,1); hold on
% idx = strcmp(event_types,'2');
% plot(EEG.times,trimmean(trimmean(data(:,:,idx),20,3),20,1),'LineWidth',1)
% idx = strcmp(event_types,'4');
% plot(EEG.times,trimmean(trimmean(data(:,:,idx),20,3),20,1),'LineWidth',1)
% idx = strcmp(event_types,'8');
% plot(EEG.times,trimmean(trimmean(data(:,:,idx),20,3),20,1),'LineWidth',1)
% legend('pleasant','neutral','unpleasant')
%
% subplot(2,1,2); hold on
% plot(EEG.times,squeeze(trimmean(betas(1,:,:),20,2)),'LineWidth',1)
% plot(EEG.times,squeeze(trimmean(betas(2,:,:),20,2)),'LineWidth',1)
% plot(EEG.times,squeeze(trimmean(betas(3,:,:),20,2)),'LineWidth',1)
% legend('pleasant','neutral','unpleasant')

% % Compare each condition directly
% plotHDI(times, squeeze(trimmean(PLEASANT,20,1)), squeeze(trimmean(BETAS(:,:,1,:),20,1)),'Trimmed mean',0.05,[],'raw','beta')
% subplot(2,1,1); title(sprintf('Pleasant (%s)',optimization))
% plotHDI(times, squeeze(trimmean(NEUTRAL,20,1)), squeeze(trimmean(BETAS(:,:,2,:),1)),'Trimmed mean',0.05,[],'raw','beta')
% subplot(2,1,1); title(sprintf('Neutral (%s)',optimization))
% plotHDI(times, squeeze(trimmean(UNPLEASANT,20,1)), squeeze(trimmean(BETAS(:,:,3,:),1)),'Trimmed mean',0.05,[],'raw','beta')
% subplot(2,1,1); title(sprintf('Unpleasant (%s)',optimization))

% Compare conditions separately to see how it impacts statistical results
elec = 64;
plotHDI(times, squeeze(BETAS(elec,:,1,:)), squeeze(BETAS(elec,:,2,:)),'Trimmed mean',0.05,[],'pleasant','neutral')
subplot(2,1,1); title(sprintf('GLM ERP (one channel - %s)',optimization))
plotHDI(times, squeeze(PLEASANT(elec,:,:)), squeeze(NEUTRAL(elec,:,:)),'Trimmed mean',0.05,[],'pleasant','neutral')
subplot(2,1,1); title('RAW ERP (one channel)')

% Run stats on each to compare at the group level
nPerm = 1000;
[~, neighbormatrix] = get_channelneighbors(chanlocs);

[tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(squeeze(BETAS(:,:,1,:)), squeeze(BETAS(:,:,2,:)), nPerm, 'trimmed mean', 'dpt');
% save(fullfile(outDir,sprintf('stats_pleasant-neutral_GLM_%s.mat',optimization)),'times','chanlocs','tvals','tvals_H0','pvals','pvals_H0','-v7.3')
% load(fullfile(outDir,sprintf('stats_pleasant-neutral_GLM_%s.mat',optimization)))
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 0, 0.05, neighbormatrix);
plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 0); title(sprintf("GLM %s",optimization))
% print(gcf, fullfile(outDir,sprintf('plot_pleasant-neutral_GLM_%s_uncorrected.png',optimization)),'-dpng','-r300');   % 300 dpi .png
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 2, 0.05, neighbormatrix);
plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 2); title(sprintf("GLM %s",optimization))
% print(gcf, fullfile(outDir,sprintf('plot_pleasant-neutral_GLM_%s_cluster.png',optimization)),'-dpng','-r300');   % 300 dpi .png
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 3, 0.05, neighbormatrix);
plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 3); title(sprintf("GLM %s",optimization))
% print(gcf, fullfile(outDir,sprintf('plot_pleasant-neutral_GLM_%s_tfce.png',optimization)),'-dpng','-r300');   % 300 dpi .png

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(PLEASANT, NEUTRAL, nPerm, 'trimmed mean', 'dpt');
% save(fullfile(outDir,'stats_pleasant-neutral_RAW.mat'),'times','chanlocs','tvals','tvals_H0','pvals','pvals_H0','-v7.3')
% % load(fullfile(outDir,'stats_pleasant-neutral_RAW.mat'))
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 0, 0.05, neighbormatrix);
% plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 0); title(sprintf("GLM %s",optimization))
% print(gcf, fullfile(outDir,'plot_pleasant-neutral_RAW_uncorrected.png'),'-dpng','-r300');   % 300 dpi .png
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 2, 0.05, neighbormatrix);
% plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 2); title("RAW ERP")
% print(gcf, fullfile(outDir,'plot_pleasant-neutral_RAW_cluster.png'),'-dpng','-r300');   % 300 dpi .png

% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(squeeze(BETAS(:,:,3,:)), squeeze(BETAS(:,:,2,:)), nPerm, 'trimmed mean', 'dpt');
% save(fullfile(outDir,sprintf('stats_unpleasant-neutral_GLM_%s.mat',optimization)),'times','chanlocs','tvals','tvals_H0','pvals','pvals_H0','-v7.3')
% % load(fullfile(outDir,sprintf('stats_unpleasant-neutral_GLM_%s.mat',optimization)))
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 0, 0.05, neighbormatrix);
% plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 0); title(sprintf("GLM %s",optimization))
% print(gcf, fullfile(outDir,sprintf('plot_unpleasant-neutral_GLM_%s_uncorrected.png',optimization)),'-dpng','-r300');   % 300 dpi .png
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 2, 0.05, neighbormatrix);
% plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 2); title(sprintf("GLM %s",optimization))
% print(gcf, fullfile(outDir,sprintf('plot_unpleasant-neutral_GLM_%s_cluster.png',optimization)),'-dpng','-r300');   % 300 dpi .png
% 
% [tvals,pvals,tvals_H0,pvals_H0] = run_stats_bootstrap(UNPLEASANT, NEUTRAL, nPerm, 'trimmed mean', 'dpt');
% save(fullfile(outDir,'stats_unpleasant-neutral_RAW.mat'),'times','chanlocs','tvals','tvals_H0','pvals','pvals_H0','-v7.3')
% % load(fullfile(outDir,'stats_unpleasant-neutral_RAW.mat'))
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 0, 0.05, neighbormatrix);
% plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 0); title(sprintf("GLM %s",optimization))
% print(gcf, fullfile(outDir,'plot_unpleasant-neutral_RAW_uncorrected.png'),'-dpng','-r300');   % 300 dpi .png
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, 2, 0.05, neighbormatrix);
% plotresults('time', times, tvals, mask, pcorr, 0.05, chanlocs, 2); title("RAW ERP")
% print(gcf, fullfile(outDir,'plot_unpleasant-neutral_RAW_cluster.png'),'-dpng','-r300');   % 300 dpi .png


gong; disp("Done!")


%% GLM method using fitglm

%
% % Parpool with max number of workers
% addons = ver;
% if any(contains({addons.Name}, 'Parallel'))
%     ps = parallel.Settings;
%     fprintf('Parallel computing set to ON. \n')
%     ps.Pool.AutoCreate = true;
%     p = gcp('nocreate');
%     % delete(gcp('nocreate')) % shut down opened parpool
%     if isempty(p) % if not already on, launch it
%         disp('Initiating parrallel computing (all cores and threads -2)...')
%         c = parcluster; % cluster profile
%         % N = feature('numcores');          % only physical cores
%         N = getenv('NUMBER_OF_PROCESSORS'); % all processor (cores + threads)
%         if ischar(N), N = str2double(N); end
%         c.NumWorkers = N-2;  % update cluster profile to include all workers
%         c.parpool();
%     end
% end
%
% filename = sprintf('%s.set',filenames{iSub});
% EEG = pop_loadset(fullfile(filepath,filename));
% data = EEG.data;
%
% % trim time window of interest
% idx = EEG.times>=tlims(1) & EEG.times<=tlims(2);
% data = data(:,idx,:);
% times = EEG.times(idx);
%
% % Remove unecessary conditions?
% event_types = {EEG.event.type};
% event_types = cellfun(@num2str, event_types, 'UniformOutput', false);
% idx = strcmp(event_types,'1');
% warning("removing %g trials from unecessary condition", sum(idx))
% event_types(idx) = [];
% data(:,:,idx) = [];
%
% % % Define one-hot encoded design matrix
% design_matrix = zeros(length(event_types), 3);
% design_matrix(strcmp(event_types,'2'), 1) = 1;  % Column 1: Pleasant stimuli
% design_matrix(strcmp(event_types,'4'), 2) = 1;  % Column 2: Neutral stimuli
% design_matrix(strcmp(event_types,'8'), 3) = 1;  % Column 3: Unpleasant stimuli
%
% % Verify data dimensions
% [nChan, nFrames, nTrials] = size(data);
% nParams = size(design_matrix, 2);
% fprintf('Data dimensions: %d channels, %d time points, %d trials\n', nChan, nFrames, nTrials);
% fprintf('Design matrix dimensions: %d trials, %d parameters \n', size(design_matrix, 1), nParams);
%
% % Initialize output variables
% betas = zeros(nChan, nFrames, nParams+1);
% residuals = zeros(nChan, nFrames, nTrials);
% R2 = zeros(nChan, nFrames);
% F = zeros(nChan, nFrames);
% pvals = zeros(nChan, nFrames,nParams);
%
% % Set convergence criteria for IRLS
% max_iter = 100;
% tol = 1e-6;
% k = 1.345; % Tuning parameter for Huber's function
%
% % Fit the GLM using IRLS
% tstart = tic;
% for iChan = 1%:nChan
%
%     disp('--------------------------------------------------------')
%     fprintf('          SUBJECT %g\n',iSub)
%     disp('--------------------------------------------------------')
%
%     parfor iFrame = 1:nFrames
%
%         % Extract the data for the current electrode and time point
%         Y = squeeze(data(iChan, iFrame, :));
%
%         h = lillietest(Y);
%         figure; histfit(Y);
%         if h, warning("Distribution is not normal. Use different distribution for GLM?"); end
%
%         % Debugging step: Check for NaN or Inf in the response variable
%         if any(isnan(Y)) || any(isinf(Y))
%             error('Response variable contains NaN or Inf at channel %d, time point %d', iChan, iFrame);
%         end
%
%         switch optimization
%
%             case 'IRLS'
%
%                 % Initialize IRLS w
%                 w = ones(nTrials, 1);
%
%                 for iter = 1:max_iter
%                     % Fit the GLM with current w
%                     mdl = fitglm(design_matrix_pca, Y, 'linear', 'Distribution', 'normal', 'w', w);
%
%                     % Compute residuals
%                     residuals_current = mdl.Residuals.Raw;
%
%                     % Update w (using Huber's weight function as an example)
%                     w_new = ones(size(residuals_current));
%                     abs_residuals = abs(residuals_current);
%                     w_new(abs_residuals > k) = k ./ abs_residuals(abs_residuals > k);
%
%                     % Check for convergence
%                     if max(abs(w_new - w)) < tol
%                         break;
%                     end
%
%                     % Update w for the next iteration
%                     w = w_new;
%                 end
%
%             case 'WLS'
%
%                 disp("Calculating weights using principal components projection...")
%                 [betas_limo, w, ~] = limo_WLS(design_matrix, Y);
%                 mdl = fitglm(design_matrix, Y, 'linear', 'Distribution', 'normal', 'w', w);
%         end
%
%         % Store the results
%         betas(iChan, iFrame, :) = mdl.Coefficients.Estimate;
%         residuals(iChan, iFrame, :) = mdl.Residuals.Raw;
%         R2(iChan, iFrame) = mdl.Rsquared.Ordinary;
%         % F(iChan, iFrame) = mdl.Coefficients.FStat(2:end); % Skip the intercept F-statistic
%         pvals(iChan, iFrame,:) = mdl.Coefficients.pValue(2:end); % Skip the intercept p-value
%     end
% end
% fprintf('Time to run GLM on one subject: %g (min)\n',round(toc(tstart)/60,1))
%
% % Save the results
% % save('GLM_IRLS_PCA_results.mat', 'betas', 'residuals', 'R2', 'F', 'pvals');

%% Same but with categorical variables instead of design matrix

% % Define the categorical design matrix
% condition_labels = cell(size(event_types)); % Create a cell array for categorical labels
% condition_labels(event_types == 2) = {'pleasant'};
% condition_labels(event_types == 4) = {'neutral'};
% condition_labels(event_types == 8) = {'unpleasant'};
% condition_labels = categorical(condition_labels);
%
% % Fit the GLM with current w
% tbl = table(condition_labels, Y, w, 'VariableNames', {'Condition', 'Response', 'w'});
% mdl = fitglm(tbl, 'Response ~ Condition', 'Distribution', 'normal', 'w', 'w');
%
% % Compute residuals
% residuals_current = mdl.Residuals.Raw;

