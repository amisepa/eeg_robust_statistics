%% Run GLM on EEG data
% 
% Vectorized approach to significantly increase
% computation cost. Vectorized code performs matrix operations directly,
% which are generally more efficient and less error-prone than iterating
% through individual elements. Matrix multiplication and other linear
% algebra operations in MATLAB are highly optimized.
%
% The core computations, such as the design matrix multiplication (XTWX, XTWY),
% fitting the GLM (pinv(XTWX) * XTWY), and computing residuals (data_reshaped - Yhat),
% are the same in both approaches.
%
% Some considerations:
% Ensure the numerical precision and stability are maintained, especially
% with operations like pseudoinverse (pinv). In some cases, the use of
% regularization or other stabilization techniques may be necessary.
% Ensure that edge cases (e.g., handling of zero or negative weights, rank
% deficiency) are managed consistently in both approaches.
% 
% Optimization options: OLS, IRLS, WLS
% 
% WLS options
%   Pernet's principal components projection (PCP): https://www.humanbrainmapping.org/files/Aperture%20Neuro/Accepted%20Works%20PDF/7_51_Cyril_ElectroEncephaloGraphy_robust_statistical_linear.pdf
%   Hubert M-estimator: Combines squared error for small residuals and absolute error for large residuals, offering a balance between efficiency and robustness.
%   Tukey's Biweight: Reduces the influence of outliers more aggressively by setting large residuals to zero weight.
% 
% EXAMPLE: 
%   [betas,rsquared,fstat,pvals] = run_glm(data,times,tlims,events,optimization,weight_method)
% 
% Copyright (C) - Cedric Cannard, May 2024

function [betas,rsquared,fstat,pvals,times] = run_glm(data,times,tlims,events,optimization,weight_method)

if strcmp(optimization,'OLS') || strcmp(optimization,'IRLS')
    weight_method = [];
end

% Add an extra dimension if needed (one channel squeezed)
if ndims(data) == 2 
    data = reshape(data, 1, size(data, 1), size(data, 2));
end

% Convergence criteria for IRLS optimization
if strcmp(optimization,'IRLS')
    max_iter = 100;
    tol = 1e-6;
    k = 1.345; % Tuning parameter for Huber's function
end

if strcmp(optimization,'WLS')
    if isempty(weight_method) || strcmpi(weight_method,'PCP')
        weight_method = 1;          % PCP by default
    elseif strcmpi(weight_method,'hubert')
        weight_method = 2;          % Hubert M-estimator
        k = 1.345; % Tuning parameter for Huber's function
    elseif strcmpi(weight_method,'tukey')
        weight_method = 3;          % Tukey's Biweight
    else
        error("You chose an unknown weighting method. Weighting method must be 'PCP', 'Hubert', or 'Tukey")
    end
end

% trim time window of interest
if ~isempty(tlims)
    idx = times>=tlims(1) & times<=tlims(2);
    data = data(:,idx,:);
    times = times(idx);
end

% Data events
event_types = {events.type};
event_types = cellfun(@num2str, event_types, 'UniformOutput', false);  % ensure events are strings
if length(event_types)~=size(data,3)
    warning("Event structure has different number of trials than EEG data. Trying to correct.")
    lats = [events.latency];
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

        % Compute betas, fitted values, residuals
        betas = pinv(XTX) * XTY;
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

            % Compute betas, fitted values and residuals
            betas = pinv(XTWX) * XTWY;
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
        
        % Using Cyril Pernet's principal components projection (PCP)
        if weight_method == 1
            disp("Calculating weights using principal components projection (PCP)...")
            [~, w, ~] = limo_WLS(design_matrix, data_reshaped);

        % Using Huber's M-estimator 
        elseif weight_method == 2
            w = compute_huber_weights(design_matrix, data_reshaped, 1.345); % Huber's M-estimator

        % Tukey's Biweight
        elseif weight_method == 3
            w = compute_tukey_weights(design_matrix, data_reshaped, 4.685); % Tukey's Biweight function
        end

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

        % Compute betas, fitted values and residuals
        betas = pinv(XTWX) * XTWY;
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

disp("GLM done.")

end



%% SUBFUNCTIONS

function w = compute_huber_weights(design_matrix, data_reshaped, k)
    % Compute initial residuals using OLS
    XTX = design_matrix' * design_matrix;
    XTY = design_matrix' * data_reshaped;
    betas = pinv(XTX) * XTY;
    Yhat = design_matrix * betas;
    residuals = data_reshaped - Yhat;

    % Compute absolute residuals
    abs_residuals = abs(residuals);

    % Initialize weights
    w = ones(size(residuals, 1), 1);

    % Update weights using Huber's function
    for i = 1:size(residuals, 1)
        if max(abs_residuals(i, :)) > k
            w(i) = k / max(abs_residuals(i, :));
        end
    end
end


function w = compute_tukey_weights(design_matrix, data_reshaped, c)
    % Compute initial residuals using OLS
    XTX = design_matrix' * design_matrix;
    XTY = design_matrix' * data_reshaped;
    betas = pinv(XTX) * XTY;
    Yhat = design_matrix * betas;
    residuals = data_reshaped - Yhat;

    % Compute absolute residuals
    abs_residuals = abs(residuals);
    mad_res = mad(residuals, 1); % Median absolute deviation

    % Initialize weights
    w = ones(size(residuals, 1), 1);

    % Update weights using Tukey's Biweight function
    for i = 1:size(residuals, 1)
        scaled_res = abs_residuals(i, :) / (c * mad_res);
        mask = scaled_res < 1;
        w(i) = sum((1 - scaled_res(mask).^2).^2);
    end
end
