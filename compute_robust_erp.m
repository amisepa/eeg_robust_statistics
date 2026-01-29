function ERP_robust = compute_robust_erp(trials, method)
% COMPUTE_ROBUST_ERP - Compute robust ERP using weighted averaging
%
% Applies robust weighting to trials before averaging to reduce impact
% of outlier trials on the subject-level ERP.
%
% SYNTAX:
%   ERP_robust = compute_robust_erp(trials, method)
%
% INPUTS:
%   trials - [nChan x nTime x nTrials] Single-subject trial data
%   method - Weighting method:
%            'mean'   - Standard arithmetic mean (no weighting)
%            'median' - Median (ultimate robustness, but lower statistical power)
%            'huber'  - Huber M-estimator weights (default; good balance, 95% efficiency, moderate outlier protection)
%            'tukey'  - Tukey bisquare weights (more aggressive, completely removes extreme outliers)
%
% OUTPUT:
%   ERP_robust - [nChan x nTime] Robust averaged ERP
%
% ALGORITHM:
%   1. For each channel-time point, compute trial weights based on
%      deviation from robust center
%   2. Down-weight outlier trials
%   3. Compute weighted average
%
% EXAMPLE:
%   % For each subject and condition
%   for iSub = 1:nSub
%       crash_robust(:,:,iSub) = compute_robust_erp(crash_trials(:,:,:,iSub), 'huber');
%       nocrash_robust(:,:,iSub) = compute_robust_erp(nocrash_trials(:,:,:,iSub), 'huber');
%   end
% 
% Cedric Cannard, Jan 2026

[nChan, nTime, nTrials] = size(trials);

% Initialize output
ERP_robust = zeros(nChan, nTime);

switch lower(method)
    case 'mean'
        % Standard mean (no robustness)
        ERP_robust = mean(trials, 3);
        
    case 'median'
        % Median (most robust, but less efficient)
        ERP_robust = median(trials, 3);
        
    case 'huber'
        % Huber M-estimator weighted mean
        k = 1.345;  % 95% efficiency
        
        for iChan = 1:nChan
            for iTime = 1:nTime
                x = squeeze(trials(iChan, iTime, :));
                
                % Robust center and scale
                center = median(x);
                mad_scale = median(abs(x - center)) / 0.6745;
                
                if mad_scale < eps
                    % No variability, just use mean
                    ERP_robust(iChan, iTime) = mean(x);
                else
                    % Standardized deviations
                    std_dev = (x - center) / mad_scale;
                    
                    % Huber weights
                    w = ones(nTrials, 1);
                    outliers = abs(std_dev) > k;
                    w(outliers) = k ./ abs(std_dev(outliers));
                    
                    % Weighted mean
                    ERP_robust(iChan, iTime) = sum(w .* x) / sum(w);
                end
            end
        end
        
    case 'tukey'
        % Tukey bisquare (more aggressive outlier rejection)
        c = 4.685;  % 95% efficiency
        
        for iChan = 1:nChan
            for iTime = 1:nTime
                x = squeeze(trials(iChan, iTime, :));
                
                % Robust center and scale
                center = median(x);
                mad_scale = median(abs(x - center)) / 0.6745;
                
                if mad_scale < eps
                    % No variability, just use mean
                    ERP_robust(iChan, iTime) = mean(x);
                else
                    % Standardized deviations
                    std_dev = (x - center) / mad_scale;
                    
                    % Tukey bisquare weights
                    w = zeros(nTrials, 1);
                    inliers = abs(std_dev) <= c;
                    w(inliers) = (1 - (std_dev(inliers) / c).^2).^2;
                    % Outliers get weight 0
                    
                    % Check if we rejected all trials
                    if sum(w) < eps
                        % Fall back to median
                        ERP_robust(iChan, iTime) = center;
                    else
                        % Weighted mean
                        ERP_robust(iChan, iTime) = sum(w .* x) / sum(w);
                    end
                end
            end
        end
        
    otherwise
        error('Unknown method: %s. Use mean, median, huber, or tukey.', method);
end

end