function w = compute_wls_weights(X, Y, weightType)
% COMPUTE_WLS_WEIGHTS - Compute robust weights for Weighted Least Squares
%
% Computes observation weights based on leverage and residuals for
% robust regression. Implements three standard weighting schemes.
%
% SYNTAX:
%   w = compute_wls_weights(X, Y, weightType)
%
% INPUTS:
%   X          - [nObs x nP] Design matrix
%   Y          - [nObs x nVox] Response matrix (vectorized data)
%   weightType - String specifying weight function:
%                'PCP'   - Principal Component Projection (default)
%                'Huber' - Huber M-estimator weights
%                'Tukey' - Tukey bisquare weights (more aggressive)
%
% OUTPUT:
%   w          - [nObs x 1] Observation weights (0 to 1)
%
% ALGORITHM:
%   1. Compute initial OLS fit
%   2. Calculate leverage (hat matrix diagonal)
%   3. Compute standardized residuals
%   4. Apply chosen weight function
%   5. Down-weight high-leverage outliers
%
% WEIGHT FUNCTIONS:
%   PCP (Principal Component Projection):
%     - Based on Mahalanobis distance in predictor space
%     - Good for multivariate outliers in X
%   
%   Huber (k=1.345):
%     - w = 1 if |r| <= k
%     - w = k/|r| if |r| > k
%     - 95% efficiency at Gaussian, moderate outlier protection
%   
%   Tukey Bisquare (c=4.685):
%     - w = [1-(r/c)^2]^2 if |r| <= c
%     - w = 0 if |r| > c
%     - More aggressive outlier rejection
%
% REFERENCES:
%   Rousseeuw, P. J., & Leroy, A. M. (1987). Robust Regression and 
%   Outlier Detection. Wiley.
%   
%   Maronna, R. A., et al. (2006). Robust Statistics: Theory and Methods. 
%   Wiley.

[nObs, nP] = size(X);
nVox = size(Y, 2);

% Initial OLS fit
B = pinv(X' * X) * (X' * Y);
Res = Y - X * B;

% Compute leverage (hat matrix diagonal)
H = X * pinv(X' * X) * X';
h = diag(H);

% Initialize weights
w = ones(nObs, 1);

switch upper(weightType)
    case 'PCP'
        % Principal Component Projection weights
        % Based on Mahalanobis distance
        
        % Robust center and covariance estimation
        X_centered = bsxfun(@minus, X, median(X, 1));
        
        % Robust covariance (using MAD)
        S = zeros(nP, nP);
        for i = 1:nP
            for j = i:nP
                mad_i = median(abs(X(:,i) - median(X(:,i)))) / 0.6745;
                mad_j = median(abs(X(:,j) - median(X(:,j)))) / 0.6745;
                if mad_i > eps && mad_j > eps
                    S(i,j) = median(X_centered(:,i) .* X_centered(:,j)) / (mad_i * mad_j);
                    S(j,i) = S(i,j);
                end
            end
        end
        
        % Add small ridge for numerical stability
        S = S + 1e-6 * eye(nP);
        
        % Mahalanobis distance
        try
            S_inv = pinv(S);
            mahal_dist = sqrt(sum((X_centered * S_inv) .* X_centered, 2));
        catch
            mahal_dist = sqrt(sum(X_centered.^2, 2));
        end
        
        % Weight based on chi-square quantile
        chi2_cutoff = sqrt(chi2inv(0.975, nP));
        w = min(chi2_cutoff ./ (mahal_dist + eps), 1);
        
    case 'HUBER'
        % Huber M-estimator weights
        k = 1.345;  % 95% efficiency
        
        % Compute robust residual weights
        w_res = zeros(nObs, 1);
        for iVox = 1:nVox
            res_col = Res(:, iVox);
            
            % Robust scale
            mad_scale = median(abs(res_col - median(res_col))) / 0.6745;
            if mad_scale < eps
                mad_scale = 1;
            end
            
            % Standardized residuals
            std_res = res_col / mad_scale;
            
            % Huber weights
            w_col = ones(nObs, 1);
            outliers = abs(std_res) > k;
            w_col(outliers) = k ./ abs(std_res(outliers));
            
            w_res = w_res + w_col;
        end
        w = w_res / nVox;
        
    case 'TUKEY'
        % Tukey bisquare weights (more aggressive)
        c = 4.685;  % 95% efficiency
        
        % Compute robust residual weights
        w_res = zeros(nObs, 1);
        for iVox = 1:nVox
            res_col = Res(:, iVox);
            
            % Robust scale
            mad_scale = median(abs(res_col - median(res_col))) / 0.6745;
            if mad_scale < eps
                mad_scale = 1;
            end
            
            % Standardized residuals
            std_res = res_col / mad_scale;
            
            % Tukey bisquare weights
            w_col = zeros(nObs, 1);
            inliers = abs(std_res) <= c;
            w_col(inliers) = (1 - (std_res(inliers) / c).^2).^2;
            % Outliers get weight 0
            
            w_res = w_res + w_col;
        end
        w = w_res / nVox;
        
    otherwise
        error('Unknown weight type: %s. Use PCP, Huber, or Tukey.', weightType);
end

% Down-weight high-leverage points
leverage_cutoff = 2 * nP / nObs;  % Standard cutoff
high_leverage = h > leverage_cutoff;
w(high_leverage) = w(high_leverage) .* (leverage_cutoff ./ h(high_leverage));

% Ensure weights are in valid range
w = max(w, 0.01);  % Minimum weight
w = min(w, 1.0);   % Maximum weight

% Normalize weights to sum to nObs (optional, for interpretability)
w = w * (nObs / sum(w));

end