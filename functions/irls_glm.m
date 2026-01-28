function [B, Res, w] = irls_glm(X, Y, w_init)
% IRLS_GLM - Iteratively Reweighted Least Squares for robust regression
%
% Implements robust regression using Huber weights with iterative reweighting.
% This is the standard IRLS algorithm used in robust statistics.
%
% SYNTAX:
%   [B, Res, w] = irls_glm(X, Y)
%   [B, Res, w] = irls_glm(X, Y, w_init)
%
% INPUTS:
%   X       - [nObs x nP] Design matrix
%   Y       - [nObs x nVox] Response matrix (vectorized data)
%   w_init  - [nObs x 1] Initial weights (optional, default: all ones)
%
% OUTPUTS:
%   B       - [nP x nVox] Regression coefficients
%   Res     - [nObs x nVox] Residuals
%   w       - [nObs x 1] Final robust weights
%
% ALGORITHM:
%   Uses Huber M-estimator with iterative reweighting until convergence.
%   Huber tuning constant k = 1.345 provides 95% efficiency at Gaussian
%   while protecting against outliers.
%
% REFERENCES:
%   Huber, P. J. (1981). Robust Statistics. Wiley.
%   Holland, P. W., & Welsch, R. E. (1977). Robust regression using 
%   iteratively reweighted least-squares. Communications in Statistics, 6(9).

[nObs, nP] = size(X);
nVox = size(Y, 2);

% Initialize weights
if nargin < 3 || isempty(w_init)
    w = ones(nObs, 1);
else
    w = w_init;
end

% Huber tuning constant (standard value for 95% efficiency)
k = 1.345;

% Convergence parameters
max_iter = 50;
tol = 1e-6;

% Initial OLS fit
B_old = pinv(X' * X) * (X' * Y);

for iter = 1:max_iter
    % Compute weighted regression
    W = spdiags(w, 0, nObs, nObs);
    B = pinv(X' * W * X) * (X' * W * Y);
    
    % Check convergence
    if iter > 1
        rel_change = max(abs(B(:) - B_old(:))) / (max(abs(B_old(:))) + eps);
        if rel_change < tol
            break;
        end
    end
    
    B_old = B;
    
    % Compute residuals
    Res = Y - X * B;
    
    % Update weights using Huber function
    % Compute MAD-based scale estimate for each column
    for iVox = 1:nVox
        res_col = Res(:, iVox);
        
        % Robust scale estimate: Median Absolute Deviation
        mad_scale = median(abs(res_col - median(res_col))) / 0.6745;
        
        % Avoid division by zero
        if mad_scale < eps
            mad_scale = 1;
        end
        
        % Standardized residuals
        std_res = res_col / mad_scale;
        
        % Huber weights: w = 1 if |r| <= k, w = k/|r| if |r| > k
        w_col = ones(nObs, 1);
        outliers = abs(std_res) > k;
        w_col(outliers) = k ./ abs(std_res(outliers));
        
        % Average weights across all voxels (for consistency)
        if iVox == 1
            w_new = w_col;
        else
            w_new = w_new + w_col;
        end
    end
    
    % Average weights across voxels
    w = w_new / nVox;
    
    % Ensure weights are in valid range
    w = max(w, 0.01);  % Minimum weight to avoid numerical issues
    w = min(w, 1.0);   % Maximum weight
end

% Final residuals
Res = Y - X * B;

end