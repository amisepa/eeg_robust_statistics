% Compute data estimator and its high density intervals (HDI) based 
% on a bayesian bootstrap. 
% 
% INPUTS 
%   Y           - 2D matrix, e.g. frames x participants
%   estimator   - 'mean', 'trimmed mean', (default) or 'median'
%   prob_cov    - the probability coverage- default 0.95%
%
% OUTPUT 
%   est         - the data estimator
%   HDI         - the high density interval (HDI)
%   bb          - a vector of (Bayes) bootstraped estimators 
% 
% USAGE 
%   [est, HDI, bb] = computeHDI(Y,estimator,prob_cov)
% 
% EXAMPLE
%   [est, HDI] = compute_HDI(data,'trimmed mean', 0.95)
% 
% Orignal R code for bayesian bootstrap: http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/
% Original R code for HDI: https://github.com/boboppie/kruschke-doing_bayesian_data_analysis/blob/master/1e/HDIofMCMC.R
%
%  Cedric Cannard, 2021 (adapted from LIMO-EEG)

function [est, HDI, bb] = compute_HDI(Y,estimator,prob_cov,method)

% number of bootstrap samples 
nboot = 1000;

if nargin == 2
    prob_cov = 0.95;
elseif nargin == 1
    prob_cov = 0.95;
    estimator = 'Trimmed mean';
    method = 'bayesian';  % 'frequentist' or 'bayesian'
end

% compute the estimator
if strcmpi(estimator,'mean')
    est = mean(Y,2,'omitnan');
elseif strcmpi(estimator,'trimmed mean')
    est = trimmean(Y,20,2);
elseif strcmpi(estimator,'median')
    est = median(Y,2,'omitnan');
end

% Sample with replacement. sampling = number of observations (e.g. participants)
n = size(Y,2); 
bb = nan(size(Y,1),nboot);
for boot = 1:nboot % bootstrap loop

    if strcmpi(method,'frequentist')
        theta    = exprnd(1,[n,1]); % generate weights using exponential distribution (related to Bayesian approach but not fully conform)
        weights  = theta ./ repmat(sum(theta,1),n,1);

    elseif strcmpi(method,'bayesian')
        theta    = gamrnd(ones(n,1),1); % generate weights using Gamma distribution
        weights  = theta / sum(theta); % Normalize to sum to 1 to align with Dirichlet distribution properties
        
    else
        error('Method to estimate confidence intervals must be defined')
    end
    
    resample = (datasample(Y',n,'Replace',true,'Weights',weights))';
    
    % compute the estimator
    if strcmpi(estimator,'mean')
        bb(:,boot) = mean(resample,2,'omitnan');
    elseif strcmpi(estimator,'trimmed Mean')
        bb(:,boot) = trimmean(resample,20,2);
    elseif strcmpi(estimator,'median')
        bb(:,boot) = median(resample,2,'omitnan');
    end
end

sorted_data   = sort(bb,2); % sort bootstrap estimates
upper_centile = floor(prob_cov*size(sorted_data,2)); % upper bound
nCIs          = size(sorted_data,2) - upper_centile;
HDI           = zeros(2,size(Y,1));

ci       = 1:nCIs;
ciWidth  = sorted_data(:,ci+upper_centile) - sorted_data(:,ci); % all centile distances
[~,J]    = min(ciWidth,[],2);
r        = size(sorted_data,1);
I        = (1:r)';
index    = I+r.*(J-1); % linear index
HDI(1,:) = sorted_data(index);
index    = I+r.*(J+upper_centile-1); % linear index
HDI(2,:) = sorted_data(index);
