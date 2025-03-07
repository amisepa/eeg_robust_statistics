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

function [est, HDI, bb] = compute_HDI(Y,estimator,prob_cov)

% number of bootstrap samples 
nboot = 1000;

if nargin == 2
    prob_cov = 0.95;
elseif nargin == 1
    prob_cov = 0.95;
    estimator = 'Trimmed mean';
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

    %     theta    = exprnd(1,[n,1]); % generate weights using exponential distribution (related to Bayesian approach but not fully conform)
    %     weights  = theta ./ repmat(sum(theta,1),n,1);
    theta    = gamrnd(ones(n,1),1); % generate weights using Gamma distribution
    weights  = theta / sum(theta); % Normalize to sum to 1 to align with Dirichlet distribution properties
            
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

% % Compute high-density intervals (HDIs)
% sorted_data   = sort(bb,2); % sort bootstrap estimates
% upper_centile = floor(prob_cov*size(sorted_data,2)); % upper bound
% nCIs          = size(sorted_data,2) - upper_centile;
% HDI           = zeros(2,size(Y,1));
% ci       = 1:nCIs;
% ciWidth  = sorted_data(:,ci+upper_centile) - sorted_data(:,ci); % all centile distances
% [~,J]    = min(ciWidth,[],2);
% r        = size(sorted_data,1);
% I        = (1:r)';
% index    = I+r.*(J-1); % linear index
% HDI(1,:) = sorted_data(index);
% index    = I+r.*(J+upper_centile-1); % linear index
% HDI(2,:) = sorted_data(index);

% % Compute quantile interval following guidelines from Etz et al. (2024)
% % HDI can lead to different conclusions based on the parameterization of 
% % the model, whereas quantile intervals are derived directly from the 
% % cumulative distribution function (from the sorted posterior bootstrap 
% % samples) and are based on probabilities rather than densities, ensuring 
% % coherence across different parameterizations.
% % REF:  Etz A, Chávez de la Peña AF, Baroja L, Medriano K, Vandekerckhove J. 
% %       The HDI + ROPE decision rule is logically incoherent but we can fix it. 
% %       Psychol Methods. 2024 May 23. doi: 10.1037/met0000660
% %       https://pubmed.ncbi.nlm.nih.gov/38780591/
HDI = zeros(2, size(Y, 1));
for i = 1:size(Y, 1)
    sampleVec = bb(i, :);
    alp = (1 - prob_cov) / 2;
    lim = quantile(sampleVec, [alp, 1 - alp]);
    HDI(:, i) = lim';
end



