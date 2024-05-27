function [r_alpha,t_alpha,p_alpha,v,vv,vvv] = MC_corrpval(n,p,method,alphav,pairs,D)

% function to compute the alpha quantile estimate of the distribution of
% minimal p-values under the null of correlations in a n*p matrix with null
% covariance but variance D (I by default)
%
% FORMAT p_alpha = MC_corrpval(n,p,method,alphav,pairs,D)
%
% INPUT n the number of observations
%       p the number of variables
%       method can be 'Pearson', 'Spearman', 'Skipped Pearson', 'Skipped Spearman'
%       alphav the alpha value
%       pairs a m*2 matrix of variables to correlate (optional)
%       D the variance of each variable (optional)
%
% p_alpha the alpha quantile estimate of the distribution of
%         minimal p-values
%
%
% Cyril Pernet v3 - Novembre 2017
% ---------------------------------------------------
%  Copyright (C) Corr_toolbox 2017

% deal with inputs
if nargin == 0
    help MC_corrpval
elseif nargin < 2
    error('at least 2 inputs requested see help MC_corrpval');
end
if ~exist('pairs','var') || isempty(pairs)
    pairs = nchoosek(1:p,2);
end
if ~exist('alphav','var') || isempty(alphav)
    alphav = 5/100;
end

% Generate the variance
SIGMA = eye(p);
if exist('D','var')
    if length(D) ~= p
        error('the vector D of variance must be of the same size as the number of variables p')
    else
        SIGMA(SIGMA==1) = D;
    end
end

% Parpool parameters
ps = parallel.Settings;
ps.Pool.AutoCreate = true;
tmp = gcp('nocreate');
% delete(gcp('nocreate')) % shut down opened parpool
if isempty(tmp) % if not already on, launch it
    disp('Initiating parrallel computing (all cores and threads -1)...')
    c = parcluster; % cluster profile
    % N = feature('numcores');          % only physical cores
    N = getenv('NUMBER_OF_PROCESSORS'); % all processor (cores + threads)
    if ischar(N), N = str2double(N); end
    c.NumWorkers = N-2;  % update cluster profile to include all workers
    c.parpool();
end

% Run the Monte Carlo simulation and keep smallest p -alues
nboots = 1000;
v = NaN(1,nboots);
vv = v; 
vvv = vv;
parfor iBoot = 1:nboots
    fprintf('Running Monte Carlo %g\n',iBoot)

    MVN = mvnrnd(zeros(1,p),SIGMA,n); % a multivariate normal distribution

    if strcmp(method,'Pearson')
        [r,t,pval] = Pearson(MVN,pairs);
    elseif strcmp(method,'Spearman')
        [r,t,pval] = Spearman(MVN,pairs);
    elseif strcmp(method,'Skipped Pearson')
        [r,t,~,pval] = skipped_Pearson(MVN,pairs);
    elseif strcmp(method,'Skipped Spearman')
        [r,t,~,pval] = skipped_Spearman(MVN,pairs);
    end
    
    v(iBoot)   = max(r);
    vv(iBoot)  = max(t);
    vvv(iBoot) = min(pval);
end

% Get the Harell-Davis estimate of the alpha quantile
n = length(v);
r_alpha = nan(length(alphav),1);
t_alpha = nan(length(alphav),1);
p_alpha = nan(length(alphav),1);
for l = 1:length(alphav)
    q       = 1-alphav(l);  % for r/t use 1-alphav
    m1      = (n+1).*q;
    m2      = (n+1).*(1-q);
    vec     = 1:n;
    w       = betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
    y       = sort(v);
    r_alpha(l,:) = sum(w(:).*y(:));
    y       = sort(vv);
    t_alpha(l,:) = sum(w(:).*y(:));
    
    q       = alphav(l);  % for p values use alphav
    m1      = (n+1).*q;
    m2      = (n+1).*(1-q);
    vec     = 1:n;
    w       = betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
    y       = sort(vvv);
    p_alpha(l,:) = sum(w(:).*y(:));
end


