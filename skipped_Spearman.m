%SKIPPED_SPEARMAN Spearman correlation after bivariate outlier removal
%
% performs a robust Spearman correlation on data cleaned up for bivariate outliers,
% that is after finding the central point in the distribution using the mid covariance
% determinant, orthogonal distances are computed to this point, and any data outside the
% bound defined by the ideal estimator of the interquartile range is removed.
%
% FORMAT: [rp,tp,CI,pval,outid,h]=skipped_Spearman(X,pairs,method,alphav,p_alpha,vis);
%
% INPUTS:  X is a matrix and correlations between all pairs (default) are computed
%          pairs (optional) is a n*2 matrix of pairs of columns to correlate
%          method (optional) is 'ECP' or 'Hochberg' (only for n>60)
%          alphav (optional, 5% by default) is the requested alpha level
%          p_alpha (optional) the critical p_value to correct for multiple
%                  comparisons (see MC_corrpval)
%          vis (optional) visualize correlation, slope, outliers, and 95% CI
%
% OUTPUTS: rs is the Spearman correlation
%          ts is the T value associated to the skipped correlation
%          CI is the robust confidence interval of r computed by bootstrapping
%             the cleaned-up data set and taking the alphav centile values
%          pval is the p value associated to t
%          outid is the index of bivariate outliers
%          h is the significance after correction for multiple comparisons
%
% This code rely on the mid covariance determinant as implemented in LIBRA
% - Verboven, S., Hubert, M. (2005), LIBRA: a MATLAB Library for Robust Analysis,
% Chemometrics and Intelligent Laboratory Systems, 75, 127-136.
% - Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
% Journal of the American Statistical Association, Vol. 79, pp. 871-881.
%
% The quantile of observations whose covariance is minimized is
% floor((n+size(X,2)*2+1)/2)),
% i.e. ((number of observations + number of variables*2)+1) / 2,
% thus for a correlation this is floor(n/2 + 5/2).
%
% The method for multiple comparisons correction is described in
% Rand R. Wilcox, Guillaume A. Rousselet & Cyril R. Pernet (2018)
% Improved methods for making inferences about multiple skipped correlations,
% Journal of Statistical Computation and Simulation, 88:16, 3116-3131,
% DOI: 10.1080/00949655.2018.1501051
%
% See also MCDCOV, IDEALF.
% 
%  Copyright (C) Corr_toolbox 2017

function [rs,ts,CI,pval,outid,h] = skipped_Spearman(varargin)

% check paths to subfunctions
tmp  = which('mcdcov.m');
if isempty(tmp)
    tmp = fileparts(which("skipped_Spearman"));
    addpath(fullfile(tmp,"LIBRA"))
    addpath(fullfile(tmp,"LIBRA"))
end

% Data
if nargin < 1
    help skipped_Spearman; return
else
    % data
    x = varargin{1};
    [n,p] = size(x);
end

% Get user inputs
for inputs = 2:nargin
    if inputs == 2
        pairs = varargin{inputs};
    elseif inputs == 3
        method = varargin{inputs};
    elseif inputs == 4
        alphav = varargin{inputs};
    elseif inputs == 5
        p_alpha = varargin{inputs};
    elseif inputs == 6
        vis = varargin{inputs};
    end
end

% Defaults and checks
if ~exist('pairs','var') || isempty(pairs) 
    pairs = nchoosek(1:p,2);
end
if size(pairs,2) ~= 2
    pairs = pairs';
end
if ~exist('method','var') || isempty(method)
    method  = 'ECP';
end
if sum(strcmpi(method,{'ECP','Hochberg'})) == 0
    error('Unknown MMC method selected, see help skipped_Spearman')
end
% if strcmp(method,'Hochberg') && n<60 || strcmp(method,'Hochberg') && n<60 && alphav == 5/100
%     error('Hochberg is only valid for n>60 and aplha 5%')
% end
if ~exist('alphav','var') || isempty(alphav)
    alphav  = 5/100;
end
if  ~exist('vis','var') || isempty(vis)
    vis = 0;  
end
nboot = 1000; % number of bootrsaps to compute 

% Create a table of resamples
if nargout > 2
    boot_index = 1;
    while boot_index <= nboot
        resample = randi(n,n,1);
        if length(unique(resample)) > 3 % at least 3 different data points
            boostrap_sampling(:,boot_index) = resample;
            boot_index = boot_index +1;
        else
            error('At least 3 different data points are needed. this can happen when matrix is not transposed the right way (col x row)')
        end
    end
    lower_bound = round((alphav*nboot)/2);
    upper_bound = nboot - lower_bound;
end

% Now for each pair to test, get the observed and boostrapped r and t
% values, then derive the p value from the bootstrap (and hboot and CI if
% requested)

% place holders
rs = nan(size(pairs,1),1);
for outputs = 2:nargout
    if outputs == 2
        ts = nan(size(pairs,1),1);
    elseif outputs == 3
        CI = nan(size(pairs,1),2);
    elseif outputs == 4
        pval = nan(size(pairs,1),1);
    elseif outputs == 5
        outid = cell(size(pairs,1),1);
    end
end

% loop for each pair to test
for row = 1:size(pairs,1)

    % select relevant columns
    X = [x(:,pairs(row,1)) x(:,pairs(row,2))];

    % get the bivariate outliers
    flag = bivariate_outliers(X);
    vec = 1:n;
    if sum(flag)==0
        outid{row}=[];
    else
        flag=(flag>=1);
        outid{row}=vec(flag);
    end
    keep=vec(~flag); % the vector of data to keep

    % Spearman correlation on cleaned data
    xrank = tiedrank(X(keep,1),0); yrank = tiedrank(X(keep,2),0);
    rs(row) = sum(detrend(xrank,'constant').*detrend(yrank,'constant')) ./ ...
        (sum(detrend(xrank,'constant').^2).*sum(detrend(yrank,'constant').^2)).^(1/2);
    ts(row) = rs(row)*sqrt((n-2)/(1-rs(row).^2));

    % Same for bootstrap samples
    if nargout > 2
        % fprintf('computing p values by bootstrapping data, pair %g %g\n',pairs(row,1),pairs(row,2))
        parfor boot = 1:nboot
            Xb = X(boostrap_sampling(:,boot),:);
            xrank = tiedrank(Xb(keep,1),0); yrank = tiedrank(Xb(keep,2),0);
            r(boot) = sum(detrend(xrank,'constant').*detrend(yrank,'constant')) ./ ...
                (sum(detrend(xrank,'constant').^2).*sum(detrend(yrank,'constant').^2)).^(1/2);
        end

        % 95% CI
        r = sort(r);
        CI(row,:) = [r(lower_bound) r(upper_bound)];

        % p-value
        Q = sum(r<0)/nboot;
        pval(row) = 2*min([Q 1-Q]);
    end
end

% Adjust for multiple comparisons (type 1 error)
if nargout == 6
    % N<60
    if strcmp(method,'ECP')
        if exist('p_alpha','var')
            h = pval < p_alpha;
        else
            disp('ECP method requested, computing p alpha ... (takes a while)')
            p_alpha = MC_corrpval(n,p,'Skipped Spearman',alphav,pairs);
            h = pval < p_alpha;
        end
        
    % N>60
    elseif strcmp(method,'Hochberg')
        [sorted_pval,index] = sort(pval,'descend');
        [~,reversed_index]=sort(index);
        k = 1; sig = 0; h = zeros(1,length(pval));
        while sig == 0
            if sorted_pval(k) <= alphav/k
                h(k:end) = 1; sig = 1;
            else
                k = k+1;
                if k == length(h)
                    break
                end
            end
        end
        h = h(reversed_index);

        % quick clean-up of individual p-values
        pval(pval==0) = 1/nboot;
    end
end

% Visualization
if vis

    x = X(:,1);
    y = X(:,2);

    % figure('Name','Skipped Spearman correlation','Color','w');

    % scatter plot
    scatter(x,y,150,'.','MarkerFaceColor',[0 0.4470 0.7410]);
    hold on

    % slope
    hh = lsline; 
    set(hh,'Color',[0.6350 0.0780 0.1840],'LineWidth',2);

    % Ellipse
    % [XEmin, YEmin] = ellipse(a{column},b{column});
    % plot(real(XEmin), real(YEmin),'LineWidth',1);
    % MM = [min(XEmin) max(XEmin) min(YEmin) max(YEmin)];

    % Outliers
    scatter(x(outid{:}),y(outid{:}),150,'.','MarkerFaceColor',[0.9290 0.6940 0.1250]);

    % scale axis
    MM2 = [min(x) max(x) min(y) max(y)];
    MM = MM2; 
    A = floor(min([MM(:,1);MM2(:,1)]) - min([MM(:,1);MM2(:,1)])*0.01);
    boot = ceil(max([MM(:,2);MM2(:,2)]) + max([MM(:,2);MM2(:,2)])*0.01);
    C = floor(min([MM(:,3);MM2(:,3)]) - min([MM(:,3);MM2(:,3)])*0.01);
    D = ceil(max([MM(:,4);MM2(:,4)]) + max([MM(:,4);MM2(:,4)])*0.01);
    axis([A boot C D]);

    title(sprintf('Spearman rho = %g, p = %g, CI = [%g %g]', ...
        round(rs,2), round(pval,3), round(CI(1),2),round(CI(2),2)),'Fontsize',12);  
    xlabel('X','Fontsize',12); ylabel('Y','Fontsize',12);
    box on; set(gca,'Fontsize',12)

    % % 95% CI
    % y1 = refline(CIpslope(1),CIpintercept(1)); set(y1,'Color','r');
    % y2 = refline(CIpslope(2),CIpintercept(2)); set(y2,'Color','r');
    % y1 = get(y1); y2 = get(y2);
    % xpoints=[y1.XData(1):y1.XData(2),y2.XData(2):-1:y2.XData(1)];
    % step1 = y1.YData(2)-y1.YData(1); step1 = step1 / (y1.XData(2)-y1.XData(1));
    % step2 = y2.YData(2)-y2.YData(1); step2 = step2 / (y2.XData(2)-y2.XData(1));
    % filled=[y1.YData(1):step1:y1.YData(2),y2.YData(2):-step2:y2.YData(1)];
    % hold on; fillhandle=fill(xpoints,filled,[1 0 0]);
    % set(fillhandle,'EdgeColor',[1 0 0],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color

end

disp('Skipped Spearman done')
