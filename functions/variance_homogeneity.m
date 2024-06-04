function [h,CI] = variance_homogeneity(x,y,condition,vis)

% Compares variances using a 95% percentile bootstrap CI.
%
% FORMAT:  [h,CI] = variance_homogeneity(x,y)
%          [h,CI] = variance_homogeneity(x,y,condition)
%
% INPUTS:  x and y are two vectors of the same length
%          condition = 0/1 (default=0) to condition x and y on each other
%
% OUTPUTS: h indicates if the data have the same variance (0) or not (1)
%          CI is the 95% confidence interval of the difference between variances
%
% see also CONDITIONAL.

% Cyril Pernet v1
% ---------------------------------
%  Copyright (C) Corr_toolbox 2012
  
if nargin == 2
    condition = 0;
    vis = 1;
end

if size(x)~=size(y)
    error('X and Y must have the same size')
end

% computes
nboot = 1000;
nm = length(x);
if nm < 40
    l=6; u=593;
elseif nm >= 40 && nm < 80
    l=7; u=592;
elseif nm >= 80 && nm < 180
    l=10; u=589;
elseif nm >= 180 && nm < 250
    l=13; u=586;
elseif nm >= 250
    l=15; u=584; 
end
 
% boostrap
table = randi(nm,nm,nboot);
Diff = nan(1, nboot);
for boot = 1:nboot
    % resample 
    a = x(table(:,boot)); 
    b = y(table(:,boot)); 
    if condition == 1
        [~,variances] = conditional(a(:),b(:));
        Diff(boot) = variances(1) - variances(2);
    else
        Diff(boot) = var(a) - var(b);
    end
end

% 95% confidence interval (CI)
Diff = sort(Diff);
CI = [Diff(l+1) Diff(u)];
if sum(isnan(Diff)) ~=0
    adj_nboot = nboot - sum(isnan(Diff));
    adj_l = (5/100*adj_nboot)/2;
    adj_u = adj_nboot - adj_l;
    CI = round([Diff(adj_l+1) Diff(adj_u)],1);
end

% Decision
if CI(1) < 0 && CI(2) > 0
    h = 0;
else
    h = 1;
end

% plot
if vis
    if condition == 1
        figure('Name','Test on conditional variances','Color','w');
    else
        figure('Name','Test on variances','Color','w');
    end

    k = round(1 + log2(nboot));
    [n,p] = hist(Diff,k); 
    bar(p,n,1,'FaceColor',[0.5 0.5 1]);
    grid on; axis tight; 
    ylabel('Frequency','Fontsize',14); hold on
    plot(repmat(CI(1), max(hist(Diff,k)),1), 1:max(hist(Diff,k)),'r','LineWidth',4);
    plot(repmat(CI(2), max(hist(Diff,k)),1), 1:max(hist(Diff,k)),'r','LineWidth',4);

    if h == 0
       title(sprintf('Data are homoscedastic (95%% CI: [%g %g])',CI(1),CI(2)),'Fontsize',12);
    else
       title(sprintf('Data are heteroscedastic (95%% CI: [%g %g])',CI(1),CI(2)),'Fontsize',12);
    end
    xlabel('Variance differences between X and Y','Fontsize',12);

    box on; set(gca,'Fontsize',12)
end

