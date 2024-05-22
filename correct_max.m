%% Correction for multiple testing using the Max method
% Calculates the probability to make at least one error with the biggest
% value (i.e. type 1 error) using null data (H0). Works for TFCE
% transformed data. 
%
% USAGE:
%       mask = correct_max(tvals, tvals_H0, p_thresh, fig)
%
% INPUT:
%       tvals     - 2D matrix of observed values 
%       tvals_H0  - 3D matrix of T^2 or F values for data bootstrapped under H0
%       pthresh   - threshold to apply e.g. 5/100
%       fig       - 1/0 to plot the maximum stat under H0
%
% OUTPUT
%       mask    - a binary matrix of the same size as M corresponding to a threshold
%                   p corrected for multiple comparisons
%       pvals   - p-values obtained via the matrix bootM (non-corrected)
%       max_th  - threshold controlling the type 1 FWER
% 
% Cedric Cannard

function [mask,pvals,max_th] = correct_max(tvals,tvals_H0,pthresh)

% check inputs 
if nargin < 3
    pthresh   = 0.05;
end

[a,b,nboot] = size(tvals_H0);
if any(size(tvals)~=[a b])
    error('Dimension error: matrices of observed and bootstrap values are different')
end

% collect highest value for each boot
for boot = 1:nboot
    data = squeeze(tvals_H0(:,:,boot));
    maxM(boot) = max(data(:)); 
end

% get threshold
maxM(maxM==Inf) = [];
sortmaxM        = sort(maxM); 
nboot           = length(sortmaxM);
U               = round((1-pthresh).*nboot);
max_th          = sortmaxM(U);
mask            = squeeze(tvals) >= max_th;
fprintf('Max threshold = %g\n', max_th)

% Get the equivalent bootstrapped p-value
smalest_pval = 1/nboot;
pvals = nan(length(a),length(b));
for row = 1:a
    for column = 1:b
        tmp = sum(tvals(row,column) >= sortmaxM) / nboot;
        pvals(row,column) = 1-tmp;
        if pvals(row,column) == 0
            pvals(row,column) = smalest_pval; 
        end
    end
end

% figure
if sum(mask(:)) == 0
    figure('Name','Results under H0 after max-correction')
    plot(sortmaxM,'LineWidth',3); grid on; hold on; 
    
    plot(find(sortmaxM==max_th, 1 ),max_th,'r*','LineWidth',5)
    txt = ['Bootstrap threshold ' num2str(max_th) '\rightarrow'];
    text(find(sortmaxM==max_th, 1 ),max_th,txt,'FontSize',12,'HorizontalAlignment','right');
    
    [~,loc] = min(abs(sortmaxM-max(tvals(:)))); 
    plot(loc,max(tvals(:)),'r*','LineWidth',5)
    txt = ['Maximum observed: ' num2str(max(tvals(:))) '\rightarrow'];
    text(loc,max(tvals(:)),txt,'FontSize',12,'HorizontalAlignment','right');
    
    title('Maxima under H0','FontSize',12)
    xlabel('Sorted bootstrap iterations','FontSize',12); 
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end
