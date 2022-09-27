function [mask,p_val,max_th] = limo_max_correction(tvals,tboot,p,fig)

% correction for multiple testing using the max stat value
% since the type 1 error is the prob to make at least one error
% we can control the prob the make an error for the biggest value
% note this works for both bootstrapped data/tfce trandformed data under H0  
%
% FORMAT mask = limo_max_correction(tvals,tboot,p,fig)
%
% INPUT
% tvals     = 2D matrix of observed values 
% tboot = 3D matrix of T^2 or F values for data bootstrapped under H0
% p     = threshold to apply e.g. 5/100
% fig   = 1/0 to plot the maximum stat under H0
%
% OUTPUT
% mask is a binary matrix of the same size as M corresponding to a threshold
%      p corrected for multiple comparisons
% p_val are the p-values obtained via the matrix bootM (non-corrected)
% max_th is the threshold controlling the type 1 FWER
% 
% Cedric Cannard

% check inputs 
% tvals       = varargin{1};
% tboot   = varargin{2};
if nargin < 3
    p   = 0.05;
    fig = 0;
end
% elseif nargin == 3
%     p   = varargin{3};
%     fig = 0;
% elseif nargin == 4
%     p   = varargin{3};
%     fig = varargin{4};
% end
% clear varargin

[a,b,nboot] = size(tboot);
if any(size(tvals)~=[a b])
    error('dimension error: matrices of observed and bootstrap values are different')
end

% collect highest value for each boot
parfor boot=1:nboot
    data = squeeze(tboot(:,:,boot));
    maxM(boot) = max(data(:)); 
end

% get threshold
maxM(maxM==Inf) = [];
sortmaxM        = sort(maxM); 
nboot           = length(sortmaxM);
U               = round((1-p).*nboot);
max_th          = sortmaxM(U);
mask            = squeeze(tvals) >= max_th;
fprintf('max threshold %g\n',max_th)

% get the equivalent bootstrapped p-value
smalest_pval = 1/nboot;
for row =1:a
    for column=1:b
        tmp = sum(tvals(row,column) >= sortmaxM) / nboot;
        p_val(row,column) = 1-tmp;
        if p_val(row,column) == 0; p_val(row,column) = smalest_pval; end
    end
end

%% figure
if sum(mask(:)) == 0
    fig = 1 ; 
end

if fig == 1
    figure('Name','Correction by max: results under H0')
    plot(sortmaxM,'LineWidth',3); grid on; hold on; 
    
    plot(min(find(sortmaxM==max_th)),max_th,'r*','LineWidth',5)
    txt = ['bootstrap threashold ' num2str(max_th) '\rightarrow'];
    text(min(find(sortmaxM==max_th)),max_th,txt,'FontSize',12,'HorizontalAlignment','right');
    
    [val,loc]=min(abs(sortmaxM-max(tvals(:)))); 
    plot(loc,max(tvals(:)),'r*','LineWidth',5)
    txt = ['maximum observed: ' num2str(max(tvals(:))) '\rightarrow'];
    text(loc,max(tvals(:)),txt,'FontSize',12,'HorizontalAlignment','right');
    
    title('Maxima under H0','FontSize',12)
    xlabel('sorted bootstrap iterations','FontSize',12); 
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end
