% Computes and plots central tendency (mean, trimmed mean, or median), and 
% the 95% high density intervals (HDI) computed with a 1000-iterations
% Bayesian bootstrap. Adapted from LIMO-EEG.
% 
% INPUTS: 
%   xAxis       - vector for x-axis (e.g. frames for ERP, frequencies for power
%               spectra).
%   data1       - 2D data for group or condition 1 (e.g. frames x participants)
%   data2       - 2D data for group or condition 2 (e.g. frames x participants)
%   estimator   - 'mean', 'trimmed Mean', 'median'
%   grp         - 'dpt' (dependent, paired), or 'idpt' (independent, unpaired)
%   a           - alpha probability coverage (default = .05)
%   h           - index of significant (true) or nonsignificant (false) values, 
%               to plot significance bars at the bottom. Should be the size of xAxis 
% 
% USAGE:
%   plotHDI(xAxis,data1,data2,estimator,a,h,data1Name,data2Name)
% 
% EXAMPLE:
%   plotHDI(freqs,data1,data2,'trimmed mean',.05,h,'condition1','condition2'); 
% 
% Cedric Cannard 2021

function plotHDI3(xAxis, data1, data2, method, a, h, data1Name, data2Name)
 
if size(xAxis,2) < size(xAxis,1)
    xAxis = xAxis';
end

% Estimator 95% high-density intervals (HDI)
fprintf('Computing estimator and high-density interval (HDI) for data 1... \n')
[est1, HDI1] = compute_HDI(data1, method, 1-a);
fprintf('Computing estimator and high-density interval (HDI) for data 2... \n')
[est2, HDI2] = compute_HDI(data2, method, 1-a);

% Difference
fprintf('Computing estimator and high-density interval (HDI) for the difference... \n')
if size(data1,2) == size(data2,2)    
    [est3, HDI3] = compute_HDI(data1-data2,method,1-a);   
else
    warning('The two datasets have a different number of participants/trials, using inpendent method')
    est3 = est1-est2;
    HDI3 = HDI1-HDI2;
end

% Colors
color3 = [0.4660, 0.6740, 0.1880];      % green

% Plot difference between data1 and data2 (mean + 95% HDI)
plot(xAxis, est3,'LineWidth',1,'Color', color3);
patch([xAxis fliplr(xAxis)], [HDI3(1,:) fliplr(HDI3(2,:))], ...
    color3,'FaceAlpha',.3,'EdgeColor',color3,'EdgeAlpha',0.9);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; box on
ylabel('Difference','FontSize',11,'FontWeight','bold')

if ~isempty(h)
    plotSigBar(h, xAxis);
end

set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');

