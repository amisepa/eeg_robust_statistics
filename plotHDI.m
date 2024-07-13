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

function plotHDI(xAxis, data1, data2, method, a, h, data1Name, data2Name)
 
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
color1 = [0, 0.4470, 0.7410];           % blue
color2 = [0.8500, 0.3250, 0.0980];      % red
color3 = [0.4660, 0.6740, 0.1880];      % green
% cb = cbrewer2('qual', 'Set3', 12, 'pchip');
% color1 = cb(5,:);  % 5=blue, 4=red, 1=green, 2=yellow
% color2 = cb(2,:);  % red
% color3 = cb(4,:);  % green

% figure('color','w'); 
% subplot(2,1,1) 
hold on

% Data1 (mean + 95% HDI)
p1 = plot(xAxis,est1,'LineWidth',2,'Color', color1);
patch([xAxis fliplr(xAxis)], [HDI1(1,:) fliplr(HDI1(2,:))], ...
    color1,'FaceAlpha',.3,'EdgeColor',color1,'EdgeAlpha',0.9);
% set(gca,'FontSize',12,'layer','top'); 

% Data2 (mean + 95% HDI)
p2 = plot(xAxis, est2,'LineWidth',2,'Color', color2);
patch([xAxis fliplr(xAxis)], [HDI2(1,:) fliplr(HDI2(2,:))], ...
    color2,'FaceAlpha',.3,'EdgeColor',color2,'EdgeAlpha',0.9);
grid off; axis tight; hold on; box on
ylabel('Potential (uV)')
title(sprintf('%s + %g%% HDI',method,(1-a)*100)); 

% % Plot difference (mean + 95% HDI)
% subplot(2,1,2)
% plot(xAxis, est3,'LineWidth',2,'Color', color3);
% patch([xAxis fliplr(xAxis)], [HDI3(1,:) fliplr(HDI3(2,:))], ...
%     color3,'FaceAlpha',.6,'EdgeColor',color3,'EdgeAlpha',0.9);
% grid off; axis tight; box on
% ylabel('Difference (uV)')

% % Add dash line to mark the null hypothesis
% hold on; plot([xAxis(1) xAxis(end)], [0 0],'k--','LineWidth',1) % thick dash line highlighting H0
% ylabel('Difference','FontSize',11,'FontWeight','bold')
% xlabel("Time (ms)",'FontSize',11,'FontWeight','bold')
% % title(sprintf('Difference (%s + 95% HDI)',method)); 

% Plot significance bar at the bottom
% if isempty(h)
%     for i = 1:length(xAxis)
%         if HDI3(1,i)<0 && HDI3(2,i)<0 || HDI3(1,i)>0 && HDI3(2,i)>0
%             h(i) = true;
%         else
%             h(i) = false;
%         end
%     end
% end
if ~isempty(h)
    plotSigBar(h, xAxis);
end

% legend([p1, p2], {data1Name,data2Name}, 'Location','SouthWest'); 
legend([p1, p2], {data1Name,data2Name}); 

set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');