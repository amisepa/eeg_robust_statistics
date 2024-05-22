% Computes and plots trimmed means and 95% High density intervals (HDI) 
% using Bayesian bootstrap. Adapted from LIMO-EEG tools. 
% 
% plotHDI2 only does upper plot, no difference plot below.
% 
% Usage:
% plotHDI2(xAxis, data1, data2, type, grp, alphaLevel, h, data1Name, data2Name, plot1Title)% 
% type: 'Mean', 'Trimmed Mean'
% alphaLevel: predetermined threshold probability (e.g. 0.05)
% h: true or false index for each xAxis value, to plot significance bars at the bottom
% 
% Cedric Cannard 2021


function plotHDI2(xAxis, data1, data2, type, alphaLevel, h, data1Name, data2Name)
 
fprintf('Computing Bayesian bootstrap for HDI estimation... \n')

% COMPUTE HDI
fprintf('Computing estimator and high-density interval (HDI) for data 1... \n')
[est1, HDI1] = compute_HDI(data1, type, 1-alphaLevel);
fprintf('Computing estimator and high-density interval (HDI) for data 2... \n')
[est2, HDI2] = compute_HDI(data2, type, 1-alphaLevel);

% Colors
color1 = [0, 0.4470, 0.7410];           %blue
color2 = [0.8500, 0.3250, 0.0980];      %red

% PLOT MEAN AND 95% HDI OF EACH CONDITION/GROUP
% figure; set(gcf,'Color','w');
p1 = plot(xAxis, est1,'LineWidth',2,'Color',color1); hold on;
fillhandle = fill([xAxis fliplr(xAxis)],[HDI1(1,:) fliplr(HDI1(2,:))], color1);
set(fillhandle,'EdgeColor', color1,'FaceAlpha',0.3,'EdgeAlpha',0.8);
p2 = plot(xAxis, est2,'LineWidth',2,'Color', color2);
fillhandle = fill([xAxis fliplr(xAxis)],[HDI2(1,:) fliplr(HDI2(2,:))], color2);
set(fillhandle,'EdgeColor', color2,'FaceAlpha',0.3,'EdgeAlpha',0.8);
set(gca,'FontSize',12,'layer','top'); 
grid on; axis tight; box on
% ylabel('PSD (dB)','FontSize',12);
% xlabel('Frequency (Hz)','FontSize',12);
% ylabel('Alpha asymmetry amplitude','FontSize',12);
% title([plot1Title ' (' type ' and 95% HDI)']); 
% title(plot1Title); 
% hPlots = flip(findall(gcf,'Type','Line')); % flipped, because the lines are found in reverse order of appearance.
% legend(hPlots, {data1Name,data2Name}, 'Orientation','vertical'); 

% Plot significance bar at the bottom
% if isempty(h)
%     diff = est1-est2;
%     for i = 1:length(xAxis)
%         if diff(1,i)<0 && diff(2,i)<0 || diff(1,i)>0 && diff(2,i)>0
% 
%             h(i) = true;
%         else
%             h(i) = false;
%         end
%     end
% end
if ~isempty(h)
    plotSigBar(h, xAxis);
end

legend([p1, p2], {data1Name,data2Name},'location','southwest'); 
    

