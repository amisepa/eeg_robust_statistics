% Plots 2 times series, their 95% CI and significance bars at the bottom
% from h vector (FDR-corrected p-values). If method is not precised,
% 10% trimmed mean is used. 
% 
% Usage:
%   - plotDiff(xAxis, data1, data2, method1, method2, h, data1Name, data2Name);
%   - plotDiff(freqs, power1, power2, 'mean', 'CI', [], 'condition 1','condition 2');
% 
% Data must be 2-D. Values in column 1 and subjects in column 2 (e.g.,freqs x subjects)
% 
% Cedric Cannard, 2021

function plotDiff(xAxis, data1, data2, method1, method2, h, data1Name, data2Name)

if size(xAxis,1) > size(xAxis,2)
    xAxis = xAxis';
end

if exist('h', 'var') && ~isempty(h)
    sigBars = true;
else
    sigBars = false;
end

color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.4660, 0.6740, 0.1880];      % green

% Variable 1
n = size(data1,2);
if strcmpi(method1, 'mean')
    data1_mean = mean(data1,2,'omitnan');
else
    data1_mean = trimmean(data1,20,2);
end
if strcmpi(method2,'SD')
    SD = std(data1,[],2,'omitnan');  % standard deviation
    data1_CI(1,:) = data1_mean + SD;
    data1_CI(2,:) = data1_mean - SD;
elseif strcmpi(method2,'SE')
    SE = std(data1,[],2,'omitnan') ./ sqrt(n)';  % standard error
    data1_CI(1,:) = data1_mean + SE;
    data1_CI(2,:) = data1_mean - SE;
elseif strcmpi(method2,'CI')
    SE = std(data1,[],2,'omitnan') ./ sqrt(n)';  % standard error
    tscore = tinv([.025 .975],n-1);  % t-score
    data1_CI = data1_mean' + (-tscore.*SE)'; % 95% confidence interval
end

% Variable 2
n = size(data2,2);
if strcmpi(method1, 'mean')
    data2_mean = mean(data2,2,'omitnan');
else
    data2_mean = trimmean(data2,20,2);
end
if strcmpi(method2,'SD')
    SD = std(data2,[],2,'omitnan');  % standard deviation
    data2_CI(1,:) = data2_mean + SD;
    data2_CI(2,:) = data2_mean - SD;
elseif strcmpi(method2,'SE')
    SE = std(data2,[],2,'omitnan') ./ sqrt(n)';  % standard error
    data2_CI(1,:) = data2_mean + SE;
    data2_CI(2,:) = data2_mean - SE;
elseif strcmpi(method2,'CI')
    SE = std(data2,[],2,'omitnan') ./ sqrt(n)';  % standard error
    tscore = tinv([.025 .975],n-1);  % t-score
    data2_CI = data2_mean' + (-tscore.*SE)'; % 95% confidence interval
end

% Difference
if strcmpi(method1, 'mean')
    data3_mean = mean(data1-data2,2,'omitnan');
else
    data3_mean = trimmean(data1-data2,20,2);
end
if strcmpi(method2,'SD')
    SD = std(data1-data2,[],2,'omitnan');  % standard deviation
    data3_CI(1,:) = data3_mean + SD;
    data3_CI(2,:) = data3_mean - SD;
elseif strcmpi(method2,'SE')
    SE = std(data1-data2,[],2,'omitnan') ./ sqrt(n)';  % standard error
    data3_CI(1,:) = data3_mean + SE;
    data3_CI(2,:) = data3_mean - SE;
elseif strcmpi(method2,'CI')
    SE = std(data1-data2,[],2,'omitnan') ./ sqrt(n)';  % standard error
    tscore = tinv([.025 .975],n-1);  % t-score
    data3_CI = data3_mean' + (-tscore.*SE)'; % 95% confidence interval
end


figure; set(gcf,'Color','w');

subplot(2,1,1); hold on;

% Plot variable 1 (mean + CI)
p1 = plot(xAxis,data1_mean,'LineWidth',2,'Color', color1);
patch([xAxis fliplr(xAxis)], [data1_CI(1,:) fliplr(data1_CI(2,:))], ...
    color1,'FaceAlpha',.3,'EdgeColor',color1,'EdgeAlpha',.9);
set(gca,'FontSize',11,'layer','top'); 


% Plot variable 2 (mean + CI)
p2 = plot(xAxis,data2_mean,'LineWidth',2,'Color', color2);
patch([xAxis fliplr(xAxis)], [data2_CI(1,:) fliplr(data2_CI(2,:))], ...
    color2,'FaceAlpha',.3,'EdgeColor',color2,'EdgeAlpha',.9);
set(gca,'FontSize',12,'layer','top'); 
grid off; axis tight; hold on; box on
ylabel('Power (db)')
title(sprintf('%s + 95%% %s',method1,method2)); 

% Plot difference (mean + CI)
subplot(2,1,2)
plot(xAxis, data3_mean,'LineWidth',2,'Color', color3);
patch([xAxis fliplr(xAxis)], [data3_CI(1,:) fliplr(data3_CI(2,:))], ...
    color3,'FaceAlpha',.6,'EdgeColor',color3,'EdgeAlpha',0.9);
grid off; axis tight; box on
ylabel('Difference (uV)')

% Add dash line to mark the null hypothesis
hold on; plot([xAxis(1) xAxis(end)], [0 0],'k--','LineWidth',1) % thick dash line highlighting H0
ylabel('Difference','FontSize',11,'FontWeight','bold')
xlabel("Frequency (Hz)",'FontSize',11,'FontWeight','bold')

% Plot significance bar at the bottom
if sigBars
    plotSigBar(h, xAxis);
end

legend([p1, p2], {data1Name,data2Name}, 'Orientation','vertical','Location','SouthWest'); 

% grid on; 
set(findall(gcf,'type','axes'),'fontSize',11,'fontweight','bold');