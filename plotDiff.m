% Plots 2 times series, their 95% CI and significance bars at the bottom
% from h vector (FDR-corrected p-values). If method is not precised,
% 10% trimmed mean is used. 
% 
% Usage:
%   - plotDiff(xAxis, data1, data2, method1, method2, h, data1Name, data2Name);
%   - plotDiff(freqs, power1, power2, 'mean', 'CI', [], 'condition 1','condition 2');
% 
% Data must be 2-D. Values in column 1 and subjects in column 2 (e.g., freqs x subjects)
% Conditions can have different numbers of subjects/trials in the last dimension.
% 
% Cedric Cannard, 2021

function [data3_mean, data3_CI] = plotDiff(xAxis, data1, data2, method1, method2, h, data1Name, data2Name)

if size(xAxis,1) > size(xAxis,2)
    xAxis = xAxis';
end

if exist('h', 'var') && ~isempty(h)
    sigBars = true;
else
    sigBars = false;
end

color1 = [0, 0.4470, 0.7410];       % blue
color2 = [0.8500, 0.3250, 0.0980];  % red
color3 = [0.4660, 0.6740, 0.1880];  % green

% Variable 1
n1 = size(data1,2);
if strcmpi(method1, 'mean')
    data1_mean = mean(data1,2,'omitnan');
else
    data1_mean = trimmean(data1,20,2);
end
if strcmpi(method2,'SD')
    SD = std(data1,[],2,'omitnan');  % standard deviation
    data1_CI(1,:) = data1_mean' + SD';
    data1_CI(2,:) = data1_mean' - SD';
elseif strcmpi(method2,'SE')
    SE = std(data1,[],2,'omitnan') ./ sqrt(n1);  % standard error
    data1_CI(1,:) = data1_mean' + SE';
    data1_CI(2,:) = data1_mean' - SE';
elseif strcmpi(method2,'CI')
    SE = std(data1,[],2,'omitnan') ./ sqrt(n1);  % standard error
    tscore = tinv([.025 .975],n1-1);  % t-score
    data1_CI(1,:) = data1_mean' + tscore(2).*SE';
    data1_CI(2,:) = data1_mean' + tscore(1).*SE';
else
    error("Method 2 not recognized. Should be 'CI', 'SD', or 'SE'.")
end

% Variable 2
n2 = size(data2,2);
if strcmpi(method1, 'mean')
    data2_mean = mean(data2,2,'omitnan');
else
    data2_mean = trimmean(data2,20,2);
end
if strcmpi(method2,'SD')
    SD = std(data2,[],2,'omitnan');  % standard deviation
    data2_CI(1,:) = data2_mean' + SD';
    data2_CI(2,:) = data2_mean' - SD';
elseif strcmpi(method2,'SE')
    SE = std(data2,[],2,'omitnan') ./ sqrt(n2);  % standard error
    data2_CI(1,:) = data2_mean' + SE';
    data2_CI(2,:) = data2_mean' - SE';
elseif strcmpi(method2,'CI')
    SE = std(data2,[],2,'omitnan') ./ sqrt(n2);  % standard error
    tscore = tinv([.025 .975],n2-1);  % t-score
    data2_CI(1,:) = data2_mean' + tscore(2).*SE';
    data2_CI(2,:) = data2_mean' + tscore(1).*SE';
end

% Difference - always compute as difference of means for unequal n
if n1 == n2
    % Paired comparison
    data_diff = data1 - data2;
    n_diff = n1;
    if strcmpi(method1, 'mean')
        data3_mean = mean(data_diff,2,'omitnan');
    else
        data3_mean = trimmean(data_diff,20,2);
    end
    
    if strcmpi(method2,'SD')
        SD = std(data_diff,[],2,'omitnan');
        data3_CI(1,:) = data3_mean' + SD';
        data3_CI(2,:) = data3_mean' - SD';
    elseif strcmpi(method2,'SE')
        SE = std(data_diff,[],2,'omitnan') ./ sqrt(n_diff);
        data3_CI(1,:) = data3_mean' + SE';
        data3_CI(2,:) = data3_mean' - SE';
    elseif strcmpi(method2,'CI')
        SE = std(data_diff,[],2,'omitnan') ./ sqrt(n_diff);
        tscore = tinv([.025 .975],n_diff-1);
        data3_CI(1,:) = data3_mean' + tscore(2).*SE';
        data3_CI(2,:) = data3_mean' + tscore(1).*SE';
    end
else
    % Independent comparison
    warning('Conditions have different sample sizes (%d vs %d). Computing unpaired difference.', n1, n2);
    data3_mean = data1_mean - data2_mean;
    
    % Welch-Satterthwaite approximation for unequal variances
    var1 = var(data1,[],2,'omitnan');
    var2 = var(data2,[],2,'omitnan');
    pooled_SE = sqrt(var1/n1 + var2/n2);
    
    if strcmpi(method2,'SD')
        pooled_SD = sqrt(var1 + var2);  % Combined SD
        data3_CI(1,:) = data3_mean' + pooled_SD';
        data3_CI(2,:) = data3_mean' - pooled_SD';
    elseif strcmpi(method2,'SE')
        data3_CI(1,:) = data3_mean' + pooled_SE';
        data3_CI(2,:) = data3_mean' - pooled_SE';
    elseif strcmpi(method2,'CI')
        % Welch-Satterthwaite degrees of freedom
        df = (var1/n1 + var2/n2).^2 ./ ((var1/n1).^2/(n1-1) + (var2/n2).^2/(n2-1));
        % Compute t-score for each frequency point
        tscore_upper = arrayfun(@(d) tinv(.975, floor(d)), df);
        tscore_lower = arrayfun(@(d) tinv(.025, floor(d)), df);
        data3_CI(1,:) = data3_mean' + tscore_upper'.*pooled_SE';
        data3_CI(2,:) = data3_mean' + tscore_lower'.*pooled_SE';
    end
end

% Create figure
figure('Color','w');
subplot(2,1,1)
hold on

% Plot variable 1 (mean + CI)
p1 = plot(xAxis,data1_mean,'LineWidth',2,'Color', color1);
patch([xAxis fliplr(xAxis)], [data1_CI(1,:) fliplr(data1_CI(2,:))], ...
    color1,'FaceAlpha',.3,'EdgeColor',color1,'EdgeAlpha',.9);

% Plot variable 2 (mean + CI)
p2 = plot(xAxis,data2_mean,'LineWidth',2,'Color', color2);
patch([xAxis fliplr(xAxis)], [data2_CI(1,:) fliplr(data2_CI(2,:))], ...
    color2,'FaceAlpha',.3,'EdgeColor',color2,'EdgeAlpha',.9);

set(gca,'FontSize',12,'layer','top');  grid off; axis tight; box on
ylims = ylim;
plot([0 0], [ylims(1) ylims(2)],'k--','LineWidth',.5) % Add dash line to mark time 0

% ylabel('Power (db)')
ylabel('Amplitude (μV)')

% Plot difference (mean + CI)
subplot(2,1,2); hold on
plot(xAxis, data3_mean,'LineWidth',2,'Color', color3);
patch([xAxis fliplr(xAxis)], [data3_CI(1,:) fliplr(data3_CI(2,:))], ...
    color3,'FaceAlpha',.3,'EdgeColor',color3,'EdgeAlpha',0.9);

% Add dash line to mark the null hypothesis
plot([xAxis(1) xAxis(end)], [0 0],'k--','LineWidth',.5)
ylims = ylim;
plot([0 0], [ylims(1) ylims(2)],'k--','LineWidth',.5) % Add dash line to mark time 0

ylabel('Difference (μV)','FontSize',11,'FontWeight','bold')
% xlabel("Frequency (Hz)",'FontSize',11,'FontWeight','bold')
xlabel("Time (ms)",'FontSize',11,'FontWeight','bold')

grid off; axis tight; box on

% Title with CT parameters
if n1 ~= n2
    title(sprintf('%s + 95%% %s (n1=%d, n2=%d)',method1,method2,n1,n2));
else
    title(sprintf('%s + 95%% %s (n=%d)',method1,method2,n1));
end

% Plot significance bar at the bottom
if sigBars
    plotSigBar(h, xAxis);
end

legend([p1, p2], {sprintf('%s (n=%d)',data1Name,n1), sprintf('%s (n=%d)',data2Name,n2)}, ...
    'Orientation','vertical','Location','best'); 

set(findall(gcf,'type','axes'),'fontSize',12,'fontweight','bold');