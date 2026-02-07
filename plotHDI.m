% Computes and plots central tendency (mean, trimmed mean, or median), and 
% the 95% high density intervals (HDI) computed with a 1000-iterations
% Bayesian bootstrap. Adapted from LIMO-EEG.
% 
% INPUTS: 
%   xAxis       - vector for x-axis (e.g. times for ERP, frequencies for power)
%   data1       - 2D data for group or condition 1 (e.g. frames x participants)
%   data2       - 2D data for group or condition 2 (e.g. frames x participants)
%   estimator   - 'mean', 'trimmed mean', 'median'
%   h           - index of significant (true) or nonsignificant (false) values, 
%                 to plot significance bars at the bottom. Should be the size of xAxis 
% 
% USAGE:
%   plotHDI(xAxis,data1,data2,estimator,a,h,data1Name,data2Name)
% 
% EXAMPLE:
%   plotHDI(times, data1, data2, 'mean', h, 'condition1', 'condition2'); 
% 
% Cedric Cannard 2021

function plotHDI(xAxis, data1, data2, estimator, h, data1Name, data2Name)
 
if size(xAxis,2) < size(xAxis,1)
    xAxis = xAxis';
end

if exist('h', 'var') && ~isempty(h)
    sigBars = true;
else
    sigBars = false;
end

a = 0.05; % alpha probability coverage (default = .05);

% Get sample sizes
n1 = size(data1,2);
n2 = size(data2,2);

% Estimator 95% high-density intervals (HDI)
fprintf('Computing estimator and quantile intervals for data 1... \n')
[est1, HDI1] = compute_HDI(data1, estimator, 1-a);
fprintf('Computing estimator and quantile intervals for data 2... \n')
[est2, HDI2] = compute_HDI(data2, estimator, 1-a);

% Difference
fprintf('Computing estimator and quantile intervals for the difference... \n')
if n1 == n2
    % Paired comparison - compute HDI of paired differences
    [est3, HDI3] = compute_HDI(data1-data2, estimator, 1-a);   
else
    % Independent comparison - use bootstrap to estimate difference distribution
    warning('The two datasets have a different number of participants/trials (%d vs %d), using independent bootstrap', n1, n2);
    
    % Point estimate of difference
    est3 = est1 - est2;
    
    % Bootstrap the difference distribution
    nboot = 1000;
    npts = size(data1, 1);  % number of time/frequency points
    boot_diff = zeros(npts, nboot);
    
    for b = 1:nboot
        % Resample each group independently
        idx1 = randi(n1, n1, 1);
        idx2 = randi(n2, n2, 1);
        
        % Compute estimator for this bootstrap sample
        if strcmpi(estimator, 'mean')
            boot1 = mean(data1(:, idx1), 2, 'omitnan');
            boot2 = mean(data2(:, idx2), 2, 'omitnan');
        elseif strcmpi(estimator, 'trimmed mean')
            boot1 = trimmean(data1(:, idx1), 20, 2);
            boot2 = trimmean(data2(:, idx2), 20, 2);
        elseif strcmpi(estimator, 'median')
            boot1 = median(data1(:, idx1), 2, 'omitnan');
            boot2 = median(data2(:, idx2), 2, 'omitnan');
        end
        
        boot_diff(:, b) = boot1 - boot2;
    end
    
    % Compute HDI from bootstrap distribution
    HDI3 = zeros(2, npts);
    for i = 1:npts
        sorted_diff = sort(boot_diff(i, :));
        lower_idx = round(nboot * a/2);
        upper_idx = round(nboot * (1 - a/2));
        HDI3(1, i) = sorted_diff(upper_idx);  % upper bound
        HDI3(2, i) = sorted_diff(lower_idx);  % lower bound
    end
end

% Colors
color1 = [0, 0.4470, 0.7410];           % blue
color2 = [0.8500, 0.3250, 0.0980];      % red
color3 = [0.4660, 0.6740, 0.1880];      % green

figure('color','w'); 
subplot(2,1,1) 
hold on

% Data1 (mean + 95% HDI)
p1 = plot(xAxis,est1,'LineWidth',2,'Color', color1);
patch([xAxis fliplr(xAxis)], [HDI1(1,:) fliplr(HDI1(2,:))], ...
    color1,'FaceAlpha',.3,'EdgeColor',color1,'EdgeAlpha',0.9);

% Data2 (mean + 95% HDI)
p2 = plot(xAxis, est2,'LineWidth',2,'Color', color2);
patch([xAxis fliplr(xAxis)], [HDI2(1,:) fliplr(HDI2(2,:))], ...
    color2,'FaceAlpha',.3,'EdgeColor',color2,'EdgeAlpha',0.9);
grid off; axis tight; hold on; box on
ylabel('Power (db)','fontsize',10,'fontweight','bold'); 


% Plot difference (mean + 95% HDI)
subplot(2,1,2)
plot(xAxis, est3,'LineWidth',2,'Color', color3);
patch([xAxis fliplr(xAxis)], [HDI3(1,:) fliplr(HDI3(2,:))], ...
    color3,'FaceAlpha',.6,'EdgeColor',color3,'EdgeAlpha',0.9);
grid off; axis tight; box on
ylabel('Difference (uV)','fontsize',10,'fontweight','bold')

% Add dash line to mark the null hypothesis
hold on; plot([xAxis(1) xAxis(end)], [0 0],'k--','LineWidth',1)
ylabel('Difference')
xlabel("Time (ms)",'fontsize',10,'fontweight','bold')
xlabel('Frequency (Hz)','fontsize',10,'fontweight','bold')

% Update title to show sample sizes
if n1 ~= n2
    title(sprintf('%s + %g%% quantile intervals (n1=%d, n2=%d)',estimator,(1-a)*100,n1,n2)); 
else
    title(sprintf('%s + %g%% quantile intervals (n=%d)',estimator,(1-a)*100,n1)); 
end

% Plot significance bar at the bottom
if sigBars
    plotSigBar(h, xAxis);
end

% Update legend to show sample sizes
legend([p1, p2], {sprintf('%s (n=%d)',data1Name,n1), sprintf('%s (n=%d)',data2Name,n2)}); 

set(findall(gcf,'type','axes'),'fontSize',10,'fontweight','bold');