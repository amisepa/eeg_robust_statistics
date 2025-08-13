%% 1. Simulate ERP trials with a nonlinear effect of trial number on P300 amplitude
% 
% NOTES:
%   - The statistic is the out-of-sample improvement in R^2 of GAM over GLM.
%   - The permutation preserves y while shuffling the continuous predictor.
%   - You can swap "trial" for any continuous covariate, or extend to multi-predictor GAMs.
%   - For mass-univariate EEG, run this per channel and time and then apply TFCE or clustering
%       to the map of delta_R2 or to an F-statistic for the smooth term.

rng(7)
fs = 250;                                % Hz
t  = -0.2:1/fs:0.8;                      % seconds
nT = numel(t);
nTrials = 200;

% Base ERP shape: a P300-like bump around 0.35 s
p300 = exp(-0.5*((t-0.35)/0.06).^2);     % Gaussian bump

% Trialwise nonlinear amplitude: smooth nonlinearity + small linear drift
trial = (1:nTrials)';                    % predictor
amp_true = 3 + 1.8*sin(2*pi*trial/220) + 0.003*trial;

% Build trialwise signals with colored noise
X = zeros(nT, nTrials);
for k = 1:nTrials
    noise = filter([1 0.6 0.2],[1 0 0],randn(1,nT))*0.25; % a bit of 1/f-ish
    X(:,k) = amp_true(k)*p300 + noise;
end

% Response variable: mean amplitude in a P300 window
win = t>=0.30 & t<=0.42;
y = mean(X(win,:),1)';                   % trialwise scalar response

%% 2. Fit GLM and GAM to predict y from trial

% 2.1 GLM with linear term
tbl = table(trial, y);
mdl_glm = fitlm(tbl,'y ~ 1 + trial');

% 2.2 GAM with smooth term
% Use k-fold CV so we can compute out-of-sample R^2 fairly
cv = cvpartition(nTrials,'KFold',10);
mdl_gam_cv = fitrgam(tbl,'y','PredictorNames','trial',...
    'CategoricalPredictors',[], 'CVPartition',cv, 'Verbose',0);

% Compute out-of-sample R^2 for GLM and GAM
% GLM CV
yhat_glm_cv = nan(nTrials,1);
for i = 1:cv.NumTestSets
    trIdx = training(cv,i); teIdx = test(cv,i);
    m = fitlm(tbl(trIdx,:),'y ~ 1 + trial');
    yhat_glm_cv(teIdx) = predict(m, tbl(teIdx,:));
end
SStot = sum( (y - mean(y)).^2 );
SSres_glm = sum( (y - yhat_glm_cv).^2 );
R2_glm = 1 - SSres_glm/SStot;

% GAM CV
yhat_gam_cv = kfoldPredict(mdl_gam_cv);
SSres_gam = sum( (y - yhat_gam_cv).^2 );
R2_gam = 1 - SSres_gam/SStot;

delta_R2_obs = R2_gam - R2_glm;

fprintf('GLM R^2 (CV): %.3f\n', R2_glm)
fprintf('GAM R^2 (CV): %.3f\n', R2_gam)
fprintf('Delta R^2 (GAM - GLM): %.3f\n', delta_R2_obs)

%% 3. Permutation test on the improvement of GAM over GLM

nPerm = 500;      % Null: the mapping between trial index and y is exchangeable.
baseSeed = 1234;   % Reproducible independent RNG per iteration

% Precompute CV masks once (as you did)
K = cv.NumTestSets;
trainMask = false(nTrials,K);
testMask  = false(nTrials,K);
for i = 1:K
    trainMask(:,i) = training(cv,i);
    testMask(:,i)  = test(cv,i);
end

% progress setup
progress = parallel.pool.DataQueue;
afterEach(progress, @(p) fprintf('Completed %d/%d permutations\n', p, nPerm));

delta_R2_null = zeros(nPerm,1);
parfor p = 1:nPerm
    s = RandStream('mrg32k3a','Seed',baseSeed);
    s.Substream = p;                      % valid for mrg32k3a
    trp = trial(randperm(s, nTrials));    % stream-aware randperm

    % GLM CV with analytic least squares
    Xp = [ones(nTrials,1) trp];
    yhat_glm = nan(nTrials,1);
    for i = 1:K
        trIdx = trainMask(:,i); teIdx = testMask(:,i);
        beta = Xp(trIdx,:) \ y(trIdx);
        yhat_glm(teIdx) = Xp(teIdx,:) * beta;
    end
    SSres_glm = sum((y - yhat_glm).^2);
    R2_glm = 1 - SSres_glm/SStot;

    % GAM CV without tables
    yhat_gam = nan(nTrials,1);
    for i = 1:K
        trIdx = trainMask(:,i); teIdx = testMask(:,i);

        % X is an n-by-1 matrix here
        Xtr = trp(trIdx);
        Xte = trp(teIdx);

        % Fit one smooth on trial, no interactions
        gamFold = fitrgam(Xtr, y(trIdx), ...
                          'Interactions', 0, ...
                          'Verbose', 0);

        yhat_gam(teIdx) = predict(gamFold, Xte);
    end
    SSres_gam = sum((y - yhat_gam).^2);
    R2_gam = 1 - SSres_gam/SStot;

    delta_R2_null(p) = R2_gam - R2_glm;

    % Update progress
    send(progress, p);
end

pval = mean(delta_R2_null >= delta_R2_obs);
fprintf('Permutation p-value for GAM improvement: %.4f\n', pval)

%% 4. Visualize fitted functions

% Fit GAM with SD so CI is available
mdl_gam = fitrgam(trial, y, ...
    'Interactions', 0, ...
    'FitStandardDeviation', true, ...
    'Verbose', 0);

% Predict fitted curve and 95% CI
tt = linspace(min(trial), max(trial), 300)';
[y_gam_line, y_gam_sd, y_gam_ci] = predict(mdl_gam, tt);  % no 'Prediction' arg
% y_gam_ci is [lower upper]; use 'Alpha',0.05 to change CI if you like

% Plot
figure('Color','w'); 
subplot(1,2,1)
plot(trial, y, '.', 'MarkerSize', 8); hold on

% GLM line
y_glm_line = predict(mdl_glm, table(tt,'VariableNames',{'trial'}));
plot(tt, y_glm_line, 'LineWidth', 1.5)

% CI patch first, then smooth
fill([tt; flipud(tt)], [y_gam_ci(:,1); flipud(y_gam_ci(:,2))], ...
     [0, 0.4470, 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.2)
plot(tt, y_gam_line, 'LineWidth', 1.5)

xlabel('Trial'); ylabel('P300 window mean (uV)')
title('GLM vs GAM fits')
legend('Trials','GLM linear','GAM 95% CI','GAM smooth','Location','best')
box off

subplot(1,2,2)
histogram(delta_R2_null, 40); hold on
xline(delta_R2_obs,'r','LineWidth',2)
xlabel('\DeltaR^2 under H0'); ylabel('Count')
title(sprintf('\\DeltaR^2 obs = %.3f, p = %.4f', delta_R2_obs, pval))
box off

% print(gcf, 'GLM_vs_GAM_fitted-functions.png','-dpng','-r300');   % 300 dpi .png



%% 5. Mass-univariate version with TFCE correction (channels × time)
% - This uses a simple 8-connected grid as the channel topology. For 
%   realistic EEG neighbors, set up a channels × channels adjacency from 
%   chanlocs and adapt TFCE to operate on a graph, or call LIMO's 
%   limo_tfce.
% 
% - The statistic is positive by construction (ΔR²), so TFCE is applied to 
%   the positive tail only.
% 
% - Start with nPermMU=200 for a quick test, then increase for final runs.
% - For your real data, you can swap the simulated Xch with your trials × 
%   time × channels array and keep the rest unchanged.


% ----- simulate multi-channel data (8x4 grid) -----
disp("Simulating 32-channel EEG data...")
nRow = 8; nCol = 4;             
nChan = nRow*nCol;

% Define conditions
condLabels = repelem([1 2], nTrials/2)';   % 1 = neutral, 2 = emotional

% Add condition-specific amplitude effect
amp_cond_effect = zeros(nTrials,1);
amp_cond_effect(condLabels == 2) = 0.8;    % ~0.8 µV larger P300 for emotional

disp("Simulating 32-channel EEG data...")
Xch = zeros(nT, nTrials, nChan);
parfor ch = 1:nChan
    Xtmp = zeros(nT, nTrials);
    topo_gain = 0.7 + 0.6*topo(ch);
    for k = 1:nTrials
        % latency-jittered P300
        p300_shift = circshift(p300, lat_jit(k));
        erp = topo_gain * (amp_base(k) + amp_cond_effect(k)) * amp_jit(k) * p300_shift;

        % noise
        n1f  = one_over_f(nT);
        nAR  = filter(1, [1 -0.6], randn(1,nT));
        slow = detrend(cumsum(0.002*randn(1,nT)));
        noise = 0.20*n1f + 0.10*nAR + 0.05*slow;

        % occasional blink artifact
        if rand < 0.04
            bl = gausswin(round(0.25*fs))';
            pos = randi([round(0.1*fs) round(0.7*fs)]);
            blink = zeros(1,nT); 
            idx = pos:min(nT, pos+numel(bl)-1);
            blink(idx) = 1.2*bl(1:numel(idx));
            noise = noise + blink;
        end

        Xtmp(:,k) = erp + noise;
    end
    Xch(:,:,ch) = Xtmp;
    % fprintf(' - channel %g/%g\n', ch, nChan)
end
disp("Done building simulated multichannel EEG data.")

% Optional spatial mixing
for k = 1:nTrials
    X2d = reshape(Xch(:,k,:), nT, nChan);
    X2d = X2d * SigmaSqrt;
    Xch(:,k,:) = reshape(X2d, nT, 1, nChan);
end

chanIdx = 19;  % example channel
erp_all = squeeze(Xch(:,:,chanIdx));

% Condition means and SE
m_neu = mean(erp_all(:, condLabels == 1), 2);
m_emo = mean(erp_all(:, condLabels == 2), 2);
se_neu = std(erp_all(:, condLabels == 1), 0, 2) / sqrt(sum(condLabels == 1));
se_emo = std(erp_all(:, condLabels == 2), 0, 2) / sqrt(sum(condLabels == 2));

% Plot
figure('Color','w'); hold on
fill([t fliplr(t)], [m_neu'+se_neu' fliplr(m_neu'-se_neu')], [0.3 0.6 1], 'FaceAlpha',0.3, 'EdgeColor','none')
plot(t, m_neu, 'Color', [0.3 0.6 1], 'LineWidth', 2)

fill([t fliplr(t)], [m_emo'+se_emo' fliplr(m_emo'-se_emo')], [1 0.4 0.4], 'FaceAlpha',0.3, 'EdgeColor','none')
plot(t, m_emo, 'Color', [1 0.4 0.4], 'LineWidth', 2)

xlabel('Time (s)')
ylabel('Amplitude (\muV)')
legend({'Neutral ±SE','Neutral mean', ...
        'Emotional ±SE','Emotional mean'}, 'Location','best')



%% ----- reuse CV masks from earlier -----

K = cv.NumTestSets;
trainMask = false(nTrials,K);
testMask  = false(nTrials,K);
for i = 1:K
    trainMask(:,i) = training(cv,i);
    testMask(:,i)  = test(cv,i);
end

% ------------------------------------------------------ %
% Compute mass-univariate ΔR^2 map on observed data
% ------------------------------------------------------ %
% disp("Computing mass-univariate ΔR^2 map on observed data...")
% deltaR2_obs_map = zeros(nChan, nT);
% parfor iChan = 1:nChan
%     Y = squeeze(Xch(:,:,iChan)).';          % trials x time
%     drow = zeros(1, nT);
%     for it = 1:nT
%         y_t = Y(:,it);
% 
%         % GLM CV (analytic)
%         Xlin = [ones(nTrials,1) trial];
%         yhat_glm = nan(nTrials,1);
%         for ii = 1:K
%             trIdx = trainMask(:,ii); teIdx = testMask(:,ii);
%             beta = Xlin(trIdx,:) \ y_t(trIdx);
%             yhat_glm(teIdx) = Xlin(teIdx,:)*beta;
%         end
%         SStot_t = sum((y_t-mean(y_t)).^2);
%         R2_glm_t = 1 - sum((y_t - yhat_glm).^2)/SStot_t;
% 
%         % GAM CV (one smooth on trial)
%         yhat_gam = nan(nTrials,1);
%         for ii = 1:K
%             trIdx = trainMask(:,ii); teIdx = testMask(:,ii);
%             gamFold = fitrgam(trial(trIdx), y_t(trIdx), ...
%                               'Interactions',0,'Verbose',0);
%             yhat_gam(teIdx) = predict(gamFold, trial(teIdx));
%         end
%         R2_gam_t = 1 - sum((y_t - yhat_gam).^2)/SStot_t;
% 
%         drow(it) = R2_gam_t - R2_glm_t;
%     end
%     deltaR2_obs_map(iChan,:) = drow;
%     % fprintf(' - channel %g/%g\n', iChan, nChan)
% end
% disp("Done computing observed ΔR² map.")

% Fast observed ΔR^2 map version (using the same spline basis and lambda)
% Reuse: trial, knots, lambda, cv, nTrials, nChan, nT, Xch
% Build CV index lists once
K = cv.NumTestSets;
trainIdx = cell(K,1); testIdx = cell(K,1);
for ii = 1:K
    trainIdx{ii} = find(training(cv,ii));
    testIdx{ii}  = find(test(cv,ii));
end

% Trials x time per channel and SStot per channel×time
Y_all = cell(nChan,1);
SStot_map = zeros(nChan, nT);
for iChan = 1:nChan
    Y = squeeze(Xch(:,:,iChan)).';        % trials x time
    Y_all{iChan} = Y;
    mu = mean(Y,1);
    SStot_map(iChan,:) = sum((Y - mu).^2, 1);
end

% Fixed designs for observed data
condVec = double(isEmotional);  % 1 = Emotional, 0 = Neutral
Xlin_obs = [ones(nTrials,1) condVec];  % GLM: intercept + condition
Xspl_obs = ncs_basis(condVec, knots);  % spline: smooth of condition if wanted
% Xlin_obs = [ones(nTrials,1) trial];       % GLM
% Xspl_obs = ncs_basis(trial, knots);       % spline basis
M = size(Xspl_obs,2);

disp("Computing observed ΔR^2 map with spline CV...")
deltaR2_obs_map = zeros(nChan, nT);
parfor iChan = 1:nChan
    Y = Y_all{iChan};                      % trials x time

    Yhat_glm = nan(nTrials, nT);
    Yhat_spl = nan(nTrials, nT);

    for ii = 1:K
        tr = trainIdx{ii}; te = testIdx{ii};

        % GLM: solve for all times jointly
        XtX = Xlin_obs(tr,:)' * Xlin_obs(tr,:);
        XtY = Xlin_obs(tr,:)' * Y(tr,:);
        beta_lin = XtX \ XtY;                          % 2 x nT
        Yhat_glm(te,:) = Xlin_obs(te,:) * beta_lin;

        % Spline ridge: solve for all times jointly
        ZtZ = Xspl_obs(tr,:)' * Xspl_obs(tr,:) + lambda * eye(M);
        ZtY = Xspl_obs(tr,:)' * Y(tr,:);
        beta_spl = ZtZ \ ZtY;                          % M x nT
        Yhat_spl(te,:) = Xspl_obs(te,:) * beta_spl;
    end

    SSres_glm = sum((Y - Yhat_glm).^2, 1);
    SSres_spl = sum((Y - Yhat_spl).^2, 1);
    R2_glm = 1 - SSres_glm ./ SStot_map(iChan,:);
    R2_spl = 1 - SSres_spl ./ SStot_map(iChan,:);

    deltaR2_obs_map(iChan,:) = R2_spl - R2_glm;

    fprintf(' - finished channel %g/%g\n', iChan, nChan)
end
disp("Done computing observed ΔR^2 map.")

%% % ------------------------------------------------------ %
% % Permutation null and TFCE max (very slow)
% % ------------------------------------------------------ %
% nPermMU    = 100;
% baseSeedMU = 5678;
% 
% maxTFCE_null = zeros(nPermMU,1);
% disp("Running permutation test for mass-univariate ΔR² map...")
% parfor p = 1:nPermMU
%     s = RandStream('mrg32k3a','Seed',baseSeedMU); s.Substream = p;
%     trp = trial(randperm(s, nTrials));
%     Xp  = [ones(nTrials,1) trp];
% 
%     delta_map = zeros(nChan, nT);
% 
%     for iChan = 1:nChan
%         Y = squeeze(Xch(:,:,iChan)).';
% 
%         for it = 1:nT
%             y_t = Y(:,it);
% 
%             % GLM CV
%             yhat_glm = nan(nTrials,1);
%             for ii = 1:K
%                 trIdx = trainMask(:,ii); teIdx = testMask(:,ii);
%                 beta = Xp(trIdx,:) \ y_t(trIdx);
%                 yhat_glm(teIdx) = Xp(teIdx,:)*beta;
%             end
%             SStot_t = sum((y_t - mean(y_t)).^2);
%             R2_glm_t = 1 - sum((y_t - yhat_glm).^2) / SStot_t;
% 
%             % GAM CV
%             yhat_gam = nan(nTrials,1);
%             for ii = 1:K
%                 trIdx = trainMask(:,ii); teIdx = testMask(:,ii);
%                 g = fitrgam(trp(trIdx), y_t(trIdx), ...
%                             'Interactions', 0, 'Verbose', 0);
%                 yhat_gam(teIdx) = predict(g, trp(teIdx));
%             end
%             R2_gam_t = 1 - sum((y_t - yhat_gam).^2) / SStot_t;
% 
%             delta_map(iChan,it) = R2_gam_t - R2_glm_t;
%         end
%     end
% 
%     tfce_perm = tfce2d(delta_map, 0.5, 2, 0.001);
%     maxTFCE_null(p) = max(tfce_perm(:));
% 
% fprintf(' - permutation %g/%g\n', p, nPermMU)
% end
% disp("Done running permutation test for mass-univariate ΔR² map.")

% ------------------------------------------------------ %
% Fast permutation null and TFCE max (vectorized, no fitrgam)
% ------------------------------------------------------ %
nPermMU    = 1000;
baseSeedMU = 5678;

% folds
K = cv.NumTestSets;
trainIdx = cell(K,1); testIdx = cell(K,1);
for ii = 1:K
    trainIdx{ii} = find(training(cv,ii));
    testIdx{ii}  = find(test(cv,ii));
end

% precompute SStot per channel×time (does not depend on permutation)
SStot_map = zeros(nChan, nT);
Y_all = cell(nChan,1); % trials x time per channel
for iChan = 1:nChan
    Y = squeeze(Xch(:,:,iChan)).';      % trials x time
    Y_all{iChan} = Y;
    mu = mean(Y,1);
    SStot_map(iChan,:) = sum((Y - mu).^2, 1);
end

% choose spline knots once (quantiles of trial)
KnotCount = 10;  % try 8-15
qs = linspace(0.05, 0.95, KnotCount);
knots = quantile(trial, qs);
lambda = 1e-2;   % ridge for spline coefficients

% Cache big, read-only data on workers
Yc      = parallel.pool.Constant(@() Y_all);           % cell: trials x time per chan
SStotc  = parallel.pool.Constant(@() SStot_map);       % [nChan x nT]
trainc  = parallel.pool.Constant(@() trainIdx);        % cell of indices
testc   = parallel.pool.Constant(@() testIdx);         % cell of indices

% Precompute spline dimension and identity once per worker
M       = size(ncs_basis(trial, knots), 2);
Ic      = parallel.pool.Constant(@() eye(M));

disp("Running permutation test for mass-univariate ΔR^2 map...")
maxTFCE_null = zeros(nPermMU,1);
parfor p = 1:nPermMU
    s = RandStream('mrg32k3a','Seed',baseSeedMU); s.Substream = p;
    trp = trial(randperm(s, nTrials));

    Xlin = [ones(nTrials,1) trp];       % GLM
    Xspl = ncs_basis(trp, knots);       % spline basis

    delta_map = zeros(nChan, nT);

    for iChan = 1:nChan
        Y = Yc.Value{iChan};            % trials x time
        Yhat_glm = nan(nTrials, nT);
        Yhat_spl = nan(nTrials, nT);

        for ii = 1:K
            tr = trainc.Value{ii}; te = testc.Value{ii};

            % GLM, all times
            XtX = Xlin(tr,:)' * Xlin(tr,:);
            XtY = Xlin(tr,:)' * Y(tr,:);
            beta_lin = XtX \ XtY;
            Yhat_glm(te,:) = Xlin(te,:) * beta_lin;

            % Spline ridge, all times
            ZtZ = Xspl(tr,:)' * Xspl(tr,:) + lambda * Ic.Value;
            ZtY = Xspl(tr,:)' * Y(tr,:);
            beta_spl = ZtZ \ ZtY;
            Yhat_spl(te,:) = Xspl(te,:) * beta_spl;
        end

        SSres_glm = sum((Y - Yhat_glm).^2, 1);
        SSres_spl = sum((Y - Yhat_spl).^2, 1);
        R2_glm = 1 - SSres_glm ./ SStotc.Value(iChan,:);
        R2_spl = 1 - SSres_spl ./ SStotc.Value(iChan,:);

        delta_map(iChan,:) = R2_spl - R2_glm;
    end

    tfce_perm = tfce2d(delta_map, 0.5, 2, 0.001);
    maxTFCE_null(p) = max(tfce_perm(:));

    fprintf(' - permutation %g/%g\n', p, nPermMU)
end
disp("Done running permutation test for mass-univariate ΔR^2 map.")


%% =======================
% TFCE on observed map and FWER p
% =======================

tfce_obs = tfce2d(deltaR2_obs_map, 0.5, 2, 0.001);
p_FWER = arrayfun(@(x) mean(maxTFCE_null >= x), tfce_obs);

alpha = 0.05;
sig_mask = p_FWER < alpha;

% =======================
% Plots
% =======================
figure('Color','w')
subplot(1,3,1)
imagesc(t, 1:nChan, deltaR2_obs_map); axis xy
title('\DeltaR^2 map (GAM - GLM)')
xlabel('Time (s)'); ylabel('Channel'); colorbar

subplot(1,3,2)
imagesc(t, 1:nChan, tfce_obs); axis xy
title('TFCE(\DeltaR^2)'); xlabel('Time (s)'); ylabel('Channel'); colorbar

subplot(1,3,3)
imagesc(t, 1:nChan, sig_mask); axis xy
title(sprintf('FWER<%.2f TFCE mask', alpha))
xlabel('Time (s)'); ylabel('Channel')
colormap(subplot(1,3,3), gray)


%% Subfunctions

function tfce = tfce2d(map, E, H, dh)
% TFCE for positive values on a 2-D array using 8-connectivity
map = max(map, 0);
mx = max(map(:));
if mx==0
    tfce = zeros(size(map));
    return
end

heights = 0:dh:mx;
tfce = zeros(size(map));

for h = heights
    bw = map > h;
    if ~any(bw(:)), continue, end
    CC = bwconncomp(bw, 8);
    for c = 1:CC.NumObjects
        idx = CC.PixelIdxList{c};
        extent = numel(idx)^E;
        tfce(idx) = tfce(idx) + (h^H) * extent * dh;
    end
end
end

% ---------- speed helpers ----------
% natural cubic spline basis (truncated power) with K internal knots
function Phi = ncs_basis(x, knots)
    x = x(:);
    K = numel(knots);
    % columns: intercept, linear x, truncated cubics for each knot
    Phi = [ones(numel(x),1) x];
    % truncated cubic helper
    function v = tc(z,k)
        v = max(z - k, 0).^3;
    end
    % use last knot for natural boundary adjustment
    kK = knots(end);
    for j = 1:K-1
        dj = (tc(x, knots(j)) - tc(x, kK)) / (kK - knots(j));
        Phi = [Phi dj]; %#ok<AGROW>
    end
end
