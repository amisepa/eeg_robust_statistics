%%  EEG simulation of GLM vs GAM (fast)
% Two-condition design: Neutral vs Emotional with nonlinear effects across trials.
% Diagnostics, visualizations, CV + permutation test for GAM>GLM.

clear; clc; close all
rng(7)

% Sampling and timeline
fs = 250;                                % Hz
t  = -0.2:1/fs:0.8;                      % seconds
nT = numel(t);
nTrials = 200;

% Base ERP shape near 350 ms
p300 = exp(-0.5*((t-0.35)/0.06).^2);

% Conditions: 1 = Neutral, 2 = Emotional
condLabels = repelem([1 2], nTrials/2)';  % equal split
isEmo = double(condLabels==2);

% Trial index and smooth practice curves
trial = (1:nTrials)';                     % continuous covariate
drift = 0.0015*trial;                     % tiny linear drift

w1 = 2*pi/120;                            % slow frequency
w2 = 2*pi/60;                             % a bit faster
nonlin_neu = 0.6*sin(w1*trial);           % neutral smooth
nonlin_emo = 0.9*sin(w1*trial + 0.7) + 0.4*sin(w2*trial);  % richer smooth

% Condition mean shift and interaction-like difference in smooths
amp_cond   = 0.4*isEmo;                   % emotional mean boost
amp_smooth = (1-isEmo).*nonlin_neu + isEmo.*nonlin_emo;

% Baseline amplitude and jitter
amp_base = 2.6 + drift + amp_cond + amp_smooth;
lat_jit = round(fs * (0.015*randn(nTrials,1)));   % ~15 ms SD
amp_jit = 1 + 0.30*randn(nTrials,1);              % 30 percent CV
one_over_f = @(n) ifft(fft(randn(1,n)).*(1:max(1,n)).^(-0.5),'symmetric');

% Simulate single channel used for y
X = zeros(nT, nTrials);
for k = 1:nTrials
    p300_shift = circshift(p300, lat_jit(k));
    signal = amp_base(k) .* amp_jit(k) .* p300_shift;
    n1f  = 0.18*one_over_f(nT);
    nAR  = 0.10*filter(1,[1 -0.6],randn(1,nT));
    slow = 0.04*detrend(cumsum(randn(1,nT)));
    X(:,k) = signal + n1f + nAR + slow;
end

% Windowed amplitude
win = t>=0.30 & t<=0.42;
y = mean(X(win,:),1)';                      % nTrials x 1

%% ERP mean ± SE per condition
neu = condLabels==1; emo = condLabels==2;
m_neu = mean(X(:,neu),2);  se_neu = std(X(:,neu),0,2)/sqrt(sum(neu));
m_emo = mean(X(:,emo),2);  se_emo = std(X(:,emo),0,2)/sqrt(sum(emo));

figure('Color','w'); hold on
fill([t fliplr(t)], [m_neu'+se_neu' fliplr(m_neu'-se_neu')], [0.3 0.6 1], 'FaceAlpha',0.3, 'EdgeColor','none')
plot(t, m_neu, 'Color',[0.3 0.6 1], 'LineWidth',2)
fill([t fliplr(t)], [m_emo'+se_emo' fliplr(m_emo'-se_emo')], [1 0.4 0.4], 'FaceAlpha',0.3, 'EdgeColor','none')
plot(t, m_emo, 'Color',[1 0.4 0.4], 'LineWidth',2)
xlabel('Time (s)'), ylabel('Amplitude (\muV)')
legend({'Neutral ±SE','Neutral mean','Emotional ±SE','Emotional mean'}, 'Location','best')
title('Condition ERPs (single channel)')
box off

%% GLM vs GAM on windowed amplitude y

% GLM: y ~ 1 + cond + trial + cond*trial  (linear)
% GAM: y ~ 1 + cond + s(trial) + cond*s(trial)  (spline with ridge)

cv = cvpartition(nTrials,'KFold',10);

% Spline basis for trial
KnotCount = 10; qs = linspace(0.05,0.95,KnotCount);
knots  = quantile(trial, qs);
B      = spline_basis_only(trial, knots);   % returns (K-1) columns
lambda = 1e-2;

% Design matrices
X_glm = [ones(nTrials,1) isEmo trial isEmo.*trial];
Z_gam = [ones(nTrials,1) isEmo B (isEmo.*B)];

% CV predictions
yhat_glm = nan(nTrials,1); yhat_gam = nan(nTrials,1);
for i = 1:cv.NumTestSets
    tr = training(cv,i); te = test(cv,i);

    % GLM
    betaL = X_glm(tr,:)\y(tr);
    yhat_glm(te) = X_glm(te,:)*betaL;

    % GAM ridge
    Ztr = Z_gam(tr,:); Zte = Z_gam(te,:);
    betaS = (Ztr'*Ztr + lambda*eye(size(Ztr,2))) \ (Ztr'*y(tr));
    yhat_gam(te) = Zte*betaS;
end

SStot = sum((y-mean(y)).^2);
R2_glm = 1 - sum((y-yhat_glm).^2)/SStot;
R2_gam = 1 - sum((y-yhat_gam).^2)/SStot;
delta_R2_obs = R2_gam - R2_glm;
fprintf('GLM R^2 (CV): %.3f\n', R2_glm)
fprintf('GAM R^2 (CV): %.3f\n', R2_gam)
fprintf('Delta R^2 (GAM - GLM): %.3f\n', delta_R2_obs)

%% Permutation test: shuffle condition labels only
nPerm = 1000; baseSeed = 1234;

% Precompute CV index lists
K = cv.NumTestSets;
trainIdx = cell(K,1); testIdx = cell(K,1);
for ii = 1:K
    trainIdx{ii} = find(training(cv,ii));
    testIdx{ii}  = find(test(cv,ii));
end

delta_R2_null = zeros(nPerm,1);
progress = parallel.pool.DataQueue;
afterEach(progress, @(p) (mod(p,25)==0) && fprintf('Completed %d/%d perms\n', p, nPerm));

parfor p = 1:nPerm
    s = RandStream('mrg32k3a','Seed',baseSeed); s.Substream = p;
    emoP = isEmo(randperm(s, nTrials));              % shuffle condition only
    Xp   = [ones(nTrials,1) emoP trial emoP.*trial];
    Zp   = [ones(nTrials,1) emoP B (emoP.*B)];

    yL = nan(nTrials,1); yS = nan(nTrials,1);
    for ii = 1:K
        tr = trainIdx{ii}; te = testIdx{ii};
        % GLM
        bL = Xp(tr,:)\y(tr);  yL(te) = Xp(te,:)*bL;
        % GAM ridge
        Ztr = Zp(tr,:); Zte = Zp(te,:);
        bS = (Ztr'*Ztr + lambda*eye(size(Ztr,2))) \ (Ztr'*y(tr));
        yS(te) = Zte*bS;
    end
    R2L = 1 - sum((y-yL).^2)/SStot;
    R2S = 1 - sum((y-yS).^2)/SStot;
    delta_R2_null(p) = R2S - R2L;

    send(progress, p);
end
pval = mean(delta_R2_null >= delta_R2_obs);
fprintf('Permutation p-value for GAM improvement: %.4f\n', pval)

%% Fitted functions plot by condition
tt  = linspace(min(trial),max(trial),300)';
Btt = spline_basis_only(tt, knots);

% Fit on full data for smooth lines
betaGLM = X_glm\y;
betaGAM = (Z_gam'*Z_gam + lambda*eye(size(Z_gam,2))) \ (Z_gam'*y);

% Lines for each condition
Xtt_neu = [ones(numel(tt),1) zeros(numel(tt),1) tt zeros(numel(tt),1)];
Xtt_emo = [ones(numel(tt),1) ones(numel(tt),1)  tt tt];
Ztt_neu = [ones(numel(tt),1) zeros(numel(tt),1) Btt zeros(size(Btt))];
Ztt_emo = [ones(numel(tt),1) ones(numel(tt),1)  Btt Btt];

y_neu_lin = Xtt_neu*betaGLM;  y_emo_lin = Xtt_emo*betaGLM;
y_neu_gam = Ztt_neu*betaGAM;  y_emo_gam = Ztt_emo*betaGAM;

figure('Color','w');
subplot(1,2,1); hold on
plot(trial(neu), y(neu), '.', 'MarkerSize',6)
plot(trial(emo), y(emo), '.', 'MarkerSize',6)
plot(tt, y_neu_lin, 'LineWidth',1.5)
plot(tt, y_emo_lin, 'LineWidth',1.5)
plot(tt, y_neu_gam, 'LineWidth',2)
plot(tt, y_emo_gam, 'LineWidth',2)
xlabel('Trial'); ylabel('P300 window mean (\muV)')
title('GLM vs GAM fits by condition')
legend('Neutral trials','Emotional trials','GLM Neutral','GLM Emotional','GAM Neutral','GAM Emotional','Location','best')
box off

subplot(1,2,2)
histogram(delta_R2_null,40); hold on
xline(delta_R2_obs,'r','LineWidth',2)
xlabel('\DeltaR^2 under H0'); ylabel('Count')
title(sprintf('\\DeltaR^2 obs = %.3f, p = %.4f', delta_R2_obs, pval))
box off





%%  Mass-univariate ΔR² map with TFCE (32 channels)
% Simple 8x4 grid, local spatial correlation, same models as above

nRow = 8; nCol = 4; nChan = nRow*nCol;

% Topography gain pattern (center-parietal)
[xg,yg] = ndgrid(linspace(-1,1,nRow), linspace(-1,1,nCol));
topo = exp(-((xg).^2 + (yg-0.4).^2)/0.5); topo = topo./max(topo(:));

% Build sparse 4-neighbor Laplacian for spatial mixing
W = spalloc(nChan,nChan,4*nChan);
node = @(r,c) (c-1)*nRow + r;
for r = 1:nRow
    for c = 1:nCol
        i = node(r,c);
        if r > 1,    W(i,node(r-1,c)) = -1; end
        if r < nRow, W(i,node(r+1,c)) = -1; end
        if c > 1,    W(i,node(r,c-1)) = -1; end
        if c < nCol, W(i,node(r,c+1)) = -1; end
        W(i,i) = -sum(W(i,:)<0);
    end
end
Sigma = expm(0.1*full(W));
[U,Sv]=svd(Sigma); SigmaSqrt = U*sqrt(Sv)*U';

% Simulate channels
disp("Simulating 32-channel EEG data...")
Xch = zeros(nT,nTrials,nChan);
parfor ch = 1:nChan
    Xtmp = zeros(nT,nTrials);
    gain = 0.7 + 0.6*topo(ch);
    for k=1:nTrials
        psh = circshift(p300, lat_jit(k));
        sig = gain*(amp_base(k)+amp_cond(k)+amp_int(k))*amp_jit(k).*psh;
        n1f = 0.18*one_over_f(nT);
        nAR = 0.10*filter(1,[1 -0.6],randn(1,nT));
        slow= 0.04*detrend(cumsum(randn(1,nT)));
        Xtmp(:,k) = sig + n1f + nAR + slow;
    end
    Xch(:,:,ch) = Xtmp;
end
% spatial mixing
for k=1:nTrials
    X2d = reshape(Xch(:,k,:),nT,nChan) * SigmaSqrt;
    Xch(:,k,:) = reshape(X2d,nT,1,nChan);
end
disp("Done building simulated multichannel EEG data.")

% Fast observed ΔR² using same designs
K = cv.NumTestSets;
trainIdx = cell(K,1); testIdx = cell(K,1);
for ii=1:K, trainIdx{ii}=find(training(cv,ii)); testIdx{ii}=find(test(cv,ii)); end
Y_all = cell(nChan,1); SStot_map = zeros(nChan,nT);
for ch=1:nChan
    Y = squeeze(Xch(:,:,ch)).'; Y_all{ch}=Y;
    mu = mean(Y,1); SStot_map(ch,:) = sum((Y-mu).^2,1);
end

X_lin = X_glm;                % reuse
Z_gam = Z_gam;                % reuse from earlier

disp("Computing observed ΔR^2 map with spline CV...")
deltaR2_obs_map = zeros(nChan,nT);
parfor ch = 1:nChan
    Y = Y_all{ch};  YhatL = nan(nTrials,nT); YhatS = nan(nTrials,nT);
    for ii=1:K
        tr = trainIdx{ii}; te = testIdx{ii};
        % GLM
        bL = (X_lin(tr,:)'*X_lin(tr,:)) \ (X_lin(tr,:)'*Y(tr,:));
        YhatL(te,:) = X_lin(te,:)*bL;
        % GAM ridge
        Ztr = Z_gam(tr,:); Zte = Z_gam(te,:);
        bS = (Ztr'*Ztr + lambda*eye(size(Ztr,2))) \ (Ztr'*Y(tr,:));
        YhatS(te,:) = Zte*bS;
    end
    SSresL = sum((Y-YhatL).^2,1);
    SSresS = sum((Y-YhatS).^2,1);
    R2L = 1 - SSresL ./ SStot_map(ch,:);
    R2S = 1 - SSresS ./ SStot_map(ch,:);
    deltaR2_obs_map(ch,:) = R2S - R2L;
end
disp("Done computing observed ΔR^2 map.")

% Fast permutation null for TFCE max (shuffle condition labels)
nPermMU = 500; baseSeedMU=5678;
maxTFCE_null = zeros(nPermMU,1);
disp("Running permutation test for mass-univariate ΔR^2 map...")
parfor p=1:nPermMU
    s = RandStream('mrg32k3a','Seed',baseSeedMU); s.Substream=p;
    emoP = isEmo(randperm(s,nTrials));
    Xp = [ones(nTrials,1) emoP trial emoP.*trial];
    Zp = [ones(nTrials,1) emoP B (emoP.*B)];
    maxv = -inf;
    for ch=1:nChan
        Y = Y_all{ch};  YhatL = nan(nTrials,nT); YhatS = nan(nTrials,nT);
        for ii=1:K
            tr=trainIdx{ii}; te=testIdx{ii};
            bL = (Xp(tr,:)'*Xp(tr,:)) \ (Xp(tr,:)'*Y(tr,:));
            YhatL(te,:) = Xp(te,:)*bL;
            Ztr = Zp(tr,:); Zte = Zp(te,:);
            bS = (Ztr'*Ztr + lambda*eye(size(Ztr,2))) \ (Ztr'*Y(tr,:));
            YhatS(te,:) = Zte*bS;
        end
        d = (1 - sum((Y-YhatS).^2,1)./SStot_map(ch,:)) ...
            - (1 - sum((Y-YhatL).^2,1)./SStot_map(ch,:));
        tf = tfce2d(d,0.5,2,0.001);
        mv = max(tf(:)); if mv>maxv, maxv=mv; end
    end
    maxTFCE_null(p)=maxv;
    if mod(p,25)==0, fprintf(' - permutation %d/%d\n',p,nPermMU); end
end
disp("Done running permutation test.")

% TFCE on observed map and FWER p
tfce_obs = tfce2d(deltaR2_obs_map,0.5,2,0.001);
p_FWER = arrayfun(@(x) mean(maxTFCE_null>=x), tfce_obs);
alpha=0.05; sig_mask = p_FWER<alpha;

figure('Color','w')
subplot(1,3,1), imagesc(t,1:nChan,deltaR2_obs_map); axis xy
xlabel('Time (s)'), ylabel('Channel'); title('\DeltaR^2 map (GAM - GLM)'); colorbar
subplot(1,3,2), imagesc(t,1:nChan,tfce_obs); axis xy
xlabel('Time (s)'), ylabel('Channel'); title('TFCE(\DeltaR^2)'); colorbar
subplot(1,3,3), imagesc(t,1:nChan,sig_mask); axis xy
xlabel('Time (s)'), ylabel('Channel'); title('FWER<0.05 TFCE mask')
colormap(subplot(1,3,3), gray)

%% ---------------- helper functions ----------------

function Phi = spline_basis_only(x, knots)
% Natural cubic spline columns excluding intercept and linear term
% Returns (K-1) columns using truncated power differences
x = x(:); K = numel(knots); Phi = zeros(numel(x), K-1);
tc = @(z,k) max(z-k,0).^3;
kK = knots(end);
for j=1:K-1
    Phi(:,j) = (tc(x,knots(j)) - tc(x,kK)) / (kK - knots(j));
end
end

function tfce = tfce2d(map,E,H,dh)
% TFCE for positive values on a 2-D array using 8-connectivity
map = max(map,0); mx = max(map(:)); if mx==0, tfce=zeros(size(map)); return; end
heights = 0:dh:mx; tfce = zeros(size(map));
for h = heights
    bw = map>h; if ~any(bw(:)), continue, end
    CC = bwconncomp(bw,8);
    for c=1:CC.NumObjects
        idx = CC.PixelIdxList{c};
        extent = numel(idx)^E;
        tfce(idx) = tfce(idx) + (h^H)*extent*dh;
    end
end
end
