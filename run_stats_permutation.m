function [tvals, pvals, tvals_H0, pvals_H0] = run_stats_permutation(data1, data2, nPerm, method, dpt)

% Add path to subfunctions
tmp = fileparts(which('run_stats_permutation'));
addpath(fullfile(tmp, 'functions'))

% Default parameters
alpha = 0.05;
if nargin < 3 || isempty(nPerm)
    nPerm = 1000;
end
if nargin < 4 || isempty(method)
    method = 'trimmed mean';
end
if nargin < 5 || isempty(dpt)
    error("Specify 'dpt' (paired) or 'idpt' (independent).")
end
if size(data1,3) == 1 && size(data1,2) > 1
    warning('Time dimension missing or collapsed — check input shape')
end

if ndims(data1) == 2
    data1 = reshape(data1, 1, size(data1,1), size(data1,2));
    data2 = reshape(data2, 1, size(data2,1), size(data2,2));
end

[nChan, nTimes, nSub] = size(data1);

tvals = nan(nChan, nTimes);
pvals = nan(nChan, nTimes);
tvals_H0 = nan(nChan, nTimes, nPerm);
pvals_H0 = nan(nChan, nTimes, nPerm);
disp('Computing observed test statistics...')
parfor iChan = 1:nChan
    for t = 1:nTimes
        tval = NaN;  % to avoid parfor warnings
        pval = NaN;  % to avoid parfor warnings
        x1 = squeeze(data1(iChan, t, :));
        x2 = squeeze(data2(iChan, t, :));
        if isequal(x1, x2)
            warning("channel %g, time %g: x1 and x2 are equal!", iChan, t)
        end
        if strcmpi(method, 'trimmed mean')
            if strcmpi(dpt, 'dpt')
                [tval, ~, ~, ~, pval] = yuend(x1, x2, 20, alpha);
            elseif strcmpi(dpt, 'idpt')
                [tval, ~, ~, ~, pval] = yuen(x1, x2, 20, alpha);
            end
        elseif strcmpi(method, 'mean')
            if strcmpi(dpt, 'dpt')
                try
                    [~, ~, ~, ~, ~, tval, pval] = limo_ttest(1, x1', x2', alpha);

                    % [~, ~, ~, ~, ~, tval, pval] = limo_ttest(1, x1, x2, alpha);
                catch ME
                    fprintf('Error at chan %d time %d: %s\n', iChan, t, ME.message)
                    tval = NaN; pval = NaN;
                end

            elseif strcmpi(dpt, 'idpt')
                [~, ~, ~, ~, ~, tval, pval] = limo_ttest(2, x1', x2', alpha);
            end
        else
            error("Unknown method: choose 'mean' or 'trimmed mean'")
        end
        tvals(iChan, t) = tval;
        pvals(iChan, t) = pval;
    end
end

disp('Running permutation test...')
parfor iChan = 1:nChan
    fprintf('Computing permutation statistics under H0 for channel %g/%g\n', iChan, nChan)
    permData1 = NaN;
    permData2 = NaN;
    tvec = NaN;
    pvec = NaN;
    for perm = 1:nPerm
        % Prepare permuted data
        if strcmpi(dpt, 'dpt')
            signs = (randi(2, [1,1,nSub]) - 1) * 2 - 1;
            permData1 = data1(iChan,:,:) + (data2(iChan,:,:) - data1(iChan,:,:)) .* signs;
            permData2 = data2(iChan,:,:) - (data2(iChan,:,:) - data1(iChan,:,:)) .* signs;
        elseif strcmpi(dpt, 'idpt')
            X = cat(3, data1(iChan,:,:), data2(iChan,:,:));
            n1 = size(data1,3);
            n2 = size(data2,3);
            idx = randperm(n1+n2);
            permData1 = X(:,:,idx(1:n1));
            permData2 = X(:,:,idx(n1+1:end));
        else
            error("Unknown dependency type: use 'dpt' or 'idpt'")
        end

        if strcmpi(method, 'mean')
            % Vectorized limo_ttest
            X1 = squeeze(permData1(1,:,:))';  % [nSub x nTimes]
            X2 = squeeze(permData2(1,:,:))';
            if strcmpi(dpt, 'dpt')
                [~, ~, ~, ~, ~, tvec, pvec] = limo_ttest(1, X1', X2', alpha);
            elseif strcmpi(dpt, 'idpt')
                [~, ~, ~, ~, ~, tvec, pvec] = limo_ttest(2, X1', X2', alpha);
            end
            tvals_H0(iChan,:,perm) = tvec;
            pvals_H0(iChan,:,perm) = pvec;

        else
            % Parallel loop for trimmed mean
            tvals_temp = nan(1, nTimes);
            pvals_temp = nan(1, nTimes);
            for t = 1:nTimes
                x1 = squeeze(permData1(1, t, :));
                x2 = squeeze(permData2(1, t, :));
                if strcmpi(dpt, 'dpt')
                    [tval_tmp, ~, ~, ~, pval_tmp] = yuend(x1, x2, 20, alpha);
                    tvals_temp(t) = tval_tmp;
                    pvals_temp(t) = pval_tmp;
                elseif strcmpi(dpt, 'idpt')
                    [tval_tmp, ~, ~, ~, pval_tmp] = yuen(x1, x2, 20, alpha);
                    tvals_temp(t) = tval_tmp;
                    pvals_temp(t) = pval_tmp;
                end
            end
            tvals_H0(iChan, :, perm) = tvals_temp;
            pvals_H0(iChan, :, perm) = pvals_temp;
        end
    end
end

disp('Permutation statistics completed.')
end
