% Assuming data1 and data2 are your two conditions with size [channels, time points, subjects]
% data1 and data2 are matrices of size [channels, time points, subjects]
%
% INPUTS:
%   data1   - condition/group 1 (channels x time/frequency x subjects)
%   data2   - condition/group 1 (channels x time/frequency x subjects)
%   a       - alpha level (e.g. 0.05)
%   nPerm   - number of permutations to perform (e.g., 1000)
%   dpt     - whether data are dependent ('dpt') or independent ('idpt')
%
% OUTPUTS:
%   tvals   - t-values for each channel and time/frequency
%   pvals   - p-values for each channel and time/frequency
%
% USAGE:
%   [tvals, pvals] = run_perm_stats(data1, data2, a, nPerm, dpt)
%
% EXAMPLE:
  % data1 = randn(1, 750, 78);
  % data2 = randn(1, 750, 78);
  % [tvals, pvals] = run_stats_perm(data1, data2, 0.05, 1000, 'dpt')
%
% Cedric Cannard, 2021

function [tvals, pvals] = run_stats_perm(data1, data2, a, nPerm, dpt)

% Ensure data has correct dimensions (add extra channel dimension if needed)
if ndims(data1) == 2
    data1 = reshape(data1, 1, size(data1, 1), size(data1, 2));
    data2 = reshape(data2, 1, size(data2, 1), size(data2, 2));
end

% Initialize variables
[nChan, nTimes, nSub] = size(data1);
tvals = zeros(nChan, nTimes);           % Observed t-values
permDistrib = zeros(nChan, nTimes, nPerm); % Permutation distribution

fprintf('Running Yuen t-tests on all electrodes against H0 estimated with %d permutations...\n', nPerm);

% Adjust permutation count if small sample size
if nSub <= 12
    nPerm = 2^nSub; % Compute all possible permutations
    exact = true;
    fprintf('Using exact permutations (%d total) due to small sample size.\n', nPerm);
else
    exact = false;
end

% Compute observed t-values
% disp("Running statistical tests on observed data (all electrodes)")
for iChan = 1:nChan
    parfor t = 1:nTimes
        if strcmpi(dpt, 'dpt') % Paired test (trimmed mean t-test)
            tvals(iChan, t) = yuend(squeeze(data1(iChan,t,:)), squeeze(data2(iChan,t,:)), 20, a);
        elseif strcmpi(dpt, 'idpt') % Independent test (trimmed mean t-test)
            tvals(iChan, t) = yuen(squeeze(data1(iChan,t,:)), squeeze(data2(iChan,t,:)), 20, a);
        else
            error("dpt variable should be 'dpt' or 'idpt' ")
        end
    end
end

% Permutation test (paired or independent)
for perm = 1:nPerm
    permData1 = data1;
    permData2 = data2;

    if strcmpi(dpt, 'dpt')  % Paired permutation (sign flipping)
        sign_flips = (randi(2, [1, 1, nSub]) - 1) * 2 - 1; % Generate +1 or -1
        permData1 = data1 + (data2 - data1) .* sign_flips;
        permData2 = data2 - (data2 - data1) .* sign_flips;
    else  % Independent samples: Shuffle labels
        % swapIdx = rand(1, nSub) > 0.5;
        % permData1(:, :, swapIdx) = data2(:, :, swapIdx);
        % permData2(:, :, swapIdx) = data1(:, :, swapIdx);
        % Independent group permutation (shuffling across subjects)
        X = cat(3, data1, data2);
        n1 = size(data1,3);
        n2 = size(data2,3);
        perm_idx = randperm(n1 + n2);
        permData1 = X(:,:,perm_idx(1:n1));
        permData2 = X(:,:,perm_idx(n1+1:end));
    end

    % Compute permuted t-values in a vectorized way
    for iChan = 1:nChan
        parfor t = 1:nTimes
            if strcmpi(dpt, 'dpt')
                permDistrib(iChan, t, perm) = yuend(squeeze(permData1(iChan,t,:)), squeeze(permData2(iChan,t,:)), 20, a);
            elseif strcmpi(dpt, 'idpt')
                permDistrib(iChan, t, perm) = yuen(squeeze(permData1(iChan,t,:)), squeeze(permData2(iChan,t,:)), 20, a);
            end
        end
    end
end

% Compute p-values (vectorized)
pvals = mean(abs(permDistrib) >= abs(tvals), 3);

end
