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
%   data1 = randn(1, 750, 78);
%   data2 = randn(1, 750, 78);
%   [tvals, pvals] = run_perm_stats(data1, data2, 0.05, 1000, 'dpt')
%
% Cedric Cannard, 2021

function [tvals, pvals] = run_stats_perm(data1, data2, a, nPerm, dpt)

% add path to subfunctions
tmp = fileparts(which('run_stats_perm'));
addpath(fullfile(tmp,'rousselet_stats'))

% Add an extra dimension if needed (one channel data squeezed)
if ndims(data1) == 2
    data1 = reshape(data1, 1, size(data1, 1), size(data1, 2));
    data2 = reshape(data2, 1, size(data2, 1), size(data2, 2));
end

% Initialize variables
nChan = size(data1, 1);
nTimes = size(data1, 2);
nSub = size(data1, 3);
tvals = zeros(nChan, nTimes);           % Observed t-values
pvals = zeros(nChan, nTimes);           % p-values from permutation test
permDistrib = zeros(nChan, nTimes, nPerm); % Permutation distribution

% Observed t-values
disp("Running statistical tests on observed data (all electrodes)")
for iChan = 1:nChan
    for t = 1:nTimes
        if strcmpi(dpt, 'dpt')
            tval = yuend(squeeze(data1(iChan,t,:)),squeeze(data2(iChan,t,:)),20,a);     % paired for 2D vector
        elseif strcmpi(dpt, 'idpt')
            tval = yuen(squeeze(data1(iChan,t,:)),squeeze(data2(iChan,t,:)),20,a);      % unpaired for 2D vector
        else
            error("dpt variable should be 'dpt' or 'idpt' ")
        end
        tvals(iChan, t) = tval;
    end
end

% Permutation test
for perm = 1:nPerm
    permData1 = data1;
    permData2 = data2;

    for subj = 1:nSub
        if rand > 0.5
            permData1(:, :, subj) = data2(:, :, subj);
            permData2(:, :, subj) = data1(:, :, subj);
        end
    end

    for iChan = 1:nChan
        for t = 1:nTimes
            if strcmpi(dpt, 'dpt')
                tval = yuend(squeeze(permData1(iChan,t,:)),squeeze(permData2(iChan,t,:)),20,a);   % paired for 2D vector
            elseif strcmpi(dpt, 'idpt')
                tval = yuen(squeeze(permData1(iChan,t,:)),squeeze(permData2(iChan,t,:)),20,a);   % paired for 2D vector
            else
                error("dpt variable should be 'dpt' or 'idpt' ")
            end
            permDistrib(iChan, t, perm) = tval;
        end
    end
end

% Calculate p-values
for iChan = 1:nChan
    for t = 1:nTimes
        pvals(iChan, t) = mean(abs(permDistrib(iChan, t, :)) >= abs(tvals(iChan, t)));
    end
end

