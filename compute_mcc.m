function [mask, pcorr, tfce_score, tfce_H0_score] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, pthresh, chanlocs)
% COMPUTE_MCC Perform multiple comparisons correction on t-values.
%
%   This function applies statistical correction methods to identify significant
%   effects in EEG/MEG t-value maps. It supports uncorrected thresholding,
%   max-t correction, cluster-based correction, and threshold-free cluster
%   enhancement (TFCE). TFCE outputs are separated for reuse in downstream
%   merging operations.
%
% USAGE:
%   [mask, pcorr, tfce_score, tfce_H0_score] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, pthresh, chanlocs)
%
% INPUTS:
%   tvals        - Observed t-values matrix [nChannels x nFreq/nTime]
%   pvals        - Uncorrected p-values matrix [nChannels x nFreq/nTime]
%   tvals_H0     - Null t-values from permutations [nChannels x nFreq/nTime x nPerms]
%   pvals_H0     - Null p-values from permutations [nChannels x nFreq/nTime x nPerms]
%   mcctype      - Type of multiple comparison correction:
%                  0 = uncorrected
%                  1 = max-t correction
%                  2 = cluster-based correction
%                  3 = TFCE correction
%   pthresh      - Significance threshold (e.g., 0.05)
%   chanlocs     - EEGLAB-style channel location structure (used to derive spatial adjacency)
%
% OUTPUTS:
%   mask           - Binary mask of significant values after correction [same size as tvals]
%   pcorr          - Corrected p-values [same size as tvals]
%   tfce_score     - TFCE-transformed observed t-values (empty unless mcctype == 3)
%   tfce_H0_score  - TFCE-transformed null t-values (empty unless mcctype == 3)
%
% NOTES:
%   - TFCE correction uses spatial adjacency estimated from channel locations.
%   - Cluster merging, effect size computation, and summary reporting are now
%     handled in a separate function (e.g., pull_clusters.m).
%
% Copyright (C) - Cedric Cannard, 2021–2025

% Initialize
tfce_score = [];
tfce_H0_score = [];

% path to functions
repo = fileparts(which("compute_mcc.m"));
addpath(genpath(repo))

% % Accept 2D [nTimes x nSub] by promoting to [1 x nTimes x nSub]
% if ndims(tvals) == 2
%     tvals = reshape(tvals, 1, size(tvals,1), size(tvals,2));
%     pvals = reshape(pvals, 1, size(pvals,1), size(pvals,2));
%     tvals_H0 = reshape(tvals_H0, 1, size(tvals_H0,1), size(tvals_H0,2), size(tvals_H0,3));
%     pvals_H0 = reshape(pvals_H0, 1, size(pvals_H0,1), size(pvals_H0,2), size(pvals_H0,3));
% end

% Get neighbors
[~, neighbormatrix] = get_channelneighbors(chanlocs);

% Ensure neighbormatrix is sanitized
neighbormatrix = logical(neighbormatrix);
neighbormatrix = neighbormatrix | neighbormatrix.';  % symmetrize
neighbormatrix(1:size(neighbormatrix,1)+1:end) = false; % zero diagonal


switch mcctype
    case 0
        % uncorrected
        pcorr = pvals;
        mask = pcorr <= pthresh;

    case 1
        % t-max correction
        [mask, pcorr] = correct_max(abs(tvals), abs(tvals_H0), pthresh);

    case 2
        % Spatiotemporal cluster-based correction
        [mask, pcorr] = correct_cluster(tvals.^2, pvals, tvals_H0.^2, pvals_H0, neighbormatrix, mcctype, pthresh);

    case 3
        % Spatiotemporal TFCE correction with FWER max control (mask + pcorr)
        % Spatiotemporal TFCE correction with FWER max control (mask + pcorr)
        ndim  = ndims(tvals);
        nPerm = size(tvals_H0, 3);

        % Neighbor matrix cleanup (once)
        if ~isempty(neighbormatrix)
            neighbormatrix = logical(neighbormatrix);
            neighbormatrix = neighbormatrix | neighbormatrix.';
            d = size(neighbormatrix,1);
            neighbormatrix(1:d+1:end) = false;
        end

        % Allocate TFCE null scores
        tfce_H0_score = nan(size(tvals_H0), 'like', tvals_H0);

        % Optional: compute TFCE on single to cut memory traffic
        use_single = ~isa(tvals_H0, 'single');
        if use_single
            tvals_H0_work = single(tvals_H0);
        else
            tvals_H0_work = tvals_H0;
        end

        disp('Computing threshold-free cluster enhancement (TFCE)...')
        [hTFCE, qTFCE] = make_parfor_waitbar(nPerm, 'TFCE perms');
        parfor b = 1:nPerm
            tfce_H0_score(:,:,b) = limo_tfce(ndim, tvals_H0_work(:,:,b), neighbormatrix, 0);
            send(qTFCE, 1);   % tick progress
        end
        close(hTFCE);

        % Cast back if we ran on single but want double outputs
        if use_single && isa(tvals_H0, 'double')
            tfce_H0_score = double(tfce_H0_score);
        end

        % Observed TFCE (quiet)
        tfce_score = limo_tfce(ndim, tvals, neighbormatrix, 0);

        % Max-TFCE FWER correction -> mask and corrected p
        [mask, pcorr] = correct_max(tfce_score, tfce_H0_score, pthresh);

    otherwise
        error('Invalid MCC type.')
end
