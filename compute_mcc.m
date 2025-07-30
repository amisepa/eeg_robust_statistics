% Perform multiple comparisons correction on mass-univariate statistical test results.
%
% INPUTS:
%   tvals         - Observed t-values [channels x time] or [space x time]
%   pvals         - Observed uncorrected p-values (same size as tvals)
%   tvals_H0      - Null distribution t-values from permutations (3D array)
%   pvals_H0      - Null distribution p-values from permutations (3D array)
%   mcctype       - Type of correction:
%                     0 = no correction (uncorrected)
%                     1 = max-t correction
%                     2 = cluster-based correction
%                     3 = TFCE (threshold-free cluster enhancement)
%   pthresh       - Significance threshold (e.g., 0.05)
%   neighbormatrix- Spatial adjacency matrix (for clustering and TFCE)
%
% OUTPUTS:
%   mask          - Binary or labeled mask of significant effects after correction
%   pcorr         - Corrected p-values
%   nClust        - Number of significant clusters
%
% NOTES:
%   - Cluster correction assumes neighborhood structure is defined in 'neighbormatrix'.
%   - TFCE assumes 'limo_tfce' and related utilities are available.
%
% Cedric Cannard, Sep 2022

% Cedric Cannard, Sep 2022

function [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, pthresh, neighbormatrix)

% add path to subfunctions
tmp = fileparts(which('compute_mcc'));
addpath(fullfile(tmp,'functions/'))

% % Parpool with max number of workers if cluster or TFCE correction selected
% if mcctype>1
%     addons = ver;
%     if any(contains({addons.Name}, 'Parallel'))
%         ps = parallel.Settings;
%         ps.Pool.AutoCreate = true;
%         p = gcp('nocreate');
%         % delete(gcp('nocreate')) % shut down opened parpool
%         if isempty(p) % if not already on, launch it
%             c = parcluster; % cluster profile
%             % N = feature('numcores');          % only physical cores
%             N = getenv('NUMBER_OF_PROCESSORS'); % all processor (cores + threads)
%             if ischar(N), N = str2double(N)-1; end
%             c.NumWorkers = N;  % update cluster profile to include all workers
%             c.parpool();
%         end
%     end
% end

mask = [];
pcorr = [];
nClust = [];
nChan = size(tvals,1);
nTimes = size(tvals,2);

% TFCE params
if mcctype==3
    if ndims(tvals)==2 && nChan==1
        warning("Data type detected: ERP/Power with one channel")
        ndim = 1;   % one channel
    elseif ndims(tvals)==2 && nChan>1 || ndims(tvals)==3 && nChan==1
        warning("Data type detected: ERP/Power with several channels or a single time-frequency map")
        ndim = 2;   % ERP/power with several channels or single time-freq map
    elseif ndims(tvals)==3 && nChan>1
        warning("Data type detected: multi-channel ERSP")
        ndim = 3;   % ERSP
    end

    % number of bootstrap iterations
    nBoot = size(tvals_H0,3);
end


switch mcctype

    case 0  % Uncorrected
        pcorr = pvals;
        mask = pcorr <= pthresh;
        nClust = length(unique(mask)) - 1; % number of significant clusters
        if nClust > 0
            fprintf('%g significant clusters (uncorrected). \n',nClust);
        end

    case 1  % Max-correction
        [mask, pcorr] = correct_max(abs(tvals), abs(tvals_H0), pthresh); % limo_max_correction(abs(M),abs(bootM),p)

    case 2 % Cluster-correction

        % % For one channel only
        % if size(tvals,1) == 1
        %     tmp = nan(1,size(tvals,2),size(tvals_H0,2));
        %     tmp(1,:,:,:) = tvals_H0; %tvals_H0 = tmp;
        %     tmp(1,:,:,:) = pvals_H0; %pvals_H0 = tmp;
        %     clear tmp
        % end

        % Get cluster mask and corrected p-values
        % [mask,pcorr] = limo_clustering(tvals.^2,pvals,tvals_H0.^2,pvals_H0,LIMO,2,0.05,0); % LIMO source code
        [mask, pcorr] = correct_cluster(tvals.^2, pvals, tvals_H0.^2, pvals_H0, neighbormatrix, mcctype, pthresh);
        nClust = length(unique(mask)) - 1; % number of significant clusters
        if nClust > 0
            fprintf('%g significant clusters (cluster-corrected).\n',nClust)
        end

    case 3  % TFCE-correction

        % Apply Threshold-free cluster enhancement (TFCE) to null data
        disp('Applying threshold-free cluster enhancement (TFCE) to null data...')
        tfce_H0_thmaps = cell(1,nBoot);
        tfce_H0_score = nan(nChan,nTimes,nBoot);
        parfor b = 1:nBoot
            fprintf('TFCE boot %g/%g \n', b,nBoot)
            [tfce_H0_score(:,:,b), tfce_H0_thmaps{b}] = limo_tfce(ndim, squeeze(tvals_H0(:,:,b)), neighbormatrix,0);
        end

        % Max-correction on TFCE data
        tfce_score = limo_tfce(ndim, tvals, neighbormatrix);  % find threshold
        [mask, pcorr] = correct_max(tfce_score,tfce_H0_score,pthresh); % adapted from limo_max_correction
        nClust = length(unique(mask))-1;
        if nClust > 0
            fprintf('%g significant clusters (TFCE-corrected).\n',nClust)
        end
end
