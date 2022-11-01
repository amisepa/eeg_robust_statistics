%% Cluster-correction of mass-univariate EEG data
% 
% Cedric Cannard, Sep 2022

function [mask, pcorr] = correct_cluster2(tvals, pvals, tvals_H0, pvals_H0, area_pairs, neighbors, neighbormatrix, mcctype, pthresh, fig)

% [mask,M] = limo_clustering(M.^2,Pval,bootM.^2,bootP,LIMO,MCC,p); % for t-test only

% fig = 0;
tvals = tvals.^2;
tvals_H0 = tvals_H0.^2;

% if 1 channel, force switch to 1D TFCE clustering
if size(tvals,1) == 1
    mcctype = 4; % TFCE 
end

% nb of boostrap performed
nboot = size(tvals_H0,3);

% Convert neighbor matrix to same dimension as connectivity data
% for iPair = 1:length(area_pairs)
%     neighbormatrix2(iPair,:) = 
% 
% end


% spatiotemporal clustering
if mcctype == 3 && size(tvals_H0,1) > 1
    minchan = 2; % the minimum number of neighbouring channels/combinations
    boot_maxclustersum = zeros(nboot,1);     % maximum cluster mass at each bootstrap
    disp('Getting spatiotemporal clusters under H0 ...');

    for boot = 1:nboot
        % 1st find the cluster, thresholding H0 pvalues <= threshold p
        for iPairs = 1:68:length(tvals)
%             [posclusterslabelmat, nposclusters] = limo_findcluster((pvals_H0(iPairs,:,boot) <= pthresh), neighbormatrix, minchan);
            onoff = (pvals_H0(iPairs,:,boot) <= pthresh);
            
            spatdimlength = size(onoff, 1);
            nfreq         = size(onoff, 2);
            ntime         = size(onoff, 3);
            
            % For each frequency, calculate number of significant neighbours 
            % If less than minchan, remove this area from onoff
%             selectmat = single(neighbormatrix);
            nremoved = 1;
            while nremoved > 0
                nsigneighb    = single(neighbormatrix*onoff);
                remove        = (onoff.*nsigneighb) < minchan;
                nremoved      = length(find(remove.*onoff));
                onoff(remove) = 0;
            end


        end

        % for each channel (combination), find the connected time-frequency clusters
        labelmat = zeros(size(onoff));
        total = 0;
    end
end


