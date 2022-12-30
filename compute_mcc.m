%% Compute multiple comparison correction (Type 1 error or FWE).
%
% Cedric Cannard, Sep 2022

function [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, pthresh, neighbormatrix)

nClust = [];

switch mcctype

    case 0  % Uncorrected
        pcorr = pvals;
        mask = pcorr <= pthresh;
        nClust = length(unique(mask)) - 1; % number of significant clusters
        if nClust > 0
            disp([num2str(nClust) ' significant clusters.'])
        end

    case 1  % Max-correction
        [mask, pcorr] = correct_max(abs(tvals),abs(tvals_H0),pthresh,1); % limo_max_correction(abs(M),abs(bootM),p)

    case 2 % Cluster-correction

        % For one channel only
        if size(tvals,1) == 1
            tmp = NaN(1,size(tvals,2),size(tvals_H0,2));
            tmp(1,:,:,:) = tvals_H0; tvals_H0 = tmp;
            tmp(1,:,:,:) = pvals_H0; pvals_H0 = tmp;
            clear tmp
        end

        % Get cluster mask and corrected p-values
        [mask, pcorr] = correct_cluster(tvals.^2, pvals, tvals_H0.^2, pvals_H0, neighbormatrix, mcctype, pthresh); % limo_clustering(M.^2,Pval,bootM.^2,bootP,LIMO,MCC,p)
        nClust = length(unique(mask)) - 1; % number of significant clusters
        if nClust > 0
            disp([num2str(nClust) ' significant clusters.'])
        end

        % References
        disp(' ');
        disp('Refs for Clustering & Bootstrap:')
        disp('Maris & Oostenveld (2007), Nonparametric statistical testing of EEG- and MEG- data.')
        disp('Journal of Neuroscience Methods.')
        disp(' ');
        disp('Pernet, Latinus, Nichols, & Rousselet, (2015).')
        disp('Cluster-based computational methods for mass univariate analyses')
        disp('of event-related brain potentials/fields: A simulation study.')
        disp('Journal of Neuroscience methods.')

    case 3  % TFCE-correction
        ndim = 2; % 1 (one channel); 2 (ERP, power, or single time-freq map); 3 (ERSP)

        % Thresholding
        tfce_score = limo_tfce(ndim, tvals, neighbormatrix);

        % Apply TFCE to null data
        disp('Applying TFCE to null data. This may take a while ...')
        nboot = size(tvals_H0,3);
        tfce_H0_thmaps = cell(1,nboot);
        tfce_H0_score = nan(size(tvals_H0,1),size(tvals_H0,2),nboot);
        progressbar('TFCE thresholding null data')
        for b = 1:nboot
            disp(['boot ' num2str(b) '/' num2str(nboot)])
            [tfce_H0_score(:,:,b), tfce_H0_thmaps{b}] = limo_tfce(ndim, squeeze(tvals_H0(:,:,b)), neighbormatrix,0);
            progressbar(b/nboot)
        end

        % Max-correction on TFCE data
        [mask, pcorr] = correct_max(tfce_score,tfce_H0_score,pthresh,1); % adapted from limo_max_correction
        nClust = length(unique(mask))-1;
        if nClust > 0
            disp([num2str(nClust) ' significant clusters.'])
        end

        % References
        disp(' ');
        disp('Ref for TFCE:')
        disp('Pernet, Latinus, Nichols, & Rousselet (2015).')
        disp('Cluster-based computational methods for mass univariate analyses')
        disp('of event-related brain potentials/fields: A simulation study.')
        disp('Journal of Neuroscience methods')
        disp(' ');

end
