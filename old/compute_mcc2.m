%% Compute multiple comparison correction (Type 1 error or FWE).
%
% Cedric Cannard, Sep 2022

function [mask, pcorr] = compute_mcc2(tvals, pvals, tvals_H0, pvals_H0, pairs, mcctype, pthresh, neighbors, neighbormatrix)

% Uncorrected
if mcctype == 1
    pcorr = pvals;
    mask = pvals <= pthresh;
    mytitle = 'Uncorrected';

    % Max-correction
elseif mcctype == 2
    try
%         [mask,M] = limo_max_correction(abs(M),abs(bootM),p);
        [mask, pcorr] = correct_max(abs(tvals),abs(tvals_H0),pthresh,1);
    catch ME
        errordlg(sprintf('error log: %s \n',ME.message),'Max-correction failure')
        return
    end

    % Cluster-correction
elseif mcctype == 3

    disp('Refs for Clustering & Bootstrap:')
    disp('Maris & Oostenveld (2007), Nonparametric statistical testing of EEG- and MEG- data.')
    disp('Journal of Neuroscience Methods.')
    disp(' ');
    disp('Pernet, Latinus, Nichols, & Rousselet, (2015).')
    disp('Cluster-based computational methods for mass univariate analyses')
    disp('of event-related brain potentials/fields: A simulation study.')
    disp('Journal of Neuroscience methods.')
    
    % for one channel only
    if size(tvals,1) == 1
        tmp = NaN(1,size(tvals,2),size(tvals_H0,2));
        tmp(1,:,:,:) = tvals_H0; tvals_H0 = tmp;
        tmp(1,:,:,:) = pvals_H0; pvals_H0 = tmp;
        clear tmp
    end

    % get cluster mask and corrected p-values
    % [mask,M] = limo_clustering(M.^2,Pval,bootM.^2,bootP,LIMO,MCC,p); % mask and cluster p values
%     [mask, pcorr] = correct_cluster(tvals, pvals, tvals_H0, pvals_H0, neighbormatrix, mcctype, pthresh, 1);
    [mask, pcorr] = correct_cluster2(tvals, pvals, tvals_H0, pvals_H0, pairs, neighbors, neighbormatrix, mcctype, pthresh, 1); %connectivity data

    % Count number of significant clusters to display in title
    Nclust = length(unique(mask)) - 1;
    if Nclust <= 1, Mclust = 'cluster'; else, Mclust = 'clusters'; end

    % TFCE-correction
elseif ~isempty(tvals) && mcctype == 4
    try
        % Thresholding
        % tfce_score = limo_tfce(2, tvals,neighbormatrix);  % 2D (power), t-values, channel neighbor matrix
        tfce_score = limo_tfce(2,tvals,neighbormatrix); 

        % Apply TFCE to null data
        disp('Applying TFCE to null data. This may take a while ...')
%         H0_tval = squeeze(results_H0(:,:,1,:));  % t-values
        nboot = size(tvals_H0,3);
        tfce_H0_thmaps = cell(1,nboot);
        tfce_H0_score = nan(size(tvals_H0,1),size(tvals_H0,2),nboot);
        progressbar('TFCE thresholding null data')
        parfor b = 1:nboot
            disp(['boot ' num2str(b) '/' num2str(nboot)])
            [tfce_H0_score(:,:,b), tfce_H0_thmaps{b}] = limo_tfce(2, squeeze(tvals_H0(:,:,b)), neighbormatrix,0);

            progressbar(b,nboot)
        end
        
        % max-correction on TFCE data
%         mask = limo_max_correction(tfce_score, tfce_H0_score, pthresh, 1);
        [mask, pcorr] = correct_max(tfce_score,tfce_H0_score,pthresh,1);
        Nclust = length(unique(mask))-1;
        if Nclust <= 1, Mclust = 'cluster'; else, Mclust = 'clusters'; end
    catch ME
        errordlg(sprintf('error log: %s \n',ME.message),'tfce correction failure')
        return
    end

    disp('Ref for TFCE:')
    disp('Pernet, Latinus, Nichols, & Rousselet (2015).')
    disp('Cluster-based computational methods for mass univariate analyses')
    disp('of event-related brain potentials/fields: A simulation study.')
    disp('Journal of Neuroscience methods')

end
