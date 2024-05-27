%% Cluster-correction of mass-univariate EEG data
% 
% Cedric Cannard, Sep 2022

function [mask, pcorr] = correct_cluster(tvals, pvals, tvals_H0, pvals_H0, neighbormatrix, mcctype, pthresh)

tvals = tvals.^2;
tvals_H0 = tvals_H0.^2;

% if 1 channel, force switch to 1D TFCE clustering
if size(tvals,1) == 1
	warning("only one EEG channel detected. Switching to TFCE correction");
    mcctype = 3; % TFCE 
end

% nb of boostrap performed
nboot = size(tvals_H0,3);      

% spatiotemporal clustering
if mcctype == 2 && size(tvals_H0,1) > 1
    minchan = 2; % the minimum number of neighbouring channels/combinations
    boot_maxclustersum = zeros(nboot,1);     % maximum cluster mass at each bootstrap
    disp('Getting spatiotemporal clusters under H0 ...');

    parfor boot = 1:nboot
        % Find the cluster, thresholding H0 pvalues <= threshold p
        sigPvals = pvals_H0(:,:,boot) <= pthresh;
        spatdimlength = size(sigPvals, 1);
        nfreq         = size(sigPvals, 2);
        ntime         = size(sigPvals, 3);

        if length(size(neighbormatrix))~=2 || ~all(size(neighbormatrix)==spatdimlength)
            error('Invalid dimension of neighbormatrix');
        end

        % Calculate number of significant neighbors for this channel, for
        % every time/frequency element. If less than minchan, the channel
        % is removed from sigPvals.
        if minchan > 0
            selectmat = single(neighbormatrix | neighbormatrix');
            nremoved = 1;
            while nremoved > 0
                nsigneighb    = reshape(selectmat*reshape(single(sigPvals),[spatdimlength (nfreq*ntime)]),[spatdimlength nfreq ntime]);
                remove        = (sigPvals.*nsigneighb)<minchan;
                nremoved      = length(find(remove.*sigPvals));
                sigPvals(remove) = 0;
            end
        end

        % for each channel (combination), find the connected time-frequency clusters
        labelmat = zeros(size(sigPvals));
        total = 0;
        
        % Label-connected components in binary image.
        if exist('bwlabeln','file') == 2
            for spatdimlev = 1:spatdimlength
                [labelmat(spatdimlev, :, :), nClusters] = bwlabeln(reshape(sigPvals(spatdimlev, :, :), nfreq, ntime), 4);
                labelmat(spatdimlev, :, :) = labelmat(spatdimlev, :, :) + (labelmat(spatdimlev, :, :)~=0)*total;
                total = total + nClusters;
            end            
        else
            error('You need the Image Processing Toolbox to do clustering');
        end

        % combine the time and frequency dimension for simplicity
        labelmat = reshape(labelmat, spatdimlength, nfreq*ntime);

        % combine clusters that are connected in neighbouring channel(s) (combinations).
        replaceby = 1:total;
        for spatdimlev = 1:spatdimlength
            neighbours = find(neighbormatrix(spatdimlev,:));
            for nbindx = neighbours
                indx = find((labelmat(spatdimlev,:)~=0) & (labelmat(nbindx,:)~=0));
                for i = 1:length(indx)
                    a = labelmat(spatdimlev, indx(i));
                    b = labelmat(nbindx, indx(i));
                    if replaceby(a)==replaceby(b)
                        % do nothing
                        continue;
                    elseif replaceby(a)<replaceby(b)
                        % replace all entries with content replaceby(b) by replaceby(a).
                        replaceby(replaceby==replaceby(b)) = replaceby(a);
                    elseif replaceby(b)<replaceby(a)
                        % replace all entries with content replaceby(a) by replaceby(b).
                        replaceby(replaceby==replaceby(a)) = replaceby(b);
                    end
                end
            end
        end

        % renumber all the clusters
        nClusters = 0;
        cluster = zeros(size(labelmat));
        uniquelabel = unique(replaceby(:))';

        for j = 1:length(uniquelabel)
            nClusters = nClusters+1;
            uniquelabel_here = find(replaceby==uniquelabel(j));
            labelmat_is_unique_here = ismember(labelmat(:),uniquelabel_here);  
%             labelmat_is_unique_here = ismembc(labelmat(:),uniquelabel_here); % axs - Using the undocumented ismembc is x3 faster than ismember
            % For ismembc to work, uniquelabel_here MUST be sorted low-high. This should be the case.
            cluster(labelmat_is_unique_here) = nClusters;
        end

        % reshape the output to the original format of the data
        cluster = reshape(cluster, spatdimlength, nfreq, ntime);

        % Compute the mass for each cluster
        bootM_b = tvals_H0(:,:,boot);
        if nClusters~=0
            tmp = zeros(1,nClusters);
            for C = 1:nClusters
                tmp(C) = sum(bootM_b(cluster==C)); % sum stat value in a cluster label
            end
            boot_maxclustersum(boot) = max(tmp(:)); % save max value only
        else
            boot_maxclustersum(boot) = 0;
        end
    end

    % find clusters in the observed data
    [posclusterslabelmat,nposclusters] = limo_findcluster(pvals<=pthresh,neighbormatrix,minchan);
    
    % bootstrap parameters
    nboot = length(boot_maxclustersum);
    sort_clustermax = sort(boot_maxclustersum);
    n = sum(isnan(sort_clustermax)); 
    
    % NaN present if there was no clusters under H0 - just warp around
    % should not happen if using limo_getclustersum as it returns 0 in that case
    if n == nboot
        error('no cluster max value found - check boot_maxclustersum parameter')
    else
        sort_clustermax(isnan(sort_clustermax))=[];
        sort_clustermax = [NaN(n,1); sort_clustermax];
    end
    max_th = sort_clustermax(round((1-pthresh)*nboot));
    fprintf('cluster mass threshold: %g\n',max_th)
    
    % compute the mask: for each cluster do the sum and set significant if > maxclustersum_th
    mask = zeros(size(tvals));
    cluster_label = 1; 
    if nposclusters~=0
        for C = nposclusters:-1:1 % compute cluster sums & compare to bootstrap threshold
            maxval(C) = sum(tvals(posclusterslabelmat==C));
            if  maxval(C)>= max_th
                mask(posclusterslabelmat==C)= cluster_label; % flag clusters above threshold
                cluster_label = cluster_label+1;
            end
        end
    end
    
    % compute corrected p-values: number of times observed mass > bootstrap
    mask2 = logical(mask);
    pcorr  = ones(size(mask));
    if any(mask2(:))
        L = posclusterslabelmat.*mask2;     % remove non significant clusters
        CL_list = setdiff(unique(L),0);     % remove label 0
        for CL = 1:length(CL_list)
            cluster_mass = sum(tvals(L==CL_list(CL)));
            if any(cluster_mass == maxval) % double checking this is in the mask
                p = 1-sum(cluster_mass >= sort_clustermax)./nboot;
                if p ==0
                    p = 1/nboot; % never 0
                end
                tmp = ones(size(mask));
                tmp(L==CL_list(CL)) = p; % set p-values for many cells
                pcorr = pcorr.*tmp;   % tmp is never at the same location so we can just add values
            else
                error('Cannot find the cluster mass for p-value? while found for the mask? something is seriously wrong')
            end
        end
    end
    pcorr(pcorr==1) = NaN;
    
    % Check again corrected p values are correct
    if any(pcorr > pthresh)
       error('some corrected p-values are above the set alpha value, which should not happen - this is a bug') 
    end
    
    % just for output
    if exist('maxval','var')
         maxval = max(maxval);   % biggest cluster
    else
         maxval = 0;
    end
end

% Temporal clustering only (1 electrode)
if mcctype == 2 && size(tvals_H0,1) == 1 || mcctype == 3
    disp('Applying temporal clustering...')
    
    % 1st get the distribution of maxima under H0
    [th, boot_maxclustersum] = limo_ecluster_make(squeeze(tvals_H0),squeeze(pvals_H0),pthresh);
    max_th = th.elec;

    % 2nd threshold observed data
    [sigcluster, pcorr, maxval] = limo_ecluster_test(squeeze(tvals),squeeze(pvals),th,pthresh, boot_maxclustersum);
    mask = sigcluster.elec_mask;
end

% Plot
if sum(mask(:)) == 0
    figure('Name','Cluster correction under H0')
    mass = sort(boot_maxclustersum);
    plot(mass,'LineWidth',3); grid on; hold on;

    plot(find(mass == max_th,1), max_th, 'r*', 'LineWidth',5)
    txt = ['bootstrap threashold ' num2str(max_th) '\rightarrow'];
    text(find(mass == max_th,1), max_th, txt, 'FontSize', 10, 'HorizontalAlignment','right');

    [~,loc] = min(abs(mass-maxval));
    plot(loc,maxval,'r*','LineWidth',5)
    txt = ['Biggest cluster mass observed: ' num2str(maxval) '\rightarrow'];
    text(loc,double(maxval),txt,'FontSize',10,'HorizontalAlignment','right');

    title('Cluster-mass Maxima under H0','FontSize',12)
    xlabel('sorted bootstrap iterations','FontSize',12);
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end
