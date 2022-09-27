%% Cluster-correction of mass-univariate EEG data
% 
% Cedric Cannard, Sep 2022

function mask = correct_cluster(tvals, pvals, tboot, pboot, neighbormatrix, mcc, p, fig)

% [mask,M] = limo_clustering(M.^2,Pval,bootM.^2,bootP,LIMO,MCC,p); % mask and cluster p values

fig = 0;
tvals = tvals.^2;
tboot = tboot.^2;

% if 1 channel, force switch to 1D TFCE clustering
if size(tvals,1) == 1
    mcc = 4; % TFCE 
end

% nb of boostrap performed
nboot = size(tboot,3);      

% spatiotemporal clustering
if mcc == 3 && size(tboot,1) > 1
    disp('Applying spatiotemporal clustering...')
    minchan = 2;
    boot_maxclustersum = zeros(nboot,1);     % maximum cluster mass at each bootstrap

    disp('getting clusters under H0 boot ...');
    parfor boot = 1:nboot
        % Find the cluster, thresholding H0 pvalues <= threshold p
        [posclusterslabelmat,nposclusters] = limo_findcluster((pboot(:,:,boot) <= p),neighbormatrix,minchan);

        % 2nd compute the mass for each cluster
        bootM_b = tboot(:,:,boot);
        if nposclusters~=0
            tmp = zeros(1,nposclusters);
            for C = 1:nposclusters
                tmp(C) = sum(bootM_b(posclusterslabelmat==C)); % sum stat value in a cluster label
            end
            boot_maxclustersum(boot) = max(tmp(:)); % save max value only
        else
            boot_maxclustersum(boot) = 0;
        end
    end

    % 3rd threshold observed cluster mass by the distribution of cluster
    % max computed in step 2
    [mask, cluster_pval, maxval, max_th] = limo_cluster_test(tvals,pvals,...
        boot_maxclustersum,neighbormatrix,minchan,p);
end

% temporal clustering
if mcc == 3 && size(tboot,1) == 1 || mcc == 4
    disp('Applying temporal clustering...')
    % 1st get the distribution of maxima under H0
    [th,boot_maxclustersum] = limo_ecluster_make(squeeze(tboot),squeeze(pboot),p);
    max_th                  = th.elec;
    % 2nd threshold observed data
    [sigcluster, cluster_pval,maxval] = limo_ecluster_test(squeeze(tvals),squeeze(pvals),th,p, boot_maxclustersum);
    mask                              = sigcluster.elec_mask;
end

%Plot
% when nothing is significant, always show why
if sum(mask(:)) == 0
    fig = 1 ;
end

if fig == 1
    figure('Name','Cluster correction under H0')
    mass = sort(boot_maxclustersum);
    plot(mass,'LineWidth',3); grid on; hold on;

    plot(min(find(mass==max_th)),max_th,'r*','LineWidth',5)
    txt = ['bootstrap threashold ' num2str(max_th) '\rightarrow'];
    text(min(find(mass==max_th)),max_th,txt,'FontSize',10,'HorizontalAlignment','right');

    [val,loc]=min(abs(mass-maxval));
    plot(loc,maxval,'r*','LineWidth',5)
    txt = ['biggest observed cluster mass: ' num2str(maxval) '\rightarrow'];
    text(loc,double(maxval),txt,'FontSize',10,'HorizontalAlignment','right');

    title('Cluster-mass Maxima under H0','FontSize',12)
    xlabel('sorted bootstrap iterations','FontSize',12);
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end
