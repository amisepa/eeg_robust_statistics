function EEGout = eeg_icainterp(EEG,urchanlocs)
% function EEGout = eeg_icainterp(EEG,urchanlocs)
%
% Interpolate missing channels of EEG.icawinv (ICA maps) and EEG data
%
% Inputs:
%
%     EEG              -     a dataset with ICA from data with rejected / missing channels
%     urchanlocs  -     array of channel locations with missing channels
%
% Outputs:
%
%     EEGout        -    dataset with ICA maps and data with rejcted /
%                               missing channels interpolated
%

% get the original channel numbers
for k = 1:EEG.nbchan
    v(k) = EEG.chanlocs(k).urchan;
end

% create fake dataset to interpolate ica maps
EEGtmp = EEG;
EEGtmp.data = EEGtmp.icawinv;
EEGtmp.trials = 1;
EEGtmp.pnts = size(EEGtmp.icawinv,2);
EEGtmp2 = eeg_interp(EEGtmp,urchanlocs);

% copy the icaact and interpolated icawinv into output
EEGout = EEG;
EEGout.nbchan = length(urchanlocs);
EEGout.icawinv = EEGtmp2.data;
EEGout.icaact = EEG.icaact;

% add zero columns to icasphere where (new) interpolated data channels are
EEGout.icasphere(:,v) = EEG.icasphere;
EEGout.icasphere(:,setdiff(1:EEGout.nbchan,v)) = 0;

% add data with interpolated channels
EEGout.data = reshape(EEGout.icawinv * double(EEGout.icaact(:,:)),EEGout.nbchan,EEGout.pnts,EEGout.trials);

% copy full chanlocs and icachansind
EEGout.chanlocs = urchanlocs;
EEGout.icachansind = 1:EEGout.nbchan;