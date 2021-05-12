if isfile(sprintf('0511_BZ_XCORR_duo_f%d.mat',fidx))
    disp('File exist'); if isunix, quit(0); else, return; end
end
disp(fidx)
fl = dir('/home/xd/data/learning/*/*/*/spike_info.hdf5');
addpath(genpath('~/buzcode'));
[spkID,spkTS,~,~,folder]=ephys.getSPKID_TS(fidx);
if isempty(spkID)
    if isunix, quit(0); else, return; end
end
mono=bz.sortSpikeIDz(spkTS,spkID); % adapted from English, Buzsaki, 2017
%bz.util.plotCCG
save(sprintf('0511_BZ_XCORR_duo_f%d.mat',fidx),'mono','fl','-v7.3')
if isunix, quit(0); else, return; end