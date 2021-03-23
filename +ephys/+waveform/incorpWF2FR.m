% incorporate wf into FR
%addpath('E:\prJ\neuropixels\learning\pixels-masterUpdate\jpsth\+ephys\+waveform')
fl = dir('E:\prJ\neuropixels\learning\dataSumLN\*\*\FR_All.hdf5');
for onefile=fl'
    disp(onefile.folder)
    FR_File = fullfile(onefile.folder,'FR_All.hdf5');
    FR_All = h5read(FR_File,'/FR_All');
    Trials = h5read(FR_File,'/Trials');
    SU_id = h5read(FR_File,'/SU_id');
    if exist(FR_File,'file')
        delete(FR_File)
    end
    fr_wf_good = ephys.waveform.goodWaveform(onefile.folder,'presel',SU_id);
    if isempty(fr_wf_good)
        disp('Empty!')
        keyboard()
        continue
    end
    h5create(FR_File,'/FR_All',size(FR_All),'Datatype','double')
    h5write(FR_File,'/FR_All',FR_All)
    h5create(FR_File,'/Trials',size(Trials),'Datatype','double')
    h5write(FR_File,'/Trials',Trials)
    h5create(FR_File,'/SU_id',size(SU_id),'Datatype','double')
    h5write(FR_File,'/SU_id',SU_id)
    h5create(FR_File,'/WF_good',size(SU_id),'Datatype','int8')
    h5write(FR_File,'/WF_good',int8(ismember(SU_id,fr_wf_good)))
end