% summarize FR data
%% 2 tracks only
sus_trans2 = h5read('F:\learning2tracks\transient6_learning_0329.hdf5','/transient6');
cluster_id2 = h5read('E:\prJ\neuropixels\learning\pixels-masterUpdate\per_sec\transient_6_2tracks.hdf5','/cluster_id');
paths2 = h5read('E:\prJ\neuropixels\learning\pixels-masterUpdate\per_sec\transient_6_2tracks.hdf5','/path');
paths2 = deblank(paths2);
path0 = '';
dataSum = cell(0,2);
for i = 1:length(cluster_id2)
    pathc = string(paths2(i));
    if ~strcmp(path0,pathc)
        FR_All = h5read(fullfile(pathc,'FR_All.hdf5'),'/FR_All');
        Trials = h5read(fullfile(pathc,'FR_All.hdf5'),'/Trials');
        SU_id = h5read(fullfile(pathc,'FR_All.hdf5'),'/SU_id');
        % fill FR data
        currClusters = cluster_id2(contains(paths2,pathc));
        [~,~,ib] = intersect(currClusters,SU_id,'stable');
        FRs = mat2cell(FR_All,size(FR_All,1),size(FR_All,2),ones(size(FR_All,3),1));
        FRs = squeeze(FRs);
        Trialist = cellfun(@(x) Trials,cell(size(FRs)),'UniformOutput',false);
        dataSum = [dataSum;FRs,Trialist];
        path0 = pathc;
    else
        continue;
    end
end
save('learningFR_2tracks.mat','dataSum','-v7.3');