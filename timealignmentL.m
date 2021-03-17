% INPUT:
%  --referenceRoute: the folder of reference imec recording system, using to
%  index events file.
%  --unCorrectRoute: the folder for file need to be corrected
%  --saveRoute: folder that save the correct file
% OUTPUT:
%  a spike_info.hdf5 file in saveRoute
%  
% DISCRIPTION OF spike_info.hdf5
% /imec*/times : time points after correct
% /imec*/clusters : cluster ID corresponding to the times in /imec*/times
%
% PRINCIPLE:
% using behavioural time point file events.hdf5 in referenceRoute as
% reference time aiming to make corrected behavioural time point equal 
% to the reference time.
%
% NOTICE: Library of npy-matlbab need to be added
% different events.hdf5 file may have differenct trial amounts.
function timealignmentL(referenceRoute, unCorrectRoute, saveRoute)
    eventReferenceRoute=fullfile(referenceRoute,'eventsRescue.hdf5');
    saveRoute=fullfile(saveRoute,'spike_info.hdf5');
    if exist(saveRoute,'file')
        delete(saveRoute);
    end
    spkTS_reference=double(readNPY(fullfile(referenceRoute,'/spike_times.npy')));
%     cluster_reference=double(readNPY(fullfile(referenceRoute,'/spike_clusters.npy')));
    imecnum=referenceRoute(1,strfind(referenceRoute,'imec')+4);
    cluster_reference=double(readNPY(fullfile(referenceRoute,'/spike_clusters.npy')))+str2num(imecnum)*10000;
    h5create(saveRoute,['/imec',imecnum,'/times'],size(spkTS_reference)); 
    h5write(saveRoute,['/imec',imecnum,'/times'],spkTS_reference);
    h5create(saveRoute,['/imec',imecnum,'/clusters'],size(cluster_reference));
    h5write(saveRoute,['/imec',imecnum,'/clusters'],cluster_reference);
    trials_reference=double(h5read(eventReferenceRoute,'/trials'));
    for i=1:size(unCorrectRoute,1)
        trials=double(h5read(fullfile(unCorrectRoute{i,1},'eventsRescue.hdf5'),'/trials'));
        spkTS=double(readNPY(fullfile(unCorrectRoute{i,1},'spike_times.npy')));
        spk=[];   % correct spike time stamp
        imecnum=unCorrectRoute{i,1}(1,strfind(unCorrectRoute{i,1},'imec')+4);
        tmax=min(size(trials,2),size(trials_reference,2));
        for t=1:tmax
            [~,index]=min(abs(trials_reference(1,:)-trials(1,t)));
            trials_reference0=trials_reference(1,index);
            if t==1
                spk=spkTS(spkTS<trials(1,t+1))-(trials(1,t)-trials_reference0);
            elseif t<tmax
                spk=[spk;spkTS(spkTS<trials(1,t+1) & spkTS>=trials(1,t))-(trials(1,t)-trials_reference0)];
            else
                spk=[spk;spkTS(spkTS>=trials(1,t))-(trials(1,t)-trials_reference0)];
            end
        end
        cluster=double(readNPY(fullfile(unCorrectRoute{i,1},'spike_clusters.npy')));
        delindex=find(spk<=0);
        spk(delindex,:)=[];
        cluster(delindex,:)=[];
        cluster=cluster+10000*str2num(imecnum);
        h5create(saveRoute,['/imec',imecnum,'/times'],size(spk));
        h5write(saveRoute,['/imec',imecnum,'/times'],spk);
        h5create(saveRoute,['/imec',imecnum,'/clusters'],size(cluster));
        h5write(saveRoute,['/imec',imecnum,'/clusters'],cluster);
        disp([saveRoute,' has finished!']);
    end
end