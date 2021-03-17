% find events with minimal trials as reference
addpath('D:\code\npy-matlab\npy-matlab');
routeinfo = dir('D:\dataSumLN\*');
SSroute=[];
for routeindex=3:size(routeinfo,1)
    route1=fullfile(routeinfo(routeindex,1).folder,routeinfo(routeindex,1).name);
    if contains(route1,'M81_6_learning_20200909_g0') % problematic session
        continue;
    end
    route1info=dir(route1);
    ii = 0;
    standardroute={};
    for i=3:size(route1info,1)
        if contains(route1info(i,1).name,'cleaned')
            ii=ii+1;
            standardroute{ii,1}=fullfile(route1info(i,1).folder,route1info(i,1).name);
        end
    end
    imecEvt = [];
    if ii >= 2
        for i = 1:ii
            if contains(standardroute{i,1},'imec')
                imecEvt = [imecEvt;str2num(standardroute{i,1}(1,strfind(standardroute{i,1},'imec')+4)),length(h5read(fullfile(standardroute{i,1},'eventsRescue.hdf5'),'/trials'))];
                if exist(fullfile(standardroute{i,1},'spike_info.hdf5'),'file')
                    delete(fullfile(standardroute{i,1},'spike_info.hdf5'));
                end
            end
        end
        [~,index]=sort(imecEvt(:,2));
        referenceRoute=standardroute{index(1),1};
        saveRoute=referenceRoute;
        unCorrectRoute=standardroute(index(2:end),1);
        timealignmentL(referenceRoute, unCorrectRoute, saveRoute);
        saveRoute=[deblank(saveRoute),'\spike_info.hdf5'];
        SSroute=[SSroute;{saveRoute}];
    else
        % No treatment for single probe session
        continue
    end
end

