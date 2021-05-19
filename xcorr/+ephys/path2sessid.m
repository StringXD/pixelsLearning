function out=path2sessid(path,opt)
arguments
    path (1,:) char
    opt.type (1,:) char {mustBeMember(opt.type,{'4n','2'})}
end
persistent map1 map2
if isempty(map1) & strcmp(opt.type,'4n')
    homedir=ephys.util.getHomedir('dtype',opt.type);
    fullpath=deblank(h5read(fullfile(homedir,'transient_6.hdf5'),'/path'));
    fullpath = cellfun(@(x) strsplit(x,'\'),fullpath,'Uniformoutput',false);
    allpath = cellfun(@(x) x(end-1),fullpath);
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map1=containers.Map('KeyType','char','ValueType','int32');
    for i=1:numel(uidx)
        map1(upath{i})=uidx(i);
    end
elseif isempty(map2) & strcmp(opt.type,'2')
    homedir=ephys.util.getHomedir('dtype',opt.type);
    fullpath=deblank(h5read(fullfile(homedir,'transient_6_2tracks0514.hdf5'),'/path'));
    fullpath = cellfun(@(x) strsplit(x,'\'),fullpath,'Uniformoutput',false);
    allpath = cellfun(@(x) x(end-1),fullpath);
    upath=unique(allpath);
    [upath,uidx]=sort(upath);
    map2=containers.Map('KeyType','char','ValueType','int32');
    for i=1:numel(uidx)
        map2(upath{i})=uidx(i);
    end
end
paths=split(path,["/","\"]);
if strcmp(opt.type,'4n')
    out=map1(char(paths(end-1)));
elseif strcmp(opt.type,'2')
    out=map2(char(paths(end-1)));
end
