function homedir=getHomedir(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'sums','raw'})} = 'sums'
    opt.dtype (1,:) char {mustBeMember(opt.dtype,{'4n','2'})}='4n'
    
end
if strcmp(opt.dtype,'4n') | strcmp(opt.dtype,'2')
    if ispc
        if strcmp(opt.type,'sums')
            homedir = fullfile('K:','code','per_sec');
        elseif strcmp(opt.type,'raw')
            homedir = fullfile('K:','neupix','SPKINFO');
        end
    elseif isunix
        if strcmp(opt.type,'sums')
            homedir = fullfile('~','data','learning');
        elseif strcmp(opt.type,'raw')
            homedir = fullfile('~','data','learning','learning2tracks');
        end
    end
else
    if ispc
        if strcmp(opt.type,'sums')
            homedir = fullfile('K:','neupix','AIOPTO','META');
        elseif strcmp(opt.type,'raw')
            homedir = fullfile('K:','neupix','AIOPTO','RECDATA');
        end
    elseif isunix
        if strcmp(opt.type,'sums')
            homedir = fullfile('~','neupix','AIOPTO','META');
        elseif strcmp(opt.type,'raw')
            homedir = fullfile('~','neupix','AIOPTO','RECDATA');
        end
    end
    
end
