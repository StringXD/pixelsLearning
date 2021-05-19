function reg_conn_bz(opt)
arguments
    opt.type (1,:) char {mustBeMember(opt.type,{'4n','2'})}='2'
    opt.data (:,1) struct =[]
    opt.prefix (1,:) char = '0511'
end

%TODO merge with load_sig_pair script
if isempty(opt.data)
    if strcmp(opt.type,'4n')
        load('sums_conn_4n.mat','sums_conn_str')
    else
        load('sums_conn_2.mat','sums_conn_str')
    end
else
    sums_conn_str=opt.data;
end

for fidx=1:length(sums_conn_str)
    disp(fidx);
    fpath=sums_conn_str(fidx).folder; %session data folder
    if strcmp(opt.type,'4n') | strcmp(opt.type,'2')
        pathsplits = strsplit(fpath,'/');
        pc_stem=fullfile(pathsplits{end-1},pathsplits{end});
        inputf=fullfile(fpath,'FR_All.hdf5');
        all_su=int32(h5read(inputf,'/SU_id'));
    else
        pc_stem=fpath;
        inputf=fullfile('K:','neupix','AIOPTO','RECDATA',fpath,'FT_SPIKE.mat');
        fstr=load(inputf);
        all_su=int32(cellfun(@(x) str2double(x),fstr.FT_SPIKE.label));
    end
    
    sig_con=int32(sums_conn_str(fidx).sig_con); % significant functional coupling
    pair_comb_one_dir=nchoosek(all_su,2); % all pairs combination
    [sig_meta,pair_meta]=bz.util.get_meta(sig_con,pair_comb_one_dir,pc_stem,'type',opt.type); % assign meta info
    
%     fields={'suid','reg','wrsp','selec','mem_type'};
    fields={'suid','reg','mem_type'};
    for fi=fields
        %TODO online genenrate session tag
%         sig.(fi{1})=cat(1,sig.(fi{1}),sig_meta.(fi{1}));
        pair_meta.(fi{1})=cat(1,pair_meta.(fi{1}),flip(pair_meta.(fi{1}),ndims(pair_meta.(fi{1}))));%uni-dir to bi-dir
%         pair.(fi{1})=cat(1,pair.(fi{1}),pair_meta.(fi{1}));
    end
    tic
    save(fullfile('xcorr_bz_2',sprintf('%s_conn_w_reg_2_%d.mat',opt.prefix,fidx)),'sig_meta','pair_meta','pc_stem','-v7','-nocompression')
    toc
end
end