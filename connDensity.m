gen_join_set = true;

for condition = {'nonsel'}
prefix = condition{1,1};

for bin = [1,5:6]
    load(sprintf('full_duo_XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
    if bin == 1
        if gen_join_set
            join_reg_list=cell(0);
            for sidx=1:length(stats)
                s=stats{sidx};
                join_reg_list{end+1}=s.reg_su1;
                join_reg_list{end+1}=s.reg_su2;
            end
            join_reg_list(cellfun(@isempty,join_reg_list))=[];
            join_reg_set=unique(join_reg_list);
        end
        reg_set=join_reg_set(~strcmp(join_reg_set,'Unlabeled'));
        reg_set=reg_set(~strcmp(reg_set,'root'));
        save('reg_keep.mat','reg_set');
    end

    if bin == 1
        continue;
    end
    for i=1:length(stats)
        s=stats{i};
        s.uid1=s.fileidx*100000+s.su1_clusterid;
        s.uid2=s.fileidx*100000+s.su2_clusterid;
        stats{i}=s;
    end

    pair_chain=zeros(0,2);
    pair_reg=zeros(0,2);
    pref_pair=zeros(0,12);

    conn_chain_S1=zeros(0,2);
    reg_chain_S1=zeros(0,2);
    pref_chain_S1=zeros(0,12);
    peaks1=zeros(0,0);
    ais1=zeros(0,0);

    conn_chain_S2=zeros(0,2);
    reg_chain_S2=zeros(0,2);
    pref_chain_S2=zeros(0,12);
    peaks2=zeros(0,0);
    ais2=zeros(0,0);

    for pidx=1:length(stats)
        s=stats{pidx};

        su1reg_idx=find(strcmp(s.reg_su1,reg_set));
        su2reg_idx=find(strcmp(s.reg_su2,reg_set));

        if isempty(su1reg_idx) || isempty(su2reg_idx)
            continue;
        end

        if s.s1_trials<20 || s.s2_trials<20
            continue
        end

        if strcmp(prefix,'sel')
            if strcmp(s.su1_sel_type,'nonsel') || strcmp(s.su2_sel_type,'nonsel')
                continue;
            end
        elseif strcmp(prefix,'nonsel')
            if ~(strcmp(s.su1_sel_type,'nonsel') && strcmp(s.su2_sel_type,'nonsel'))
                continue;
            end
        end

        pair_chain(end+1,:)=[s.uid1,s.uid2];
        pair_reg(end+1,:)=[su1reg_idx,su2reg_idx];
        pref_pair(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];

%            if isempty(su1reg_idx) || isempty(su2reg_idx)
%                fprintf('%s, %s\n', s.reg_su1, s.reg_su2);
%                continue
%            end
        if s.totalcount<250
            continue
        end

        %%S1
        if isfield(s,'s1_peak_significant') && s.s1_peak_significant
            if isfield(s,'AIs1') && s.AIs1>0.4
                conn_chain_S1(end+1,:)=[s.uid1,s.uid2];
                reg_chain_S1(end+1,:)=[su1reg_idx,su2reg_idx];
                pref_chain_S1(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                %peaks1(end+1)=s.peaks1;
                ais1(end+1)=s.AIs1;
            elseif isfield(s,'AIs1') && s.AIs1<-0.4
                conn_chain_S1(end+1,:)=[s.uid2,s.uid1];
                reg_chain_S1(end+1,:)=[su2reg_idx,su1reg_idx];
                pref_chain_S1(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                %peaks1(end+1)=s.peaks1;
                ais1(end+1)=-s.AIs1;
            end
        end
        if isfield(s,'s2_peak_significant') && s.s2_peak_significant
            if isfield(s,'AIs2') && s.AIs2>0.4
                conn_chain_S2(end+1,:)=[s.uid1,s.uid2];
                reg_chain_S2(end+1,:)=[su1reg_idx,su2reg_idx];
                pref_chain_S2(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
                %peaks2(end+1)=s.peaks2;
                ais2(end+1)=s.AIs2;
            elseif isfield(s,'AIs2') && s.AIs2<-0.4
                conn_chain_S2(end+1,:)=[s.uid2,s.uid1];
                reg_chain_S2(end+1,:)=[su2reg_idx,su1reg_idx];
                pref_chain_S2(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
                %peaks2(end+1)=s.peaks2;
                ais2(end+1)=-s.AIs2;
            end
        end
%            if s.s1_peak_significant && s.s2_peak_significant 
%               if  s.AIs1>0.4 & s.AIs2>0.4
%                    conn_chain_both(end+1,:)=[s.uid1,s.uid2];
%                    reg_chain_both(end+1,:)=[su1reg_idx,su2reg_idx];
%                    pref_chain_both(end+1,:)=[s.prefered_sample_su1(2:end),s.prefered_sample_su2(2:end)];
%               end
%
%               if  s.AIs1<-0.4 & s.AIs2<-0.4
%                    conn_chain_both(end+1,:)=[s.uid2,s.uid1];
%                    reg_chain_both(end+1,:)=[su2reg_idx,su1reg_idx];
%                    pref_chain_both(end+1,:)=[s.prefered_sample_su2(2:end),s.prefered_sample_su1(2:end)];
%               end
%            end
    end
    disp('check file name')
    disp(bin);
    save(sprintf('%s_conn_chain_duo_6s_%d_%d.mat',prefix,bin,bin+1),'conn_chain_S1','reg_chain_S1','pref_chain_S1','conn_chain_S2','reg_chain_S2','pref_chain_S2','pair_chain','pair_reg','pref_pair','ais1','ais2','peaks1','peaks2');
end
end

%% 3 conditions
conncount=[];
% overall density selec
for bin=1:6
    load(sprintf('sel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    if bin==1
        paircount=nnz(diff(pair_reg,1,2));
    end
    conncount=[conncount,nnz(diff(reg_chain_S1,1,2)),nnz(diff(reg_chain_S2,1,2))];    
end
mmm=mean(conncount)/paircount;
stdm=std(conncount/paircount);

% overall density nonsel
conncount=[];
for bin=1:6
    load(sprintf('nonsel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    if bin==1
        paircount=nnz(diff(pair_reg,1,2));
    end
    conncount=[conncount,nnz(diff(reg_chain_S1,1,2)),nnz(diff(reg_chain_S2,1,2))];
%     conncount=[conncount,length(reg_chain_S1),length(reg_chain_S2)];
end
mmn=mean(conncount)/paircount;
stdn=std(conncount/paircount);


% congruent pairs
conncount=[];
for bin=1:6
    load(sprintf('sel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
%     if bin==1
%         paircount=nnz(diff(pair_reg,1,2));
%     end
    congruS1=nnz(reg_chain_S1(:,1)~=reg_chain_S1(:,2) & pref_chain_S1(:,bin)==pref_chain_S1(:,bin+6) & pref_chain_S1(:,bin)>0);
    congruS2=nnz(reg_chain_S2(:,1)~=reg_chain_S2(:,2) & pref_chain_S2(:,bin)==pref_chain_S2(:,bin+6) & pref_chain_S2(:,bin)>0);
    congruPair=nnz(pair_reg(:,1)~=pair_reg(:,2) & pref_pair(:,bin)==pref_pair(:,bin+6) & pref_pair(:,bin)>0);
    conncount=[conncount,congruS1/congruPair,congruS2/congruPair];
    
end
mmc=mean(conncount);
stdc=std(conncount);

close all
figure('Color','w','Position',[100,100,230,130])
hold on
bar([mmn,mmm,mmc],'FaceColor','k')
errorbar(1:3,[mmn,mmm,mmc],[stdn/sqrt(12),stdm/sqrt(12),stdc/sqrt(12)],'k.','CapSize',16,'Color',[0.5,0.5,0.5])
xlim([0.5,3.5])
ylim([0,0.25]);
set(gca,'XTick',1:3,'XTickLabel',{'non-memory','memory units','congruent units'},'XTickLabelRotation',15)
ylabel('connection density')
saveas(gcf,'connDensity.fig');
exportgraphics(gcf,'pairwise-conn-dens.pdf','ContentType','vector');


