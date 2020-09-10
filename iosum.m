load('full_duo_XCORR_stats_delay_6_1_2_2msbin.mat');
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

plot_entire=true;
type = 'selec';
plot_per_bin=false;

ioselstats=cell(1,6);
for bin=1:6
    if strcmp(type,'error')
        load(sprintf('correct_error\\0820_error_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    elseif strcmp(type,'correct_resample') && exist('rpt','var')
        load(sprintf('correct_error\\0820_correct_resample_%03d__conn_chain_duo_6s_%d_%d.mat',rpt,bin,bin+1));
    elseif strcmp(type,'baseline')
        if bin>1
            break;
        end
        load('0826_selec_conn_chain_duo_6s_-2_-1.mat');
    else
        load(sprintf('sel_conn_chain_duo_6s_%d_%d.mat',bin,bin+1));
    end
    in_out_sel=nan(length(reg_set),12);
    for reg_idx=1:length(reg_set)
        pair_fw= nnz((pair_reg(:,1)~=reg_idx) & (pair_reg(:,2)==reg_idx));
        pair_rev= nnz((pair_reg(:,1)==reg_idx) & (pair_reg(:,2)~=reg_idx));
        pair_count=pair_fw+pair_rev;
        
        in_conn_S1= nnz((reg_chain_S1(:,1)~=reg_idx) & (reg_chain_S1(:,2)==reg_idx));
        in_sel_S1=nnz((reg_chain_S1(:,1)~=reg_idx) & (reg_chain_S1(:,2)==reg_idx) & pref_chain_S1(:,bin)>0);
        
        out_conn_S1= nnz((reg_chain_S1(:,2)~=reg_idx) & (reg_chain_S1(:,1)==reg_idx));
        out_sel_S1=nnz((reg_chain_S1(:,2)~=reg_idx) & (reg_chain_S1(:,1)==reg_idx) & pref_chain_S1(:,bin+6)>0);
        
        auto_pair=nnz((pair_reg(:,1)==reg_idx) & (pair_reg(:,2)==reg_idx))*2;
        auto_conn_S1= nnz((reg_chain_S1(:,1)==reg_idx) & (reg_chain_S1(:,2)==reg_idx));
        
        in_out_sel(reg_idx,:)=[pair_count,in_conn_S1,in_conn_S1/pair_count, ...%1 2 3
            in_sel_S1,in_sel_S1/pair_count,...% 4 5
            out_conn_S1,out_conn_S1/pair_count,...% 6 7
            out_sel_S1,out_sel_S1/pair_count,...% 8 9
            auto_pair,auto_conn_S1,auto_conn_S1/auto_pair]; % 10 11 12
    end
    
    ioselstats{bin}=in_out_sel;
    if plot_per_bin
        plotOne(in_out_sel,reg_set,sprintf('bin %d',bin),sprintf('0714_io_selec_bin%d.png',bin));
    end
end

%% entire delay
io_entire_delay=ioselstats{1};
for bin=2:6
    io_entire_delay=io_entire_delay+ioselstats{bin};
end
io_entire_delay(:,[3 5 7 9])=io_entire_delay(:,[2 4 6 8])./io_entire_delay(:,1);
io_entire_delay(:,12)=io_entire_delay(:,11)./io_entire_delay(:,10);
if plot_entire
    %     plotOne(io_entire_delay,reg_set,'sum of all bins','0714_io_congruent_bin_sum.png',[4,5,8,9]);
    % plotOne(io_entire_delay,reg_set,'sum of all bins','0714_io_congruent_bin_sum.png',[2,3,11,12]);
    plotOne(io_entire_delay,reg_set,'sum of all bins','learning_io_sel_in_out_bin_sum_.png',[2,3,6,7]);
end


function plotOne(in_out_sel,reg_set,str_title,fname,v_idx)
    if ~exist('v_idx','var')
        v_idx=[4,5,8,9];
    end
    t=num2cell(v_idx);
    [inCount,inFrac,outCount,outFrac]=t{:};
    
    greymatter=cellfun(@(x) ~isempty(regexp(x,'[A-Z]','match','once')), reg_set);
    countsel=in_out_sel(:,1)>=100 & greymatter; %126 115 104
    iosel=in_out_sel(countsel,:);
    sink_sel=false(size(iosel,1),1);
    source_sel=false(size(iosel,1),1);
    for i=1:size(iosel,1)
        currcount=iosel(i,1);
        [~,~,p]=crosstab(1:2*currcount>currcount,[1:currcount>iosel(i,inCount),1:currcount>iosel(i,outCount)]);
        %     chisqsel(i)=p*numel(chisqsel)<0.05; %bonferroni adjust
        if p<0.05 && iosel(i,outCount)>iosel(i,inCount)
            source_sel(i)=true;
        elseif p<0.05 && iosel(i,outCount)<iosel(i,inCount)
            sink_sel(i)=true;
        end
    end
    
    regsel=reg_set(countsel);
    %     regid=zeros(size(regsel));
    %     for i=1:length(regsel)
    %         regid(i)=regclass.id(strcmp(regclass.reg,regsel{i})); %CTX=0,TH=1,STR=2,ELSE=3
    %     end
    %
    fh=figure('Color','w','Position',[100,100,280,280]);
    hold on
    otherh=scatter(iosel(~(sink_sel|source_sel),inFrac),iosel(~(sink_sel|source_sel),outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.5,0.5,0.5],'MarkerFaceAlpha',0.5);
    gainh=scatter(iosel(source_sel,inFrac),iosel(source_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerFaceAlpha',0.5);
    damph=scatter(iosel(sink_sel,inFrac),iosel(sink_sel,outFrac),'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerFaceAlpha',0.5);
    
    % disp('source')
    % disp(regsel(source_sel));
    % disp('sink')
    % disp(regsel(sink_sel));
    %     for i=reshape(find(source_sel | sink_sel),1,[])
    %         text(iosel(i,inFrac),iosel(i,outFrac),regsel(i),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
    %     end
    
    plot([0,1],[0,1],'k:');
    xlim([0,0.2]);
    ylim([0,0.2]);
    xlabel('in-density');
    ylabel('out-density');
    title(str_title);
    legend([gainh,damph,otherh],{'gain','damping','others'})
    
    [r,p]=corr(iosel(:,inFrac),iosel(:,outFrac));
    saveas(gcf,'iodensity.fig');
    save('corrstat.mat','r','p');
    
    %keyboard
    %     print(fh,fname,'-dpng','-r300');
    exportgraphics(fh,replace(fname,'.png','.pdf'),'ContentType','vector');
    
    end