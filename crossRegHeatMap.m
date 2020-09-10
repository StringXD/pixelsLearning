% Cross-region connectivity density and difference between learning and welltrain
prefix='welltrain';
% load('full_duo_XCORR_stats_delay_6_1_2_2msbin.mat');
% gen_join_set = true;
% if gen_join_set
%     join_reg_list=cell(0);
%     for sidx=1:length(stats)
%         s=stats{sidx};
%         join_reg_list{end+1}=s.reg_su1;
%         join_reg_list{end+1}=s.reg_su2;
%     end
%     join_reg_list(cellfun(@isempty,join_reg_list))=[];
%     join_reg_set=unique(join_reg_list);
% end
% reg_set=join_reg_set(~strcmp(join_reg_set,'Unlabeled'));
% reg_set=reg_set(~strcmp(reg_set,'root'));
% save('reg_keep.mat','reg_set');
load('reg_keep.mat');

for bin = 1:6
    load(sprintf('full_duo_XCORR_stats_delay_6_%d_%d_2msbin',bin,bin+1));
    gen_pair_mat=true;
    if gen_pair_mat
        pair_mat=zeros(length(reg_set),length(reg_set));
        pair_sel_mat=zeros(length(reg_set),length(reg_set));
        for sidx=1:length(stats)
            s=stats{sidx};
            
            if s.s1_trials<20 || s.s2_trials<20 
                continue
            end

            su1reg_idx=find(strcmp(s.reg_su1,reg_set));
            su2reg_idx=find(strcmp(s.reg_su2,reg_set));
            if isempty(su1reg_idx) || isempty(su2reg_idx)
                fprintf('%s, %s\n',s.reg_su1,s.reg_su2);
                %keyboard
                continue
            end
            
            pair_mat(su1reg_idx,su2reg_idx)=pair_mat(su1reg_idx,su2reg_idx)+1;
            pair_mat(su2reg_idx,su1reg_idx)=pair_mat(su2reg_idx,su1reg_idx)+1;

            if bin>0 && s.prefered_sample_su1(bin+1) && s.prefered_sample_su2(bin+1) && s.prefered_sample_su1(bin+1)==s.prefered_sample_su2(bin+1)  
                pair_sel_mat(su1reg_idx,su2reg_idx)=pair_sel_mat(su1reg_idx,su2reg_idx)+1;
                pair_sel_mat(su2reg_idx,su1reg_idx)=pair_sel_mat(su2reg_idx,su1reg_idx)+1;
            end
        end
        save(sprintf('%s_pair_mat_duo_6s_%d_%d.mat',prefix,bin_range(1),bin_range(2)),'pair_mat','pair_sel_mat');
    end

    gen_conn_mat=true;
    if gen_conn_mat
        %conn_mat_all=cell(0);
        
%         load(sprintf('XCORR_stats_delay_6_%d_%d_2msbin.mat',bin,bin+1));
        conn_mat_S1=zeros(length(reg_set),length(reg_set));
        conn_sel_mat_S1=zeros(length(reg_set),length(reg_set));
        conn_mat_S2=zeros(length(reg_set),length(reg_set));
        conn_sel_mat_S2=zeros(length(reg_set),length(reg_set));
        for pidx=1:length(stats)
            s=stats{pidx};
            if s.totalcount<250 || s.s1_trials<20 || s.s2_trials<20 
                continue
            end
            
            su1reg_idx=find(strcmp(s.reg_su1,reg_set));
            su2reg_idx=find(strcmp(s.reg_su2,reg_set));
            if isempty(su1reg_idx) || isempty(su2reg_idx)
                fprintf('%s, %s\n', s.reg_su1, s.reg_su2);
                continue
            end
            
            
            if bin>0 && s.prefered_sample_su1(bin+1) && s.prefered_sample_su2(bin+1) && s.prefered_sample_su1(bin+1)==s.prefered_sample_su2(bin+1)  
                sel_flag=true;
            else
                sel_flag=false;
            end

            if s.s1_peak_significant
                if s.AIs1>0.4 % su1 to su2
                    conn_mat_S1(su1reg_idx,su2reg_idx)=conn_mat_S1(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat_S1(su1reg_idx,su2reg_idx)=conn_sel_mat_S1(su1reg_idx,su2reg_idx)+1;
                    end                    
                elseif s.AIs1<-0.4  %=0 is possible
                    conn_mat_S1(su2reg_idx,su1reg_idx)=conn_mat_S1(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat_S1(su2reg_idx,su1reg_idx)=conn_sel_mat_S1(su2reg_idx,su1reg_idx)+1;
                    end
                end
            end    
            if s.s2_peak_significant
                if s.AIs2>0.4 % su1 to su2
                    conn_mat_S2(su1reg_idx,su2reg_idx)=conn_mat_S2(su1reg_idx,su2reg_idx)+1;
                    if sel_flag
                        conn_sel_mat_S2(su1reg_idx,su2reg_idx)=conn_sel_mat_S2(su1reg_idx,su2reg_idx)+1;
                    end
                elseif s.AIs2<-0.4  %=0 is possible
                    conn_mat_S2(su2reg_idx,su1reg_idx)=conn_mat_S2(su2reg_idx,su1reg_idx)+1;
                    if sel_flag
                        conn_sel_mat_S2(su2reg_idx,su1reg_idx)=conn_sel_mat_S2(su2reg_idx,su1reg_idx)+1;
                    end
                end
            end
                
        end
    disp(sprintf('%s_conn_mat_duo_6s_%d_%d.mat',prefix,bin_range(1),bin_range(2)));
    save(sprintf('%s_conn_mat_duo_6s_%d_%d.mat',prefix,bin_range(1),bin_range(2)),'conn_mat_S1','conn_sel_mat_S1','conn_mat_S2','conn_sel_mat_S2')    
    % return
    end    
end


