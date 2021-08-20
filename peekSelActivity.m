% peek sample, test, pair selective neuron activity
load('I:\FC data\FR_modulated_lnwt.mat');
allpath = deblank(h5read('I:\FC data\transient_6_complete.hdf5','/path'));
wtpaths = deblank(h5read('I:\FC data\transient_6.hdf5','/path'));
sel_test_or_reward = load('E:\prJ\neuropixels\learning\sel_after_6.mat').selectivityAfterDelay;
pairSel = load('E:\prJ\neuropixels\learning\pairSel.mat').pairSel;
% sample selective
wrs_p = h5read('I:\FC data\transient_6_complete.hdf5','/wrs_p');
selectivity = h5read('I:\FC data\transient_6_complete.hdf5','/selectivity');
[mem_type,per_bin] = get_mem_type(wrs_p,selectivity);
% S1T1 S1T2 S2T1 S2T2 raw spk 
for learningphase = 1:2 % 1 for learning, 2 for welltrain
    switch learningphase
        case 1
            phase = 'learning';
            phasesel = ~contains(allpath,wtpaths);
        case 2
            phase = 'welltrain';
            phasesel = contains(allpath,wtpaths);
    end
    for celltype = 1:12
        switch celltype
            case 1
                suffix = 'S1 sus';
                susel = find(mem_type'==1&phasesel);
            case 2
                suffix = 'S2 sus';
                susel = find(mem_type'==3&phasesel);
            case 3
                suffix = 'S1 trans';
                susel = find(mem_type'==2&phasesel);
            case 4
                suffix = 'S2 trans';
                susel = find(mem_type'==4&phasesel);
            case 5
                suffix = 'T1';
                susel = find(sel_test_or_reward(:,1)==1&phasesel);
            case 6
                suffix = 'T2';
                susel = find(sel_test_or_reward(:,1)==2&phasesel);
            case 7
                suffix = 'S1T1';
                susel = find(pairSel(:,1)==1&phasesel);
            case 8 
                suffix = 'S1T2';
                susel = find(pairSel(:,2)==1&phasesel);
            case 9
                suffix = 'S2T1';
                susel = find(pairSel(:,3)==1&phasesel);
            case 10
                suffix = 'S2T2';
                susel = find(pairSel(:,4)==1&phasesel);
            case 11
                suffix = 'Go';
                susel = find(pairSel(:,5)==1&phasesel);
            case 12
                suffix = 'Nogo';
                susel = find(pairSel(:,6)==1&phasesel);
        end
        savepath = fullfile('E:\prJ\neuropixels\learning','raw',suffix,phase);
        for suidx = susel'
            for trialtype = 1:4
                subplot(4,1,trialtype);
                switch trialtype
                    case 1
                        trialsel = FR_modulated{suidx,2}(:,5)==4 & FR_modulated{suidx,2}(:,6)==8 & FR_modulated{suidx,2}(:,8)==6;
                    case 2
                        trialsel = FR_modulated{suidx,2}(:,5)==4 & FR_modulated{suidx,2}(:,6)==4 & FR_modulated{suidx,2}(:,8)==6;
                    case 3
                        trialsel = FR_modulated{suidx,2}(:,5)==8 & FR_modulated{suidx,2}(:,6)==8 & FR_modulated{suidx,2}(:,8)==6;
                    case 4
                        trialsel = FR_modulated{suidx,2}(:,5)==8 & FR_modulated{suidx,2}(:,6)==4 & FR_modulated{suidx,2}(:,8)==6;
                end
                % plot
                hold on;
                imagesc(FR_modulated{suidx,1}(trialsel,:));
                caxis([0 3]);
                line([60.5,60.5],[0,sum(trialsel)],'Color','k','LineStyle','--');
                line([80.5,80.5],[0,sum(trialsel)],'Color','k','LineStyle','--');
                line([200.5,200.5],[0,sum(trialsel)],'Color','k','LineStyle','--');
                line([220.5,220.5],[0,sum(trialsel)],'Color','k','LineStyle','--');
            end
            saveas(gcf,fullfile(savepath,sprintf('%d.png',suidx)));
            close all;
        end              
    end
end






% functions
function [mem_type,per_bin]=get_mem_type(wrs_p,selec)
% 0=NM,1=S1 sust, 2=S1 trans, 3=S2 sust, 4=S2 trans,-1=switched
arguments
    wrs_p (14,:) double % assuming 1 sec bin
    selec (14,:) double % assuming 1 sec bin
end
mem_type=nan(1,size(selec,2));
per_bin=nan(6,size(selec,2));
for i=1:size(selec,2)
    sel_bin=find(wrs_p(5:10,i)<0.05);
    if isempty(sel_bin)
        mem_type(i)=0;
        per_bin(:,i)=0;
    else
        ssign=sign(selec(sel_bin+4,i));
        if numel(unique(ssign))==1 % non-switched 
            if numel(ssign)==6  % sustained
                if ssign(1)==1
                    mem_type(i)=1;
                    per_bin(:,i)=1;
                else
%                     if ssign(1)~=-1,keyboard;end % one time check
                    mem_type(i)=3;
                    per_bin(:,i)=2;
                end
            else % transient
                per_bin(:,i)=0;
                if ssign(1)==1
                    mem_type(i)=2;
                    per_bin(sel_bin,i)=1;
                else
%                     if ssign(1)~=-1,keyboard;end % one time check
                    mem_type(i)=4;
                    per_bin(sel_bin,i)=2;
                end
            end
        else %switched
            mem_type(i)=-1;
            per_bin(:,i)=-1;
        end
    end
end
mem_type=int32(mem_type);
per_bin=int32(per_bin);
end