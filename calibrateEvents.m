% interactively calibrate events
% Sptrialtype 58, SpITI 59, SpSess 61,SpTrain 62, SpLaser 65, Splaserswitch 79: 1 start 0 end.
% Sample odor 10; Test odor 9
addpath('E:\matlab scripts\lickAnalysis');
dataSumFolder = 'E:\prJ\neuropixels\learning\dataSumLN';
serFolder = 'E:\prJ\neuropixels\learning\sync\Learning';
flEvt=dir([dataSumFolder '\*\*\events.hdf5']);
for oneEvt=flEvt'
    rootpath=oneEvt.folder;
    trials = h5read(fullfile(rootpath,'events.hdf5'),'/trials');
    events = h5read(fullfile(rootpath,'events.hdf5'),'/events');
    % extract mice id and date
    Fm = regexp(rootpath,'M\d*','match');
    miceid = Fm{1};
    midSer = (miceid - 'M75')*[100;10;1];
    Fd = regexp(rootpath,'2020\d*','match');
    date = Fd{1};
    serFl = dir(fullfile(serFolder,['M' num2str(midSer) '*' date '*.ser']));
    if ~isempty(serFl)
        serData = ser2mat(fullfile(serFl(1).folder,serFl(1).name));
    else
        disp('Can not find serial file!')
        disp(['Check path: ',rootpath])
        continue;
    end
    % serData as ground
    refTime = serData(:,1);
    refStatus = serData(:,3);
    refBefAllSess_serial = find(refStatus == 51);
    refAftAllSess_serial = find(refStatus == 62);
    if isempty(refAftAllSess_serial)
        refAftAllSess_serial = find(refStatus == 61,1,'last');
    end
    refSession = refStatus(refBefAllSess_serial:refAftAllSess_serial);
    refSess_time = refTime(refBefAllSess_serial:refAftAllSess_serial);
    refTrial_begin_serial = find(refSession == 58);
    refTrial_end_serial = find(refSession == 59);
    trialNum = length(refTrial_end_serial);
    serTable = zeros(4,trialNum); % Sample, Test, Lick, Delay
    serTime = zeros(2,trialNum);
    uncertainTrials = [];
    for n = 1:trialNum
        current_trial = refSession(refTrial_begin_serial(n):refTrial_end_serial(n));
        current_trial_period = refSess_time(refTrial_begin_serial(n):refTrial_end_serial(n));
        current_trial_period = double(current_trial_period);
        odorTag = current_trial == 9 | current_trial == 10;
        sertimes = current_trial_period(odorTag);
        serTime(1,n) = sertimes(1);
        serTime(2,n) = sertimes(3);
        current_trial_period = (current_trial_period - current_trial_period(1))/1000;
        try
            serTable(1:2,n) = 4+[4 0 0 0;0 0 4 0]*(current_trial(odorTag) == 10);
        catch
            temp = current_trial(odorTag);
            serTable(1,n) = 4 + 4*(temp(1)==10);
            serTable(2,n) = 4 + 4*(temp(3)==10);
        end
        serTable(3,n) = 2*any(current_trial == 4 | current_trial == 7) - 1;
        delay = round(max(diff(current_trial_period(odorTag))));
        if delay == 6 || delay == 3
            serTable(4,n) = delay;
        else
            timelists = round(diff(current_trial_period(odorTag)));
            if any(timelists == 6) && ~any(timelists == 3)
                serTable(4,n) = 6;
            elseif any(timelists == 3) 
                serTable(4,n) = 3;
            else
                uncertainTrials = [uncertainTrials;n];
                disp('Uncertain Delay')
                disp(['Check path: ',rootpath])
                %pause;
            end
        end
    end
    % sequence matching
    trials = double(trials);
    flag = 0;
    errList = [];
    if length(serTable) == length(trials) && sum(abs(trials(7:8,:) - serTable(3:4,:)),'all') == 0
        trials(5:6,:) = serTable(1:2,:);
        trials(end+1,:) = ones(1,length(trials));
    elseif length(serTable) < length(trials)
        offset = 0;
        while length(serTable) + offset <= length(trials)
            err = sum(trials(7:8,1+offset:length(serTable)+offset) - serTable(3:4,:)~=0,'all');
            errList = [errList;offset,err];
            if err == 0
                trials(5:6,1+offset:length(serTable)+offset) = serTable(1:2,:);
                trials = trials(:,1+offset:length(serTable)+offset);
                flag = 1;
                break;
            end
            offset = offset + 1;
        end
        if flag ~= 1
            % remove uncertain trials
            % do not use for trial-history based analysis
            if min(errList(:,2)) < 10
                [~,I] = min(min(errList(:,2)));
                offset = errList(I,1);
                [~,tid] = find(trials(7:8,1+offset:length(serTable)+offset) - serTable(3:4,:));
                trials(:,tid) = [];
            elseif contains(rootpath,'M81_6_learning_20200909_g0')
                continue;
            elseif contains(rootpath,'M79_4_learning_20200903_g1')
                continue;
            else
                disp('Uncertain Events')
                pause;
            end
        end
        trials(end+1,:) = ones(1,length(trials)); 
    elseif length(serTable) > length(trials)
        disp('Fewer trials than ser!')
        disp(['Check path: ',rootpath])
        offset = 0;
        hastrim = 0;
        while length(trials) + offset <= length(serTable) && flag == 0
            err = sum(trials(7:8,:) - serTable(3:4,1+offset:length(trials)+offset)~=0,'all');
            errList = [errList;offset,err];
            if err == 0
                flag = 1;
                % make sure there are spikes in missing trials or discard
                % them
                spikeTS = readNPY(fullfile(rootpath,'spike_times.npy'));
                % fill the blank data
                % linear regression: t(pixel_t2) - t(pixel_t1) = 30[t(ser_t2) - t(ser_t1)]
                tl = length(trials);
                trials = [zeros(size(trials,1),offset),trials,zeros(size(trials,1),length(serTable)-tl-offset)];
                regressTimeS1 = trials(1,1) + 30*(serTime(1,:) - serTime(1,offset+1));
                regressTimeS2 = trials(2,1) + 30*(serTime(2,:) - serTime(2,offset+1));
                trials(end+1,:) = logical(trials(1,:));
                trials(1,:) = regressTimeS1.* ~trials(1,:) + trials(1,:);
                trials(2,:) = regressTimeS2.* ~trials(2,:) + trials(2,:);
                trials(3:4,:) = round(trials(1:2,:)/30000);
                trials(5:8,:) = serTable;
                % To avoid spikeGLX recording missing
                missProneMask = zeros(1,length(trials));
                for tidx = 1:length(trials)
                    missProneMask(tidx) = length(spikeTS(spikeTS>trials(1,tidx) & spikeTS <= trials(2,tidx))) < 500;
                end
                trials = trials(:,~missProneMask);
                break;
            end
            offset = offset + 1;
            if flag == 0 && length(trials) + offset == length(serTable) && hastrim == 0
                trials = trials(:,2:end-1);
                offset = 0;
                hastrim = 1;
            end
        end
        if flag == 0 && length(trials) + offset == length(serTable) && hastrim == 1
            % remove uncertain trials
            % do not use for trial-history based analysis
            if min(errList(:,2)) < 10
                [~,I] = min(min(errList(:,2)));
                offset = errList(I,1);
                [~,tid] = find(trials(7:8,:) - serTable(3:4,1+offset:length(trials)+offset));
                trials(:,tid) = [];
            else
                pause;
            end
        end
    end
    h5create(fullfile(rootpath,'eventsRescue.hdf5'),'/trials',size(trials),'Datatype','int32')   
    h5write(fullfile(rootpath,'eventsRescue.hdf5'),'/trials',int32(trials))
    h5create(fullfile(rootpath,'eventsRescue.hdf5'),'/events',size(events),'Datatype','int32')   
    h5write(fullfile(rootpath,'eventsRescue.hdf5'),'/events',int32(events))
end
disp('Done!')
    
%% peek events
fl=dir('E:\prJ\neuropixels\learning\dataSumLN\*\*\events.hdf5');
eventList = cell(0,3);
for onefile=fl'
    rootpath=onefile.folder;
    trials = h5read(fullfile(rootpath,'events.hdf5'),'/trials');
    eventList{end+1,1} = rootpath;
    eventList{end,2} = trials;
    try
        trials1 = h5read(fullfile(rootpath,'eventsRescue.hdf5'),'/trials');
        eventList{end,3} = trials1;
    catch
    end
end