function analyzeSessions_optoPair(sessionpath,task,stimPattern,options)

arguments
    sessionpath string
    task string
    stimPattern cell
    options.redo logical = true % Recalculate trial table and all preprocessing
    options.round logical = false % Round reward/airpuff/tone to get duration data
    options.performing logical = false; % Only plot traces where the animal performs
    
    options.plotPhotometry logical = true % Plot photometry summary plot
    options.plotLicks logical = true % Plot lick raster summary plot
end

%% Notes
% Shun_analyzeBehavior_optoPair
% Shun Li, 11/8/2022
% 02/14/2023: tidied up code, renamed to analyzeBehavior_optoPair
% 2023/07/28: packaged trial table into a function
% 2023/09/02: added camera plotting
% 2023/09/05: changed baselineIdx to selecting baseline licks 

%% Load data

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
             
% 1. Select session via uigetdir
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; 
if ispc; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
elseif isunix; session.projectPath = strcat('/',fullfile(dirsplit{2:end-1}));
end
clear dirsplit

disp(strcat('********** ',sessionName,'**********'));
load(strcat(sessionpath,filesep,'sync_',sessionName,'.mat'));
if ~isfield(session,'name'); session.name = sessionName; end
disp(['Session ',sessionName,' loaded']);

% Store opto stim params to sessions
% Opto stim params
params.stim.pulseFreq = str2double(stimPattern{1}); 
params.stim.pulseDuration = str2double(stimPattern{2}); 
params.stim.stimDuration = str2double(stimPattern{3});
params.stim.nPulsesPerStim = (params.stim.pulseDuration/1000) * params.stim.pulseFreq;
params.stim.pulseInterval = (1/params.stim.pulseFreq) - params.stim.pulseDuration;

%% Preprocess outcome and opto data

disp('Ongoing: preprocess outcome and opto data');

% Reward/punishment params
rewardUnit = 0.012; % 8ms opening to dispense 1ul water
rewardList = [0 2 5]; % in ul
% punishList = [0 0.1];
% toneList = [0 0.5 1]; % in sec
if options.round
    if ~exist('rightSolenoid_rounded','var') || options.redo
        % Round reward and tone
        rightSolenoid = rightSolenoid ./ rewardUnit;
        rightSolenoid_rounded = roundToTarget(rightSolenoid, rewardList); disp('Finished rounding: rightSolenoid');
        airpuff_rounded = roundToTarget(airpuff,punishList); disp('Finished rounding: airpuff');
        
        % Save rounded data
        save(strcat(sessionpath,filesep,'sync_',session.name),'rightSolenoid_rounded','airpuff_rounded','-append');
        disp('Finished: rounding cue/outcome data');
    end
else
    rightSolenoid_rounded = rightSolenoid;
    airpuff_rounded = airpuff;
end


% Find start of opto
if ~exist('firstPulse','var') || options.redo
    if ~isempty(find(redLaser, 1))
        % Find the first pulse of each stim pattern if nPulsePerPattern>1
        if params.stim.nPulsesPerStim > 1
            allPulses = find(redLaser);
            intervalThreshold = 10000;
            temp_interval = [100000,diff(allPulses)];
            firstPulse = allPulses(temp_interval > intervalThreshold);
    
            % Save first pulse data
            save(strcat(sessionpath,filesep,'sync_',session.name),"firstPulse",'-append');
        else
            firstPulse = find(redLaser);
            save(strcat(sessionpath,filesep,'sync_',session.name),"firstPulse",'-append');
        end
        disp('Finished: first pulse data');
    else 
        firstPulse = [];
    end
    save(strcat(sessionpath,filesep,'sync_',session.name),"firstPulse",'-append');
end


% Find start of lick bout
rightLickON = find(rightLick);
if ~exist('lickBout','var') || options.redo
    % Get lick bout start time (ILI < 0.5s)
    lickBout = getLickBout(rightLickON);
    save(strcat(sessionpath,filesep,'sync_',session.name),"lickBout",'-append');
end


% Combine stim&tone to form trial start
waterIdx = find(rightSolenoid_rounded);  
airpuffIdx = find(airpuff_rounded);

if ~exist('trials','var') || options.redo
    if strcmp(task,'random')
        [allTrials,~] = getTrials(find(leftTone),firstPulse,...
                             waterIdx,airpuffIdx);
    elseif contains(task,'punish')
        [allTrials,~] = getTrials(find(leftTone),firstPulse,waterIdx);
    else
        [allTrials,~] = getTrials(find(leftTone),firstPulse);
    end
    save(strcat(sessionpath,filesep,'sync_',session.name),"allTrials",'-append');
end


disp('Finished: preprocess outcome and opto data');

%% Generate trial table

if (~exist('trials','var') || options.redo)

    disp('Ongoing: making trial table');
    % Find digital events
    events{1} = allTrials;
    events{2} = airpuffIdx;
    events{3} = waterIdx;
    events{4} = rightLickON;
    events{5} = find(leftTone);
    events{6} = firstPulse;

    trials = getTrialTable(task,events,rightSolenoid_rounded,airpuff_rounded);

    % Calculate performance cutoff
    if contains(task,'reward'); [trials,cutoff_sample] = getSessionCutoff(trials,"->reward");
    elseif contains(task,'punish'); [trials,cutoff_sample] = getSessionCutoff(trials,"->punish");
    else; [trials,cutoff_sample] = getSessionCutoff(trials,"random"); 
    end
    params.analysis.cutoff_sample = cutoff_sample;
    disp(['Session cutoff calculated: ',num2str(cutoff_sample)]);

    % Save to sync.mat
    save(strcat(sessionpath,filesep,'sync_',session.name),'trials','-append');
    disp('Finished: trial table saved');
end


%% Task specific params

if strcmp(task,'random')
    
    % Select event idx
    toneIdx = find(leftTone);
    stimIdx = firstPulse;

    % Select baseline idx (align to each baseline lick)
    baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{:,'BaselineLicks'}),'BaselineLicks'});
    if isempty(baselineLicks); baselineIdx = [];
    else; baselineIdx = baselineLicks(:,1); end
    
    % baselineIdx = trials{:,"ENL"} - 5*params.sync.behaviorFs;
    
    % Create task legend
    events = {waterIdx,airpuffIdx,toneIdx,stimIdx,baselineIdx};
    labels = {'Water','Airpuff','Tone','Stim','Baseline'};
    taskLegend = getLegend(events,labels);

else
    % Select baseline idx (align to each baseline lick)
    baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{:,'BaselineLicks'}),'BaselineLicks'});
    if isempty(baselineLicks); baselineIdx = round((trials{2:end,"CueTime"} - trials{2:end,"ENL"}) - 5*params.sync.behaviorFs);
    else; baselineIdx = baselineLicks(:,1); end

    if contains(task,'reward')
        if options.performing
            stimIdx = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isReward == 1 & trials.performing == 1, "CueTime"};
            stimOmissionIdx = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isReward == 0 & trials.performing == 1, "CueTime"};
            pairIdx = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isReward == 1 & trials.performing == 1, "CueTime"};
            pairOmissionIdx = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isReward == 0 & trials.performing == 1, "CueTime"};
            toneIdx = trials{trials.isTone == 1 & trials.isStim == 0 ...
                & trials.performing == 1,"CueTime"};
        else
            stimIdx = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isReward == 1,"CueTime"};
            stimOmissionIdx = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isReward == 0,"CueTime"};
            pairIdx = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isReward == 1,"CueTime"};
            pairOmissionIdx = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isReward == 0,"CueTime"};
            toneIdx = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isReward == 1,"CueTime"};
            toneOmissionIdx = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isReward == 0,"CueTime"};
            % if isempty(toneIdx) && 
        end
    elseif contains(task,'punish')
        if options.performing
            stimIdx = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isPunishment == 1 & trials.performing == 1, "CueTime"};
            stimOmissionIdx = trials{trials.isTone == 0 & trials.isStim == 1 ...
                & trials.isPunishment == 0 & trials.performing == 1, "CueTime"};
            pairIdx = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isPunishment == 1 & trials.performing == 1, "CueTime"};
            pairOmissionIdx = trials{trials.isTone == 1 & trials.isStim == 1 ...
                & trials.isPunishment == 0 & trials.performing == 1, "CueTime"};
            toneIdx = trials{trials.isTone == 1 & trials.isStim == 0 ...
                & trials.performing == 1,"CueTime"};
        else
            stimIdx = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isPunishment == 1,"CueTime"};
            stimOmissionIdx = trials{trials.isTone == 0 & trials.isStim == 1 & trials.isPunishment == 0,"CueTime"};
            pairIdx = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isPunishment == 1,"CueTime"};
            pairOmissionIdx = trials{trials.isTone == 1 & trials.isStim == 1 & trials.isPunishment == 0,"CueTime"};
            toneIdx = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isPunishment == 1,"CueTime"};
            toneOmissionIdx = trials{trials.isTone == 1 & trials.isStim == 0 & trials.isPunishment == 0,"CueTime"};
        end
    end

    events = {waterIdx,toneIdx,stimIdx,pairIdx,airpuffIdx,baselineIdx};
    labels = {'Water','Tone only','Stim only','Pair','Airpuff','Baseline'};
    taskLegend = getLegend(events,labels);
end

%% Plot photometry summary plots

if options.plotPhotometry
    
    %% Set photometry signal to use (esp for LJ)
    % photometryNI = photometryNI;
    if session.withPhotometry; photometryLJ = rollingGreen; end

    %% Plot pre-processing steps

    if session.withPhotometryNI && session.withPhotometry
        initializeFig(.5,.5);
        histogram(normrnd(0,1,size(photometryLJ)),200); hold on
        histogram(photometryLJ,200); hold on
        histogram(photometryNI,200); hold on
        skew_lj = skewness(photometryLJ); kur_lj = kurtosis(photometryLJ);
        skew_ni = skewness(photometryNI); kur_ni = kurtosis(photometryNI);
        xlabel('z-score'); ylabel('Count'); legend({'Normal distribution','LJ Photometry','NI Photometry'});
        dim = [0.8205 0.58 0.55 0.27];
        str = {strcat("NI Skewness: ",num2str(skew_ni)),strcat("NI Kurtosis: ",num2str(kur_ni)),...
            strcat("LJ Skewness: ",num2str(skew_lj)),strcat("LJ Kurtosis: ",num2str(kur_lj))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        title('Histogram of photometry traces');
    
    elseif session.withPhotometryNI && ~session.withPhotometry
        initializeFig(.5,.5);
        histogram(normrnd(0,1,size(photometryNI)),200); hold on
        histogram(photometryNI,200); hold on
        xlabel('z-score'); ylabel('Count'); legend({'Normal distribution','Photometry'});
        dim = [0.8205 0.001 0.25 0.27];
        str = {strcat("Skewness: ",num2str(skewness(photometryNI))),strcat("Kurtosis: ",num2str(kurtosis(photometryNI)))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        title('Histogram of z-scored photometry');
    
    else
        initializeFig(.5,.5);
        histogram(normrnd(0,1,size(rollingGreen)),200); hold on
        histogram(rollingGreen,200); hold on
        skew = skewness(rollingGreen); kur = kurtosis(rollingGreen);
        xlabel('z-score'); ylabel('Count'); legend({'Normal distribution','Photometry'});
        dim = [0.8205 0.6 0.55 0.27];
        str = {strcat("Skewness: ",num2str(skew)),strcat("Kurtosis: ",num2str(kur))};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        title(['Histogram of ',getVarName(rollingGreen)]);
    end
    
    % Save figure
    saveas(gcf,strcat(sessionpath,filesep,'Summary_signal_',session.name,'.png'));


    %% (LJ) Plot combined photometry PSTHs
    
    if session.withPhotometry
        timeRange = [-1,5]; lick_binSize = 0.2;
        
        % 2. Plot traces
        % Determine figure layout 
        if params.session.withCamera && params.session.withEyeTracking
            initializeFig(0.5, 1);
            tiledlayout(4,1);
        else
            initializeFig(0.5,0.5);
            tiledlayout(2,1);
        end
    
        if strcmp(task,'random')
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterIdx,timeRange,photometryLJ,bluePurpleRed(1,:),params);
            [~,~] = plotTraces(airpuffIdx,timeRange,photometryLJ,[0.2, 0.2, 0.2],params);
            [~,~] = plotTraces(toneIdx,timeRange,photometryLJ,bluePurpleRed(350,:),params);
            [~,~] = plotTraces(stimIdx,timeRange,photometryLJ,bluePurpleRed(end,:),params);
            [~,~] = plotTraces(baselineIdx,timeRange,photometryLJ,[.75 .75 .75],params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('z-score');
            legend(taskLegend,'Location','northeast');
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
            plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
            plotLicks(baselineIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('Licks/s'); 
            legend(taskLegend,'Location','best');

            if params.session.withCamera && params.session.withEyeTracking
                % 2.3 Plot eye area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,eyeArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,eyeArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,eyeArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,eyeArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (z-score)');
                legend(taskLegend,'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,pupilArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,pupilArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (z-score)');
                legend(taskLegend,'Location','northeast');
            end
    
        else
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterIdx,timeRange,photometryLJ,bluePurpleRed(1,:),params);
            [~,~] = plotTraces(toneIdx,timeRange,photometryLJ,bluePurpleRed(350,:),params);
            [~,~] = plotTraces(stimIdx,timeRange,photometryLJ,bluePurpleRed(end,:),params);
            [~,~] = plotTraces(pairIdx,timeRange,photometryLJ,bluePurpleRed(150,:),params);
            [~,~] = plotTraces(airpuffIdx,timeRange,photometryLJ,[0.2, 0.2, 0.2],params);
            [~,~] = plotTraces(baselineIdx,timeRange,photometryLJ,[.75 .75 .75],params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('z-score');
            legend(taskLegend,'Location','northeast');
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
            plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
            plotLicks(pairIdx,timeRange,lick_binSize,bluePurpleRed(150,:),[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(baselineIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('Licks/s'); 
            legend(taskLegend,'Location','best');

            if params.session.withCamera && params.session.withEyeTracking
                % 2.3 Plot eye area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,eyeArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,eyeArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,eyeArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(pairIdx,timeRange,eyeArea_detrend,bluePurpleRed(150,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,eyeArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (z-score)');
                legend(taskLegend,'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,pupilArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,pupilArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(pairIdx,timeRange,pupilArea_detrend,bluePurpleRed(150,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (z-score)');
                legend(taskLegend,'Location','northeast');
            end
        end 
        
        saveas(gcf,strcat(sessionpath,filesep,'Summary_events_LJ_',session.name,'.png'));
    end
    
    %% (LJ) Plot single stimulus PSTH
    
    if session.withPhotometry
    
        if strcmp(task,'random')
            eventIdxes = {stimIdx,waterIdx,toneIdx,airpuffIdx};
            labels = {'Stim','Water','Tone','Airpuff'};
            eventDurations = [0.5,0,0.5,0.02];
            groupSizes = [20,30,10,30];
            longTimeRange = [-5,10];
            shortTimeRange = [-1,5]; 
            
            for event = 1:length(eventIdxes)
                eventIdx = eventIdxes{event};
                if isempty(eventIdx); continue; end
                label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.67,1);
                tiledlayout(4,4);
                
                % Plot short timescale
                nexttile(3,[1 2]);
                [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryLJ,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(baselineIdx,shortTimeRange,photometryLJ,[.75 .75 .75],params);
                plotEvent(label,eventDuration); 
                xlabel('Time (s)'); ylabel('z-score'); 
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Baseline (n=',num2str(length(baselineIdx)),')']},...
                    'Location','northeast');
                
                nexttile(7,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);


                % Plot heatmap
                nexttile(1,[4 2]);
                imagesc(t,1:length(eventIdx),traces); 
                set(gca,'YDir','normal');
                colorbar; box off
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');

                % Plot long timescale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,photometryLJ,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(baselineIdx,longTimeRange,photometryLJ,[.75 .75 .75],params);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Baseline (n=',num2str(length(baselineIdx)),')']},...
                    'Location','northeast');
                
                nexttile(15,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);
                
                saveas(gcf,strcat(sessionpath,filesep,'Events_LJ_',label,'_',session.name,'.png'));
            end
        else
            eventIdxes = {stimIdx,pairIdx,toneIdx,waterIdx,airpuffIdx};
            omissionIdxes = {stimOmissionIdx, pairOmissionIdx,toneOmissionIdx,[],[]};
            labels = {'Stim','Pair','Tone','Water','Airpuff'};
            eventDurations = [0.5,0.5,0.5,0,0.02];
            groupSizes = [10,10,10,30,30];
            longTimeRange = [-5,10];
            shortTimeRange = [-1,5]; 
            
            for event = 1:length(eventIdxes)
                eventIdx = eventIdxes{event};
                omissionIdx = omissionIdxes{event};
                if isempty(eventIdx); continue; end
                label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.67,1);
                tiledlayout(4,4);

                % Plot short time scale
                nexttile(3,[1 2]);
                [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryLJ,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(baselineIdx,shortTimeRange,photometryLJ,[.75 .75 .75],params);
                [~,~] = plotTraces(omissionIdx,shortTimeRange,photometryLJ,[0.3, 0.3, 0.3],params);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score'); 
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                        ['Baseline (n=',num2str(length(baselineIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdxes)),')']},...
                        'Location','northeast');
                
                nexttile(7,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);

                % Plot heatmap
                nexttile(1,[4 2]);
                imagesc(t,1:length(eventIdx),traces);
                set(gca,'YDir','normal');
                colorbar; box off
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');

                % Plot long time scale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,photometryLJ,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(baselineIdx,longTimeRange,photometryLJ,[.75 .75 .75],params);
                [~,~] = plotTraces(omissionIdx,longTimeRange,photometryLJ,[0.2, 0.2, 0.2],params);
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                        ['Baseline (n=',num2str(length(baselineIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdxes)),')']},...
                        'Location','northeast');
                
                nexttile(15,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);
                
                saveas(gcf,strcat(sessionpath,filesep,'Events_LJ_',label,'_',session.name,'.png'));
            end
        end
    end
    
    %% (NI) Plot combined PSTH
    
    if session.withPhotometryNI
        timeRange = [-1,5]; lick_binSize = 0.2;
        
        % 2. Plot traces
        % Determine figure layout 
        if params.session.withCamera && params.session.withEyeTracking
            initializeFig(0.5, 1);
            tiledlayout(4,1);
        else
            initializeFig(0.5,0.5);
            tiledlayout(2,1);
        end
    
        if strcmp(task,'random')
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterIdx,timeRange,photometryNI,bluePurpleRed(1,:),params,signalSystem='ni');
            [~,~] = plotTraces(airpuffIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,signalSystem='ni');
            [~,~] = plotTraces(toneIdx,timeRange,photometryNI,bluePurpleRed(350,:),params,signalSystem='ni');
            [~,~] = plotTraces(stimIdx,timeRange,photometryNI,bluePurpleRed(end,:),params,signalSystem='ni');
            [~,~] = plotTraces(baselineIdx,timeRange,photometryNI,[.75 .75 .75],params,signalSystem='ni');
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('z-score'); 
            legend(taskLegend,'Location','best'); 
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
            plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
            plotLicks(baselineIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('Licks/s'); 
            legend(taskLegend,'Location','best');

            if params.session.withCamera && params.session.withEyeTracking
                % 2.3 Plot eye area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,eyeArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,eyeArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,eyeArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,eyeArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (z-score)');
                legend(taskLegend,'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,pupilArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,pupilArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (z-score)');
                legend(taskLegend,'Location','northeast');
            end

        else
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterIdx,timeRange,photometryNI,bluePurpleRed(1,:),params,signalSystem='ni');
            [~,~] = plotTraces(toneIdx,timeRange,photometryNI,bluePurpleRed(350,:),params,signalSystem='ni');
            [~,~] = plotTraces(stimIdx,timeRange,photometryNI,bluePurpleRed(end,:),params,signalSystem='ni');
            [~,~] = plotTraces(pairIdx,timeRange,photometryNI,bluePurpleRed(150,:),params,signalSystem='ni');
            [~,~] = plotTraces(airpuffIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,signalSystem='ni');
            [~,~] = plotTraces(baselineIdx,timeRange,photometryNI,[.75 .75 .75],params,signalSystem='ni');
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('z-score'); 
            legend(taskLegend,'Location','best'); 
            
            % 2.2 Plot lick traces
            nexttile
            plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
            plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
            plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
            plotLicks(pairIdx,timeRange,lick_binSize,bluePurpleRed(150,:),[],rightLick,params);
            plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
            plotLicks(baselineIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
            plotEvent('',0);
            xlabel('Time (s)'); ylabel('Licks/s'); 
            legend(taskLegend,'Location','best');

            if params.session.withCamera && params.session.withEyeTracking
                % 2.3 Plot eye area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,eyeArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,eyeArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,eyeArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(pairIdx,timeRange,eyeArea_detrend,bluePurpleRed(150,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,eyeArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,eyeArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Eye area (z-score)');
                legend(taskLegend,'Location','northeast');
    
                % 2.4 Plot pupil area traces
                nexttile
                [~,~] = plotTraces(waterIdx,timeRange,pupilArea_detrend,bluePurpleRed(1,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(toneIdx,timeRange,pupilArea_detrend,bluePurpleRed(350,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(stimIdx,timeRange,pupilArea_detrend,bluePurpleRed(end,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(pairIdx,timeRange,pupilArea_detrend,bluePurpleRed(150,:),params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(airpuffIdx,timeRange,pupilArea_detrend,[0.2, 0.2, 0.2],params,signalSystem='camera',smooth=15);
                [~,~] = plotTraces(baselineIdx,timeRange,pupilArea_detrend,[.75 .75 .75],params,signalSystem='camera',smooth=15);
                plotEvent('',0);
                xlabel('Time (s)'); ylabel('Pupil area (z-score)');
                legend(taskLegend,'Location','northeast');
            end
        end
     
        saveas(gcf,strcat(sessionpath,filesep,'Summary_events_NI_',session.name,'.png'));
    end
    
    %% (NI) Plot single stimulus PSTH
    
    if session.withPhotometryNI
    
        if strcmp(task,'random')
            eventIdxes = {stimIdx,waterIdx,toneIdx,airpuffIdx};
            labels = {'Stim','Water','Tone','Airpuff'};
            eventDurations = [0.5,0,0.5,0.02];
            groupSizes = [20,30,10,30];
            longTimeRange = [-5,10];
            shortTimeRange = [-1,5]; 
            
            for event = 1:length(eventIdxes)
                eventIdx = eventIdxes{event};
                if isempty(eventIdx); continue; end
                label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.67,1);
                tiledlayout(4,4);

                % Plot short time scale
                nexttile(3,[1 2]);
                [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryNI,bluePurpleRed(1,:),params,signalSystem='ni');
                [~,~] = plotTraces(baselineIdx,shortTimeRange,photometryNI,[.75 .75 .75],params,signalSystem='ni');
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Baseline (n=',num2str(length(baselineIdx)),')']},...
                    'Location','northeast'); 
                
                nexttile(7,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);

                % Plot heatmap
                nexttile(1,[4 2]);
                imagesc(t,1:length(eventIdx),traces);
                set(gca,'YDir','normal');
                colorbar; box off
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');


                % Plot long time scale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,photometryNI,bluePurpleRed(1,:),params,signalSystem='ni');
                [~,~] = plotTraces(baselineIdx,longTimeRange,photometryNI,[.75 .75 .75],params,signalSystem='ni');
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Baseline (n=',num2str(length(baselineIdx)),')']},...
                    'Location','northeast'); 
                
                nexttile(15,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);
                
                
                saveas(gcf,strcat(sessionpath,filesep,'Events_NI_',label,'_',session.name,'.png'));
            end
        else
            eventIdxes = {stimIdx,pairIdx,toneIdx,waterIdx,airpuffIdx};
            omissionIdxes = {stimOmissionIdx, pairOmissionIdx,toneOmissionIdx,[],[]};
            labels = {'Stim','Pair','Tone','Water','Airpuff'};
            eventDurations = [0.5,0.5,0.5,0,0.02];
            groupSizes = [10,10,10,30,30];
            longTimeRange = [-5,10];
            shortTimeRange = [-1,5]; 
            
            for event = 1:length(eventIdxes)
                eventIdx = eventIdxes{event};
                omissionIdx = omissionIdxes{event};
                if isempty(eventIdx); continue; end
                label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.67,1);
                tiledlayout(4,4);

                % Plot short time scale
                nexttile(3,[1 2]);
                [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryNI,bluePurpleRed(1,:),params,signalSystem='ni');
                [~,~] = plotTraces(baselineIdx,shortTimeRange,photometryNI,[.75 .75 .75],params,signalSystem='ni');
                [~,~] = plotTraces(omissionIdx,shortTimeRange,photometryNI,[0.2, 0.2, 0.2],params,signalSystem='ni');
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                        ['Baseline (n=',num2str(length(baselineIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdxes)),')']},...
                        'Location','northeast');
                
                nexttile(7,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);

                % Plot heatmap
                nexttile(1,[4 2]);
                imagesc(t,1:length(eventIdx),traces);
                set(gca,'YDir','normal');
                colorbar; box off
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('Trials');

                % Plot long time scale
                nexttile(11,[1 2]);
                [traces,t] = plotTraces(eventIdx,longTimeRange,photometryNI,bluePurpleRed(1,:),params,signalSystem='ni');
                [~,~] = plotTraces(baselineIdx,longTimeRange,photometryNI,[.75 .75 .75],params,signalSystem='ni');
                [~,~] = plotTraces(omissionIdx,longTimeRange,photometryNI,[0.2, 0.2, 0.2],params,signalSystem='ni');
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                        ['Baseline (n=',num2str(length(baselineIdx)),')'],...
                        [label,' omission (n=',num2str(length(omissionIdxes)),')']},...
                        'Location','northeast');
                
                nexttile(15,[1 2]);
                nLines = ceil(size(traces,1)/groupSize);
                legendList = cell(nLines,1);
                nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
                for i = 1:nLines
                    startTrial = (i-1)*groupSize+1; 
                    if i == nLines; endTrial = size(traces,1);
                    else; endTrial = i*groupSize; end
                    plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
                    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
                end
                plotEvent(label,eventDuration);
                xlabel('Time (s)'); ylabel('z-score');
                legend(legendList);
                
                saveas(gcf,strcat(sessionpath,filesep,'Events_NI_',label,'_',session.name,'.png'));
            end
        end    
    end


end

%% Plot behavior related plots

% Plot lick bout distribution
if options.plotLicks
    initializeFig(0.5,0.5); tiledlayout(2,2);

    ENLinSec = trials{:,"ENL"} / params.sync.behaviorFs;
    ITIextra = trials{:,'ITI'} - ENLinSec;

    % Plot lick bout count distribution
    nexttile;
    histogram(lickBout(:,2),30); 
    xlabel('Licks per lick bout'); ylabel('Occurance'); box off

    % Plot ITI-ENL distribution
    nexttile;
    histogram(trials{:,'ITI'},30); 
    xlabel('ITI (s)'); ylabel('Occurance'); box off

    % Plot lick per trial vs trials
    nexttile;
    plot(trials{:,'TrialNumber'},trials{:,'nLicks'},'Color',bluePurpleRed(1,:),LineWidth=2);
    xlabel('Trials'); ylabel('Licks per trial'); box off

    % Plot ITI-ENL vs trials
    nexttile;
    plot(trials{:,'TrialNumber'},ENLinSec,'Color',bluePurpleRed(500,:),LineWidth=2); hold on
    plot(trials{:,'TrialNumber'},ITIextra,'Color',bluePurpleRed(1,:),LineWidth=2);
    xlabel('Trials'); ylabel('Time (s)'); ylim([0,Inf]); box off
    legend({'ENL','ITI-ENL'});

    saveas(gcf,strcat(sessionpath,filesep,'Summary_ITI&Bout_',session.name,'.png'));

end

if options.plotLicks && contains(task,'pairing')
        
    %% Plot lick raster plot and eye traces
    timeRange = [-5,10]; cameraTimeRange = [-1,5];
    markerSize = 20; lick_binSize = 0.2;

    % Get event time and number by trial type
    stimIdx = trials{trials.isTone == 0 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
    toneIdx = trials{trials.isTone == 1 & trials.isStim == 0,["TrialNumber","CueTime","OutcomeTime","ENL"]};
    pairIdx = trials{trials.isTone == 1 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","ENL"]};
    stimIdx(:,3) = stimIdx(:,3)./params.sync.behaviorFs;
    toneIdx(:,3) = toneIdx(:,3)./params.sync.behaviorFs;
    pairIdx(:,3) = pairIdx(:,3)./params.sync.behaviorFs;

    % getLicks by trial type
    [stimLickRate,~,stimLicks] = getLicks(timeRange,stimIdx(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [pairLickRate,~,pairLicks] = getLicks(timeRange,pairIdx(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [toneLickRate,~,toneLicks] = getLicks(timeRange,toneIdx(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);


    % Plot overall raster plot (color coded by trial type)
    if params.session.withCamera && params.session.withEyeTracking
        initializeFig(1,1); nCols = 5;
        % get camera trace by trial type
        [stimEyeArea,t_cam] = plotTraces(stimIdx(:,2),cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [pairEyeArea,~] = plotTraces(pairIdx(:,2),cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [toneEyeArea,~] = plotTraces(toneIdx(:,2),cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [stimPupilArea,~] = plotTraces(stimIdx(:,2),cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [pairPupilArea,~] = plotTraces(pairIdx(:,2),cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [tonePupilArea,~] = plotTraces(toneIdx(:,2),cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
    else
        initializeFig(0.5,0.67); nCols = 3;   
    end
    tiledlayout(3,nCols); 

    % Plot lick raster plot for session
    nexttile([3,2]);
    for i = 1:size(stimLicks,1)
        scatter(stimLicks{i},stimIdx(i,1),markerSize,'filled','MarkerFaceColor',bluePurpleRed(end,:)); hold on
        scatter(stimIdx(i,3),stimIdx(i,1),markerSize+10,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    for i = 1:size(toneLicks,1)
        scatter(toneLicks{i},toneIdx(i,1),markerSize,'filled','MarkerFaceColor',bluePurpleRed(350,:)); hold on
        scatter(toneIdx(i,3),toneIdx(i,1),markerSize+10,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    for i = 1:size(pairLicks,1)
        scatter(pairLicks{i},pairIdx(i,1),markerSize,'filled','MarkerFaceColor',bluePurpleRed(150,:)); hold on
        scatter(pairIdx(i,3),pairIdx(i,1),markerSize+10,bluePurpleRed(1,:),'pentagram','filled'); hold on
    end
    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
    ylabel('Trial'); ylim([0,size(trials,1)]);
    plotEvent("",0.5);

    % Plot lick traces across session
    traces = {stimLickRate,pairLickRate,toneLickRate};
    labels = {'Stim','Pair','Tone'};
    groupSizes = [20, 20, 10];
    for event = 1:length(traces)
        trace = traces{event};
        label = labels{event};
        groupSize = groupSizes(event);

        nexttile(3 + nCols*(event-1));
        t = linspace(timeRange(1),timeRange(2),size(trace,2));
        nLines = ceil(size(trace,1)/groupSize); legendList = cell(nLines,1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
        for i = 1:nLines
            startTrial = (i-1)*groupSize+1; 
            if i == nLines; endTrial = size(trace,1);
            else; endTrial = i*groupSize; end
            plotSEM(t,trace(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
            legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
        end
        plotEvent(label,0.5);
        xlabel('Time (s)'); ylabel('Licks/s');
        legend(legendList);
    end


    if params.session.withCamera && params.session.withEyeTracking
        % Plot eye area across session
        traces = {stimEyeArea,pairEyeArea,toneEyeArea};
        labels = {'Stim','Pair','Tone'};
        groupSizes = [20, 20, 10];
        for event = 1:length(traces)
            trace = traces{event};
            label = labels{event};
            groupSize = groupSizes(event);
    
            nexttile(4 + nCols*(event-1));
            nLines = ceil(size(trace,1)/groupSize); legendList = cell(nLines,1);
            nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
            for i = 1:nLines
                startTrial = (i-1)*groupSize+1; 
                if i == nLines; endTrial = size(trace,1);
                else; endTrial = i*groupSize; end
                plotSEM(t_cam,trace(startTrial:endTrial,:),bluePurpleRed(nColors(i),:),smooth=15);
                legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
            end
            plotEvent(label,0.5);
            xlabel('Time (s)'); ylabel('Eye area (z-score)');
            legend(legendList);
        end
    
    
        % Plot pupil area across session
        traces = {stimPupilArea,pairPupilArea,tonePupilArea};
        labels = {'Stim','Pair','Tone'};
        groupSizes = [20, 20, 10];
        for event = 1:length(traces)
            trace = traces{event};
            label = labels{event};
            groupSize = groupSizes(event);
    
            nexttile(5 + nCols*(event-1));
            nLines = ceil(size(trace,1)/groupSize); legendList = cell(nLines,1);
            nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
            for i = 1:nLines
                startTrial = (i-1)*groupSize+1; 
                if i == nLines; endTrial = size(trace,1);
                else; endTrial = i*groupSize; end
                plotSEM(t_cam,trace(startTrial:endTrial,:),bluePurpleRed(nColors(i),:),smooth=15);
                legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
            end
            plotEvent(label,0.5);
            xlabel('Time (s)'); ylabel('Pupil area (z-score)');
            legend(legendList);
        end
    end

    saveas(gcf,strcat(sessionpath,filesep,'Summary_behavior_',session.name,'.png'));

    %% Plot ENL analysis (whether animal licks to reward specifically or randomlly licking)
    initializeFig(1,1);
    tiledlayout(5,3);

    % 1. Bin trials based on ENL length, plot lick trace aligning to ENL start
    % should not smear if animal lick specifically to the cue
    timeRangeENL = [0,10];
    [stimLickRateENL,~,~] = getLicks(timeRangeENL,stimIdx(:,2)-stimIdx(:,4),lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI);
    [pairLickRateENL,~,~] = getLicks(timeRangeENL,pairIdx(:,2)-pairIdx(:,4),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [toneLickRateENL,~,~] = getLicks(timeRangeENL,toneIdx(:,2)-toneIdx(:,4),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);

    nENLBins = [2,3,4];
    for i = 1:length(nENLBins)

        % 1.1 Sort trials based on ENL duration
        % For stim only
        nexttile;
        trialIdx = stimIdx; traces = stimLickRateENL;
        [~,sortIdx] = sortrows(trialIdx,4); 
        nModTrials = mod(height(trialIdx),nENLBins(i));
        ENLsortIdx = mat2cell(reshape(sortIdx(1:end-nModTrials),[],nENLBins(i))',ones(1,nENLBins(i)));
        ENLsortIdx{end} = [ENLsortIdx{end}, sortIdx(end-nModTrials+1:end)'];

        t = linspace(timeRangeENL(1),timeRangeENL(2),size(traces,2));
        legendList = cell(nENLBins(i),1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nENLBins(i)));
        % 1.2 Plot lick trace
        for j = 1:nENLBins(i)
            plotSEM(t,traces(ENLsortIdx{j},:),bluePurpleRed(nColors(j),:));
            legendList{j} = ['ENL bin ', num2str(j),' (n=',num2str(length(ENLsortIdx{j})),')'];
        end
        xlabel('Time from ENL start (s)'); ylabel('Licks/s');
        legend(legendList);

        % For pair
        nexttile;
        trialIdx = pairIdx; traces = pairLickRateENL;
        [~,sortIdx] = sortrows(trialIdx,4); 
        nModTrials = mod(height(trialIdx),nENLBins(i));
        ENLsortIdx = mat2cell(reshape(sortIdx(1:end-nModTrials),[],nENLBins(i))',ones(1,nENLBins(i)));
        ENLsortIdx{end} = [ENLsortIdx{end}, sortIdx(end-nModTrials+1:end)'];

        t = linspace(timeRangeENL(1),timeRangeENL(2),size(traces,2));
        legendList = cell(nENLBins(i),1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nENLBins(i)));
        % 1.2 Plot lick trace
        for j = 1:nENLBins(i)
            plotSEM(t,traces(ENLsortIdx{j},:),bluePurpleRed(nColors(j),:));
            legendList{j} = ['ENL bin ', num2str(j),' (n=',num2str(length(ENLsortIdx{j})),')'];
        end
        xlabel('Time from ENL start (s)'); ylabel('Licks/s');
        legend(legendList);

        % For tone only
        nexttile;
        trialIdx = toneIdx; traces = toneLickRateENL;
        [~,sortIdx] = sortrows(trialIdx,4); 
        nModTrials = mod(height(trialIdx),nENLBins(i));
        ENLsortIdx = mat2cell(reshape(sortIdx(1:end-nModTrials),[],nENLBins(i))',ones(1,nENLBins(i)));
        ENLsortIdx{end} = [ENLsortIdx{end}, sortIdx(end-nModTrials+1:end)'];

        t = linspace(timeRangeENL(1),timeRangeENL(2),size(traces,2));
        legendList = cell(nENLBins(i),1);
        nColors = round(linspace(1,size(bluePurpleRed,1),nENLBins(i)));
        % 1.2 Plot lick trace
        for j = 1:nENLBins(i)
            plotSEM(t,traces(ENLsortIdx{j},:),bluePurpleRed(nColors(j),:));
            legendList{j} = ['ENL bin ', num2str(j),' (n=',num2str(length(ENLsortIdx{j})),')'];
        end
        xlabel('Time from ENL start (s)'); ylabel('Licks/s');
        legend(legendList);
        
    end


    % 2. Find out baseline period (eg 5s-10s after reward), bootstrap and
    % compare with stim onset
    % (i.e.) how likely I get 2+ licks within a 2s window randomly selected
    % from baseline (a window of 10 secs)?
    timeRangeBaseline = [-10,0]; 
    nboot = 100000; nBins = 50;
    [~,allLicksBaseline,~] = getLicks(timeRangeBaseline,trials{:,'CueTime'}-trials{:,"ENL"},...
        lick_binSize,[],rightLick,params.sync.behaviorFs,params.sync.timeNI);
    % Bootstrap baseline licking success rate
    hitPercentBaseline = zeros(nboot,1);
    for i = 1:nboot
        % Randomly select columns for 2 sec (2/0.2 = 10)
        windowInBin = 2 / lick_binSize;
        baselineSamples = datasample(allLicksBaseline,windowInBin,2,'Replace',false);
        hitPercentBaseline(i) = sum(sum(baselineSamples,2) >= 2)/size(baselineSamples,1);
    end


    % For stim only
    % 2.1 Calculate cue window success rate
    hitPercentCue = length(find(trials.isTone == 0 & trials.isStim == 1 & trials.Outcome == "H"))/size(stimIdx,1);
    % 2.2 Plot distribution
    nexttile([2 1]);
    h = histogram(hitPercentBaseline,nBins);
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:);
    hold on
    xline(hitPercentCue,'-r',{'Hit rate of','stim only trials'},...
        'LineWidth',3,...
        'LabelOrientation','horizontal');
    xlim([0,1]); xlabel('Hit rate'); ylabel('Occurance'); box off
   

    % For pair
    % 2.1 Calculate cue window success rate
    hitPercentCue = length(find(trials.isTone == 1 & trials.isStim == 1 & trials.Outcome == "H"))/size(pairIdx,1);
    % 2.2 Plot distribution
    nexttile([2 1]);
    h = histogram(hitPercentBaseline,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:);
    hold on
    xline(hitPercentCue,'-r',{'Hit rate of','pair trials'},...
        'LineWidth',3,...
        'LabelOrientation','horizontal');
    xlim([0,1]); xlabel('Hit rate'); ylabel ('Occurance'); box off


    % For tone only
    % 2.1 Calculate cue window success rate
    hitPercentCue = length(find(trials.isTone == 1 & trials.isStim == 0 & trials.Outcome == "H"))/size(toneIdx,1);
    % 2.2 Plot distribution
    nexttile([2 1]);
    h = histogram(hitPercentBaseline,nBins);
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:);
    hold on
    xline(hitPercentCue,'-r',{'Hit rate of','tone only trials'},...
        'LineWidth',3,...
        'LabelOrientation','horizontal');
    xlim([0,1]); xlabel('Hit rate'); ylabel('Occurance'); box off

    saveas(gcf,strcat(sessionpath,filesep,'Summary_ENL_',session.name,'.png'));

end

% save(strcat(session.path,filesep,'sync_',sessionName),'params','-append');
return

% %% test
% 
% waterIdx = find(rightSolenoid); 
% timeRange = [-1, 5];
% 
% initializeFig(0.67,0.4); tiledlayout(1,4);
% nexttile;
% [~,~] = plotTraces(waterIdx,timeRange,detrendGreen,bluePurpleRed(1,:),params);
% nexttile;
% [~,~] = plotTraces(waterIdx,timeRange,demodGreen,bluePurpleRed(150,:),params);
% nexttile;
% [~,~] = plotTraces(waterIdx,timeRange,rollingGreen,bluePurpleRed(350,:),params);
% nexttile;
% [~,~] = plotTraces(waterIdx,timeRange,rollingGreenLP,bluePurpleRed(500,:),params);

return
end