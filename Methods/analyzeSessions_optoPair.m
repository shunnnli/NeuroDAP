function analyzeSessions_optoPair(sessionpath,task,stimPattern,options)

arguments
    sessionpath string
    task string
    stimPattern cell

    options.pavlovian logical = true
    options.reactionTime double = 1.5

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
% 2023/10/23: changed how to plot photometry signal, assume everything
% recorded in labjack

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
params.stim.nPulsesPerStim = (params.stim.stimDuration/1000) * params.stim.pulseFreq;
params.stim.pulseInterval = (1/params.stim.pulseFreq) - params.stim.pulseDuration;

% Invert ShutterRed and ShutterBlue if neccessary
% params.stim.invertedSignal = str2double(stimPattern{4});
% if params.stim.invertedSignal
%     redLaser = ~redLaser;
%     blueLaser = ~blueLaser;
% end

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


% Find start of opto cue
if ~exist('optoCue','var') || options.redo
    if ~isempty(find(redLaser, 1))
        % Find the first pulse of each stim pattern if nPulsePerPattern>1
        if params.stim.nPulsesPerStim > 1
            allPulses = find(redLaser);
            intervalThreshold = 10000;
            temp_interval = [100000,diff(allPulses)];
            optoCue = allPulses(temp_interval > intervalThreshold);
    
            % Save first pulse data
            save(strcat(sessionpath,filesep,'sync_',session.name),"optoCue",'-append');
        else
            optoCue = find(redLaser);
            save(strcat(sessionpath,filesep,'sync_',session.name),"optoCue",'-append');
        end
        disp('Finished: saved optoCue data');
    else 
        optoCue = [];
    end
    save(strcat(sessionpath,filesep,'sync_',session.name),"optoCue",'-append');
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
        [allTrials,~] = getTrials(find(leftTone),optoCue,...
                             waterIdx,airpuffIdx);
    elseif contains(task,'punish')
        [allTrials,~] = getTrials(find(leftTone),optoCue,waterIdx);
    else
        [allTrials,~] = getTrials(find(leftTone),optoCue);
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
    events{6} = optoCue;

    trials = getTrialTable(task,events,rightSolenoid_rounded,airpuff_rounded,...
                pavlovian=options.pavlovian,reactionTime=options.reactionTime);

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
    stimIdx = optoCue;

    % Select baseline idx (align to each baseline lick)
    baselineLicks = cell2mat(trials{~cellfun(@isempty,trials{:,'BaselineLicks'}),'BaselineLicks'});
    if isempty(baselineLicks); baselineIdx = [];
    else; baselineIdx = baselineLicks(:,1); end
    
    % baselineIdx = trials{:,"ENL"} - 5*params.sync.behaviorFs;
    
    % Create task legend
    events = {waterIdx,airpuffIdx,toneIdx,stimIdx,baselineIdx};
    labels = {'Water','Airpuff','Tone','Stim','Baseline'};
    taskLegend = getLegend(events,labels);
    
    for i = 1:length(events)
        disp(['Total ',labels{i},': ',num2str(length(events{i}))]);
    end

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
            if isempty(toneIdx) && ~isempty(toneOmissionIdx)
                toneIdx = toneOmissionIdx;
            end
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
            if isempty(toneIdx) && ~isempty(toneOmissionIdx)
                toneIdx = toneOmissionIdx;
            end
        end
    end

    events = {waterIdx,toneIdx,stimIdx,pairIdx,airpuffIdx,baselineIdx};
    labels = {'Water','Tone only','Stim only','Pair','Airpuff','Baseline'};
    taskLegend = getLegend(events,labels);

    for i = 1:length(events)
        disp(['Total ',labels{i},': ',num2str(length(events{i}))]);
    end

end

%% Plot photometry summary plots

if options.plotPhotometry
    %% Plot pre-processing steps
    nSignals = size(processed,2);

    initializeFig(.67,.67); tiledlayout('flow');
    
    for i = 1:nSignals
        nexttile;
        histogram(normrnd(0,1,size(processed(i).signal)),200); hold on
        histogram(processed(i).signal,200); hold on
        box off

        skew_lj = skewness(processed(i).signal); 
        kur_lj = kurtosis(processed(i).signal);
        xlabel('z-score'); ylabel('Count'); legend({'Normal distribution',processed(i).name});
        
        title(processed(i).name);
        subtitle(strcat("Skewness: ",num2str(skew_lj),", Kurtosis: ",num2str(kur_lj)));
    end
    
    % Save figure
    saveas(gcf,strcat(sessionpath,filesep,'Summary_photometry_distribution.png'));

    %% Loop through processed struct
    for photometry = 1:nSignals

        % Load signal of interest
        signal = processed(photometry).signal;
        finalFs = processed(photometry).finalFs;
        system = processed(photometry).system;
        
        %% Plot combined PSTH
        timeRange = [-1,5]; lick_binSize = 0.2;
        % 2. Plot traces
        % Determine figure layout 
        if params.session.withCamera && params.session.withEyeTracking
            initializeFig(0.5, 1); tiledlayout(4,1);
        else
            initializeFig(0.5,0.5); tiledlayout(2,1);
        end

        if strcmp(task,'random')
            % 2.1 Plot photometry traces
            nexttile
            [~,~] = plotTraces(waterIdx,timeRange,signal,bluePurpleRed(1,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(airpuffIdx,timeRange,signal,[0.2, 0.2, 0.2],params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(toneIdx,timeRange,signal,bluePurpleRed(350,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(stimIdx,timeRange,signal,bluePurpleRed(end,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(baselineIdx,timeRange,signal,[.75 .75 .75],params,...
                        signalFs=finalFs,signalSystem=system);
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
            [~,~] = plotTraces(waterIdx,timeRange,signal,bluePurpleRed(1,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(toneIdx,timeRange,signal,bluePurpleRed(350,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(stimIdx,timeRange,signal,bluePurpleRed(end,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(pairIdx,timeRange,signal,bluePurpleRed(150,:),params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(airpuffIdx,timeRange,signal,[0.2, 0.2, 0.2],params,...
                        signalFs=finalFs,signalSystem=system);
            [~,~] = plotTraces(baselineIdx,timeRange,signal,[.75 .75 .75],params,...
                        signalFs=finalFs,signalSystem=system);
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
        saveas(gcf,strcat(sessionpath,filesep,'Summary_events_',processed(photometry).name,'.png'));

        %% Plot single stimulus PSTH
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
                label = labels{event}; eventDuration = eventDurations(event); 
                groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.67,1); tiledlayout(4,4);
                
                % Plot short timescale
                nexttile(3,[1 2]);
                [traces,t] = plotTraces(eventIdx,shortTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system);
                [~,~] = plotTraces(baselineIdx,shortTimeRange,signal,[.75 .75 .75],params,...
                                signalFs=finalFs,signalSystem=system);
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
                [traces,t] = plotTraces(eventIdx,longTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system);
                [~,~] = plotTraces(baselineIdx,longTimeRange,signal,[.75 .75 .75],params,...
                                signalFs=finalFs,signalSystem=system);
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
                
                saveas(gcf,strcat(sessionpath,filesep,'Events_',processed(photometry).name,'_',label,'.png'));
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
                [traces,t] = plotTraces(eventIdx,shortTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system);
                [~,~] = plotTraces(baselineIdx,shortTimeRange,signal,[.75 .75 .75],params,...
                                signalFs=finalFs,signalSystem=system);
                [~,~] = plotTraces(omissionIdx,shortTimeRange,signal,[0.3, 0.3, 0.3],params,...
                                signalFs=finalFs,signalSystem=system);
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
                [traces,t] = plotTraces(eventIdx,longTimeRange,signal,bluePurpleRed(1,:),params,...
                                signalFs=finalFs,signalSystem=system);
                [~,~] = plotTraces(baselineIdx,longTimeRange,signal,[.75 .75 .75],params,...
                                signalFs=finalFs,signalSystem=system);
                [~,~] = plotTraces(omissionIdx,longTimeRange,signal,[0.2, 0.2, 0.2],params,...
                                signalFs=finalFs,signalSystem=system);
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

                saveas(gcf,strcat(sessionpath,filesep,'Events_',processed(photometry).name,'_',label,'.png'));
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
    xlabel('Licks per lick bout'); ylabel('Count'); box off

    % Plot ITI-ENL distribution
    nexttile;
    histogram(trials{:,'ITI'},30); 
    xlabel('ITI (s)'); ylabel('Count'); box off

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

    saveas(gcf,strcat(sessionpath,filesep,'Behavior_ITI&LickBout.png'));
end

if options.plotLicks && contains(task,'pairing')
        
    %% Plot session overview for licking
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

    % 1. Plot lick raster and trace
    initializeFig(0.67,0.67); tiledlayout(3,2);
    % 1.1 Plot lick raster plot for session
    nexttile([3,1]);
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
    % 1.2 Plot lick traces across session
    traces = {stimLickRate,pairLickRate,toneLickRate};
    labels = {'Stim','Pair','Tone'};
    groupSizes = [20, 20, 10];
    for event = 1:length(traces)
        trace = traces{event};
        label = labels{event};
        groupSize = groupSizes(event);

        nexttile;
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
    saveas(gcf,strcat(sessionpath,filesep,'Behavior_LickOverview.png'));

    
    %% Plot session overview for eye
    if params.session.withCamera && params.session.withEyeTracking
        initializeFig(0.67,0.67); tiledlayout(3,2); 
        % get camera trace by trial type
        [stimEyeArea,t_cam] = plotTraces(stimIdx(:,2),cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [pairEyeArea,~] = plotTraces(pairIdx(:,2),cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [toneEyeArea,~] = plotTraces(toneIdx(:,2),cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [stimPupilArea,~] = plotTraces(stimIdx(:,2),cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [pairPupilArea,~] = plotTraces(pairIdx(:,2),cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
        [tonePupilArea,~] = plotTraces(toneIdx(:,2),cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');

        % Plot eye area across session
        traces = {stimEyeArea,pairEyeArea,toneEyeArea};
        labels = {'Stim','Pair','Tone'};
        groupSizes = [20, 20, 10];
        for event = 1:length(traces)
            trace = traces{event};
            label = labels{event};
            groupSize = groupSizes(event);
    
            nexttile(1+ 2*(event-1));
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
    
            nexttile(2+ 2*(event-1));
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
        saveas(gcf,strcat(sessionpath,filesep,'Behavior_EyeOverview.png'));
    end

    %% Plot ENL aligned lick trace
    % Bin trials based on ENL length, plot lick trace aligning to ENL start
    % should not smear if animal lick specifically to the cue

    initializeFig(0.67,0.67); tiledlayout(3,3);
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
    saveas(gcf,strcat(sessionpath,filesep,'Behavior_ENLAlignedLick.png'));

    %% Plot distributions
    initializeFig(1,1);
    tiledlayout(4,4);

    % 1. Distribution of baseline licks
    % Find out baseline period (eg 5s-10s after reward), bootstrap and
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
        windowInBin = round(options.reactionTime / lick_binSize);
        baselineSamples = datasample(allLicksBaseline,windowInBin,2,'Replace',true);
        hitPercentBaseline(i) = sum(sum(baselineSamples,2) >= 2)/size(baselineSamples,1);
    end

    % For pair
    % Calculate cue window success rate
    hitPercentCue = length(find(trials.isTone == 1 & trials.isStim == 1 & trials.nAnticipatoryLicks >= 3))/size(pairIdx,1);
    % Plot distribution
    nexttile; title('Distribution of baseline licks for pair trials (bootstrapped)');
    h = histogram(hitPercentBaseline,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:);
    hold on
    xline(hitPercentCue,'-r',{'Hit rate of','pair trials'},...
        'LineWidth',3,...
        'LabelOrientation','horizontal');
    xlim([0,1]); xlabel('Hit rate'); ylabel ('Count'); box off
    title("Baseline licks (pair trials)");

    % For stim only
    % Calculate cue window success rate
    hitPercentCue = length(find(trials.isTone == 0 & trials.isStim == 1 & trials.nAnticipatoryLicks >= 3))/size(stimIdx,1);
    % Plot distribution
    nexttile; title('Distribution of baseline licks for stim only trials (bootstrapped)');
    h = histogram(hitPercentBaseline,nBins);
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:);
    hold on
    xline(hitPercentCue,'-r',{'Hit rate of','stim only trials'},...
        'LineWidth',3,...
        'LabelOrientation','horizontal');
    xlim([0,1]); xlabel('Hit rate'); ylabel('Count'); box off
    title("Baseline licks (stim only trials)");

    % For tone only
    % Calculate cue window success rate
    hitPercentCue = length(find(trials.isTone == 1 & trials.isStim == 0 & trials.nAnticipatoryLicks >= 3))/size(toneIdx,1);
    % Plot distribution
    nexttile; title('Distribution of baseline licks for tone only trials (bootstrapped)');
    h = histogram(hitPercentBaseline,nBins);
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:);
    hold on
    xline(hitPercentCue,'-r',{'Hit rate of','tone only trials'},...
        'LineWidth',3,...
        'LabelOrientation','horizontal');
    xlim([0,1]); xlabel('Hit rate'); ylabel('Count'); box off
    title("Baseline licks (tone only trials)");

    % Plot trend
    nexttile;
    plot(trials{:,"TrialNumber"},trials{:,"nBaselineLicks"},Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairIdx(:,1),trials{pairIdx(:,1),"nBaselineLicks"},100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimIdx(:,1),trials{stimIdx(:,1),"nBaselineLicks"},100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneIdx(:,1),trials{toneIdx(:,1),"nBaselineLicks"},100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Baseline licks"); box off
    title("Baseline licks (all trials)");


    stageCutoff = linspace(1,size(trials,1),4);
    stageCutoff = stageCutoff(2:3);
    
    % 2. Distribution of first lick reaction time
    rt_pair = trials{pairIdx(:,1),["TrialNumber","ReactionTime"]};
    rt_stim = trials{stimIdx(:,1),["TrialNumber","ReactionTime"]};
    rt_tone = trials{toneIdx(:,1),["TrialNumber","ReactionTime"]};
    % For early stage of session
    % Calculate bootstrap distribution
    rt_pair_early = rt_pair(rt_pair(:,1)<stageCutoff(1),2);
    rt_stim_early = rt_pair(rt_stim(:,1)<stageCutoff(1),2);
    rt_tone_early = rt_pair(rt_tone(:,1)<stageCutoff(1),2);
    [~,bootsam] = bootstrp(nboot,[],rt_pair_early);
    rt_pair_early_bs = rt_pair_early(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],rt_stim_early);
    rt_stim_early_bs = rt_stim_early(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],rt_tone_early);
    rt_tone_early_bs = rt_tone_early(bootsam) / params.sync.behaviorFs;
    % Plot distribution
    nexttile; title("Cue reaction time distribution during early stage of session (bootstrapped)");
    h = histogram(rt_pair_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(rt_stim_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(rt_tone_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (early stage)");

    % For middle stage of session
    % Calculate bootstrap distribution
    rt_pair_mid = rt_pair(rt_pair(:,1)>=stageCutoff(1) & rt_pair(:,1)<stageCutoff(2),2);
    rt_stim_mid = rt_pair(rt_stim(:,1)>=stageCutoff(1) & rt_stim(:,1)<stageCutoff(2),2);
    rt_tone_mid = rt_pair(rt_tone(:,1)>=stageCutoff(1) & rt_tone(:,1)<stageCutoff(2),2);
    [~,bootsam] = bootstrp(nboot,[],rt_pair_mid);
    rt_pair_mid_bs = rt_pair_mid(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],rt_stim_mid);
    rt_stim_mid_bs = rt_stim_mid(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],rt_tone_mid);
    rt_tone_mid_bs = rt_tone_mid(bootsam) / params.sync.behaviorFs;
    % Plot distribution
    nexttile; title("Cue reaction time distribution during middle stage of session (bootstrapped)");
    h = histogram(rt_pair_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(rt_stim_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(rt_tone_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (middle stage)");

    % For late stage of session
    % Calculate bootstrap distribution
    rt_pair_late = rt_pair(rt_pair(:,1)<=stageCutoff(2),2);
    rt_stim_late = rt_pair(rt_stim(:,1)<=stageCutoff(2),2);
    rt_tone_late = rt_pair(rt_tone(:,1)<=stageCutoff(2),2);
    [~,bootsam] = bootstrp(nboot,[],rt_pair_late);
    rt_pair_late_bs = rt_pair_late(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],rt_stim_late);
    rt_stim_late_bs = rt_stim_late(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],rt_tone_late);
    rt_tone_late_bs = rt_tone_late(bootsam) / params.sync.behaviorFs;
    % Plot distribution
    nexttile; title("Cue reaction time distribution during late stage of session (bootstrapped)");
    h = histogram(rt_pair_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(rt_stim_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(rt_tone_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('First lick reaction time (s)'); ylabel('Count'); box off
    title("First lick reaction time (late stage)");

    % Plot trend
    nexttile;
    plot(trials{:,"TrialNumber"},trials{:,"ReactionTime"}/params.sync.behaviorFs,Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairIdx(:,1),trials{pairIdx(:,1),"ReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimIdx(:,1),trials{stimIdx(:,1),"ReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneIdx(:,1),trials{toneIdx(:,1),"ReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Reaction time (s)"); box off
    title("First lick reaction time (all trials)");


    % 3. Distribution of anticipatory licks
    al_pair = trials{pairIdx(:,1),["TrialNumber","nAnticipatoryLicks"]};
    al_stim = trials{stimIdx(:,1),["TrialNumber","nAnticipatoryLicks"]};
    al_tone = trials{toneIdx(:,1),["TrialNumber","nAnticipatoryLicks"]};
    % For early stage of session
    % Calculate bootstrap distribution
    al_pair_early = al_pair(al_pair(:,1)<stageCutoff(1),2);
    al_stim_early = al_pair(al_stim(:,1)<stageCutoff(1),2);
    al_tone_early = al_pair(al_tone(:,1)<stageCutoff(1),2);
    [~,bootsam] = bootstrp(nboot,[],al_pair_early);
    al_pair_early_bs = al_pair_early(bootsam);
    [~,bootsam] = bootstrp(nboot,[],al_stim_early);
    al_stim_early_bs = al_stim_early(bootsam);
    [~,bootsam] = bootstrp(nboot,[],al_tone_early);
    al_tone_early_bs = al_tone_early(bootsam);
    % Plot distribution
    nexttile; title("Cue reaction time distribution during early stage of session (bootstrapped)");
    h = histogram(al_pair_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(al_stim_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(al_tone_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('Anticipatory licks'); ylabel('Count'); box off
    title("Anticipatory licks (early stage)");

    % For middle stage of session
    % Calculate bootstrap distribution
    al_pair_mid = al_pair(al_pair(:,1)>=stageCutoff(1) & al_pair(:,1)<stageCutoff(2),2);
    al_stim_mid = al_pair(al_stim(:,1)>=stageCutoff(1) & al_stim(:,1)<stageCutoff(2),2);
    al_tone_mid = al_pair(al_tone(:,1)>=stageCutoff(1) & al_tone(:,1)<stageCutoff(2),2);
    [~,bootsam] = bootstrp(nboot,[],al_pair_mid);
    al_pair_mid_bs = al_pair_mid(bootsam);
    [~,bootsam] = bootstrp(nboot,[],al_stim_mid);
    al_stim_mid_bs = al_stim_mid(bootsam);
    [~,bootsam] = bootstrp(nboot,[],al_tone_mid);
    al_tone_mid_bs = al_tone_mid(bootsam);
    % Plot distribution
    nexttile; title("Cue reaction time distribution during middle stage of session (bootstrapped)");
    h = histogram(al_pair_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(al_stim_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(al_tone_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('Anticipatory licks'); ylabel('Count'); box off
    title("Anticipatory licks (middle stage)");

    % For late stage of session
    % Calculate bootstrap distribution
    al_pair_late = al_pair(al_pair(:,1)<=stageCutoff(2),2);
    al_stim_late = al_pair(al_stim(:,1)<=stageCutoff(2),2);
    al_tone_late = al_pair(al_tone(:,1)<=stageCutoff(2),2);
    [~,bootsam] = bootstrp(nboot,[],al_pair_late);
    al_pair_late_bs = al_pair_late(bootsam);
    [~,bootsam] = bootstrp(nboot,[],al_stim_late);
    al_stim_late_bs = al_stim_late(bootsam);
    [~,bootsam] = bootstrp(nboot,[],al_tone_late);
    al_tone_late_bs = al_tone_late(bootsam);
    % Plot distribution
    nexttile; title("Cue reaction time distribution during late stage of session (bootstrapped)");
    h = histogram(al_pair_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(al_stim_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(al_tone_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('Anticipatory licks'); ylabel('Count'); box off
    title("Anticipatory licks (late stage)");

    % Plot trend
    nexttile;
    plot(trials{:,"TrialNumber"},trials{:,"nAnticipatoryLicks"},Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairIdx(:,1),trials{pairIdx(:,1),"nAnticipatoryLicks"},100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimIdx(:,1),trials{stimIdx(:,1),"nAnticipatoryLicks"},100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneIdx(:,1),trials{toneIdx(:,1),"nAnticipatoryLicks"},100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Anticipatory licks"); box off
    title("Anticipatory licks (all trials)");


    % 4. Distribution of reward rection time
    ort_pair = trials{pairIdx(:,1),["TrialNumber","OutcomeReactionTime"]};
    ort_stim = trials{stimIdx(:,1),["TrialNumber","OutcomeReactionTime"]};
    ort_tone = trials{toneIdx(:,1),["TrialNumber","OutcomeReactionTime"]};
    % For early stage of session
    % Calculate bootstrap distribution
    ort_pair_early = ort_pair(ort_pair(:,1)<stageCutoff(1),2);
    ort_stim_early = ort_pair(ort_stim(:,1)<stageCutoff(1),2);
    ort_tone_early = ort_pair(ort_tone(:,1)<stageCutoff(1),2);
    [~,bootsam] = bootstrp(nboot,[],ort_pair_early);
    ort_pair_early_bs = ort_pair_early(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],ort_stim_early);
    ort_stim_early_bs = ort_stim_early(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],ort_tone_early);
    ort_tone_early_bs = ort_tone_early(bootsam) / params.sync.behaviorFs;
    % Plot distribution
    nexttile; title("Cue reaction time distribution during early stage of session (bootstrapped)");
    h = histogram(ort_pair_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(ort_stim_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(ort_tone_early_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (early stage)");

    % For middle stage of session
    % Calculate bootstrap distribution
    ort_pair_mid = ort_pair(ort_pair(:,1)>=stageCutoff(1) & ort_pair(:,1)<stageCutoff(2),2);
    ort_stim_mid = ort_pair(ort_stim(:,1)>=stageCutoff(1) & ort_stim(:,1)<stageCutoff(2),2);
    ort_tone_mid = ort_pair(ort_tone(:,1)>=stageCutoff(1) & ort_tone(:,1)<stageCutoff(2),2);
    [~,bootsam] = bootstrp(nboot,[],ort_pair_mid);
    ort_pair_mid_bs = ort_pair_mid(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],ort_stim_mid);
    ort_stim_mid_bs = ort_stim_mid(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],ort_tone_mid);
    ort_tone_mid_bs = ort_tone_mid(bootsam) / params.sync.behaviorFs;
    % Plot distribution
    nexttile; title("Cue reaction time distribution during middle stage of session (bootstrapped)");
    h = histogram(ort_pair_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(ort_stim_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(ort_tone_mid_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (middle stage)");

    % For late stage of session
    % Calculate bootstrap distribution
    ort_pair_late = ort_pair(ort_pair(:,1)<=stageCutoff(2),2);
    ort_stim_late = ort_pair(ort_stim(:,1)<=stageCutoff(2),2);
    ort_tone_late = ort_pair(ort_tone(:,1)<=stageCutoff(2),2);
    [~,bootsam] = bootstrp(nboot,[],ort_pair_late);
    ort_pair_late_bs = ort_pair_late(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],ort_stim_late);
    ort_stim_late_bs = ort_stim_late(bootsam) / params.sync.behaviorFs;
    [~,bootsam] = bootstrp(nboot,[],ort_tone_late);
    ort_tone_late_bs = ort_tone_late(bootsam) / params.sync.behaviorFs;
    % Plot distribution
    nexttile; 
    h = histogram(ort_pair_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(150,:); h.EdgeColor = bluePurpleRed(150,:); hold on
    h = histogram(ort_stim_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(end,:); h.EdgeColor = bluePurpleRed(end,:); hold on
    h = histogram(ort_tone_late_bs,nBins); 
    h.FaceColor = bluePurpleRed(350,:); h.EdgeColor = bluePurpleRed(350,:); hold on
    xlabel('Outcome reaction time (s)'); ylabel('Count'); box off
    title("Outcome reaction time (late stage)");

    % Plot trend
    nexttile;
    plot(trials{:,"TrialNumber"},trials{:,"OutcomeReactionTime"}/params.sync.behaviorFs,Color=[0.75,0.75,0.75],LineWidth=2); hold on
    scatter(pairIdx(:,1),trials{pairIdx(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(150,:),'filled'); hold on
    scatter(stimIdx(:,1),trials{stimIdx(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(end,:),'filled'); hold on
    scatter(toneIdx(:,1),trials{toneIdx(:,1),"OutcomeReactionTime"}/params.sync.behaviorFs,100,bluePurpleRed(350,:),'filled'); hold on
    xlabel("Trials"); ylabel("Outcome reaction time (s)"); box off
    title("Outcome reaction time (all trials)");

    
    saveas(gcf,strcat(sessionpath,filesep,'Summary_Distributions.png'));

end

% save(strcat(session.path,filesep,'sync_',sessionName),'params','-append');
return

end