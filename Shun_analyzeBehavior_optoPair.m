% Opto pairing behavior analysis

clear; close all;
% addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed,blueGreenPurple] = loadColors;
             
% Manually enter session
%session.name = "20221206-SL046-D2_g0";
% session.name = "20221201-M38-473_g0";
%sessionpath = strcat('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Project valence\Recordings\',session.name);
% sessionpath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Test photometry\';

% Select session via uigetdir
sessionpath = uigetdir('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
clear dirsplit
load(strcat(sessionpath,'\','sync_',sessionName,'.mat'));
disp(['Session ',sessionName,' loaded']);

timeRange = [-0.5,3]; minLicks = 1; reactionTimeSamp = 2 * nidq.Fs;
blockLength = 20; toneDuration = 0.5;

lick_binSize = 0.1; blink_thresh = 2.5; % in turns of z score

%% Photometry signal summary

fig_signal_summary = initializeFig(1,1);

subplot(3,2,1)
plot(demodGreen);xlabel('Time'); ylabel('z-score');
title('detrend->demod');

subplot(3,2,3)
plot(rollingGreen);xlabel('Time'); ylabel('z-score');
title('detrend->demod->rolling');

subplot(3,2,5)
plot(rollingGreenLP);xlabel('Time'); ylabel('z-score');
title('detrend->demod->LP->rolling');

% Create the uitable
subplot(3,2,[2 4 6])
histogram(normrnd(0,1,size(rollingGreenLP)),200); hold on
histogram(rollingGreenLP,200); hold on
skew = skewness(rollingGreenLP); kur = kurtosis(rollingGreenLP);
xlabel('z-score'); ylabel('Count'); legend({'Normal distribution','Photometry'});
dim = [0.8205 0.6 0.55 0.27];
str = {strcat("Skewness: ",num2str(skew)),strcat("Kurtosis: ",num2str(kur))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title(['Histogram of ',getVarName(rollingGreenLP)]);

saveas(fig_signal_summary,strcat(sessionpath,'\signal_summary_',session.name,'.png'));

%% Generate trial table

if ~exist('trials','var')
    % Extract data
    allTrialsTime = find(redLaser);
    rightSolenoidON = find(rightSolenoid);
    rightLickON = find(rightLick);
    airpuffON = find(airpuff);
    toneON = find(leftTone_rounded);
    
    firstPulse = find(redLaser);
    stimON = firstPulse;
    
    % Initialize trial table
    varTypes = {'double','logical','logical','logical','logical',...
                'double','double','double',...
                'double','double','double'};
    varNames = {'TrialNumber','isReward','isPunish','isTone','isStim',...
                'CueTime','OutcomeTime','NextCue',...
                'nLicks','nAnticipatoryLicks','FirstLickTime'};
    trials = table('Size',[length(allTrialsTime) length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    
    % Loop over session
    for i=1:length(allTrialsTime)
        cur_cue = allTrialsTime(i);
        if i == length(allTrialsTime); next_cue = length(redLaser);
        else; next_cue = allTrialsTime(i+1); end
    
        % Trial cue/stim
        toneTime = toneON(toneON >= cur_cue-gracePeriod & toneON < next_cue-gracePeriod) - cur_cue;
        isTone = ~isempty(toneTime);
        stimTime = stimON(stimON >= cur_cue-gracePeriod & stimON < next_cue-gracePeriod) - cur_cue;
        isStim = ~isempty(stimTime);
    
        rewardTime = 0; punishTime = 0; toneTime = 0;
        gracePeriod = 0.1 * params.finalFs;
        rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue-gracePeriod & rightSolenoidON<next_cue-gracePeriod) - cur_cue;
        punishTime = airpuffON(airpuffON>=cur_cue-gracePeriod & airpuffON<next_cue-gracePeriod) - cur_cue;
        toneTime = toneON(toneON >= cur_cue-gracePeriod & toneON < next_cue-gracePeriod) - cur_cue;
        isReward = ~isempty(rewardTime);
        isPunish = ~isempty(punishTime);
        isTone = ~isempty(toneTime);
        if isReward && isPunish; outcomeTime = min([rewardTime,punishTime]);
        elseif ~isReward && ~isPunish; outcomeTime = min([rewardTime,punishTime,0]);
        else; outcomeTime = max([rewardTime,punishTime]); 
        end
    
        trialLicksTime = rightLickON(rightLickON > cur_cue & rightLickON < next_cue)  - cur_cue;
        if isempty(trialLicksTime); firstLick = nan;
        else; firstLick = trialLicksTime(1); end
        anticipatoryLicks = trialLicksTime(trialLicksTime <= outcomeTime);
    
        % Trial table
        trials(i,:) = {i,isReward,isPunish,isTone,isStim,...
            cur_cue,outcomeTime,next_cue,...
            length(trialLicksTime),length(anticipatoryLicks),firstLick};
    
        % Save to sync.mat
        save(strcat(sessionpath,'\','sync_',session.name),'trials','-append');
        disp('Finished: trial table saved');
    end
end

%% Calculate PSTHs

binSize = params.finalTimeStep; 
timeRange = [-1,5]; lick_binSize = 0.1;



% For pairing or random stim
waterIdx = find(rightSolenoid);  
% stimIdx = find(redLaser); % for random stim
% stimIdx = trials{trials.isReward == 0,"CueTime"}; % for pairing

% For flexible learning
stimIdx = trials{trials.isTone == 0,"CueTime"};
toneIdx = trials{trials.isTone == 1,"CueTime"};
% hitIdx = trials{,"CueTime"};

% 1. Calculate photometry PSTHs
baselineIdx = find(blueLaser,500);
waterIdxinLJ = findCorrespondingTime(waterIdx,timeNI,timePhotometry);
baselineInLJ = findCorrespondingTime(baselineIdx,timeNI,timePhotometry);
stimIdxinLJ = findCorrespondingTime(stimIdx,timeNI,timePhotometry);

[waterTraces,t] = getTraces(waterIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
[baselineTraces,~] = getTraces(baselineInLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
[stimTraces,~] = getTraces(stimIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);

% 2. Calculate lick PSTHs
[~,lickTraces_water] = getLicks(timeRange,waterIdx,lick_binSize,[],rightLick,session.behaviorFs,timeNI);
[~,lickTraces_stim] = getLicks(timeRange,stimIdx,lick_binSize,[],rightLick,session.behaviorFs,timeNI);
[~,lickTraces_shutter] = getLicks(timeRange,baselineIdx,lick_binSize,[],rightLick,session.behaviorFs,timeNI);
lickRate_water = lickTraces_water{2}/lick_binSize;
lickRate_stim = lickTraces_stim{2}/lick_binSize;
lickRate_shutter = lickTraces_shutter{2}/lick_binSize;

fig = initializeFig(0.5,0.5);
subplot(2,1,1)
plotCI(t,baselineTraces(1:end,:),[.75 .75 .75]);
plotCI(t,waterTraces(1:end,:),blueWhiteRed(1,:));
plotCI(t,stimTraces(1:end,:),blueWhiteRed(end,:));
plotEvent('',0,'r');
xlabel('Time (s)'); ylabel('z-score'); %legend('Shutter','Water','Stim');
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    ['Water (n=',num2str(length(waterIdx)),')'],...
    ['Stim (n=',num2str(length(stimIdx)),')']},...
    'Location','best'); 

subplot(2,1,2)
t = linspace(timeRange(1),timeRange(2),size(lickRate_shutter,2));
plotCI(t,lickRate_shutter,[.75 .75 .75]); hold on
plotCI(t,lickRate_water,blueWhiteRed(1,:)); hold on
plotCI(t,lickRate_stim,blueWhiteRed(end,:)); hold on
plotEvent('',0,'r');
xlabel('Time (s)'); ylabel('Licks/s'); %legend('Shutter','Water','Stim'); 
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    ['Water (n=',num2str(length(waterIdx)),')'],...
    ['Stim (n=',num2str(length(stimIdx)),')']},...
    'Location','best'); 
saveas(fig,strcat(sessionpath,'\psth_combined_',session.name,'.png'));

%% Plot photometry PSTHs

eventIdx = stimIdx; eventInLJ = zeros(size(eventIdx));
label = 'Stim'; eventDuration = 0.5;
% eventIdx = waterIdx; eventInLJ = zeros(size(eventIdx));
% label = 'Water'; eventDuration = 0;

longTimeRange = [-5,10];
shortTimeRange = [-1,5];

baselineIdx = find(blueLaser,length(eventIdx)); baselineInLJ = zeros(size(eventIdx));
binSize = params.finalTimeStep; groupSize = 20; % num of trials to plot in one line

for i = 1:length(eventIdx)
    [~, eventTime_lj] = min(abs(timePhotometry-timeNI(eventIdx(i))));
    [~, baselineTime_lj] = min(abs(timePhotometry-timeNI(baselineIdx(i))));
    eventInLJ(i) = eventTime_lj;
    baselineInLJ(i) = baselineTime_lj;
end

fig = initializeFig(0.5,1);
subplot(4,1,1)
[traces,t] = getTraces(eventInLJ/params.finalFs,rollingGreenLP,longTimeRange,binSize);
[baseline,~] = getTraces(baselineInLJ/params.finalFs,rollingGreenLP,longTimeRange,binSize);
plotCI(t,baseline(1:end,:),[.75 .75 .75]);
plotCI(t,traces(1:end,:),blueWhiteRed(1,:));
plotEvent(label,eventDuration,blueWhiteRed(end,:));
xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label); 
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    [label,' (n=',num2str(length(eventIdx)),')']},...
    'Location','northeast'); 

subplot(4,1,3)
nLines = ceil(size(traces,1)/groupSize);
legendList = cell(nLines,1);
nColors = round(linspace(1,size(blueWhiteRed,1),nLines));
for i = 1:nLines
    startTrial = (i-1)*groupSize+1; 
    if i == nLines; endTrial = size(traces,1);
    else; endTrial = i*groupSize; end
    plotCI(t,traces(startTrial:endTrial,:),blueWhiteRed(nColors(i),:));
    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
end
plotEvent(label,eventDuration,blueWhiteRed(end,:));
legend(legendList);

subplot(4,1,2)
[traces,t] = getTraces(eventInLJ/params.finalFs,rollingGreenLP,shortTimeRange,binSize);
[baseline,~] = getTraces(baselineInLJ/params.finalFs,rollingGreenLP,shortTimeRange,binSize);
plotCI(t,baseline(1:end,:),[.75 .75 .75]);
plotCI(t,traces(1:end,:),blueWhiteRed(1,:));
plotEvent(label,eventDuration,blueWhiteRed(end,:)); 
xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label);
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    [label,' (n=',num2str(length(eventIdx)),')']},...
    'Location','northeast'); 

subplot(4,1,4)
nLines = ceil(size(traces,1)/groupSize);
legendList = cell(nLines,1);
nColors = round(linspace(1,size(blueWhiteRed,1),nLines));
for i = 1:nLines
    startTrial = (i-1)*groupSize+1; 
    if i == nLines; endTrial = size(traces,1);
    else; endTrial = i*groupSize; end
    plotCI(t,traces(startTrial:endTrial,:),blueWhiteRed(nColors(i),:));
    legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
end
plotEvent(label,eventDuration,blueWhiteRed(end,:));
legend(legendList);

saveas(gcf,strcat(sessionpath,'\psth_',label,'_',session.name,'.png'));
return

%% Trial history

lick_summary_fig = figure('Position', get(0,'Screensize'));
eventIdx = find(allTrials); 
drawTrials_optoPair(timeRange,eventIdx,[],rightLick,nidq,timeNI,[],rightSolenoid,airpuff);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
        'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off'); 
box off

saveas(lick_summary_fig,strcat(sessionpath,'\summary_lick_',session.name,'.png'));

%% Calculate trial PSTHs (lick rate, eye)

eventIdx = find(redLaser);
eyeTraces_reward = getEye(timeRange,eventIdx,eye_pixel_detrend,sync,timeNI,timeCamera);
[lickSummary_reward,lickTraces_reward] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);


lickRate_reward = lickTraces_reward{2}/lick_binSize;
lickRate_punish = lickTraces_punish{2}/lick_binSize;
lickRate_hit = lickTraces_hit{2}/lick_binSize;
lickRate_miss = lickTraces_miss{2}/lick_binSize;
lickRate_fa = lickTraces_fa{2}/lick_binSize;
lickRate_cr = lickTraces_cr{2}/lick_binSize;

%% Plot trial PSTHs (lick rate, eye)

psth_summary_fig = figure('Position', get(0,'Screensize'));

subplot(2,2,1)
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotCI(t,lickRate_reward,colors(1)); hold on
plotCI(t,lickRate_punish,colors(2)); hold on
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); legend(["Reward","Punish"]);
xlim([timeRange(1),timeRange(2)]);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
    'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,2)
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotCI(t,eyeTraces_reward,colors(1)); hold on
plotCI(t,eyeTraces_punish,colors(2)); hold on
xlabel("Time (s)"); ylabel("Eye closing (a.u.)"); legend(["Reward","Punish"]);
xlim([timeRange(1),timeRange(2)]);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
    'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,3)
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotCI(t,lickRate_hit,blueWhiteRed(1,:)); hold on
plotCI(t,lickRate_miss,blueWhiteRed(150,:)); hold on
plotCI(t,lickRate_fa,blueWhiteRed(350,:)); hold on
plotCI(t,lickRate_cr,blueWhiteRed(500,:)); hold on
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); legend(["Hit","Miss","FA","CR"]);
xlim([timeRange(1),timeRange(2)]);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
    'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,4)
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotCI(t,eyeTraces_hit,blueWhiteRed(1,:)); hold on
plotCI(t,eyeTraces_miss,blueWhiteRed(150,:)); hold on
plotCI(t,eyeTraces_fa,blueWhiteRed(350,:)); hold on
plotCI(t,eyeTraces_cr,blueWhiteRed(500,:)); hold on
xlabel("Time (s)"); ylabel("Eye closing (a.u.)"); legend(["Hit","Miss","FA","CR"]);
xlim([timeRange(1),timeRange(2)]);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
    'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

saveas(psth_summary_fig,strcat(sessionpath,'\summary_psth_',session.name,'.png'));

%% Calculate switch PSTHs

switchTrials = trials{trials.TrialInBlock == 1,"TrialNumber"}; 
switchTrials = switchTrials(2:end-1); % Remove the first and last switch

choiceLicks_rewardToPunish = []; choiceLicks_punishToReward = [];
meanEye_rewardToPunish = []; meanEye_punishToReward = [];

for i = 1:length(switchTrials)
    cur_trial = switchTrials(i);
    trialRange = cur_trial-(blockLength/2):1:cur_trial+(blockLength/2);

    % Collect behavior variables
    choiceLicks_switch = trials{trialRange,"nAnticipatoryLicks"}';
    meanEye_switch = trials{trialRange,"MeanEyeIntensity"}';

    if trials{cur_trial,"TrialType"} == 1 % switching to punish from reward
        choiceLicks_rewardToPunish = [choiceLicks_rewardToPunish; choiceLicks_switch];
        meanEye_rewardToPunish = [meanEye_rewardToPunish; meanEye_switch];
    else
        choiceLicks_punishToReward = [choiceLicks_punishToReward; choiceLicks_switch];
        meanEye_punishToReward = [meanEye_punishToReward; meanEye_switch];
    end
end

lickRate_rewardToPunish = choiceLicks_rewardToPunish/0.5;
lickRate_punishToReward = choiceLicks_punishToReward/0.5;

%% Calculate first reward/punish PSTHs (trial switch from animal's perspective)

%% Calculate outcome composition during switch

switchTrials = trials{trials.TrialInBlock == 1,"TrialNumber"}; 
switchTrials = switchTrials(2:end-1); % Remove the first and last switch

outcome_rewardToPunish = []; outcome_punishToReward = [];

for i = 1:length(switchTrials)
    cur_trial = switchTrials(i);
    trialRange = cur_trial-(blockLength/2):1:cur_trial+(blockLength/2);

    % Collect behavior variables
    outcome_switch = trials{trialRange,"Outcome"}';

    if trials{cur_trial,"TrialType"} == 1 % switching to punish from reward
        outcome_rewardToPunish = [outcome_rewardToPunish; outcome_switch];
    else
        outcome_punishToReward = [outcome_punishToReward; outcome_switch];
    end
end

hit_rewardToPunish = zeros(1,size(outcome_rewardToPunish,2)); 
hit_punishToReward = zeros(1,size(outcome_punishToReward,2));
miss_rewardToPunish = zeros(1,size(outcome_rewardToPunish,2)); 
miss_punishToReward = zeros(1,size(outcome_punishToReward,2));
fa_rewardToPunish = zeros(1,size(outcome_rewardToPunish,2)); 
fa_punishToReward = zeros(1,size(outcome_punishToReward,2));
cr_rewardToPunish = zeros(1,size(outcome_rewardToPunish,2)); 
cr_punishToReward = zeros(1,size(outcome_punishToReward,2));

for i = 1:size(outcome_rewardToPunish,2)
    col = outcome_rewardToPunish(:,i);
    [~,labels,perc] = groupcounts(col);
    
    for j = 1:length(labels)
        label = labels(j);
        if label == "H"; hit_rewardToPunish(i) = perc(j); end
        if label == "M"; miss_rewardToPunish(i) = perc(j); end
        if label == "FA"; fa_rewardToPunish(i) = perc(j); end
        if label == "CR"; cr_rewardToPunish(i) = perc(j); end
    end
end

for i = 1:size(outcome_punishToReward,2)
    col = outcome_punishToReward(:,i);
    [~,labels,perc] = groupcounts(col);
    
    for j = 1:length(labels)
        label = labels(j);
        if label == "H"; hit_punishToReward(i) = perc(j); end
        if label == "M"; miss_punishToReward(i) = perc(j); end
        if label == "FA"; fa_punishToReward(i) = perc(j); end
        if label == "CR"; cr_punishToReward(i) = perc(j); end
    end
end

%% Plot switch PSTHs

switch_summary_fig = figure('Position', get(0,'Screensize').*[1,1,1,1]);
trialRange = -10:1:10;

subplot(2,2,1);
plotCI(trialRange, lickRate_rewardToPunish,colors(2));
plotCI(trialRange, lickRate_punishToReward,colors(1));
xlabel("Distance from switch trial"); ylabel("Anticipatory lick rate (Hz)"); 
legend(["Reward to punish","Punish to reward"]); ylim([0 inf]);
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,2);
plotCI(trialRange, meanEye_rewardToPunish,colors(2));
plotCI(trialRange, meanEye_punishToReward,colors(1));
xlabel("Distance from switch trial"); ylabel("Mean eye intensity (a.u.)"); 
legend(["Reward to punish","Punish to reward"]);
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,3);
plotCI(trialRange,hit_rewardToPunish,blueWhiteRed(1,:)); hold on
plotCI(trialRange,miss_rewardToPunish,blueWhiteRed(150,:)); hold on
plotCI(trialRange,fa_rewardToPunish,blueWhiteRed(350,:)); hold on
plotCI(trialRange,cr_rewardToPunish,blueWhiteRed(500,:)); hold on
xlabel("Distance from switch trial"); ylabel("Percentage (%)"); 
legend(["Hit","Miss","FA","CR"]); ylim([0 inf]); title("Reward \rightarrow Punish");
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,4);
plotCI(trialRange,hit_punishToReward,blueWhiteRed(1,:)); hold on
plotCI(trialRange,miss_punishToReward,blueWhiteRed(150,:)); hold on
plotCI(trialRange,fa_punishToReward,blueWhiteRed(350,:)); hold on
plotCI(trialRange,cr_punishToReward,blueWhiteRed(500,:)); hold on
xlabel("Distance from switch trial"); ylabel("Percentage (%)"); 
legend(["Hit","Miss","FA","CR"]); ylim([0 inf]); title("Punish \rightarrow Reward");
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

saveas(switch_summary_fig,strcat(sessionpath,'\summary_switch_',session.name,'.png'));

%% Test photometry data (BS version)
eventIdx = find(processed.behavior.downSampled(7,:));
timeRange = [-5,10]; binSize = params.finalTimeStep;

figure;
[traces,t] = getTraces(eventIdx/params.finalFs,processed.photometry.signals{1},timeRange,binSize);
plotCI(t,traces,'g');