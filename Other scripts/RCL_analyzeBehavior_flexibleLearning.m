% Shun_analyzeBehavior_flexibleLearning
% Shun Li, 11/8/2022

%% Load data

clear; close all;
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
[~,~,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed] = loadColors;
             
% 1. Select session via uigetdir
sessionpath = uigetdir('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
clear dirsplit

% 2. Decide which task it is
answer = questdlg("Which task is this session running?",...
    "Select task","Flexible learning","OptoPair","Random",'Random');
switch answer
    case 'Flexible learning'; task = 'flexible learning';
    case 'OptoPair'; task = 'optopair';
    case 'Random'; task = 'random';
end

disp(strcat('********** ',sessionName,'**********'));
load(strcat(sessionpath,'\','sync_',sessionName,'.mat'));
if ~isfield(session,'name'); session.name = sessionName; end
disp(['Session ',sessionName,' loaded']);

% Set other trial related variables
timeRange = [-0.5,3]; minLicks = 1; reactionTimeSamp = 2 * nidq.Fs;
blockLength = 2000; toneDuration = 0.5;

lick_binSize = 0.1; blink_thresh = 2.5; % in turns of z score

%% Preprocess outcome and opto data

% Reward/punishment params
rewardUnit = 0.012; % 8ms opening to dispense 1ul water
rewardList = [0 2 5]; % in ul
punishList = [0 0.05];
toneList = [0 0.5 1]; % in sec 

% Opto stim params
nPulsePerPattern = 25; % PulseNum in arduino
pulseDuration = 0.005; pulseFreq = 50;
pulseInterval = (1/pulseFreq) - pulseDuration;

disp('Ongoing: preprocess outcome and opto data');

if ~exist('rightSolenoid_rounded','var')
    % leftSolenoid = leftSolenoid ./ rewardUnit;
    rightSolenoid = rightSolenoid ./ rewardUnit;
    
    % Round reward and tone
    leftTone_rounded = roundToTarget(leftTone, toneList); disp('Finished rounding: leftTone');
    % leftSolenoid_rounded = roundToTarget(leftSolenoid, rewardList); disp('Finished rounding: leftSolenoid');
    rightSolenoid_rounded = roundToTarget(rightSolenoid, rewardList); disp('Finished rounding: rightSolenoid');
    airpuff_rounded = roundToTarget(airpuff,punishList); disp('Finished rounding: airpuff');

    % Analyze rounded data instead
    % leftTone = leftTone_rounded; airpuff = airpuff_rounded;
    % leftSolenoid = leftSolenoid_rounded; 
    % rightSolenoid = rightSolenoid_rounded;
    
    % Save rounded data
    save(strcat(sessionpath,'\','sync_',session.name),...
        'leftTone_rounded','rightSolenoid_rounded','airpuff_rounded','-append');
    disp('Finished: rounding cue/outcome data');
end

if ~exist('firstPulse','var') 
    if ~isempty(find(redLaser, 1))
        % Find the first pulse of each stim pattern if nPulsePerPattern>1
        if nPulsePerPattern > 1
            allPulses = find(redLaser);
            intervalThreshold = 10000;
            temp_interval = [100000,diff(allPulses)];
            firstPulse = allPulses(temp_interval > intervalThreshold);
    
            % Save first pulse data
            save(strcat(sessionpath,'\','sync_',session.name),"firstPulse",'-append');
        else
            firstPulse = find(redLaser);
            save(strcat(sessionpath,'\','sync_',session.name),"firstPulse",'-append');
        end
        disp('Finished: first pulse data');
    else 
        firstPulse = [];
    end
end

%% Combine stim&tone to form trial start

% toneON = find(leftTone_rounded);
% 
% % Add stim and tone together
% combined = sort([toneON, firstPulse]);
% diff_combined = [5*params.sync.behaviorFs, diff(combined)];
% 
% % Combine stim and tone that are sufficiently close to each other
% seperationThreshold = 2.5; % in sec
% allTrials = combined(diff_combined > seperationThreshold*params.sync.behaviorFs);

if ~exist('trials','var') && strcmp(task,'random')
    [allTrials,~] = getTrials(find(leftTone_rounded),firstPulse,...
                          find(rightSolenoid),find(airpuff));
else
    [allTrials,~] = getTrials(find(leftTone_rounded),firstPulse);
end


%% Generate trial table (flexible learning)

if ~exist('trials','var') && strcmp(task,'flexible learning')
    disp('Ongoing: making trial table for flexible learning');
    % Find digital events
    allTrialsTime = allTrials;
    airpuffON = find(airpuff_rounded);
    % leftSolenoidON = find(leftSolenoid_rounded);
    rightSolenoidON = find(rightSolenoid_rounded);
    leftLickON = find(leftLick);
    rightLickON = find(rightLick);
    toneON = find(leftTone_rounded);
    stimON = firstPulse;
    
    % Initialize trial table
    % Selection time: time of last choice lick (before outcome)
    % Reaction time: time of first lick
    varTypes = {'double','double','double','double','double','string',...
                'logical','logical','logical','logical',...
                'double','double','double','double','double',...
                'double','double','double','double','double','double'};
    varNames = {'Block','TrialNumber','TrialInBlock','TrialType','Choice','Outcome',...
                'isTone','isStim','isReward','isPunishment',...
                'CueTime','ReactionTime','OutcomeTime','LastLickTime','NextCue',...
                'RewardSize','PunishSize',...
                'nLicks','nAnticipatoryLicks','MeanEyeIntensity','PeakEyeIntensity'};
    
    % Initialize trial subtypes
    trials = table('Size',[length(allTrialsTime) length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    go = []; nogo = [];
    hit = []; miss = []; fa = []; cr = [];
    
    gracePeriod = floor(0.2 * params.finalFs);
    
    % Loop over all trials
    for i=1:length(allTrialsTime)
        cur_cue = allTrialsTime(i);
        if i == length(allTrialsTime); next_cue = length(allTrials);
        else; next_cue = allTrialsTime(i+1); end
    
        % Trial related
        TrialInBlock = mod(i,blockLength);
        if TrialInBlock == 0; TrialInBlock = blockLength; end
        block = mod(floor((i-1)/blockLength),2); % 0: reward, 1: punish
        trialType = sum(leftTone_rounded(cur_cue-gracePeriod:cur_cue+gracePeriod));
    
        % Trial cue/stim
        toneTime = toneON(toneON >= cur_cue-gracePeriod & toneON < next_cue-gracePeriod) - cur_cue;
        isTone = ~isempty(toneTime);
        stimTime = stimON(stimON >= cur_cue-gracePeriod & stimON < next_cue-gracePeriod) - cur_cue;
        isStim = ~isempty(stimTime);
    
        % Trial outcome
        rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue-gracePeriod & rightSolenoidON<next_cue-gracePeriod) - cur_cue;
        isReward = ~isempty(rewardTime);
        if isReward; rewardSize = sum(rightSolenoid_rounded(rewardTime + cur_cue));
        else; rewardSize = 0; rewardTime = reactionTimeSamp; end
        punishTime = airpuffON(airpuffON>=cur_cue-gracePeriod & airpuffON<next_cue-gracePeriod) - cur_cue;
        isPunishment = ~isempty(punishTime);
        if isPunishment; punishSize = sum(airpuff_rounded(punishTime + cur_cue));
        else; punishSize = 0; punishTime = reactionTimeSamp; end
        outcomeTime = min([reactionTimeSamp, rewardTime, punishTime]);
    
        % Trial licks
        trial_l_licks = (leftLickON(leftLickON > cur_cue & leftLickON < next_cue) - cur_cue)';
        trial_l_licks = [trial_l_licks, ones(length(trial_l_licks),1)];
        trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
        trial_r_licks = [trial_r_licks, 2*ones(length(trial_r_licks),1)];
        % First col is lick time in sample, second col is side of lick (1->left, 2->right)
        trialLicks = sortrows([trial_l_licks;trial_r_licks]);
        nLicks = size(trialLicks,1);
        % Initialize trial table value
        reactionTime = 0; lastLickTime = 0;
        if ~isempty(trialLicks)
            if block == 0
                reactionTime = trialLicks(1,1); 
                lastLickTime = trialLicks(end,1); 
            else
                reactionTime = min([trialLicks(1,1),outcomeTime]); 
                lastLickTime = trialLicks(end,1);
            end
        end
    
        % Separate choice licks (i.e. anticipatory licks: licks before outcome)
        if isempty(trialLicks); choiceLicks = [];
        else
            choiceLicks = trialLicks(trialLicks(:,1)<=outcomeTime,:); 
            consecLicks = getConsecutive(choiceLicks(:,2)); % Get consecutive lick number
        end
        nAnticipatoryLicks = size(choiceLicks,1);
    
        % Trial choice (1: go; 0: nogo)
        if isempty(choiceLicks); nogo = [nogo;cur_cue]; choice = 0;
        else
            if consecLicks(end) >= minLicks && choiceLicks(end,2) == 1
                go = [go;cur_cue]; choice = 1;
            elseif consecLicks(end) >= minLicks && choiceLicks(end,2) == 2
                go = [go;cur_cue]; choice = 1;
            else; nogo = [nogo;cur_cue]; choice = 0;
            end
        end
    
        % Update outcomeTime for omission trials (no outcome but licked)
        if ~(isReward || isPunishment) && (choice ~= 0)
            choiceLicks_idx = find(consecLicks>=minLicks,minLicks);
            choiceLicks = choiceLicks(1:choiceLicks_idx(end),:);
            outcomeTime = min(outcomeTime,choiceLicks(end,1));
        end
    
        % Trial outcome
        if block == 0 % Reward block
            if choice == 1
                hit = [hit; cur_cue]; outcome = 'H';
            elseif choice == 0
                miss = [miss; cur_cue]; outcome = 'M';
            end
        else % Punish block
            if choice == 1 && isPunishment
                fa = [fa; cur_cue]; outcome = 'FA';
            elseif choice == 0 && ~isPunishment
                cr = [cr; cur_cue]; outcome = 'CR';
            elseif choice == 1 && isReward
                hit = [hit; cur_cue]; outcome = 'H';
            elseif choice == 0 && isReward
                miss = [miss; cur_cue]; outcome = 'M';
            end
        end
    
        % Extract eye intensity
        if session.withCamera == 1
            cur_eye = getEye(timeRange,cur_cue,eye_pixel_detrend,session,timeNI,timeCamera);
            % Average eye intensity
            mean_eye = mean(cur_eye);
            % Peak eye intensity
            peak_eye = max(cur_eye);
        else 
            mean_eye = 0;
            peak_eye = 0;
        end
    
        % Trial table
        trials(i,:) = {block,i,TrialInBlock,trialType,choice,outcome,...
            isTone,isStim,isReward,isPunishment,...
            cur_cue,reactionTime,outcomeTime,lastLickTime,next_cue,...
            rewardSize,punishSize,nLicks,nAnticipatoryLicks,mean_eye,peak_eye};
    end
    
    % Sanity check
    disp(['Total hit = ',num2str(length(hit))]);
    disp(['Total miss = ',num2str(length(miss))]);
    disp(['Total FA = ',num2str(length(fa))]);
    disp(['Total CR = ',num2str(length(cr))]);
    
    % Save to sync.mat
    save(strcat(sessionpath,'\','sync_',session.name),'trials','-append');
    disp('Finished: trial table saved');
end

%% Generate trial table (optoPair)

if ~exist('trials','var') && strcmp(task,'optopair')
    disp('Ongoing: making trial table for optopair experiments');
    % Find digital events
    allTrialsTime = allTrials;
    airpuffON = find(airpuff_rounded);
    % leftSolenoidON = find(leftSolenoid_rounded);
    rightSolenoidON = find(rightSolenoid_rounded);
    leftLickON = find(leftLick);
    rightLickON = find(rightLick);
    toneON = find(leftTone_rounded);
    stimON = firstPulse;
    
    % Initialize trial table
    % Selection time: time of last choice lick (before outcome)
    % Reaction time: time of first lick
    varTypes = {'double','double','double','double','double','string',...
                'logical','logical','logical','logical',...
                'double','double','double','double','double',...
                'double','double','double','double','double','double'};
    varNames = {'Block','TrialNumber','TrialInBlock','TrialType','Choice','Outcome',...
                'isTone','isStim','isReward','isPunishment',...
                'CueTime','ReactionTime','OutcomeTime','LastLickTime','NextCue',...
                'RewardSize','PunishSize',...
                'nLicks','nAnticipatoryLicks','MeanEyeIntensity','PeakEyeIntensity'};
    
    % Initialize trial subtypes
    trials = table('Size',[length(allTrialsTime) length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    go = []; nogo = [];
    hit = []; miss = []; fa = []; cr = [];
    
    gracePeriod = floor(0.2 * params.finalFs);
    
    % Loop over all trials
    for i=1:length(allTrialsTime)
        cur_cue = allTrialsTime(i);
        if i == length(allTrialsTime); next_cue = length(allTrials);
        else; next_cue = allTrialsTime(i+1); end
    
        % Trial related
        TrialInBlock = mod(i,blockLength);
        if TrialInBlock == 0; TrialInBlock = blockLength; end
        block = mod(floor((i-1)/blockLength),2); % 0: reward, 1: punish
        trialType = sum(leftTone_rounded(cur_cue-gracePeriod:cur_cue+gracePeriod));
    
        % Trial cue/stim
        toneTime = toneON(toneON >= cur_cue-gracePeriod & toneON < next_cue-gracePeriod) - cur_cue;
        isTone = ~isempty(toneTime);
        stimTime = stimON(stimON >= cur_cue-gracePeriod & stimON < next_cue-gracePeriod) - cur_cue;
        isStim = ~isempty(stimTime);
    
        % Trial outcome
        rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue-gracePeriod & rightSolenoidON<next_cue-gracePeriod) - cur_cue;
        isReward = ~isempty(rewardTime);
        if isReward; rewardSize = sum(rightSolenoid_rounded(rewardTime + cur_cue));
        else; rewardSize = 0; rewardTime = reactionTimeSamp; end
        punishTime = airpuffON(airpuffON>=cur_cue-gracePeriod & airpuffON<next_cue-gracePeriod) - cur_cue;
        isPunishment = ~isempty(punishTime);
        if isPunishment; punishSize = sum(airpuff_rounded(punishTime + cur_cue));
        else; punishSize = 0; punishTime = reactionTimeSamp; end
        outcomeTime = min([reactionTimeSamp, rewardTime, punishTime]);
    
        % Trial licks
        trial_l_licks = (leftLickON(leftLickON > cur_cue & leftLickON < next_cue) - cur_cue)';
        trial_l_licks = [trial_l_licks, ones(length(trial_l_licks),1)];
        trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
        trial_r_licks = [trial_r_licks, 2*ones(length(trial_r_licks),1)];
        % First col is lick time in sample, second col is side of lick (1->left, 2->right)
        trialLicks = sortrows([trial_l_licks;trial_r_licks]);
        nLicks = size(trialLicks,1);
        % Initialize trial table value
        reactionTime = 0; lastLickTime = 0;
        if ~isempty(trialLicks)
            if block == 0
                reactionTime = trialLicks(1,1); 
                lastLickTime = trialLicks(end,1); 
            else
                reactionTime = min([trialLicks(1,1),outcomeTime]); 
                lastLickTime = trialLicks(end,1);
            end
        end
    
        % Separate choice licks (i.e. anticipatory licks: licks before outcome)
        if isempty(trialLicks); choiceLicks = [];
        else
            choiceLicks = trialLicks(trialLicks(:,1)<=outcomeTime,:); 
            consecLicks = getConsecutive(choiceLicks(:,2)); % Get consecutive lick number
        end
        nAnticipatoryLicks = size(choiceLicks,1);
    
        % Trial choice (1: go; 0: nogo)
        if isempty(choiceLicks); nogo = [nogo;cur_cue]; choice = 0;
        else
            if consecLicks(end) >= minLicks && choiceLicks(end,2) == 1
                go = [go;cur_cue]; choice = 1;
            elseif consecLicks(end) >= minLicks && choiceLicks(end,2) == 2
                go = [go;cur_cue]; choice = 1;
            else; nogo = [nogo;cur_cue]; choice = 0;
            end
        end
    
        % Update outcomeTime for omission trials (no outcome but licked)
        if ~(isReward || isPunishment) && (choice ~= 0)
            choiceLicks_idx = find(consecLicks>=minLicks,minLicks);
            choiceLicks = choiceLicks(1:choiceLicks_idx(end),:);
            outcomeTime = min(outcomeTime,choiceLicks(end,1));
        end
    
        % Trial outcome
        if block == 0 % Reward block
            if choice == 1
                hit = [hit; cur_cue]; outcome = 'H';
            elseif choice == 0
                miss = [miss; cur_cue]; outcome = 'M';
            end
        else % Punish block
            if choice == 1 && isPunishment
                fa = [fa; cur_cue]; outcome = 'FA';
            elseif choice == 0 && ~isPunishment
                cr = [cr; cur_cue]; outcome = 'CR';
            elseif choice == 1 && isReward
                hit = [hit; cur_cue]; outcome = 'H';
            elseif choice == 0 && isReward
                miss = [miss; cur_cue]; outcome = 'M';
            end
        end
    
        % Extract eye intensity
        if session.withCamera == 1
            cur_eye = getEye(timeRange,cur_cue,eye_pixel_detrend,session,timeNI,timeCamera);
            % Average eye intensity
            mean_eye = mean(cur_eye);
            % Peak eye intensity
            peak_eye = max(cur_eye);
        else 
            mean_eye = 0;
            peak_eye = 0;
        end
    
        % Trial table
        trials(i,:) = {block,i,TrialInBlock,trialType,choice,outcome,...
            isTone,isStim,isReward,isPunishment,...
            cur_cue,reactionTime,outcomeTime,lastLickTime,next_cue,...
            rewardSize,punishSize,nLicks,nAnticipatoryLicks,mean_eye,peak_eye};
    end
    
    % Sanity check
    disp(['Total hit = ',num2str(length(hit))]);
    disp(['Total miss = ',num2str(length(miss))]);
    disp(['Total FA = ',num2str(length(fa))]);
    disp(['Total CR = ',num2str(length(cr))]);
    
    % Save to sync.mat
    save(strcat(sessionpath,'\','sync_',session.name),'trials','-append');
    disp('Finished: trial table saved');
end
%% Generate trial table (random)
if ~exist('trials','var') && strcmp(task,'random')
    disp('Ongoing: making trial table for random outcome');
    % Find digital events
    allTrialsTime = allTrials;
    airpuffON = find(airpuff_rounded);
    rightSolenoidON = find(rightSolenoid_rounded);
    leftLickON = find(leftLick);
    rightLickON = find(rightLick);
    toneON = find(leftTone_rounded);
    stimON = firstPulse;
    
    % Initialize trial table
    % Selection time: time of last choice lick (before outcome)
    % Reaction time: time of first lick
    varTypes = {'double',...
                'logical','logical','logical','logical',...
                'double','double','double','double','double',...
                'double','double','double','double','double','double'};
    varNames = {'TrialNumber',...
                'isTone','isStim','isReward','isPunishment',...
                'CueTime','ReactionTime','OutcomeTime','LastLickTime','NextCue',...
                'RewardSize','PunishSize',...
                'nLicks','nAnticipatoryLicks','MeanEyeIntensity','PeakEyeIntensity'};
    
    % Initialize trial subtypes
    trials = table('Size',[length(allTrialsTime) length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);

    gracePeriod = floor(0.2 * params.finalFs);
    
    % Loop over all trials
    for i=1:length(allTrialsTime)
        cur_cue = allTrialsTime(i);
        if i == length(allTrialsTime); next_cue = length(allTrials);
        else; next_cue = allTrialsTime(i+1); end
    
        % Trial cue/stim
        toneTime = toneON(toneON >= cur_cue-gracePeriod & toneON < next_cue-gracePeriod) - cur_cue;
        isTone = ~isempty(toneTime);
        stimTime = stimON(stimON >= cur_cue-gracePeriod & stimON < next_cue-gracePeriod) - cur_cue;
        isStim = ~isempty(stimTime);
    
        % Trial outcome
        rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue-gracePeriod & rightSolenoidON<next_cue-gracePeriod) - cur_cue;
        isReward = ~isempty(rewardTime);
        if isReward; rewardSize = sum(rightSolenoid_rounded(rewardTime + cur_cue));
        else; rewardSize = 0; rewardTime = reactionTimeSamp; end
        punishTime = airpuffON(airpuffON>=cur_cue-gracePeriod & airpuffON<next_cue-gracePeriod) - cur_cue;
        isPunishment = ~isempty(punishTime);
        if isPunishment; punishSize = sum(airpuff_rounded(punishTime + cur_cue));
        else; punishSize = 0; punishTime = reactionTimeSamp; end
        outcomeTime = min([reactionTimeSamp, rewardTime, punishTime]);
    
        % Trial licks
        trial_l_licks = (leftLickON(leftLickON > cur_cue & leftLickON < next_cue) - cur_cue)';
        trial_l_licks = [trial_l_licks, ones(length(trial_l_licks),1)];
        trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
        trial_r_licks = [trial_r_licks, 2*ones(length(trial_r_licks),1)];
        % First col is lick time in sample, second col is side of lick (1->left, 2->right)
        trialLicks = sortrows([trial_l_licks;trial_r_licks]);
        nLicks = size(trialLicks,1);
        % Initialize trial table value
        reactionTime = 0; lastLickTime = 0;
        if ~isempty(trialLicks)
            reactionTime = trialLicks(1,1); 
            lastLickTime = trialLicks(end,1); 
        end
    
        % Separate choice licks (i.e. anticipatory licks: licks before outcome)
        if isempty(trialLicks); choiceLicks = [];
        else
            choiceLicks = trialLicks(trialLicks(:,1)<=outcomeTime,:); 
            consecLicks = getConsecutive(choiceLicks(:,2)); % Get consecutive lick number
        end
        nAnticipatoryLicks = size(choiceLicks,1);
    
        % Extract eye intensity
        if session.withCamera == 1
            cur_eye = getEye(timeRange,cur_cue,eye_pixel_detrend,session,timeNI,timeCamera);
            % Average eye intensity
            mean_eye = mean(cur_eye);
            % Peak eye intensity
            peak_eye = max(cur_eye);
        else 
            mean_eye = 0;
            peak_eye = 0;
        end
    
        % Trial table
        trials(i,:) = {i,...
            isTone,isStim,isReward,isPunishment,...
            cur_cue,reactionTime,outcomeTime,lastLickTime,next_cue,...
            rewardSize,punishSize,nLicks,nAnticipatoryLicks,mean_eye,peak_eye};
    end
    
    % Save to sync.mat
    save(strcat(sessionpath,'\','sync_',session.name),'trials','-append');
    disp('Finished: trial table saved');
end

%% Plot pre-processing steps

if session.ni_photometryON && session.withPhotometry
    initializeFig(1,1);
    
    % Raw photometry
    subplot(3,3,1);
    plot(photometry_raw);
    xlabel(['Time (',num2str(1000/nidq.Fs),' ms)']); ylabel('Signal (V)');
    title('NIDAQ photometry: raw');
    
    % After detrend (1min)
    subplot(3,3,4);
    plot(photometry_detrended); hold on
    xlabel(['Time (',num2str(1000/nidq.Fs),' ms)']); ylabel('z-score');
    title('NIDAQ photometry: after detrend (60s window)');
    
    % After downsample to 50Hz
    subplot(3,3,7);
    plot(photometryNI);
    xlabel(['Time (',num2str(1000/50),' ms)']); ylabel('z-score');
    title('NIDAQ photometry: Downsampled 100Hz -> rolling -> Downsampled to 50Hz');
    
    subplot(3,3,2)
    plot(demodGreen);xlabel('Time'); ylabel('z-score');
    title('Labjack photometry: detrend->demod');
    
    subplot(3,3,5)
    plot(rollingGreen);xlabel('Time'); ylabel('z-score');
    title('Labjack photometry: detrend->demod->rolling');
    
    subplot(3,3,8)
    plot(rollingGreenLP);xlabel('Time'); ylabel('z-score');
    title('Labjack photometry: detrend->demod->LP->rolling');
    
    % Create the uitable
    subplot(3,3,[3 6 9])
    histogram(normrnd(0,1,size(rollingGreenLP)),200); hold on
    histogram(rollingGreenLP,200); hold on
    histogram(photometryNI,200); hold on
    skew_lj = skewness(rollingGreenLP); kur_lj = kurtosis(rollingGreenLP);
    skew_ni = skewness(photometryNI); kur_ni = kurtosis(photometryNI);
    xlabel('z-score'); ylabel('Count'); legend({'Normal distribution','LJ Photometry','NI Photometry'});
    dim = [0.8205 0.58 0.55 0.27];
    str = {strcat("NI Skewness: ",num2str(skew_ni)),strcat("NI Kurtosis: ",num2str(kur_ni)),...
        strcat("LJ Skewness: ",num2str(skew_lj)),strcat("LJ Kurtosis: ",num2str(kur_lj))};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    title('Histogram of photometry traces');

elseif session.ni_photometryON && ~session.withPhotometry
    initializeFig(1,1);

    % Raw photometry
    subplot(3,2,1);
    plot(photometry_raw);
    xlabel(['Time (',num2str(1000/nidq.Fs),' ms)']); ylabel('Signal (V)');
    title('Raw photometry');
    
    % After detrend (1min)
    subplot(3,2,3);
    plot(photometry_detrended); hold on
    xlabel(['Time (',num2str(1000/nidq.Fs),' ms)']); ylabel('z-score');
    title('After detrend (60s window)');
    
    % After downsample to 50Hz
    subplot(3,2,5);
    plot(photometryNI);
    xlabel(['Time (',num2str(1000/50),' ms)']); ylabel('z-score');
    title('Downsampled to 50Hz');
    
    % Create the uitable
    subplot(3,2,[2 4 6]);
    histogram(normrnd(0,1,size(photometryNI)),200); hold on
    histogram(photometryNI,200); hold on
    xlabel('z-score'); ylabel('Count'); legend({'Normal distribution','Photometry'});
    dim = [0.8205 0.001 0.25 0.27];
    str = {strcat("Skewness: ",num2str(skewness(photometryNI))),strcat("Kurtosis: ",num2str(kurtosis(photometryNI)))};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    title('Histogram of z-scored photometry');

else
    initializeFig(1,1);

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
end

% Save figure
saveas(gcf,strcat(sessionpath,'\signal_summary_',session.name,'.png'));


%% Task specific params

if strcmp(task,'random')
%     waterIdx = find(rightSolenoid);  
%     airpuffIdx = find(airpuff_rounded);
%     stimIdx = firstPulse;
%     toneIdx = find(leftTone_rounded);
%     % pairIdx= [];
% %     taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
% %         ['Stim only (n=',num2str(length(stimIdx)),')'],...
% %         ['Airpuff (n=',num2str(length(airpuffIdx)),')']};

    waterIdx = find(rightSolenoid);  
    airpuffIdx = trials{trials.isPunishment == 1 & trials.isStim == 0,"CueTime"};
    toneIdx = trials{trials.isTone == 1 & trials.isStim == 0,"CueTime"};
    pairIdx = trials{trials.isPunishment == 1 & trials.isStim == 1,"CueTime"};
    stimIdx = trials{trials.isPunishment == 0 & trials.isStim == 1,"CueTime"};

    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
            ['Tone only (n=',num2str(length(toneIdx)),')'],...
            ['Airpuff (n=',num2str(length(airpuffIdx)),')'],...
            ['Stim only (n=',num2str(length(stimIdx)),')'],...
            ['Tone + stim (n=',num2str(length(airpuffIdx)),')'],...
            ['Airpuff (n=',num2str(length(airpuffIdx)),')']};

%     taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
%         ['Tone (n=',num2str(length(toneIdx)),')'],...
%         ['Airpuff (n=',num2str(length(airpuffIdx)),')']};
else
    waterIdx = find(rightSolenoid);  
    airpuffIdx = find(airpuff_rounded);
    stimIdx = trials{trials.isTone == 0 & trials.isStim == 1,"CueTime"};
    toneIdx = trials{trials.isTone == 1 & trials.isStim == 0,"CueTime"};
    pairIdx = trials{trials.isTone == 1 & trials.isStim == 1,"CueTime"};
    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
        ['Tone only (n=',num2str(length(toneIdx)),')'],...
        ['Stim only (n=',num2str(length(stimIdx)),')'],...
        ['Pair (n=',num2str(length(pairIdx)),')'],...
        ['Airpuff (n=',num2str(length(airpuffIdx)),')']};
end

if ~params.session.ni_photometryON; baselineIdx = find(blueLaser,500);
else; baselineIdx = []; end

%% (LJ) Plot combined photometry PSTHs

if session.withPhotometry
    timeRange = [-1,4]; lick_binSize = 0.2;
    
    % 2. Plot traces
    initializeFig(0.5,0.5);
    
    % 2.1 Plot photometry traces
    subplot(2,1,1)

    % Uncomment the bottom two lines when analyzing opto punish pair
    % [testTraces,t] = plotTraces(stimIdx,timeRange,photometry_50,bluePurpleRed(1,:),params,photometrySystem='ni',photometryFs=ni_photometryFs);
    % [~,~] = plotTraces(baselineIdx,timeRange,rollingGreenLP,[.75 .75 .75],params);
    [~,~] = plotTraces(waterIdx,timeRange,rollingGreenLP,bluePurpleRed(1,:),params);
    [~,~] = plotTraces(toneIdx,timeRange,rollingGreenLP,bluePurpleRed(350,:),params);
    %[~,~] = plotTraces(stimIdx,timeRange,rollingGreenLP,bluePurpleRed(end,:),params);
    %[~,~] = plotTraces(pairIdx,timeRange,rollingGreenLP,bluePurpleRed(150,:),params);
    [~,~] = plotTraces(airpuffIdx,timeRange,rollingGreenLP,[0.2, 0.2, 0.2],params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score'); %legend('Shutter','Water','Stim');
    legend(taskLegend,'Location','northeast'); ylim([-2,2.3]);
    
    % 2.2 Plot lick traces
    subplot(2,1,2)
    plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
    plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);

    %Uncomment the bottom two lines when analyzing opto punish pair
    %plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
    %plotLicks(pairIdx,timeRange,lick_binSize,bluePurpleRed(150,:),[],rightLick,params);
    plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('Licks/s'); %legend('Shutter','Water','Stim'); 
    legend(taskLegend,'Location','best'); 
    
    saveas(gcf,strcat(sessionpath,'\psth_lj_combined_',session.name,'.png'));
end

%% (LJ) Plot single stimulus PSTH

eventIdxes = {stimIdx,waterIdx,toneIdx,pairIdx,airpuffIdx};
labels = {'Stim','Water','Tone','Pair','Airpuff'};
eventDurations = [0.5,0,0.5,0.5,0.01];
groupSizes = [10,30,10,30,30];
longTimeRange = [-5,10];
shortTimeRange = [-0.5,3]; 

if session.withPhotometry
    for event = 1:length(eventIdxes)
        eventIdx = eventIdxes{event};
        if isempty(eventIdx); continue; end
        label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
        
        % binSize = params.finalTimeStep; 
        % baselineInLJ = findCorrespondingTime(baselineIdx,timeNI,timePhotometry);
        % eventInLJ = findCorrespondingTime(eventIdx,timeNI,timePhotometry);
        
        initializeFig(0.5,1);
        subplot(4,1,1)
        [~,~] = plotTraces(baselineIdx,longTimeRange,rollingGreenLP,[.75 .75 .75],params);
        [traces,t] = plotTraces(eventIdx,longTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
        plotEvent(label,eventDuration,bluePurpleRed(end,:));
        xlabel('Time (s)'); ylabel('z-score');
        legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
            [label,' (n=',num2str(length(eventIdx)),')']},...
            'Location','northeast'); 
        
        subplot(4,1,3)
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
        plotEvent(label,eventDuration,bluePurpleRed(end,:));
        legend(legendList);
        
        subplot(4,1,2)
        [~,~] = plotTraces(baselineIdx,shortTimeRange,rollingGreenLP,[.75 .75 .75],params);
        [traces,t] = plotTraces(eventIdx,shortTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
        plotEvent(label,eventDuration,bluePurpleRed(end,:)); 
        xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label);
        legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
            [label,' (n=',num2str(length(eventIdx)),')']},...
            'Location','northeast'); 
        
        subplot(4,1,4)
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
        plotEvent(label,eventDuration,bluePurpleRed(end,:));
        legend(legendList);
        
        saveas(gcf,strcat(sessionpath,'\psth_lj_',label,'_',session.name,'.png'));
    end
end

%% (NI) Plot combined PSTH

if session.ni_photometryON
    timeRange = [-1,5]; lick_binSize = 0.2;
    
    % 2. Plot traces
    initializeFig(0.5,0.5);
    
    % 2.1 Plot photometry traces
    subplot(2,1,1)
    [~,~] = plotTraces(waterIdx,timeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
    [~,~] = plotTraces(toneIdx,timeRange,photometryNI,bluePurpleRed(350,:),params,photometrySystem='ni');

    %Uncomment the bottom two lines when analyzing opto punish pair
%     [~,~] = plotTraces(stimIdx,timeRange,photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni');
%     [~,~] = plotTraces(pairIdx,timeRange,photometryNI,bluePurpleRed(150,:),params,photometrySystem='ni');

    [~,~] = plotTraces(airpuffIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni');
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score'); 
    legend(taskLegend,'Location','best'); 
    
    % 2.2 Plot lick traces
    subplot(2,1,2)
    plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
    plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);

    %Uncomment the bottom two lines when analyzing opto punish pair
%     plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
%     plotLicks(pairIdx,timeRange,lick_binSize,bluePurpleRed(150,:),[],rightLick,params);

    plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('Licks/s'); 
    legend(taskLegend,'Location','best');
    
    saveas(gcf,strcat(sessionpath,'\psth_ni_combined_',session.name,'.png'));
end

%% (NI) Plot single stimulus PSTH

eventIdxes = {stimIdx,waterIdx,toneIdx,pairIdx,airpuffIdx};
labels = {'Stim','Water','Tone','Pair','Airpuff'};
eventDurations = [0.5,0,0.5,0.5,0.01];
groupSizes = [10,30,10,30,30];
longTimeRange = [-5,10];
shortTimeRange = [-0.5,3]; 

binSize = 1/50; %binSize = params.finalTimeStep; 
Fs = params.sync.behaviorFs;
timestamp = timeNI;

if session.ni_photometryON
    for event = 1:length(eventIdxes)
        eventIdx = eventIdxes{event};
        if isempty(eventIdx); continue; end
        label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
        
        initializeFig(0.5,1);
        subplot(4,1,1)
        [~,~] = plotTraces(baselineIdx,longTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
        [traces,t] = plotTraces(eventIdx,longTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
        plotEvent(label,eventDuration,bluePurpleRed(end,:));
        xlabel('Time (s)'); ylabel('z-score');
        legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
            [label,' (n=',num2str(length(eventIdx)),')']},...
            'Location','northeast'); 
        
        subplot(4,1,3)
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
        plotEvent(label,eventDuration,bluePurpleRed(end,:));
        legend(legendList);
        
        subplot(4,1,2)
        [~,~] = plotTraces(baselineIdx,shortTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
        [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
        plotEvent(label,eventDuration,bluePurpleRed(end,:)); 
        xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label);
        legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
            [label,' (n=',num2str(length(eventIdx)),')']},...
            'Location','northeast'); 
        
        subplot(4,1,4)
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
        plotEvent(label,eventDuration,bluePurpleRed(end,:));
        legend(legendList);
        
        saveas(gcf,strcat(sessionpath,'\psth_ni_',label,'_',session.name,'.png'));
    end
end

return

%% Optosweeps for Rico


initializeFig(0.5,0.5);
timeRange = [-0.5,3];

% plotTraces(firstPulse(1:30),timeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
% plotTraces(firstPulse(31:60),timeRange,photometryNI,bluePurpleRed(100,:),params,photometrySystem='ni');
% plotTraces(firstPulse(61:90),timeRange,photometryNI,bluePurpleRed(200,:),params,photometrySystem='ni');
% plotTraces(firstPulse(91:120),timeRange,photometryNI,bluePurpleRed(300,:),params,photometrySystem='ni');
% plotTraces(firstPulse(121:150),timeRange,photometryNI,bluePurpleRed(400,:),params,photometrySystem='ni');
% plotTraces(firstPulse(151:180),timeRange,photometryNI,bluePurpleRed(500,:),params,photometrySystem='ni');
% plotTraces(firstPulse(181:210),timeRange,photometryNI,[.2,.8,.3],params,photometrySystem='ni');


plotTraces(airpuffIdx,timeRange,photometryNI,[.2,.2,.2],params,photometrySystem='ni');
plotTraces(stimIdx,timeRange,photometryNI,[178,102,255],params,photometrySystem='ni');
plotTraces(pairIdx,timeRange,photometryNI,bluePurpleRed(300,:),params,photometrySystem='ni');
plotTraces(toneIdx,timeRange,photometryNI,bluePurpleRed(400,:),params,photometrySystem='ni');

% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% plotTraces(waterIdx,timeRange,photometryNI,[.2,.8,.3],params,photometrySystem='ni');

% legend({'15z 250ms (n=30)',...
%         '20Hz 250ms (n=30)' ...
%         '25Hz 250ms (n=30)' ...
%         '30Hz 250ms (n=30)' ...
%         '35Hz 250ms (n=30)' ...
%         ['Airpuff (n=',num2str(length(airpuffIdx)),')'] ...
%         ['Reward (n=',num2str(length(waterIdx)),')']},...
%         'Location','northeast'); 

legend({
        ['Airpuff (n=',num2str(length(airpuffIdx)),')'] ...
        ['Tone (n=',num2str(length(toneIdx)),')'] ...
        ['Stim (n=',num2str(length(stimIdx)),')'] ...
        ['Airpuff + Stim (n=',num2str(length(pairIdx)),')']},...
        'Location','northeast'); 

% saveFigures(gcf,'Summary',sessionpath);

%% Optosweeps 2 (Rico)

initializeFig(0.5,0.5);

% plotTraces(firstPulse(1:30),timeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
plotTraces(firstPulse(31:60),timeRange,photometryNI,bluePurpleRed(100,:),params,photometrySystem='ni');
% plotTraces(firstPulse(61:90),timeRange,photometryNI,bluePurpleRed(200,:),params,photometrySystem='ni');
% plotTraces(firstPulse(91:120),timeRange,photometryNI,bluePurpleRed(300,:),params,photometrySystem='ni');
% plotTraces(firstPulse(121:150),timeRange,photometryNI,bluePurpleRed(400,:),params,photometrySystem='ni');
% plotTraces(firstPulse(151:180),timeRange,photometryNI,bluePurpleRed(500,:),params,photometrySystem='ni');

plotTraces(airpuffIdx,timeRange,photometryNI,[.2,.2,.2],params,photometrySystem='ni');
plotTraces(waterIdx,timeRange,photometryNI,[.2,.8,.3],params,photometrySystem='ni');

legend({'50Hz 200ms (n=30)',...
        ['Airpuff (n=',num2str(length(airpuffIdx)),')'] ...
        ['Reward (n=',num2str(length(waterIdx)),')']},...
        'Location','northeast'); 



%% Correlation analysis of GABASnFR trace



%% (NI + LJ) Plot airpuf/water on both sides

waterIdx = find(rightSolenoid);  
airpuffIdx = find(airpuff_rounded);

initializeFig(0.5,0.5);
timeRange = [-0.5,2];
[waterLJ,~] = plotTraces(waterIdx,timeRange,rollingGreenLP,blueWhiteRed(1,:),params);
[waterNI,~] = plotTraces(waterIdx,timeRange,photometryNI,blueWhiteRed(100,:),params,photometrySystem='ni');

[airpuffLJ,~] = plotTraces(airpuffIdx,timeRange,rollingGreenLP,blueWhiteRed(500,:),params);
[airpuffNI,~] = plotTraces(airpuffIdx,timeRange,photometryNI,blueWhiteRed(400,:),params,photometrySystem='ni');

plotEvent('',0,'r');
xlabel('Time (s)'); ylabel('z-score'); 
legend({['Water-LJ (n=',num2str(length(waterIdx)),')'],...
    ['Water-NI (n=',num2str(length(waterIdx)),')'],...
    ['Airpuff-LJ (n=',num2str(length(airpuffIdx)),')'],...
    ['Airpuff-NI (n=',num2str(length(airpuffIdx)),')']},...
    'Location','best'); 

% saveas(gcf,strcat(sessionpath,'\psth_combined_reward&airpuff_',session.name,'.png'));

% cross correlation
initializeFig(0.5,0.5);
[r_water, l_water] = xcorr(mean(waterLJ), mean(waterNI));
[r_airpuff, l_airpuff] = xcorr(mean(airpuffLJ), mean(airpuffNI));
plot(l_water,r_water); hold on
plot(l_airpuff,r_airpuff);

%% Trial type PSTH

%hitIdx = trials{trials.Outcome == 'H',"CueTime"};
%missIdx = trials{trials.Outcome == 'M',"CueTime"};
%faIdx = trials{trials.Outcome == 'FA',"CueTime"};
%crIdx = trials{trials.Outcome == 'CR',"CueTime"};

% bigRewardIdx = trials{trials.RewardSize == 5,"CueTime"};
% smallRewardIdx = trials{trials.RewardSize == 2,"CueTime"};

baseTrialIdx = trials{trials.TrialType == 0.5,"CueTime"};
surpriseTrialIdx = trials{trials.TrialType == 1,"CueTime"};

%hitIdxinLJ = findCorrespondingTime(hitIdx,timeNI,timePhotometry);
%missIdxinLJ = findCorrespondingTime(missIdx,timeNI,timePhotometry);

baseTrialIdxinLJ = findCorrespondingTime(baseTrialIdx,timeNI,timePhotometry);
surpriseTrialIdxinLJ = findCorrespondingTime(surpriseTrialIdx,timeNI,timePhotometry);

[baseTrialTraces,~] = getTraces(baseTrialIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
[surpriseTrialTraces,~] = getTraces(surpriseTrialIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);

[~,lickTraces_baseTrial] = getLicks(timeRange,baseTrialIdx,lick_binSize,[],rightLick,session.behaviorFs,timeNI);
[~,lickTraces_surpriseTrial] = getLicks(timeRange,surpriseTrialIdx,lick_binSize,[],rightLick,session.behaviorFs,timeNI);
lickRate_baseTrial = lickTraces_baseTrial{2}/lick_binSize;
lickRate_surpriseTrial = lickTraces_surpriseTrial{2}/lick_binSize;

fig = initializeFig(0.5,0.5);
subplot(2,1,1)
plotSEM(t,baselineTraces(1:end,:),[.75 .75 .75]);
plotSEM(t,baseTrialTraces(1:end,:),bluePurpleRed(1,:));
plotSEM(t,surpriseTrialTraces(1:end,:),bluePurpleRed(end,:));
plotEvent('',0,'r');
xlabel('Time (s)'); ylabel('z-score'); %legend('Shutter','Water','Stim');
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    ['Baseline trials (n=',num2str(length(baseTrialIdx)),')'],...
    ['Surprise trials (n=',num2str(length(surpriseTrialIdx)),')']},...
    'Location','best'); 

subplot(2,1,2)
t = linspace(timeRange(1),timeRange(2),size(lickRate_shutter,2));
plotSEM(t,lickRate_shutter,[.75 .75 .75],LineWidth=5); hold on
plotSEM(t,lickRate_baseTrial,bluePurpleRed(1,:)); hold on
plotSEM(t,lickRate_surpriseTrial,bluePurpleRed(end,:)); hold on
plotEvent('',0,'r');
xlabel('Time (s)'); ylabel('Licks/s'); %legend('Shutter','Water','Stim'); 
legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
    ['Baseline trials (n=',num2str(length(baseTrialIdx)),')'],...
    ['Surprise trials (n=',num2str(length(surpriseTrialIdx)),')']},...
    'Location','best'); 
saveas(fig,strcat(sessionpath,'\psth_trialType_',session.name,'.png'));
return

%% Trial history

lick_summary_fig = figure('Position', get(0,'Screensize'));
eventIdx = find(allTrials); 
drawTrials_flexibleLearning(timeRange,eventIdx,trials,[],rightLick,nidq,timeNI,[],rightSolenoid,airpuff);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
        'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off'); 
box off

saveas(lick_summary_fig,strcat(sessionpath,'\summary_lick_',session.name,'.png'));

%% Calculate trial PSTHs (lick rate, eye)

eventIdx = trials{trials.TrialType == 0,"CueTime"};
eyeTraces_reward = getEye(timeRange,eventIdx,eye_pixel_detrend,session,timeNI,timeCamera);
[lickSummary_reward,lickTraces_reward] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

eventIdx = trials{trials.TrialType == 1,"CueTime"};
eyeTraces_punish = getEye(timeRange,eventIdx,eye_pixel_detrend,session,timeNI,timeCamera);
[lickSummary_punish,lickTraces_punish] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

eventIdx = trials{trials.Outcome == 'H',"CueTime"};
eyeTraces_hit = getEye(timeRange,eventIdx,eye_pixel_detrend,session,timeNI,timeCamera);
[lickSummary_hit,lickTraces_hit] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

eventIdx = trials{trials.Outcome == 'M',"CueTime"};
eyeTraces_miss = getEye(timeRange,eventIdx,eye_pixel_detrend,session,timeNI,timeCamera);
[lickSummary_miss,lickTraces_miss] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

eventIdx = trials{trials.Outcome == 'FA',"CueTime"};
eyeTraces_fa = getEye(timeRange,eventIdx,eye_pixel_detrend,session,timeNI,timeCamera);
[lickSummary_fa,lickTraces_fa] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

eventIdx = trials{trials.Outcome == 'CR',"CueTime"};
eyeTraces_cr = getEye(timeRange,eventIdx,eye_pixel_detrend,session,timeNI,timeCamera);
[lickSummary_cr,lickTraces_cr] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

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
plotSEM(t,lickRate_reward,colors(1)); hold on
plotSEM(t,lickRate_punish,colors(2)); hold on
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); legend(["Reward","Punish"]);
xlim([timeRange(1),timeRange(2)]);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
    'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,2)
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotSEM(t,eyeTraces_reward,colors(1)); hold on
plotSEM(t,eyeTraces_punish,colors(2)); hold on
xlabel("Time (s)"); ylabel("Eye closing (a.u.)"); legend(["Reward","Punish"]);
xlim([timeRange(1),timeRange(2)]);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
    'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,3)
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotSEM(t,lickRate_hit,bluePurpleRed(1,:)); hold on
plotSEM(t,lickRate_miss,bluePurpleRed(150,:)); hold on
plotSEM(t,lickRate_fa,bluePurpleRed(350,:)); hold on
plotSEM(t,lickRate_cr,bluePurpleRed(500,:)); hold on
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); legend(["Hit","Miss","FA","CR"]);
xlim([timeRange(1),timeRange(2)]);
ylimit = ylim;
patch([0 toneDuration toneDuration 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
    'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,4)
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotSEM(t,eyeTraces_hit,bluePurpleRed(1,:)); hold on
plotSEM(t,eyeTraces_miss,bluePurpleRed(150,:)); hold on
plotSEM(t,eyeTraces_fa,bluePurpleRed(350,:)); hold on
plotSEM(t,eyeTraces_cr,bluePurpleRed(500,:)); hold on
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
plotSEM(trialRange, lickRate_rewardToPunish,colors(2));
plotSEM(trialRange, lickRate_punishToReward,colors(1));
xlabel("Distance from switch trial"); ylabel("Anticipatory lick rate (Hz)"); 
legend(["Reward to punish","Punish to reward"]); ylim([0 inf]);
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,2);
plotSEM(trialRange, meanEye_rewardToPunish,colors(2));
plotSEM(trialRange, meanEye_punishToReward,colors(1));
xlabel("Distance from switch trial"); ylabel("Mean eye intensity (a.u.)"); 
legend(["Reward to punish","Punish to reward"]);
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,3);
plotSEM(trialRange,hit_rewardToPunish,bluePurpleRed(1,:)); hold on
plotSEM(trialRange,miss_rewardToPunish,bluePurpleRed(150,:)); hold on
plotSEM(trialRange,fa_rewardToPunish,bluePurpleRed(350,:)); hold on
plotSEM(trialRange,cr_rewardToPunish,bluePurpleRed(500,:)); hold on
xlabel("Distance from switch trial"); ylabel("Percentage (%)"); 
legend(["Hit","Miss","FA","CR"]); ylim([0 inf]); title("Reward \rightarrow Punish");
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

subplot(2,2,4);
plotSEM(trialRange,hit_punishToReward,bluePurpleRed(1,:)); hold on
plotSEM(trialRange,miss_punishToReward,bluePurpleRed(150,:)); hold on
plotSEM(trialRange,fa_punishToReward,bluePurpleRed(350,:)); hold on
plotSEM(trialRange,cr_punishToReward,bluePurpleRed(500,:)); hold on
xlabel("Distance from switch trial"); ylabel("Percentage (%)"); 
legend(["Hit","Miss","FA","CR"]); ylim([0 inf]); title("Punish \rightarrow Reward");
xline(0,'-','Switch','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
box off

saveas(switch_summary_fig,strcat(sessionpath,'\summary_switch_',session.name,'.png'));

%% (Ver 2) Calculate photometry PSTHs

% binSize = params.finalTimeStep; 
% timeRange = [-1,5]; lick_binSize = 0.1;
% 
% waterIdx = find(rightSolenoid);  
% airpuffIdx = find(airpuff_rounded);
% stimIdx = trials{trials.isTone == 0,"CueTime"};
% toneIdx = trials{trials.isTone == 1,"CueTime"};
% 
% % 1. Calculate photometry PSTHs
% baselineIdx = find(blueLaser,500);
% waterIdxinLJ = findCorrespondingTime(waterIdx,timeNI,timePhotometry);
% toneIdxinLJ = findCorrespondingTime(toneIdx,timeNI,timePhotometry);
% baselineInLJ = findCorrespondingTime(baselineIdx,timeNI,timePhotometry);
% stimIdxinLJ = findCorrespondingTime(stimIdx,timeNI,timePhotometry);
% airpuffIdxinLJ = findCorrespondingTime(airpuffIdx,timeNI,timePhotometry);
% 
% [waterTraces,t] = getTraces(waterIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
% [toneTraces,~] = getTraces(toneIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
% [baselineTraces,~] = getTraces(baselineInLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
% [stimTraces,~] = getTraces(stimIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
% [airpuffTraces,~] = getTraces(airpuffIdxinLJ/params.finalFs,rollingGreenLP,timeRange,binSize);
% 
% % 2. Calculate lick PSTHs
% [lickRate_water,~,~] = getLicks(timeRange,waterIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_tone,~,~] = getLicks(timeRange,toneIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_stim,~,~] = getLicks(timeRange,stimIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_shutter,~,~] = getLicks(timeRange,baselineIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_airpuff,~,~] = getLicks(timeRange,baselineIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% lickRate_water = lickRate_water{2};
% lickRate_tone = lickRate_tone{2};
% lickRate_stim = lickRate_stim{2};
% lickRate_shutter = lickRate_shutter{2};
% lickRate_airpuff = lickRate_airpuff{2};
% 
% fig = initializeFig(0.5,0.5);
% subplot(2,1,1)
% plotSEM(t,baselineTraces,[.75 .75 .75]);
% plotSEM(t,waterTraces,bluePurpleRed(1,:));
% plotSEM(t,toneTraces,bluePurpleRed(350,:));
% plotSEM(t,stimTraces,bluePurpleRed(end,:));
% plotSEM(t,airpuffTraces,bluePurpleRed(150,:));
% plotEvent('',0,'r');
% xlabel('Time (s)'); ylabel('z-score'); %legend('Shutter','Water','Stim');
% legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
%     ['Water (n=',num2str(length(waterIdx)),')'],...
%     ['Tone (n=',num2str(length(toneIdx)),')'],...
%     ['Stim (n=',num2str(length(stimIdx)),')'],...
%     ['Airpuff (n=',num2str(length(airpuffIdx)),')']},...
%     'Location','best'); 
% 
% subplot(2,1,2)
% t = linspace(timeRange(1),timeRange(2),size(lickRate_shutter,2));
% plotSEM(t,lickRate_shutter,[.75 .75 .75]); 
% plotSEM(t,lickRate_water,bluePurpleRed(1,:));
% plotSEM(t,lickRate_tone,bluePurpleRed(350,:)); 
% plotSEM(t,lickRate_stim,bluePurpleRed(end,:)); 
% plotSEM(t,lickRate_airpuff,bluePurpleRed(150,:));
% plotEvent('',0,'r');
% xlabel('Time (s)'); ylabel('Licks/s'); %legend('Shutter','Water','Stim'); 
% legend({['Shutter (n=',num2str(length(baselineIdx)),')'],...
%     ['Water (n=',num2str(length(waterIdx)),')'],...
%     ['Tone (n=',num2str(length(toneIdx)),')'],...
%     ['Stim (n=',num2str(length(stimIdx)),')'],...
%     ['Airpuff (n=',num2str(length(airpuffIdx)),')']},...
%     'Location','best'); 
% saveas(fig,strcat(sessionpath,'\psth_combined_',session.name,'.png'));

%% (Ver 2) Plot licks
% [lickRate_water,~,~] = getLicks(timeRange,waterIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_tone,~,~] = getLicks(timeRange,toneIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_stim,~,~] = getLicks(timeRange,stimIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_shutter,~,~] = getLicks(timeRange,baselineIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_airpuff,~,~] = getLicks(timeRange,baselineIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% [lickRate_pair,~,~] = getLicks(timeRange,pairIdx,lick_binSize,[],rightLick,params.sync.behaviorFs,timeNI);
% 
% t = linspace(timeRange(1),timeRange(2),size(lickRate_shutter,2));
% % plotSEM(t,lickRate_shutter,[.75 .75 .75]); 
% plotSEM(t,lickRate_water,bluePurpleRed(1,:));
% plotSEM(t,lickRate_tone,bluePurpleRed(350,:)); 
% plotSEM(t,lickRate_stim,bluePurpleRed(end,:)); 
% plotSEM(t,lickRate_pair,bluePurpleRed(150,:));
% plotSEM(t,lickRate_airpuff,[0.2, 0.2, 0.2]);