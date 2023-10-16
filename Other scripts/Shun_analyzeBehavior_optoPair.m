% Shun_analyzeBehavior_optoPair
% Shun Li, 11/8/2022
% 02/14/2023: tidied up code, renamed to analyzeBehavior_optoPair

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
    "Select task","OptoPair reward","OptoPair punish","Random",'Random');
switch answer
    case 'OptoPair reward'; task = 'Optopair reward';
    case 'OptoPair punish'; task = 'Optopair punish';
    case 'Random'; task = 'random';
end

% Opto stim params
stim_prompt = {'Enter # pulses per pattern: ','Enter pulse duration: ',...
               'Enter pulse frequency: '};
definput = {'1','0.5','1'};
answer = inputdlg(stim_prompt,'Opto stim params',[1,35],definput);
nPulsePerPattern = str2double(answer{1}); pulseDuration = str2double(answer{2}); pulseFreq = str2double(answer{3});
pulseInterval = (1/pulseFreq) - pulseDuration;

disp(strcat('********** ',sessionName,'**********'));
load(strcat(sessionpath,'\','sync_',sessionName,'.mat'));
if ~isfield(session,'name'); session.name = sessionName; end
disp(['Session ',sessionName,' loaded']);

% Set other trial related variables
timeRange = [-0.5,3]; minLicks = 1; reactionTimeSamp = 2 * nidq.Fs;
toneDuration = 0.5;

lick_binSize = 0.2; blink_thresh = 2.5; % in turns of z score

%% Preprocess outcome and opto data

% Reward/punishment params
rewardUnit = 0.012; % 8ms opening to dispense 1ul water
rewardList = [0 2 5]; % in ul
punishList = [0 0.1];
toneList = [0 0.5 1]; % in sec

disp('Ongoing: preprocess outcome and opto data');

if ~exist('rightSolenoid_rounded','var')
    rightSolenoid = rightSolenoid ./ rewardUnit;
    
    % Round reward and tone
    leftTone_rounded = roundToTarget(leftTone, toneList); disp('Finished rounding: leftTone');
    rightSolenoid_rounded = roundToTarget(rightSolenoid, rewardList); disp('Finished rounding: rightSolenoid');
    airpuff_rounded = roundToTarget(airpuff,punishList); disp('Finished rounding: airpuff');
    
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

if ~exist('trials','var') && strcmp(task,'random')
    [allTrials,~] = getTrials(find(leftTone_rounded),firstPulse,...
                         find(rightSolenoid_rounded),find(airpuff_rounded));
%     [allTrials,~] = getTrials(find(leftTone_rounded),...
%                          find(rightSolenoid_rounded),find(airpuff_rounded));
else
    [allTrials,~] = getTrials(find(leftTone_rounded),firstPulse);
end

%% Generate trial table (optopair reward)

if ~exist('trials','var') && strcmp(task,'Optopair reward')
    disp('Ongoing: making trial table for optopair reward');
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
    varTypes = {'double','double','double','string',...
                'logical','logical','logical','logical',...
                'double','double','double','double','double',...
                'double','double','double','double','double','double'};
    varNames = {'TrialNumber','TrialType','Choice','Outcome',...
                'isTone','isStim','isReward','isPunishment',...
                'CueTime','ReactionTime','OutcomeTime','LastLickTime','NextCue',...
                'RewardSize','PunishSize',...
                'nLicks','nAnticipatoryLicks','MeanEyeIntensity','PeakEyeIntensity'};
    
    % Initialize trial subtypes
    trials = table('Size',[length(allTrialsTime) length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    go = []; nogo = [];
    hit = []; miss = []; fa = []; cr = [];
    
    gracePeriod = floor(0.2 * params.sync.behaviorFs);
    
    % Loop over all trials
    for i=1:length(allTrialsTime)
        cur_cue = allTrialsTime(i);
        if i == length(allTrialsTime); next_cue = length(allTrials);
        else; next_cue = allTrialsTime(i+1); end
    
        % Trial related
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
        trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
        trialLicks = sortrows([trial_r_licks, 2*ones(length(trial_r_licks),1)]);
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
    
        % Trial choice (1: go; 0: nogo)
        if isempty(choiceLicks); nogo = [nogo;cur_cue]; choice = 0;
        else
            if consecLicks(end) >= minLicks; go = [go;cur_cue]; choice = 1;
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
        if choice == 1; hit = [hit; cur_cue]; outcome = 'H';
        elseif choice == 0; miss = [miss; cur_cue]; outcome = 'M';
        end
    
        % Extract eye intensity
        if session.withCamera == 1
            cur_eye = getEye(timeRange,cur_cue,eye_pixel_detrend,session,timeNI,timeCamera);   
            mean_eye = mean(cur_eye); % Average eye intensity
            peak_eye = max(cur_eye); % Peak eye intensity
        else 
            mean_eye = 0; peak_eye = 0;
        end
    
        % Trial table
        trials(i,:) = {i,trialType,choice,outcome,...
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

%% Generate trial table (optoPair punish)

if ~exist('trials','var') && strcmp(task,'Optopair punish')
    disp('Ongoing: making trial table for optopair punish');
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
    varTypes = {'double','double','double',...
                'logical','logical','logical','logical',...
                'double','double','double','double','double',...
                'double','double','double','double','double','double'};
    varNames = {'TrialNumber','TrialType','Choice',...
                'isTone','isStim','isReward','isPunishment',...
                'CueTime','ReactionTime','OutcomeTime','LastLickTime','NextCue',...
                'RewardSize','PunishSize',...
                'nLicks','nAnticipatoryLicks','MeanEyeIntensity','PeakEyeIntensity'};
    
    % Initialize trial subtypes
    trials = table('Size',[length(allTrialsTime) length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    go = []; nogo = [];
    hit = []; miss = []; fa = []; cr = [];
    
    gracePeriod = floor(0.2 * params.sync.behaviorFs);
    
    % Loop over all trials
    for i=1:length(allTrialsTime)
        cur_cue = allTrialsTime(i);
        if i == length(allTrialsTime); next_cue = length(allTrials);
        else; next_cue = allTrialsTime(i+1); end
    
        % Trial related
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
        trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
        trialLicks = sortrows([trial_r_licks, 2*ones(length(trial_r_licks),1)]);
        nLicks = size(trialLicks,1);
        % Initialize trial table value
        reactionTime = 0; lastLickTime = 0;
        if ~isempty(trialLicks)
            reactionTime = min([trialLicks(1,1),outcomeTime]); 
            lastLickTime = trialLicks(end,1);
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
            if consecLicks(end) >= minLicks && choiceLicks(end,2) == 2
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
%         if choice == 1 && isPunishment
%             fa = [fa; cur_cue]; outcome = 'FA';
%         elseif choice == 0 && ~isPunishment
%             cr = [cr; cur_cue]; outcome = 'CR';
%         elseif choice == 1 && isReward
%             hit = [hit; cur_cue]; outcome = 'H';
%         elseif choice == 0 && isReward
%             miss = [miss; cur_cue]; outcome = 'M';
%         end
    
        % Extract eye intensity
        if session.withCamera == 1
            cur_eye = getEye(timeRange,cur_cue,eye_pixel_detrend,session,timeNI,timeCamera);
            mean_eye = mean(cur_eye); % Average eye intensity
            peak_eye = max(cur_eye); % Peak eye intensity
        else 
            mean_eye = 0; peak_eye = 0;
        end
    
        % Trial table
        trials(i,:) = {i,trialType,choice,...
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

    gracePeriod = floor(0.2 * params.sync.behaviorFs);
    
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
        trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
        trialLicks = sortrows([trial_r_licks, 2*ones(length(trial_r_licks),1)]);
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
    plot(modgreen);xlabel('Time'); ylabel('z-score');
    title('Labjack photometry: green modulation');
    
    subplot(3,3,5)
    plot(green);xlabel('Time'); ylabel('z-score');
    title('Labjack photometry: detrend');
    
    subplot(3,3,8)
    plot(rollingGreenLP);xlabel('Time'); ylabel('z-score');
    title('Labjack photometry: detrend->(demod)->LP->rolling');
    
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
    
    % Select event idx
    waterIdx = find(rightSolenoid);  
    airpuffIdx = find(airpuff_rounded);
    toneIdx = find(leftTone);
    stimIdx = firstPulse;
    randShutterIdx = round(rand([500,1])*length(blueLaser));
    
    % Create task legend
    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
            ['Airpuff (n=',num2str(length(airpuffIdx)),')'],...
            ['Tone only (n=',num2str(length(toneIdx)),')'],...
            ['Stim only (n=',num2str(length(stimIdx)),')'],...
            ['Random shutter (n=',num2str(length(randShutterIdx)),')']};

else
    waterIdx = find(rightSolenoid);  
    airpuffIdx = find(airpuff_rounded);
    stimIdx = trials{trials.isTone == 0 & trials.isStim == 1,"CueTime"};
    toneIdx = trials{trials.isTone == 1 & trials.isStim == 0,"CueTime"};
    pairIdx = trials{trials.isTone == 1 & trials.isStim == 1,"CueTime"};
    randShutterIdx = round(rand([500,1])*length(blueLaser));

    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
    ['Tone only (n=',num2str(length(toneIdx)),')'],...
    ['Stim only (n=',num2str(length(stimIdx)),')'],...
    ['Pair (n=',num2str(length(pairIdx)),')'],...
    ['Airpuff (n=',num2str(length(airpuffIdx)),')'],...
    ['Random shutter (n=',num2str(length(randShutterIdx)),')']};
end

%% (LJ) Plot combined photometry PSTHs

if session.withPhotometry
    timeRange = [-1,3]; lick_binSize = 0.2;
    
    % 2. Plot traces
    initializeFig(0.5,0.5);

    if strcmp(task,'random')
        % 2.1 Plot photometry traces
        subplot(2,1,1)
        [~,~] = plotTraces(waterIdx,timeRange,rollingGreenLP,bluePurpleRed(1,:),params);
        [~,~] = plotTraces(airpuffIdx,timeRange,rollingGreenLP,[0.2, 0.2, 0.2],params);
        [~,~] = plotTraces(toneIdx,timeRange,rollingGreenLP,bluePurpleRed(350,:),params);
        [~,~] = plotTraces(stimIdx,timeRange,rollingGreenLP,bluePurpleRed(end,:),params);
        [~,~] = plotTraces(randShutterIdx,timeRange,rollingGreenLP,[.75 .75 .75],params);
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('z-score');
        legend(taskLegend,'Location','northeast');
        
        % 2.2 Plot lick traces
        subplot(2,1,2)
        plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
        plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
        plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
        plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
        plotLicks(randShutterIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('Licks/s');  %legend('Shutter','Water','Stim'); 
        legend(taskLegend,'Location','best');

    else
        % 2.1 Plot photometry traces
        subplot(2,1,1)
        [~,~] = plotTraces(waterIdx,timeRange,rollingGreenLP,bluePurpleRed(1,:),params);
        [~,~] = plotTraces(toneIdx,timeRange,rollingGreenLP,bluePurpleRed(350,:),params);
        [~,~] = plotTraces(stimIdx,timeRange,rollingGreenLP,bluePurpleRed(end,:),params);
        [~,~] = plotTraces(pairIdx,timeRange,rollingGreenLP,bluePurpleRed(150,:),params);
        [~,~] = plotTraces(airpuffIdx,timeRange,rollingGreenLP,[0.2, 0.2, 0.2],params);
        [~,~] = plotTraces(randShutterIdx,timeRange,rollingGreenLP,[.75 .75 .75],params);
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('z-score');
        legend(taskLegend,'Location','northeast');
        
        % 2.2 Plot lick traces
        subplot(2,1,2)
        plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
        plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
        plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
        plotLicks(pairIdx,timeRange,lick_binSize,bluePurpleRed(150,:),[],rightLick,params);
        plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
        plotLicks(randShutterIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('Licks/s'); %legend('Shutter','Water','Stim'); 
        legend(taskLegend,'Location','best');
    end 
    
    saveas(gcf,strcat(sessionpath,'\psth_lj_combined_',session.name,'.png'));
end

%% (LJ) Plot single stimulus PSTH

if session.withPhotometry

    if strcmp(task,'random')
        eventIdxes = {stimIdx,waterIdx,toneIdx,airpuffIdx};
        labels = {'Stim','Water','Tone','Airpuff'};
        eventDurations = [0.5,0,0.5,0.01];
        groupSizes = [10,30,10,30];
        longTimeRange = [-5,10];
        shortTimeRange = [-0.5,3]; 
        
        if session.withPhotometry
            for event = 1:length(eventIdxes)
                eventIdx = eventIdxes{event};
                if isempty(eventIdx); continue; end
                label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
                
                initializeFig(0.5,1);
                subplot(4,1,1)
                [traces,t] = plotTraces(eventIdx,longTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(randShutterIdx,longTimeRange,rollingGreenLP,[.75 .75 .75],params);
                plotEvent(label,eventDuration,bluePurpleRed(end,:));
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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
                [traces,t] = plotTraces(eventIdx,shortTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(randShutterIdx,shortTimeRange,rollingGreenLP,[.75 .75 .75],params);
                plotEvent(label,eventDuration,bluePurpleRed(end,:)); 
                xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label);
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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
    else
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
                
                initializeFig(0.5,1);
                subplot(4,1,1)
                [traces,t] = plotTraces(eventIdx,longTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(randShutterIdx,longTimeRange,rollingGreenLP,[.75 .75 .75],params);
                plotEvent(label,eventDuration,bluePurpleRed(end,:));
                xlabel('Time (s)'); ylabel('z-score');
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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
                [traces,t] = plotTraces(eventIdx,shortTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
                [~,~] = plotTraces(randShutterIdx,shortTimeRange,rollingGreenLP,[.75 .75 .75],params);
                plotEvent(label,eventDuration,bluePurpleRed(end,:)); 
                xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label);
                legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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
    end
end

%% (NI) Plot combined PSTH

if session.ni_photometryON
    timeRange = [-1,3]; lick_binSize = 0.2;
    
    % 2. Plot traces
    initializeFig(0.5,0.5);

    if strcmp(task,'random')
        % 2.1 Plot photometry traces
        subplot(2,1,1)
        [~,~] = plotTraces(waterIdx,timeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
        [~,~] = plotTraces(airpuffIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni');
        [~,~] = plotTraces(toneIdx,timeRange,photometryNI,bluePurpleRed(350,:),params,photometrySystem='ni');
        [~,~] = plotTraces(stimIdx,timeRange,photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni');
        [~,~] = plotTraces(randShutterIdx,timeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('z-score'); 
        legend(taskLegend,'Location','best'); 
        
        % 2.2 Plot lick traces
        subplot(2,1,2)
        plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
        plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
        plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
        plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
        plotLicks(randShutterIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('Licks/s'); 
        legend(taskLegend,'Location','best');
    else
        % 2.1 Plot photometry traces
        subplot(2,1,1)
        [~,~] = plotTraces(waterIdx,timeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
        [~,~] = plotTraces(toneIdx,timeRange,photometryNI,bluePurpleRed(350,:),params,photometrySystem='ni');
        [~,~] = plotTraces(stimIdx,timeRange,photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni');
        [~,~] = plotTraces(pairIdx,timeRange,photometryNI,bluePurpleRed(150,:),params,photometrySystem='ni');
        [~,~] = plotTraces(airpuffIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni');
        [~,~] = plotTraces(randShutterIdx,timeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('z-score'); 
        legend(taskLegend,'Location','best'); 
        
        % 2.2 Plot lick traces
        subplot(2,1,2)
        plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
        plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
        plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
        plotLicks(pairIdx,timeRange,lick_binSize,bluePurpleRed(150,:),[],rightLick,params);
        plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
        plotLicks(randShutterIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
        plotEvent('',0,'r');
        xlabel('Time (s)'); ylabel('Licks/s'); 
        legend(taskLegend,'Location','best');
    end
 
    saveas(gcf,strcat(sessionpath,'\psth_ni_combined_',session.name,'.png'));
end

%% (NI) Plot single stimulus PSTH

if session.ni_photometryON

    if strcmp(task,'random')
        eventIdxes = {stimIdx,waterIdx,toneIdx,airpuffIdx};
        labels = {'Stim','Water','Tone','Airpuff'};
        eventDurations = [0.5,0,0.5,0.01];
        groupSizes = [10,30,10,30];
        longTimeRange = [-5,10];
        shortTimeRange = [-0.5,3]; 
        
        binSize = 1/50; %binSize = params.finalTimeStep; 
        Fs = params.sync.behaviorFs;
        timestamp = timeNI;
        
        for event = 1:length(eventIdxes)
            eventIdx = eventIdxes{event};
            if isempty(eventIdx); continue; end
            label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
            
            initializeFig(0.5,1);
            subplot(4,1,1)
            [traces,t] = plotTraces(eventIdx,longTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
            [~,~] = plotTraces(randShutterIdx,longTimeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
            plotEvent(label,eventDuration,bluePurpleRed(end,:));
            xlabel('Time (s)'); ylabel('z-score');
            legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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
            [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
            [~,~] = plotTraces(randShutterIdx,shortTimeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
            plotEvent(label,eventDuration,bluePurpleRed(end,:));
            xlabel('Time (s)'); ylabel('z-score');
            legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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
    else
        eventIdxes = {stimIdx,waterIdx,toneIdx,pairIdx,airpuffIdx};
        labels = {'Stim','Water','Tone','Pair','Airpuff'};
        eventDurations = [0.5,0,0.5,0.5,0.01];
        groupSizes = [10,30,30,30,30];
        longTimeRange = [-5,10];
        shortTimeRange = [-0.5,3]; 
        
        binSize = 1/50; %binSize = params.finalTimeStep; 
        Fs = params.sync.behaviorFs;
        timestamp = timeNI;
        
        for event = 1:length(eventIdxes)
            eventIdx = eventIdxes{event};
            if isempty(eventIdx); continue; end
            label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
            
            initializeFig(0.5,1);
            subplot(4,1,1)
            [traces,t] = plotTraces(eventIdx,longTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
            [~,~] = plotTraces(randShutterIdx,longTimeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
            plotEvent(label,eventDuration,bluePurpleRed(end,:));
            xlabel('Time (s)'); ylabel('z-score');
            legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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
            [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
            [~,~] = plotTraces(randShutterIdx,shortTimeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
            plotEvent(label,eventDuration,bluePurpleRed(end,:));
            xlabel('Time (s)'); ylabel('z-score');
            legend({[label,' (n=',num2str(length(eventIdx)),')'],...
                    ['Shutter (n=',num2str(length(randShutterIdx)),')']},...
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

    
end

%% Plot lick raster plot

timeRange = [-1,5]; markerSize = 20;

if strcmp(task,'Optopair reward')
    initializeFig(1,.67);
    
    % Get event time and number by trial type
    stimIdx = trials{trials.isTone == 0 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime"]};
    toneIdx = trials{trials.isTone == 1 & trials.isStim == 0,["TrialNumber","CueTime","OutcomeTime"]};
    pairIdx = trials{trials.isTone == 1 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime"]};
    stimIdx(:,3) = stimIdx(:,3)./params.sync.behaviorFs;
    toneIdx(:,3) = toneIdx(:,3)./params.sync.behaviorFs;
    pairIdx(:,3) = pairIdx(:,3)./params.sync.behaviorFs;

    % getLicks by trial type
    [stimLickRate,~,stimLicks] = getLicks(timeRange,stimIdx(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [toneLickRate,~,toneLicks] = getLicks(timeRange,toneIdx(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);
    [pairLickRate,~,pairLicks] = getLicks(timeRange,pairIdx(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI);


    % Plot overall raster plot (color coded by trial type)
    tiledlayout(3,2); nexttile([3,1]);
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
    plotEvent("",0.5,'r');

    
    % Plot lick traces across session (stim only)
    nexttile;
    traces = stimLickRate; groupSize = 10;
    t = linspace(timeRange(1),timeRange(2),size(traces,2));
    nLines = ceil(size(traces,1)/groupSize); legendList = cell(nLines,1);
    nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
    for i = 1:nLines
        startTrial = (i-1)*groupSize+1; 
        if i == nLines; endTrial = size(traces,1);
        else; endTrial = i*groupSize; end
        plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
        legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
    end
    plotEvent('Stim',0.5,bluePurpleRed(end,:));
    xlabel('Time (s)'); ylabel('Licks/s');
    legend(legendList);

    % Plot lick traces across session (pair)
    nexttile;
    traces = pairLickRate; groupSize = 30;
    t = linspace(timeRange(1),timeRange(2),size(traces,2));
    nLines = ceil(size(traces,1)/groupSize); legendList = cell(nLines,1);
    nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
    for i = 1:nLines
        startTrial = (i-1)*groupSize+1; 
        if i == nLines; endTrial = size(traces,1);
        else; endTrial = i*groupSize; end
        plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
        legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
    end
    plotEvent('Stim & tone',0.5,bluePurpleRed(end,:));
    xlabel('Time (s)'); ylabel('Licks/s');
    legend(legendList);

    % Plot lick traces across session (stim only)
    nexttile;
    traces = toneLickRate; groupSize = 10;
    t = linspace(timeRange(1),timeRange(2),size(traces,2));
    nLines = ceil(size(traces,1)/groupSize); legendList = cell(nLines,1);
    nColors = round(linspace(1,size(bluePurpleRed,1),nLines));
    for i = 1:nLines
        startTrial = (i-1)*groupSize+1; 
        if i == nLines; endTrial = size(traces,1);
        else; endTrial = i*groupSize; end
        plotSEM(t,traces(startTrial:endTrial,:),bluePurpleRed(nColors(i),:));
        legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial)];
    end
    plotEvent('Tone',0.5,bluePurpleRed(end,:));
    xlabel('Time (s)'); ylabel('Licks/s');
    legend(legendList);
    
    saveas(gcf,strcat(sessionpath,'\lick_summary_',session.name,'.png'));
    
end

return

%% Plot AUC during stim vs performance params

%% Test photometryNI