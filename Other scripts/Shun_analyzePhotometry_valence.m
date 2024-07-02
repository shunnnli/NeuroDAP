% Plot photometry signal aligned to behavioral event

clear; close all;
% addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Neuropixel analysis\Methods'));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed,blueGreenPurple] = loadColors;
             
             
sessionName = '20221201-SL043-D2-laser_g0';
% sessionName = "20221004-TP173-D9_g0";
% sessionName = "20220920-SL032-test-photometry2_g0";
session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Project valence\Recordings\';
% session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Test photometry\';
load(strcat(session.path,sessionName,'\','sync_',sessionName,'.mat'));

%% Old stim pattern log
% Find first laser sample time of each strength/pattern

% 20220518: 15, 10, 5, 0 mW
% npulse_per_stim = 20; stim_per_pattern = 10; nPatterns = 4;

% 20220527
% npulse_per_stim = 10; stim_per_pattern = 20; nPatterns = 8;

% 20220602
% npulse_per_stim = 10; stim_per_pattern = 20; nPatterns = 8;
% npulse_per_stim = 20; stim_per_pattern = 20; nPatterns = 2;

% 20220603
% titles = {'10 5mW 2ms pulse @ 20Hz','10 4mW 2ms pulse @ 20Hz',...
%           '10 3mW 2ms pulse @ 20Hz','10 2.5mW 2ms pulse @ 20Hz',...
%           '10 2mW 2ms pulse @ 20Hz','10 1mW 2ms pulse @ 20Hz',...
%           '10 0mW 2ms pulse @ 20Hz'};
% stim_per_pattern = [20 20 20 20 20 20 20];
% npulse_per_stim = [10 10 10 10 10 10 10];
% nPatterns = length(titles);
% 20220604
% titles = {'20Hz 4mW 2ms pulse (10 pulse)','20Hz 3mW 2ms pulse (10 pulse)',...
%           '20Hz 2mW 2ms pulse (10 pulse)','20Hz 1mW 2ms pulse (10 pulse)',...
%           '20Hz 4.5mW 5ms pulse (5 pulse)','20Hz 3mW 5ms pulse (5 pulse)',...
%           '20Hz 2.5mW 5ms pulse (5 pulse)','20Hz 2mW 5ms pulse (5 pulse)'};
% stim_per_pattern = [20 20 20 20 40 20 20 20];
% npulse_per_stim = [10 10 10 10 5 5 5 5];
% nPatterns = length(titles);

% 20220606
% redLaserON(501:520) = [];
% redLaserON_original = redLaserON;
% redLaserON(1:200) = redLaserON_original(101:300);
% redLaserON(201:300) = redLaserON_original(1:100);
% titles = {'20Hz 5mW 5ms pulse (5 pulse)','20Hz 4.5mW 5ms pulse (5 pulse)',...
%           '20Hz 4mW 5ms pulse (5 pulse)','20Hz 2.5mW 5ms pulse (5 pulse)',...
%           '20Hz 2mW 5ms pulse (5 pulse)','20Hz 5mW 2ms pulse (5 pulse)',...
%           '20Hz 4.5mW 2ms pulse (5 pulse)','20Hz 0mW 2ms pulse (5 pulse)'};
% stim_per_pattern = [40 60 40 40 40 40 40 40];
% npulse_per_stim = [5 5 5 5 5 5 5 5];
% nPatterns = length(titles);

%20220908
% titles = {'20Hz 1mW 5ms pulse (20 pulse)','20Hz 2mW 5ms pulse (20 pulse)',...
%           '20Hz 5mW 5ms pulse (20 pulse)','20Hz 8mW 5ms pulse (20 pulse)',...
%           '20Hz 0mW 5ms pulse (20 pulse)'};
% stim_per_pattern = [20 20 20 20 20];
% npulse_per_stim = [20 20 20 20 20];
% nPatterns = length(titles);

%% Set up analysis parameters
%redLaserON = find(redLaser);
timeRange = [-1,5]; % in sec, the first number is considered as baseline period
timebin = 1/nidq.Fs;
reactionTime = 1;

% Package into one struct
params.timeRange = timeRange;
params.timebin = timebin;
    
%% Preprocess photometry data

photometry_raw = photometry;
rollingSize = 60; % rolling window in sec

% Rolling z score to detrend (60s window)
rollingmean = movmean(photometry_raw,rollingSize*nidq.Fs);
rollingstd = movstd(photometry_raw,rollingSize*nidq.Fs);
photometry_detrended = (photometry_raw - rollingmean)./rollingstd;

% Downsample to 100Hz
downsample_freq = 100; nSampPerBin = (1/downsample_freq)*nidq.Fs;
photometry_100 = downsamplePhotometry(photometry_detrended,nSampPerBin);

% Rolling z score (60s window)
rollingmean = movmean(photometry_100,rollingSize*100);
rollingstd = movstd(photometry_100,rollingSize*100);
photometry_zscore = (photometry_100 - rollingmean)./rollingstd;

% (optional?) Downsample to 50Hz
downsample_freq = 50; nSampPerBin = 1/downsample_freq*100;
photometry_50 = downsamplePhotometry(photometry_zscore,nSampPerBin);

skewness = skewness(photometry_50);
kurtosis = kurtosis(photometry_50);
disp('Preprocessing finished');

%% (Summary) Plot pre-processing steps
preprocess_summary_fig = figure('Position', get(0, 'Screensize'));

% Raw photometry
subplot(3,2,1);
plot(photometry_raw);
xlabel(['Time (',num2str(1000/nidq.Fs),' ms)']); ylabel('Signal (V)');
title('Raw photometry');

% % After downsample
% subplot(3,2,2);
% plot(photometry_200);
% xlabel(['Time (',num2str(1000/200),' ms)']); ylabel('Signal (a.u.)');
% title('Downsampled to 200Hz');

% After detrend (1min)
subplot(3,2,2);
plot(photometry_detrended); hold on
xlabel(['Time (',num2str(1000/nidq.Fs),' ms)']); ylabel('z-score');
title('After detrend (60s window)');

% After downsample to 100Hz
subplot(3,2,3);
plot(photometry_100);
xlabel(['Time (',num2str(1000/100),' ms)']); ylabel('z-score');
title('Downsampled to 100Hz');

% After rolling z score (1min)
subplot(3,2,4);
plot(photometry_zscore); hold on
xlabel(['Time (',num2str(1000/100),' ms)']); ylabel('z-score');
title('After z-score (60s window)');

% After downsample to 50Hz
subplot(3,2,5);
plot(photometry_50);
xlabel(['Time (',num2str(1000/50),' ms)']); ylabel('z-score');
title('Downsampled to 50Hz');

% Create the uitable
subplot(3,2,6);
histogram(normrnd(0,1,size(photometry_50)),200); hold on
histogram(photometry_50,200); hold on
xlabel('z-score'); ylabel('Count'); legend({'Normal distribution','Photometry'});
dim = [0.8205 0.001 0.25 0.27];
str = {strcat("Skewness: ",num2str(skewness)),strcat("Kurtosis: ",num2str(kurtosis))};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
title('Histogram of z-scored photometry');

% Save figure
saveas(preprocess_summary_fig,strcat(session.path,sessionName,'\summary_preprocess_',sessionName,'.png'));
return

%% (Without stim) Plot raw and downsampled data

figure;
subplot(1,2,1);
eventTime = find(allTrials == 1); 
eventInSec = eventTime/nidq.Fs; binSize = 1/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
t = timeRange(1):binSize:timeRange(2);
plotCI(t,traces,colors(1)); hold on
eventTime = find(allTrials == 2); binSize = 1/nidq.Fs;
eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,colors(2));
xlabel('Time (s)'); ylabel('V');
xline(0,'-','Tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
legend({'Go trials','No-go trials'});

subplot(1,2,2);
eventTime = find(allTrials == 1); 
eventInSec = eventTime/nidq.Fs; binSize = downsample_bin;
traces = getTraces(eventInSec,photometry_downsample,timeRange,binSize);
t = timeRange(1):binSize:timeRange(2);
plotCI(t,traces,colors(1)); hold on
eventTime = find(allTrials == 2); binSize = downsample_bin;
eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry_downsample,timeRange,binSize);
plotCI(t,traces,colors(2));
xlabel('Time (s)'); ylabel('\DeltaF/F');
xline(0,'-','Tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
legend({'Go trials','No-go trials'});

%% (With red stim) Take out optical crosstalk from the laser, 0-5ms after laser stim onset
photometry_raw = photometry;
for i=1:length(redLaserON)
    idx_start=round(redLaserON(i)+0.001*nidq.Fs);
    idx_end=round(redLaserON(i)+0.004*nidq.Fs);
    
    photometry(idx_start:idx_end) = mean(photometry_raw([idx_start,idx_end]));
end

figure(10);
plot(photometry_raw(redLaserON(1):redLaserON(1)+0.005*nidq.Fs)); hold on
plot(photometry(redLaserON(1):redLaserON(1)+0.005*nidq.Fs));
legend({'Raw','Removed crosstalk'});


%% (With red stim) Find stimulation event time

nSampPerBinNI = length(allTrials)/nDownsampleBins;

% Find the first pulse of each stimulation of each pattern
stimOnset = NaN(nPatterns,max(stim_per_pattern));
offset = 0;
for i = 1:nPatterns
    for j = 1:stim_per_pattern(i)
        % disp(num2str((i-1)*npulse_per_stim+1));
        % redStimOnset(i) = redLaserON(1600+(i-1)*npulse_per_stim+1);
        stimOnset(i,j) = redLaserON(offset+(j-1)*npulse_per_stim(i)+1);
    end
    offset = offset + npulse_per_stim(i)*stim_per_pattern(i);
end

%% (Go/Nogo) Find behavior event subtypes (hit, miss, FA, CR)

% Set cutoff
%cutoff = max(max(find(allTrials==1,210),find(allTrials==2,210)));
cutoff = length(allTrials);
allTrials = allTrials(1:cutoff);

% Find tones
leftToneTime = find(allTrials==1);
rightToneTime = find(allTrials==2);
reactionTimeSamp = reactionTime * nidq.Fs;
rightSolenoidON = find(rightSolenoid);
rightLickON = find(rightLick);

% Initialize trial subtypes
hit = []; miss = [];
false_alarm = []; correct_reject = [];
reward_omission = []; free_reward = []; free_hit = [];

% temp data structs
hitLicks = {};

allTrialsTime = find(allTrials);
for i=1:length(allTrialsTime)
    cur_cue = allTrialsTime(i);
    if i == length(allTrialsTime); next_cue = length(allTrials);
    else; next_cue = allTrialsTime(i+1); end

    % Left cue: hit/miss
    if allTrials(cur_cue) == 1
        rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue & rightSolenoidON<next_cue);
        isReward = ~isempty(rewardTime);
        trialLicksTime = find(rightLickON > cur_cue & rightLickON < next_cue);
        trialLicks = rightLickON(trialLicksTime) - cur_cue;
        firstLick = find(rightLickON > cur_cue,3);
        if (rewardTime-cur_cue) < 0.1*nidq.Fs; free_hit = [free_hit; cur_cue];
        elseif ~isempty(firstLick) && rightLickON(firstLick(end)) <= cur_cue + reactionTimeSamp
            if isReward; hit = [hit;cur_cue]; hitLicks{end+1} = trialLicks;
            else; reward_omission = [reward_omission; cur_cue]; end
        else
            miss = [miss;cur_cue];
        end
    elseif allTrials(cur_cue) == 2
        %isReward = logical(sum(find(rightSolenoidON >= cur_cue & rightSolenoidON < next_cue)));
        rewardTime = rightSolenoid(rightSolenoidON>=cur_cue & rightSolenoidON<next_cue);
        isReward = ~isempty(rewardTime);
        allLicks = find(rightLickON > cur_cue & rightLickON < next_cue);
        firstLick = find(rightLickON > cur_cue,3);
        if ~isempty(firstLick) && rightLickON(firstLick(end)) <= cur_cue + reactionTimeSamp
            false_alarm = [false_alarm;cur_cue];
        else
            if isReward; free_reward = [free_reward;cur_cue];
            else; correct_reject = [correct_reject;cur_cue]; end
        end
    end
end

% Sanity check
disp(['Hit = ',num2str(length(hit))]);
disp(['Miss = ',num2str(length(miss))]);
disp(['FA = ',num2str(length(false_alarm))]);
disp(['CR = ',num2str(length(correct_reject))]);
disp(['Reward omission = ',num2str(length(reward_omission))]);
disp(['Free reward = ',num2str(length(free_reward))]);
disp(['Free hit = ',num2str(length(free_hit))]);


%% (Anything) Find lick subtypes

% Find licks
leftLickTime = find(leftLick==1);
rightLickTime = find(rightLick==1);

% Find lick subtypes
% generate a list of timestamps where corresponding events happened

thresholdSamp = 0.2 * nidq.Fs;
reactionTimeSamp = reactionTime * nidq.Fs;

allTrialsTime = find(allTrials);
choice_lick = []; % First lick after the cue
rewarded_lick = []; % Second lick in solenoid on session
unrewarded_lick = []; % Second lick in solenoid off session
spontaneous_lick = []; % Spontaneous licks after bout
omitted_lick = [];
miss_trial = []; % Trials where animal have choice lick

correct_lick = [];
incorrect_lick = [];

total_analyzed = 0;
for i = 1:length(allTrialsTime)
    % Sanity check
    %total_analyzed = total_analyzed + length(correct_lick)+length(incorrect_lick);
    %disp('----------')
    %disp(i-1);
    %disp(num2str(total_analyzed));
    %disp(num2str(length(rewarded_lick)+length(unrewarded_lick)+...
    %    length(choice_lick)+length(spontaneous_lick)));
    %disp('----------');
    
    cur_cue = allTrialsTime(i);
    
    % Licks before first cue as spontaneous licks
    if i == 1
        spon_l = leftLickTime((1<leftLickTime)&(leftLickTime<cur_cue));
        spon_r = rightLickTime((1<rightLickTime)&(rightLickTime<cur_cue));
        spontaneous_lick = [spontaneous_lick,spon_l,spon_r];
    end
    
    % Left cue
    if allTrials(cur_cue) == 1
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick for the current trial
        correct_lick = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        incorrect_lick = rightLickTime((cur_cue<rightLickTime)&(rightLickTime<next_cue));
        
        % Assign licks to subtypes
        if isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = incorrect_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                unrewarded_lick = [unrewarded_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && isempty(incorrect_lick)
            choice_lick_time = correct_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                rewarded_lick = [rewarded_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = min(correct_lick(1),incorrect_lick(1));
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                if correct_lick(1) < incorrect_lick(1)
                    rewarded_lick = [rewarded_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                else
                    unrewarded_lick = [unrewarded_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                end
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        else; miss_trial = [miss_trial, cur_cue]; 
        end

    % Right cue
    elseif allTrials(cur_cue) == 2
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick for the current trial
        incorrect_lick = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        correct_lick = rightLickTime((cur_cue<rightLickTime)&(rightLickTime<next_cue));
        
        % Assign licks to subtypes
        if isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = incorrect_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign lick
                unrewarded_lick = [unrewarded_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && isempty(incorrect_lick)
            choice_lick_time = correct_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                rewarded_lick = [rewarded_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = min(correct_lick(1),incorrect_lick(1));
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                if correct_lick(1) < incorrect_lick(1)
                    rewarded_lick = [rewarded_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                else
                    unrewarded_lick = [unrewarded_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                end
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        else
            miss_trial = [miss_trial, cur_cue];
        end
    
    % Left omission
    elseif allTrials(cur_cue) == 3
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick for the current trial
        correct_lick = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        incorrect_lick = rightLickTime((cur_cue<rightLickTime)&(rightLickTime<next_cue));
        
        % Assign licks to subtypes
        if isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = incorrect_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign lick
                unrewarded_lick = [unrewarded_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && isempty(incorrect_lick)
            choice_lick_time = correct_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                omitted_lick = [omitted_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = min(correct_lick(1),incorrect_lick(1));
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                if correct_lick(1) < incorrect_lick(1)
                    omitted_lick = [omitted_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                else
                    unrewarded_lick = [unrewarded_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                end
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        else; miss_trial = [miss_trial, cur_cue]; 
        end
    
    % Right omission trial
    elseif allTrials(cur_cue) == 4
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick for the current trial
        incorrect_lick = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        correct_lick = rightLickTime((cur_cue<rightLickTime)&(rightLickTime<next_cue));
        
        % Assign licks to subtypes
        if isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = incorrect_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign lick
                unrewarded_lick = [unrewarded_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && isempty(incorrect_lick)
            choice_lick_time = correct_lick(1);
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                omitted_lick = [omitted_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        elseif ~isempty(correct_lick) && ~isempty(incorrect_lick)
            choice_lick_time = min(correct_lick(1),incorrect_lick(1));
            if choice_lick_time < cur_cue + reactionTimeSamp
                choice_lick = [choice_lick, choice_lick_time];
                [bout,spontaneous] = findLickBout(thresholdSamp,choice_lick_time,correct_lick,incorrect_lick);
                % Assign licks
                if correct_lick(1) < incorrect_lick(1)
                    omitted_lick = [omitted_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                else
                    unrewarded_lick = [unrewarded_lick,bout];
                    spontaneous_lick = [spontaneous_lick,spontaneous];
                end
            else
                % All remaining licks are spontaneous since they're outside
                % reaction time and no reward
                spontaneous_lick = [spontaneous_lick,correct_lick,incorrect_lick];
                continue;
            end
        else; miss_trial = [miss_trial, cur_cue]; 
        end
        
    else
        continue
    end
    
end

% Sanity check
disp(num2str(length(leftLickTime)+length(rightLickTime)));
disp(num2str(length(rewarded_lick)+length(unrewarded_lick)+...
             length(choice_lick)+length(spontaneous_lick)));

%% (With red stim) Plot downsampled photometry trace vs stimulation

timestamp = -baseline_duration:downsample_bin:duration-baseline_duration;
traces_stim = NaN(max(stim_per_pattern),length(timestamp),nPatterns);

for i = 1:nPatterns
    eventTimeInSec = stimOnset(i,:)/nidq.Fs;
    traces_stim(:,:,i) = getTraces(eventTimeInSec,timestamp,photometry_downsample,params);
end

% Plot mean response overlaied
figure(16);
for i = 1:nPatterns
    traces = traces_stim(:,:,i);
    plotCI(timestamp,traces,blueRedYellow(i)); hold on
end
xlabel('Time (s)');
ylabel('\DeltaF/F');
legend(titles(1:nPatterns));
box off


%% (Normal Go/Nogo) Plot daily summary plot

daily_summary_fig = figure('Position', get(0,'Screensize'));
photometry = photometry_50; binSize = 1/50;

% Go vs nogo
subplot(2,2,1); 
eventTime = find(allTrials == 1); eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
t = timeRange(1):binSize:timeRange(2);
plotCI(t,traces,colors(1)); hold on

eventTime = find(allTrials == 2); eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,colors(2));

eventTime = find(airpuff); eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,colors(3));

xlabel('Time (s)'); ylabel('z-score');
xline(0,'-','Tone/Airpuff','Color','r','LineWidth',1.5,'HandleVisibility','off');
legend({['Go trials (n=',num2str(length(hit)+length(miss)),')'],...
    ['No-go trials (n=',num2str(length(false_alarm)+length(correct_reject)),')'],...
    ['Airpuffs (n=',num2str(length(eventTime)),')']},...
    'Location','northwest'); 
box off

% Hit/miss/FA/CR
subplot(2,2,2); 
eventTime = hit; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
t = timeRange(1):binSize:timeRange(2);
plotCI(t,traces,blueWhiteRed(1,:));

eventTime = miss; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueWhiteRed(150,:));

eventTime = false_alarm; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueWhiteRed(350,:));

eventTime = correct_reject; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueWhiteRed(500,:));

eventTime = reward_omission; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueGreenYellow(12));

eventTime = free_reward; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueGreenYellow(6));

xlabel('Time (s)'); ylabel('z-score'); %ylim([-3,4]);
xline(0,'-','Tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
legend({['Hit (n=',num2str(length(hit)),')'],...
    ['Miss (n=',num2str(length(miss)),')'],...
    ['FA (n=',num2str(length(false_alarm)),')'],...
    ['CR (n=',num2str(length(correct_reject)),')'],...
    ['Reward omission (n=',num2str(length(reward_omission)),')'],...
    ['Free reward (n=',num2str(length(free_reward)),')']},...
    'Location','northwest'); 
box off

% Lick rasters for trial subtypes
subplot(2,2,3); 
eventIdx = find(allTrials == 1); drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI);
xline(0,'-','Go tone','Color','r','LineWidth',1.5,'HandleVisibility','off'); 
xlabel('Time (s)'); ylabel('Trials'); box off

subplot(2,2,4); 
eventIdx = find(allTrials == 2); drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI);
xline(0,'-','No-go tone','Color','r','LineWidth',1.5,'HandleVisibility','off'); 
xlabel('Time (s)'); ylabel('Trials'); box off


% Save figure
saveas(daily_summary_fig,strcat(session.path,sessionName,'\summary_photometry_airpuff_',sessionName,'.png'));
% saveas(daily_summary_fig,strcat(session.path,sessionName,'\summary_photometry_',sessionName,'.png'));

%% (Catch) Plot catch progress by trials

trial_summary_fig = figure('Position', get(0,'Screensize'));
photometry = photometry_50; binSize = 1/50;
timeRange = [-1,2];
t = timeRange(1):binSize:timeRange(2);

% Plot go trials
subplot(3,2,1); nTrialsPerGroup = 10;
eventTime = find(allTrials == 1); eventInSec = eventTime/nidq.Fs;
go_traces = getTraces(eventInSec,photometry,timeRange,binSize);
nTraces = length(1:nTrialsPerGroup:size(go_traces,1));
paletteLoc = floor(linspace(1,500,nTraces)); counter = 0;
for i = 1:nTrialsPerGroup:size(go_traces,1)
    counter = counter + 1;
    trialGroup = go_traces(i:min(i+nTrialsPerGroup,size(go_traces,1)),:);
    plotCI(t,trialGroup,blueGreenPurple(paletteLoc(counter),:));
    % legendLabel{counter} = ['Trials ',num2str(i),'-',num2str(min(i+nTrialsPerGroup,size(go_traces,1)))];
end
xlabel('Time (s)'); ylabel('z-score'); %legend(legendLabel); 
xline(0,'-','Go tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
box off; 


% Plot nogo trials
subplot(3,2,2); nTrialsPerGroup = 10;
eventTime = find(allTrials == 2); eventInSec = eventTime/nidq.Fs;
nogo_traces = getTraces(eventInSec,photometry,timeRange,binSize);
nTraces = length(1:nTrialsPerGroup:size(nogo_traces,1));
paletteLoc = floor(linspace(1,500,nTraces)); counter = 0;
for i = 1:nTrialsPerGroup:size(nogo_traces,1)
    counter = counter + 1;
    trialGroup = nogo_traces(i:min(i+nTrialsPerGroup,size(nogo_traces,1)),:);
    plotCI(t,trialGroup,blueGreenPurple(paletteLoc(counter),:));
    % legendLabel{counter} = ['Trials ',num2str(i),'-',num2str(min(i+nTrialsPerGroup,size(go_traces,1)))];
end
xlabel('Time (s)'); ylabel('z-score'); %legend(legendLabel); 
xline(0,'-','No-go tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
box off


% Plot reward omission trials
subplot(3,2,3); nTrialsPerGroup = 10;
eventTime = reward_omission; eventInSec = eventTime/nidq.Fs;
ro_traces = getTraces(eventInSec,photometry,timeRange,binSize);
nTraces = length(1:nTrialsPerGroup:size(ro_traces,1));
paletteLoc = floor(linspace(1,500,nTraces)); counter = 0;
for i = 1:nTrialsPerGroup:size(ro_traces,1)
    counter = counter + 1;
    trialGroup = ro_traces(i:min(i+nTrialsPerGroup,size(ro_traces,1)),:);
    plotCI(t,trialGroup,blueGreenPurple(paletteLoc(counter),:));
    legendLabel{counter} = ['Trials ',num2str(i),'-',num2str(min(i+nTrialsPerGroup,size(go_traces,1)))];
end
xlabel('Time (s)'); ylabel('z-score'); legend(legendLabel,"Location",'northwest'); 
xline(0,'-','Go tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
title('Reward omission trials'); box off


% Plot free reward trials
subplot(3,2,4); nTrialsPerGroup = 10;
eventTime = free_reward; eventInSec = eventTime/nidq.Fs;
fr_traces = getTraces(eventInSec,photometry,timeRange,binSize);
nTraces = length(1:nTrialsPerGroup:size(fr_traces,1));
paletteLoc = floor(linspace(1,500,nTraces)); counter = 0;
for i = 1:nTrialsPerGroup:size(fr_traces,1)
    counter = counter + 1;
    trialGroup = fr_traces(i:min(i+nTrialsPerGroup,size(fr_traces,1)),:);
    plotCI(t,trialGroup,blueGreenPurple(paletteLoc(counter),:));
    legendLabel{counter} = ['Trials ',num2str(i),'-',num2str(min(i+nTrialsPerGroup,size(go_traces,1)))];
end
xlabel('Time (s)'); ylabel('z-score'); legend(legendLabel,"Location",'northwest'); 
xline(0,'-','No-go tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
title('Free reward trials'); box off


% Plot reward omission lick raster
subplot(3,2,5); 
eventIdx = reward_omission; drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI,leftSolenoid,rightSolenoid,airpuff);
xline(0,'-','Go tone','Color','r','LineWidth',1.5,'HandleVisibility','off'); 
xlabel('Time (s)'); ylabel('Trials'); title('Reward omission trials'); box off

% Plot free reward lick raster
subplot(3,2,6); 
eventIdx = free_reward; drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI,leftSolenoid,rightSolenoid,airpuff);
xline(0,'-','No-go tone','Color','r','LineWidth',1.5,'HandleVisibility','off'); 
xlabel('Time (s)'); ylabel('Trials'); title('Free reward trials'); box off

% Save figure
saveas(trial_summary_fig,strcat(session.path,sessionName,'\summary_catch_',sessionName,'.png'));

%% (Switch Go/Nogo) Plot daily summary fig

daily_summary_fig = figure('Position', get(0,'Screensize'));
photometry = photometry_50; binSize = 1/50;
timeRange = [-1,5];
t = timeRange(1):binSize:timeRange(2);

% Go vs nogo
subplot(3,2,1); 
eventTime = find(allTrials == 1); eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,colors(1)); hold on

eventTime = find(allTrials == 2); eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,colors(2));

eventTime = find(airpuff); eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,colors(3));

xlabel('Time (s)'); ylabel('z-score');
xline(0,'-','Tone/Airpuff','Color','r','LineWidth',1.5,'HandleVisibility','off');
legend({['Go trials (n=',num2str(length(hit)+length(miss)),')'],...
    ['No-go trials (n=',num2str(length(false_alarm)+length(correct_reject)),')'],...
    ['Airpuffs (n=',num2str(length(eventTime)),')']},...
    'Location','northeast'); 
box off

% Hit/miss/FA/CR
subplot(3,2,2); 
eventTime = hit; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
t = timeRange(1):binSize:timeRange(2);
plotCI(t,traces,blueWhiteRed(1,:));

eventTime = miss; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueWhiteRed(150,:));

eventTime = false_alarm; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueWhiteRed(350,:));

eventTime = correct_reject; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueWhiteRed(500,:));

eventTime = free_hit; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,blueGreenYellow(6));

xlabel('Time (s)'); ylabel('z-score'); ylim([-3,3]);
xline(0,'-','Tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
legend({['Hit (n=',num2str(length(hit)),')'],...
    ['Miss (n=',num2str(length(miss)),')'],...
    ['FA (n=',num2str(length(false_alarm)),')'],...
    ['CR (n=',num2str(length(correct_reject)),')'],...
    ['Free reward (n=',num2str(length(free_hit)),')']},...
    'Location','northeast'); 
box off

% Plot go trials
subplot(3,2,3);
eventTime = find(allTrials == 1); eventInSec = eventTime/nidq.Fs;
go_traces = getTraces(eventInSec,photometry,timeRange,binSize);
nTrialsPerGroup = 10; totalTrialsPlotted = 50;%size(go_traces,1);
nTraces = length(1:nTrialsPerGroup:totalTrialsPlotted);
paletteLoc = floor(linspace(1,500,nTraces)); counter = 0;
for i = 1:nTrialsPerGroup:totalTrialsPlotted
    counter = counter + 1;
    trialGroup = go_traces(i:min(i+nTrialsPerGroup,size(go_traces,1)),:);
    plotCI(t,trialGroup,blueGreenPurple(paletteLoc(counter),:));
    legendLabel{counter} = ['Trials ',num2str(i),'-',num2str(min(i+nTrialsPerGroup,size(go_traces,1)))];
end
xlabel('Time (s)'); ylabel('z-score'); legend(legendLabel); 
xline(0,'-','Go tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
box off; 


% Plot nogo trials
subplot(3,2,4);
eventTime = find(allTrials == 2); eventInSec = eventTime/nidq.Fs;
nogo_traces = getTraces(eventInSec,photometry,timeRange,binSize);
nTrialsPerGroup = 10; totalTrialsPlotted = 50;%size(nogo_traces,1);
nTraces = length(1:nTrialsPerGroup:totalTrialsPlotted);
paletteLoc = floor(linspace(1,500,nTraces)); counter = 0;
for i = 1:nTrialsPerGroup:totalTrialsPlotted
    counter = counter + 1;
    trialGroup = nogo_traces(i:min(i+nTrialsPerGroup,size(nogo_traces,1)),:);
    plotCI(t,trialGroup,blueGreenPurple(paletteLoc(counter),:));
    legendLabel{counter} = ['Trials ',num2str(i),'-',num2str(min(i+nTrialsPerGroup,size(nogo_traces,1)))];
end
xlabel('Time (s)'); ylabel('z-score'); legend(legendLabel); ylim([-3,4]);
xline(0,'-','No-go tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
box off

% Lick rasters for trial subtypes
subplot(3,2,5); 
eventIdx = find(allTrials == 1); drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI,leftSolenoid,rightSolenoid,airpuff);
xline(0,'-','Go tone','Color','r','LineWidth',1.5,'HandleVisibility','off'); 
xlabel('Time (s)'); ylabel('Trials'); box off

subplot(3,2,6); 
eventIdx = find(allTrials == 2); drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI,leftSolenoid,rightSolenoid,airpuff);
xline(0,'-','No-go tone','Color','r','LineWidth',1.5,'HandleVisibility','off'); 
xlabel('Time (s)'); ylabel('Trials'); box off

% Save figure
saveas(daily_summary_fig,strcat(session.path,sessionName,'\summary_photometry_airpuff_',sessionName,'.png'));

%% Second bump: does the magnitude of second peak correlates with end of lick bout?

photometry = photometry_50; binSize = 1/50;
eventTime = hit; eventInSec = eventTime/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
t = timeRange(1):binSize:timeRange(2);

%correlation_summary_fig = figure('Position', get(0,'Screensize'));

% Plot single trial traces
subplot(1,3,1);
i = 8;
[p,S] = polyfit(t,traces(i,:),10); yfit = polyval(p,t);
ysmooth = smooth(traces(i,:),10);

secondPeakCutoff = -timeRange(1) + 1/binSize; % second bump happens after 1sec
[max1,max_fit] = max(yfit(secondPeakCutoff:end));
[max2,max_smooth] = max(ysmooth(secondPeakCutoff:end));

plot(t,traces(i,:)); hold on
plot(t,yfit); hold on
plot(t,ysmooth); hold on
scatter(max_fit*binSize,max1); hold on
scatter(max_smooth*binSize,max2); hold on
legend(["Raw trace","10 degree polyfit","Smooth to 5Hz"],"Location","northwest"); box off


% Find peak of second bump for hit trials
peakLocFit = []; peakLocSmooth = []; peakValFit = []; peakValSmooth = [];
lastLicks = [];
% Get peak of second bump of hit trials
for i = 1:size(traces,1)
    [p,S] = polyfit(t,traces(i,:),10); yfit = polyval(p,t);
    ysmooth = smooth(traces(i,:),10);

    secondPeakCutoff = (1-timeRange(1))/binSize -1; % second bump happens after 1sec
    [maxFit,maxFitLoc] = max(yfit(secondPeakCutoff:end));
    [maxSmooth,maxSmoothLoc] = max(ysmooth(secondPeakCutoff:end));

    peakValFit = [peakValFit; maxFit];
    peakValSmooth = [peakValSmooth; maxSmooth];
    peakLocFit = [peakLocFit;maxFitLoc + secondPeakCutoff];
    peakLocSmooth = [peakLocSmooth;maxSmoothLoc + secondPeakCutoff];

end

peakLocFitTime = peakLocFit .* binSize;
peakLocSmoothTime = peakLocSmooth .* binSize;

% Get last lick of hit trials
for i = 1:size(hitLicks,2)
    lastLicks = [lastLicks; hitLicks{i}(end)/nidq.Fs];
end

% Scatter plot
subplot(1,3,2); 
scatter(lastLicks,peakLocSmoothTime-1); hold on
%scatter(lastLicks,peakLocFitTime); hold on
[fit1,S1]= polyfit(lastLicks,peakLocSmoothTime-1,1); plot(lastLicks,polyval(fit1,lastLicks)); hold on
xlabel("Last lick time (s)"); ylabel("Second peak time (s)"); box off

subplot(1,3,3); 
scatter(lastLicks,peakValSmooth-1); hold on
%scatter(lastLicks,peakLocFitTime); hold on
[fit2,S2] = polyfit(lastLicks,peakValSmooth,1); plot(lastLicks,polyval(fit2,lastLicks)); hold on
xlabel("Last lick time (s)"); ylabel("Second peak value (a.u.)");  box off


%% Plot downsample photometry vs lick subtypes
timestamp = -baseline_duration:downsample_bin:duration-baseline_duration;
leftLickInSec = leftLickTime/nidq.Fs; % event time in sec
rightLickInSec = rightLickTime/nidq.Fs; % event time in sec
choiceLickInSec = choice_lick/nidq.Fs;
rewardLickInSec = rewarded_lick/nidq.Fs;
unrewardLickInSec = unrewarded_lick/nidq.Fs;
sponLickInSec = spontaneous_lick/nidq.Fs;

traces_leftLick = getTraces(leftLickInSec,timestamp,photometry_downsample,params);
traces_rightLick = getTraces(rightLickInSec,timestamp,photometry_downsample,params);
traces_rewardLick = getTraces(rewardLickInSec,timestamp,photometry_downsample,params);
traces_unrewardLick = getTraces(unrewardLickInSec,timestamp,photometry_downsample,params);
traces_sponLick = getTraces(sponLickInSec,timestamp,photometry_downsample,params);

figure(15);
plotCI(timestamp,traces_leftLick,colors(1)); hold on
plotCI(timestamp,traces_rightLick,colors(2)); hold on
plotCI(timestamp,traces_rewardLick,colors(3)); hold on
plotCI(timestamp,traces_unrewardLick,colors(4)); hold on
plotCI(timestamp,traces_sponLick,colors(5));
% ste=std(traces,0,1)/sqrt(size(traces,1));
% confplot(timestamp(1:l),mean(temp,1),ste(1:l),ste(1:l));
xlabel('Time (s)');
ylabel('\DeltaF/F');
legend({'Left licks','Right licks',...
    'Rewarded licks','Unrewarded licks','Spontaneous licks'});  box off;
% return;

%% Plot overlay
timestamp = -baseline_duration:downsample_bin:duration-baseline_duration;

% Chrimson stim and left tone
figure; pattern_num = 7;
plotCI(timestamp,traces_leftTone,colors(1)); hold on
plotCI(timestamp,traces_stim(:,:,pattern_num),colors(2));
xlabel('Time (s)'); ylabel('\DeltaF/F');
legend({'Left cue',titles{pattern_num}});  box off

% Chrimson stim and right tone
figure; pattern_num = 7;
plotCI(timestamp,traces_rightTone,colors(1)); hold on
plotCI(timestamp,traces_stim(:,:,pattern_num),colors(2));
xlabel('Time (s)'); ylabel('\DeltaF/F');
legend({'Right cue',titles{pattern_num}});  box off

%% Plot downsample photometry vs some event

timestamp = -baseline_duration:downsample_bin:duration-baseline_duration;
rightSolenoidON = find(rightSolenoid);
eventTime = rightSolenoidON([2,3,4]);

traces = getTraces(eventTime/nidq.Fs,timestamp,photometry_downsample,params);

figure(14);
plotCI(timestamp,traces,colors(3));
xlabel('Time (s)');
ylabel('\DeltaF/F');
% legend({'Left cue','Right cue','Miss trials'});  
box off;
% return;

%% (With stim) Plot raw photometry trace
timestamp = -baseline_duration:timebin:duration-baseline_duration;
eventTime = redLaserON;

counter = 0;
for i = 1:nPatterns
    cmap = hsv(stim_per_pattern(i));
    
    for j = 1:stim_per_pattern(i)
        % Ignore the last pulse of the last pattern
        if i == nPatterns && j == stim_per_pattern(i); continue; end
        
        counter = counter+1;
        if counter <= length(eventTime)    
            firstSamp = eventTime(counter)-floor(baseline_duration*nidq.Fs);
            lastSamp = eventTime(counter)+floor((duration-baseline_duration)*nidq.Fs);
        
            baseline = mean(photometry(firstSamp:eventTime(counter)-1));
            dff = (photometry(firstSamp:lastSamp) - baseline)/baseline;

            l = min(length(dff),length(timestamp));
            figure(i+1);
            plot(timestamp(1:l),dff(1:l),'color',cmap(j,:));
            xlabel('Time (s)');
            ylabel('\DeltaF/F');
            hold on;
        end
    end
    autoArrangeFigures();
end

%% (Without stim) Plot raw photometry trace

figure(1);
plot(photometry); hold on
xlabel('Time (s)');
ylabel('V');

figure(2);
eventTime = find(allTrials == 1); 
eventInSec = eventTime/nidq.Fs; binSize = 1/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
t = timeRange(1):binSize:timeRange(2);
plotCI(t,traces,colors(1)); hold on

eventTime = find(allTrials == 2); 
eventInSec = eventTime/nidq.Fs; binSize = 1/nidq.Fs;
traces = getTraces(eventInSec,photometry,timeRange,binSize);
plotCI(t,traces,colors(2));
xlabel('Time (s)'); ylabel('\DeltaF/F');
xline(0,'-','Tone','Color','r','LineWidth',1.5,'HandleVisibility','off');
legend({'Go trials','No-go trials'});

autoArrangeFigures();

%% Test photometry trace

fig = initializeFig(0.5,0.33);
eventTime = find(rightSolenoid); 
eventInSec = eventTime/nidq.Fs; binSize = 1/50;
[traces,t] = getTraces(eventInSec,photometry_50,timeRange,binSize);
[dff,t] = getdff(eventInSec,photometry_50,timeRange,binSize);
% plotCI(t,traces,colors(1)); hold on
plotCI(t,dff,colors(2)); hold on
plotEvent('Water',0,blueWhiteRed(end,:));
xlabel('Time (s)'); ylabel('z-score');

saveas(fig,strcat(session.path,sessionName,'\psth_Water_',sessionName,'.png'));