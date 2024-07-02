% Plot photometry signal aligned to behavioral event

clear; close all;
addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed] = loadColors;
             
             
% sessionName ='20220512 SJ512-L_g0';
% sessionName = '20220518-SJ512-R-photometry only_g0';
% sessionName = '20220527-SJ512-L_g0';
% sessionName = '20220530-SJ512-L_g0';
% sessionName = '20220602-SJ513-L_g0';
% sessionName = '20220603-SJ513-L_g0';
% sessionName = '20220604-SJ513-L_g0';
% sessionName = '20220606-SJ513-R_g0';
sessionName = '20220908-SJ522-DLS_g0';
load(strcat('sync_',sessionName,'.mat'),'photometry','nidq','redLaser',...
    'blueLaser','allTrials','leftLick','rightLick','leftSolenoid','rightSolenoid');

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

%% Set up analysis parameters
redLaserON = find(redLaser);
duration = 10; % duration of plotting in sec
baseline_duration = 2; % in sec
timebin = 1/nidq.Fs;

% Stimulation trials params
%20220908
titles = {'20Hz 1mW 5ms pulse (20 pulse)','20Hz 2mW 5ms pulse (20 pulse)',...
          '20Hz 5mW 5ms pulse (20 pulse)','20Hz 8mW 5ms pulse (20 pulse)',...
          '20Hz 0mW 5ms pulse (20 pulse)'};
stim_per_pattern = [20 20 20 20 20];
npulse_per_stim = [20 20 20 20 20];
nPatterns = length(titles);

% Downsample
downsample_bin = 0.1; % photometry timebin: 10 ms
nSampPerBin = downsample_bin*nidq.Fs;
nDownsampleBins = floor(length(photometry)/(downsample_bin*nidq.Fs));
windowSize = 2; % moving average filter

% Package into one struct
params.duration = duration; params.baseline_duration = baseline_duration;
params.timebin = timebin; params.downsample_bin = downsample_bin;
params.nSampPerBin = nSampPerBin; params.nDownsampleBins = nDownsampleBins;
params.windowSize = windowSize;

%% Plot raw photometry trace
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
    
%% Take out optical crosstalk from the laser, 0-5ms after laser stim onset
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

%% Smooth photometry data

% Downsample
nSampPerBin = downsample_bin*nidq.Fs;
nDownsampleBins = floor(length(photometry)/(downsample_bin*nidq.Fs));

% Putting raw traces into a user defined timebin
photometry_downsample = zeros(1,floor(length(photometry)/(downsample_bin*nidq.Fs)));
for k = 1:length(photometry_downsample)
    firstBinTime = floor(nSampPerBin*(k-1)+1);
    lastBinTime = floor(nSampPerBin*k);
    photometry_downsample(k) = sum(photometry(floor(firstBinTime:lastBinTime)));
    %data(k) = median(photometry(floor(downsample_bin*nidq.Fs*(k-1)+1):floor(downsample_bin*nidq.Fs*k)));
end
clear nSampPerBin firstBinTime lastBinTime

% Moving average filter
% windowDuration = 0.2; % Moving window: 200 ms
% windowSize = windowDuration * nidq.Fs;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
photometry_downsample_unfiltered = photometry_downsample;
photometry_downsample = filter(b,a,photometry_downsample);

figure(11); 
plot(photometry_downsample_unfiltered); hold on
plot(photometry_downsample); hold on
xlabel('Time (s)');
ylabel('\DeltaF/F');
legend({'Downsample & filtered','Downsample & unfiltered'});

figure(12);
plot(photometry_downsample_unfiltered(25332-100:25332+100)); hold on
plot(photometry_downsample(25332-100:25332+100)); hold on
xlabel('Time (s)');
ylabel('\DeltaF/F');
legend({'Downsample & filtered','Downsample & unfiltered'});

autoArrangeFigures();

%% Find stimulation event time

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

%% Find behavior event times

% Find tones
leftToneTime = find(allTrials==1);
rightToneTime = find(allTrials==2);

% Find licks
leftLickTime = find(leftLick==1);
rightLickTime = find(rightLick==1);

% Find lick subtypes
% generate a list of timestamps where corresponding events happened

thresholdSamp = 0.2 * nidq.Fs;
reactionTimeSamp = 0.5 * nidq.Fs;

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

%% Plot downsampled photometry trace vs stimulation

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

%% Plot downsample photometry vs tone 
timestamp = -baseline_duration:downsample_bin:duration-baseline_duration;
leftToneInSec = setdiff(leftToneTime,miss_trial)/nidq.Fs; %event time in sec
rightToneInSec = setdiff(rightToneTime,miss_trial)/nidq.Fs; %event time in sec
missTrialInSec = miss_trial/nidq.Fs;

traces_leftTone = getTraces(leftToneInSec,timestamp,photometry_downsample,params);
traces_rightTone = getTraces(rightToneInSec,timestamp,photometry_downsample,params);
traces_missTrial = getTraces(missTrialInSec,timestamp,photometry_downsample,params);

figure(14);
plotCI(timestamp,traces_leftTone,colors(1)); hold on
plotCI(timestamp,traces_rightTone,colors(2)); hold on
plotCI(timestamp,traces_missTrial,colors(3));
xlabel('Time (s)');
ylabel('\DeltaF/F');
legend({'Left cue','Right cue','Miss trials'});  box off;
% return;

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