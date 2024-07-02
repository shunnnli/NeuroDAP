%% (Required) Load sync data 
clear; close all;
% addpath(genpath('/Users/shunli/Downloads/Sabatini lab/Methods'))
addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));

sessionName = '20220918-SL029-D7_g0';
session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Project valence switching\Recordings\';

load(strcat(session.path,sessionName,'\','sync_',sessionName,'.mat'));
eye = readmatrix(strcat(session.path,sessionName, '\','eye-pixel-intensity.csv'));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed] = loadColors;

timeRange = [-1,5]; lick_binSize = 0.2; blink_thresh = 2.5; % in turns of z score

disp(strcat('Session'," ", sessionName,' loaded'));

%% (Training day) Get PSTHs of everything
eventIdx = find(allTrials==1);
eyeTraces_reward = getEye(timeRange,eventIdx,eye,sync,timeNI,timeCamera);
[lickSummary_reward,lickTraces_reward] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

eventIdx = find(allTrials==2);
eyeTraces_punish = getEye(timeRange,eventIdx,eye,sync,timeNI,timeCamera);
[lickSummary_punish,lickTraces_punish] = getLicks(timeRange,eventIdx,lick_binSize,[],rightLick,nidq,timeNI);

lickRate_reward = lickTraces_reward{2}/lick_binSize;
lickRate_punish = lickTraces_punish{2}/lick_binSize;

blink_occurance_reward = getBlinks(eyeTraces_reward,blink_thresh);
blink_occurance_punish = getBlinks(eyeTraces_punish,blink_thresh);

blink_rate_reward = blink_occurance_reward / size(blink_occurance_reward,1);
blink_rate_punish = blink_occurance_punish / size(blink_occurance_punish,1);

%% (Switching day) Get PSTHs of everything

switchAfter = 75;
allTrials_time = find(allTrials > 0);
allTrials_type = allTrials(allTrials_time);
beforeSwitch = allTrials_type(1:switchAfter);
afterSwitch = allTrials_type(switchAfter+1:end);
rewardTrials = [allTrials_time(beforeSwitch == 1), allTrials_time(afterSwitch == 2)];
punishTrials = [allTrials_time(beforeSwitch == 2), allTrials_time(afterSwitch == 1)];

beforeSwitch_reward = allTrials_time(beforeSwitch == 1);
beforeSwitch_punish = allTrials_time(beforeSwitch == 2);
afterSwitch_reward = allTrials_time(afterSwitch == 2);
afterSwitch_punish = allTrials_time(afterSwitch == 1);

% All reward trials
eyeTraces_reward = getEye(timeRange,rewardTrials,eye,sync,timeNI,timeCamera);
[lickSummary_reward,lickTraces_reward] = getLicks(timeRange,rewardTrials,lick_binSize,[],rightLick,nidq,timeNI);
% Before switch reward
eyeTraces_reward_beforeSwitch = getEye(timeRange,beforeSwitch_reward,eye,sync,timeNI,timeCamera);
[lickSummary_reward_beforeSwitch,lickTraces_reward_beforeSwitch] = ...
    getLicks(timeRange,beforeSwitch_reward,lick_binSize,[],rightLick,nidq,timeNI);
% After switch reward
eyeTraces_reward_afterSwitch = getEye(timeRange,afterSwitch_reward,eye,sync,timeNI,timeCamera);
[lickSummary_reward_afterSwitch,lickTraces_reward_afterSwitch] = ...
    getLicks(timeRange,afterSwitch_reward,lick_binSize,[],rightLick,nidq,timeNI);

% All punish trials
eyeTraces_punish = getEye(timeRange,punishTrials,eye,sync,timeNI,timeCamera);
[lickSummary_punish,lickTraces_punish] = getLicks(timeRange,punishTrials,lick_binSize,[],rightLick,nidq,timeNI);
% Before switch punish
eyeTraces_punish_beforeSwitch = getEye(timeRange,beforeSwitch_punish,eye,sync,timeNI,timeCamera);
[lickSummary_punish_beforeSwitch,lickTraces_punish_beforeSwitch] = ...
    getLicks(timeRange,beforeSwitch_punish,lick_binSize,[],rightLick,nidq,timeNI);
% After switch punish
eyeTraces_punish_afterSwitch = getEye(timeRange,afterSwitch_punish,eye,sync,timeNI,timeCamera);
[lickSummary_punish_afterSwitch,lickTraces_punish_afterSwitch] = ...
    getLicks(timeRange,afterSwitch_punish,lick_binSize,[],rightLick,nidq,timeNI);

lickRate_reward = lickTraces_reward{2}/lick_binSize;
lickRate_punish = lickTraces_punish{2}/lick_binSize;
lickRate_reward_beforeSwitch = lickTraces_reward_beforeSwitch{2}/lick_binSize;
lickRate_punish_beforeSwitch = lickTraces_punish_beforeSwitch{2}/lick_binSize;
lickRate_reward_afterSwitch = lickTraces_reward_afterSwitch{2}/lick_binSize;
lickRate_punish_afterSwitch = lickTraces_punish_afterSwitch{2}/lick_binSize;

% blink_occurance_reward = getBlinks(eyeTraces_reward,blink_thresh);
% blink_occurance_punish = getBlinks(eyeTraces_punish,blink_thresh);
% blink_occurance_reward_beforeSwitch = getBlinks(eyeTraces_reward_beforeSwitch,blink_thresh);
% blink_occurance_punish_beforeSwitch = getBlinks(eyeTraces_punish_beforeSwitch,blink_thresh);
% 
% blink_rate_reward = blink_occurance_reward / size(blink_occurance_reward,1);
% blink_rate_punish = blink_occurance_punish / size(blink_occurance_punish,1);

%% Plot avg lick rates for reward vs punish

t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotCI(t,lickRate_reward,colors(1)); hold on
plotCI(t,lickRate_punish,colors(2));
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); legend(["Reward","punish"]);
xlim([timeRange(1),timeRange(2)]); drawCSUS(2);

%% Plot avg eye pixel intensity
figure; t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotCI(t,eyeTraces_reward,colors(1)); hold on
plotCI(t,eyeTraces_punish,colors(2));
drawCSUS(1);

%% Plot trials vs licks
figure; plot(smooth(lickSummary_reward(:,1))); hold on
plot(smooth(lickSummary_punish(:,1)));
box off

%% Plot trials vs avg eye pixel intensity (over 10 trials)
figure; t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
maxValue = max(max(eyeTraces_punish));
counter = 0;
for i = 1:10:size(eyeTraces_punish,1)
    counter = counter + 1;
    if i+10 < size(eyeTraces_punish,1)
        plotCI(t,eyeTraces_punish(i:i+10,:),blueGreenYellow(counter)); hold on
        %plotCI(t,counter + eyeTraces_punish(i:i+10,:)/maxValue,colors(2)); hold on
    else
        %plotCI(t,eyeTraces_punish(i:end,:),blueGreenYellow(counter)); hold on
        %plotCI(t,counter + eyeTraces_punish(i:end,:)/maxValue,colors(2)); hold on
    end
end
box off

%% Plot raster plot of licking

timeRange = [-1,5];

% Plot left-cue-triggered licks
eventIdx = find(allTrials==1);
figure; drawLicks('Left cue',timeRange,eventIdx,...
                        leftLick,rightLick,nidq,timeNI)

% Plot right-cue-triggered licks
eventIdx = find(allTrials==2);
figure; drawLicks('Right cue',timeRange,eventIdx,...
                        leftLick,rightLick,nidq,timeNI)



%% (Training day) Plot daily summary plots
daily_summary_fig = figure('Position', get(0, 'Screensize'));

subplot(2,3,1); 
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotCI(t,lickRate_reward,colors(1)); hold on
plotCI(t,lickRate_punish,colors(2));
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); legend(["Reward","punish"]); 
xlim([timeRange(1),timeRange(2)]); drawCSUS(2);

subplot(2,3,2);
% Plot left-cue-triggered licks
eventIdx = find(allTrials==1);
drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI);
drawCSUS(0);

subplot(2,3,3);
% Plot right-cue-triggered licks
eventIdx = find(allTrials==2);
drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI);
drawCSUS(1);

subplot(2,3,4);
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotCI(t,eyeTraces_reward,colors(1)); hold on
plotCI(t,eyeTraces_punish,colors(2));
xlabel("Time (s)"); ylabel("Eye pixel intensity"); legend(["Reward","punish"]);
% plotCI(t,blink_rate_reward * 100,colors(1)); hold on
% plotCI(t,blink_rate_punish * 100,colors(2));
% xlabel("Time (s)"); ylabel("Blink occurance (%)"); legend(["Reward","punish"]);
drawCSUS(1);

subplot(2,3,5);
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
counter = 0;
for i = 1:10:size(eyeTraces_punish,1)
    counter = counter + 1;
    if i+10 < size(eyeTraces_punish,1)
        plotCI(t,eyeTraces_punish(i:i+10,:),blueGreenYellow(counter)); hold on
    else
        %plotCI(t,eyeTraces_punish(i:end,:),blueGreenYellow(counter)); hold on
    end
end
xlabel("Time (s)"); ylabel("Eye pixel intensity"); drawCSUS(1);

subplot(2,3,6);
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
maxValue = max(max(eyeTraces_punish)); counter = 0;
for i = 1:10:size(eyeTraces_punish,1)
    counter = counter + 1;
    if i+10 < size(eyeTraces_punish,1)
        plotCI(t,counter + eyeTraces_punish(i:i+10,:)/maxValue,colors(2)); hold on
    else
        plotCI(t,counter + eyeTraces_punish(i:end,:)/maxValue,colors(2)); hold on
    end
end
xlabel("Time (s)"); ylabel("Trials (groups of 10)"); drawCSUS(1);

% Save figure
saveas(daily_summary_fig,strcat(session.path,sessionName,'\summary_',sessionName,'.png'));

%% (Switching day) Plot daily summary plots
daily_summary_fig = figure('Position', get(0, 'Screensize'));
nTrials_per_group = 8; nGroups = 5;

subplot(3,3,1); 
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotCI(t,lickRate_reward,colors(1)); hold on
plotCI(t,lickRate_punish,colors(2));
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); legend(["Reward","punish"]); 
xlim([timeRange(1),timeRange(2)]); drawCSUS(2);

subplot(3,3,2);
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotCI(t,lickRate_reward_beforeSwitch,colors(1)); hold on
counter = 0; 
% palatte = round(linspace(1,size(blueGreenYellow,1),5));
for i = 1:nTrials_per_group:nGroups*nTrials_per_group
    counter = counter + 1;
    if i+nTrials_per_group < size(lickRate_reward_afterSwitch,1)
        %plotCI(t,lickRate_reward_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(palatte(counter),:)); hold on
        plotCI(t,lickRate_reward_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(counter)); hold on
    else
        % plotCI(t,lickRate_reward_afterSwitch(i:end,:),blueGreenYellow(palatte(counter),:)); hold on
    end
end
xlabel("Time (s)"); ylabel("Lick rate (Hz)");
xlim([timeRange(1),timeRange(2)]); drawCSUS(0);

subplot(3,3,3);
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward{1,2},2));
plotCI(t,lickRate_punish_beforeSwitch,colors(1)); hold on
counter = 0;
% palatte = round(linspace(1,size(blueGreenYellow,1),5));
for i = 1:nTrials_per_group:nGroups*nTrials_per_group
    counter = counter + 1;
    if i+nTrials_per_group < size(lickRate_punish_afterSwitch,1)
        %plotCI(t,lickRate_punish_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(palatte(counter),:)); hold on
        plotCI(t,lickRate_punish_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(counter)); hold on
    else
        % plotCI(t,lickRate_punish_afterSwitch(i:end,:),blueGreenYellow(palatte(counter),:)); hold on
    end
end
xlabel("Time (s)"); ylabel("Lick rate (Hz)");
xlim([timeRange(1),timeRange(2)]); drawCSUS(1);


subplot(3,3,4);
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotCI(t,eyeTraces_reward,colors(1)); hold on
plotCI(t,eyeTraces_punish,colors(2));
xlabel("Time (s)"); ylabel("Eye pixel intensity"); legend(["Reward","punish"]);
% plotCI(t,blink_rate_reward * 100,colors(1)); hold on
% plotCI(t,blink_rate_punish * 100,colors(2));
% xlabel("Time (s)"); ylabel("Blink occurance (%)"); legend(["Reward","punish"]);
drawCSUS(1);

subplot(3,3,5);
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_reward,2));
plotCI(t,eyeTraces_reward_beforeSwitch,colors(1)); hold on
counter = 0;
% palatte = round(linspace(1,size(blueGreenYellow,1),5));
for i = 1:nTrials_per_group:nGroups*nTrials_per_group
    counter = counter + 1;
    if i+nTrials_per_group < size(eyeTraces_reward_afterSwitch,1)
        %plotCI(t,eyeTraces_reward_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(palatte(counter),:)); hold on
        plotCI(t,eyeTraces_reward_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(counter)); hold on
    else
        % plotCI(t,eyeTraces_reward_afterSwitch(i:end,:),blueGreenYellow(palatte(counter),:)); hold on
    end
end
xlabel("Time (s)"); ylabel("Lick rate (Hz)");
xlim([timeRange(1),timeRange(2)]); drawCSUS(0);

subplot(3,3,6);
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_punish,2));
plotCI(t,eyeTraces_punish_beforeSwitch,colors(1)); hold on
counter = 0; 
% palatte = round(linspace(1,size(blueGreenYellow,1),5));
for i = 1:nTrials_per_group:nGroups*nTrials_per_group
    counter = counter + 1;
    if i+nTrials_per_group < size(eyeTraces_punish_afterSwitch,1)
        %plotCI(t,eyeTraces_punish_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(palatte(counter),:)); hold on
        plotCI(t,eyeTraces_punish_afterSwitch(i:i+nTrials_per_group,:),blueGreenYellow(counter)); hold on
    else
        % plotCI(t,eyeTraces_punish_afterSwitch(i:end,:),blueGreenYellow(palatte(counter),:)); hold on
    end
end
xlabel("Time (s)"); ylabel("Lick rate (Hz)");
xlim([timeRange(1),timeRange(2)]); drawCSUS(1);

subplot(3,3,7);
% Plot left-cue-triggered licks
eventIdx = find(allTrials==1);
drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI);
drawCSUS(0);

subplot(3,3,8);
% Plot right-cue-triggered licks
eventIdx = find(allTrials==2);
drawLicks(timeRange,eventIdx,[],rightLick,nidq,timeNI);
drawCSUS(1);

% Save figure
saveas(daily_summary_fig,strcat(session.path,sessionName,'\summary_',sessionName,'.png'));

%% Store cross-session data

% Package summary data into struct
D7.lickSummary_reward = lickSummary_reward;
D7.lickSummary_punish = lickSummary_punish;
D7.lickTraces_reward = lickTraces_reward;
D7.lickTraces_punish = lickTraces_punish;
D7.eyeTraces_reward = eyeTraces_reward;
D7.eyeTraces_punish = eyeTraces_punish;

% Create new data
% save(strcat('summary_SL029'),'D1');
save(strcat('summary_SL029'),'D7','-append');
% save(strcat('summary_SJ519'),'D10','-append');

%% (Multi-session) Load session across animals to common matrix
clear; close all;
% addpath(genpath('/Users/shunli/Downloads/Sabatini lab/Methods'))
addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));

animals = {'SL028';'SL029';'SL011';'SJ519'};

lickTraces_reward_D1 = []; lickTraces_punish_D1 = [];
lickTraces_reward_D2 = []; lickTraces_punish_D2 = [];
lickTraces_reward_D3 = []; lickTraces_punish_D3 = [];
lickTraces_reward_D4 = []; lickTraces_punish_D4 = [];
lickTraces_reward_D5 = []; lickTraces_punish_D5 = [];

eyeTraces_reward_D1 = []; eyeTraces_punish_D1 = [];
eyeTraces_reward_D2 = []; eyeTraces_punish_D2 = [];
eyeTraces_reward_D3 = []; eyeTraces_punish_D3 = [];
eyeTraces_reward_D4 = []; eyeTraces_punish_D4 = [];
eyeTraces_reward_D5 = []; eyeTraces_punish_D5 = [];

for i = 1:length(animals)
    load(strcat('summary_',animals{i},'.mat'));
    disp(strcat('Animal'," ", animals{i},' loaded')); 

    lickTraces_reward_D1 = [lickTraces_reward_D1; D1.lickTraces_reward{1,2}];
    lickTraces_punish_D1 = [lickTraces_punish_D1; D1.lickTraces_punish{1,2}]; 
    lickTraces_reward_D2 = [lickTraces_reward_D2; D2.lickTraces_reward{1,2}];
    lickTraces_punish_D2 = [lickTraces_punish_D2; D2.lickTraces_punish{1,2}];
    lickTraces_reward_D3 = [lickTraces_reward_D3; D3.lickTraces_reward{1,2}];
    lickTraces_punish_D3 = [lickTraces_punish_D3; D3.lickTraces_punish{1,2}];
    lickTraces_reward_D4 = [lickTraces_reward_D4; D4.lickTraces_reward{1,2}];
    lickTraces_punish_D4 = [lickTraces_punish_D4; D4.lickTraces_punish{1,2}];
    %lickTraces_reward_D5 = [lickTraces_reward_D5; D5.lickTraces_reward{1,2}];
    %lickTraces_punish_D5 = [lickTraces_punish_D5; D5.lickTraces_punish{1,2}];

    eyeTraces_reward_D1 = [eyeTraces_reward_D1; D1.eyeTraces_reward];
    eyeTraces_punish_D1 = [eyeTraces_punish_D1; D1.eyeTraces_punish]; 
    eyeTraces_reward_D2 = [eyeTraces_reward_D2; D2.eyeTraces_reward];
    eyeTraces_punish_D2 = [eyeTraces_punish_D2; D2.eyeTraces_punish];
    eyeTraces_reward_D3 = [eyeTraces_reward_D3; D3.eyeTraces_reward];
    eyeTraces_punish_D3 = [eyeTraces_punish_D3; D3.eyeTraces_punish];
    eyeTraces_reward_D4 = [eyeTraces_reward_D4; D4.eyeTraces_reward];
    eyeTraces_punish_D4 = [eyeTraces_punish_D4; D4.eyeTraces_punish];
    %eyeTraces_reward_D5 = [eyeTraces_reward_D5; D5.eyeTraces_reward];
    %eyeTraces_punish_D5 = [eyeTraces_punish_D5; D5.eyeTraces_punish];
end

[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed] = loadColors;
timeRange = [-1,5]; lick_binSize = 0.2;

disp(strcat('All animals loaded')); 

%% Plot cross-session summary

session_summary_fig = figure('Position', get(0, 'Screensize'));

% Plot cross-session summary of lick rate
subplot(1,2,1);
t = linspace(timeRange(1),timeRange(2),size(lickTraces_reward_D1,2));
plotCI(t,lickTraces_reward_D1/lick_binSize,blueWhiteRed(60,:)); hold on
plotCI(t,lickTraces_reward_D2/lick_binSize,blueWhiteRed(40,:)); hold on
plotCI(t,lickTraces_reward_D3/lick_binSize,blueWhiteRed(20,:)); hold on
plotCI(t,lickTraces_reward_D4/lick_binSize,blueWhiteRed(1,:)); hold on
plotCI(t,lickTraces_punish_D1/lick_binSize,blueWhiteRed(140,:)); hold on
plotCI(t,lickTraces_punish_D2/lick_binSize,blueWhiteRed(160,:)); hold on
plotCI(t,lickTraces_punish_D3/lick_binSize,blueWhiteRed(180,:)); hold on
plotCI(t,lickTraces_punish_D4/lick_binSize,blueWhiteRed(200,:)); hold on
xlabel("Time (s)"); ylabel("Lick rate (Hz)"); 
legend(["D1-Reward","D2-Reward","D3-Reward","D4-Reward",...
    "D1-punish","D2-punish","D3-punish","D4-punish"]); 
xlim([timeRange(1),timeRange(2)]); drawCSUS(2);


% Plot cross-session summary of eye blink
subplot(1,2,2);
t = linspace(timeRange(1),timeRange(2),size(eyeTraces_reward_D1,2));
plotCI(t,eyeTraces_reward_D1,blueWhiteRed(60,:)); hold on
plotCI(t,eyeTraces_reward_D2,blueWhiteRed(40,:)); hold on
plotCI(t,eyeTraces_reward_D3,blueWhiteRed(20,:)); hold on
plotCI(t,eyeTraces_reward_D4,blueWhiteRed(1,:)); hold on
plotCI(t,eyeTraces_punish_D1,blueWhiteRed(140,:)); hold on
plotCI(t,eyeTraces_punish_D2,blueWhiteRed(160,:)); hold on
plotCI(t,eyeTraces_punish_D3,blueWhiteRed(180,:)); hold on
plotCI(t,eyeTraces_punish_D4,blueWhiteRed(200,:)); hold on
xlabel("Time (s)"); ylabel("Eye pixel intensity (a.u.)"); 
legend(["D1-Reward","D2-Reward","D3-Reward","D4-Reward",...
    "D1-punish","D2-punish","D3-punish","D4-punish"]); 
xlim([timeRange(1),timeRange(2)]); drawCSUS(2);