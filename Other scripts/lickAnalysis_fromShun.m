%% For Paolo from Shun

clear; close all;

%% Set up
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
[~,colors,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed] = loadColors;

resultspath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Paolo\Licking_Beh\ToRaster';

%% Read matlab file
% move file to workspace

%% Extract time points

leftLick = RawData{RawData.Events == 'LickLeft',"Times"};
rightLick = RawData{RawData.Events == 'LickRight',"Times"};

incorrectLick = RawData{RawData.Events == 'IncorrectLick',"Times"};
missTrials = RawData{RawData.Events == 'MissTrial',"Times"};

ITI = RawData{RawData.Events == 'ITI',"Times"};
leftTrials = RawData{RawData.Events == 'LEFTTrial',"Times"};
rightTrials = RawData{RawData.Events == 'RIGHTTrial',"Times"};

noStim = RawData{RawData.Events == 'NoStim',"Times"};
laserStim = RawData{RawData.Events == 'LaserStim',"Times"};
water = RawData{RawData.Events == 'RewardDelivery',"Times"};

%% Determine whether trial is opto stim or not

gracePeriod = 30; % in ms

% Left trials
leftTrials(:,2) = zeros(length(leftTrials),1);
for i = 1:size(leftTrials,1)
    % Find the closest stim/no stim time
    [noStimTime,~] = min(abs(noStim-leftTrials(i,1)));
    [stimTime,~] = min(abs(laserStim-leftTrials(i,1)));
    
    if stimTime <= gracePeriod && stimTime > 0
        leftTrials(i,2) = 1; % stim trials
    elseif noStimTime <= gracePeriod && noStimTime > 0
        leftTrials(i,2) = 0; % no stim trials
    else
        disp(['No appropriate stim/no stim trials! Trial number is ',num2str(i)]); 
    end
end

% Right trials
rightTrials(:,2) = zeros(length(rightTrials),1);
for i = 1:size(rightTrials,1)
    % Find the closest stim/no stim time
    [noStimTime,~] = min(abs(noStim-rightTrials(i,1)));
    [stimTime,~] = min(abs(laserStim-rightTrials(i,1)));
    
    if stimTime <= gracePeriod && stimTime > 0
        rightTrials(i,2) = 1; % stim trials
    elseif noStimTime <= gracePeriod && noStimTime > 0
        rightTrials(i,2) = 0; % no stim trials
    else
        disp(['No appropriate stim/no stim trials! Trial number is ',num2str(i)]); 
    end
end

clearvars noStimTime stimTime

%% Plot lick raster

timeRange = [-0.01,1]; markerSize = 20;
timeNI = linspace(1,max(RawData.Times),max(RawData.Times));

initializeFig(1,.67);
tiledlayout(1,2); 
noStimAlpha = 0.5;

% Plot left trial licks
nexttile; trials = leftTrials;
[~,~,leftTrialLicks] = getLicks(timeRange,trials,0.1,leftLick,rightLick,1000,timeNI,side=[1,1],inputLickIdx=true,getRate=false);
for i = 1:size(trials,1)
    % determine whether trial have stim
    if trials(i,2) == 1
        %scatter(0,trials(i,1),markerSize,'filled','MarkerFaceColor',colors{1}); hold on
        scatter(leftTrialLicks{i,1},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(1,:)); hold on
        scatter(leftTrialLicks{i,2},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(end,:)); hold on
    else
        %scatter(0,trials(i,1),markerSize,'filled','MarkerFaceColor',colors{2}); hold on
        %scatter(leftTrialLicks{i,1},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(1,:),'MarkerFaceAlpha',noStimAlpha); hold on
        %scatter(leftTrialLicks{i,2},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(end,:),'MarkerFaceAlpha',noStimAlpha); hold on
    end
end
xlabel('Time (ms)'); xlim([timeRange(1)*1000,timeRange(2)*1000]);
ylabel('Trial'); ylim([0,size(trials,1)]);
plotEvent("",0.5,'r');
% title("Left trials");

% Plot right trial licks
nexttile; trials = rightTrials;
[~,~,rightTrialLicks] = getLicks(timeRange,trials,0.1,leftLick,rightLick,1000,timeNI,side=[1,1],inputLickIdx=true,getRate=false);
for i = 1:size(trials,1)
    % determine whether trial have stim
    if trials(i,2) == 1
        %scatter(0,trials(i,1),markerSize,'filled','MarkerFaceColor',colors{1}); hold on
        scatter(rightTrialLicks{i,1},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(1,:)); hold on
        scatter(rightTrialLicks{i,2},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(end,:)); hold on
    else
        %scatter(0,trials(i,1),markerSize,'filled','MarkerFaceColor',colors{2}); hold on
        %scatter(rightTrialLicks{i,1},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(1,:),'MarkerFaceAlpha',noStimAlpha); hold on
        %scatter(rightTrialLicks{i,2},i,markerSize,'filled','MarkerFaceColor',bluePurpleRed(end,:),'MarkerFaceAlpha',noStimAlpha); hold on
    end
end
xlabel('Time (ms)'); xlim([timeRange(1)*1000,timeRange(2)*1000]);
ylabel('Trial'); ylim([0,size(trials,1)]);
plotEvent("",0.5,'r');
% title("Right trials");

% Save 
saveFigures(gcf,'Adora_DAT_6_1_coupling_crop2_onlyOpto',resultspath);



