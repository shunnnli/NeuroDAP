function analyzeSessions_forEmily(sessionpath,options)

arguments
    sessionpath string
    options.plot logical = true % Plot session figures
end

% Shun_analyzeBehavior_optoPair
% Shun Li, 6/23/2023
% branched of from Shun_analyzeSessions_optoPair, rewrite to fit emily's
% go/nogo task

%% Load data

[~,~,~,~,blueWhiteRed,~,~] = loadColors;
             
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

% Set other trial related variables
minLicks = 3; reactionTimeSamp = 1 * params.sync.behaviorFs;

lick_binSize = 0.2; blink_thresh = 2.5; % in turns of z score

%% Generate trial table 

% Combine stim&tone to form trial start
goCue = find(leftTone); nogoCue = find(rightTone);
allTrials = sortrows([[goCue';nogoCue'], [ones(length(goCue),1); zeros(length(nogoCue),1)]],1);

disp('Ongoing: making trial table for random outcome');
% Find digital events
rightSolenoidON = find(rightSolenoid);
rightLickON = find(rightLick);
toneON = find(allTrials);

% Initialize trial outcome for easy access
hit = []; miss = []; fa = []; cr = [];

% Initialize trial table
% Selection time: time of last choice lick (before outcome)
% Reaction time: time of first lick
varTypes = {'double','double','double','logical','logical','string',...
            'double','double','double','double',...
            'double','double','double','double'};
varNames = {'TrialNumber','CueTime','go/nogo','isReward','choice','outcome',...
            'ReactionTime','OutcomeTime','LastLickTime','NextCue',...
            'nLicks','nAnticipatoryLicks','MeanEyeIntensity','PeakEyeIntensity'};

% Initialize trial subtypes
trials = table('Size',[length(allTrials) length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

% Loop over all trials
for i=1:length(allTrials)
    cur_cue = allTrials(i);
    if i == length(allTrials); next_cue = length(timeNI);
    else; next_cue = allTrials(i+1); end
    
    % Trial outcome
    rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue & rightSolenoidON<next_cue) - cur_cue;
    isReward = ~isempty(rewardTime);
    if ~isReward; rewardTime = reactionTimeSamp; end
    outcomeTime = min([reactionTimeSamp, rewardTime]);

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
    if isempty(choiceLicks); choice = 0; % nogo = [nogo;cur_cue];
    else
        if consecLicks(end) >= minLicks && allTrials(i,2) == 1
            choice = 1; % go = [go;cur_cue];
        elseif ~isempty(choiceLicks) && allTrials(i,2) == 0
            choice = 1;
        else
            choice = 0; % nogo = [nogo;cur_cue];
        end
    end
    
    % % Update outcomeTime for omission trials (no outcome but licked)
    % if ~isReward && choice ~= 0
    %     choiceLicks_idx = find(consecLicks>=minLicks,minLicks);
    %     choiceLicks = choiceLicks(1:choiceLicks_idx(end),:);
    %     outcomeTime = min(outcomeTime,choiceLicks(end,1));
    % end

    % Trial outcome
    if choice == 1 && allTrials(i,2) == 1 && isReward; hit = [hit; cur_cue]; outcome = 'H';
    elseif choice == 0 && allTrials(i,2) == 1; miss = [miss; cur_cue]; outcome = 'M';
    elseif choice == 1 && allTrials(i,2) == 0 && ~isReward; fa = [fa; cur_cue]; outcome = "FA";
    elseif choice == 0 && allTrials(i,2) == 0 && ~isReward; cr = [cr; cur_cue]; outcome = "CR";

    elseif allTrials(i,2) == 1 && rewardTime < 0.05*params.sync.behaviorFs; outcome = "FH";
    elseif allTrials(i,2) == 1 && choice == 1 && ~isReward; outcome = "C";
    elseif allTrials(i,2) == 0 && choice == 0 && isReward; outcome = "FR";
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
    trials(i,:) = {i,cur_cue,allTrials(i,2),isReward,choice,outcome,...
        reactionTime,rewardTime,lastLickTime,next_cue,...
        nLicks,nAnticipatoryLicks,mean_eye,peak_eye};
end

% Save to sync.mat
save(strcat(sessionpath,filesep,'sync_',session.name),'trials','-append');
disp('Finished: trial table saved');

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
saveas(gcf,strcat(sessionpath,filesep,'signal_summary_',session.name,'.png'));

%% Task specific params

% Select event idx
waterIdx = find(rightSolenoid);  
hitIdx = trials{strcmp('H',trials.outcome),'CueTime'};
missIdx = trials{strcmp('M',trials.outcome),'CueTime'};
faIdx = trials{strcmp('FA',trials.outcome),'CueTime'};
crIdx = trials{strcmp('CR',trials.outcome),'CueTime'};
goIdx = trials{trials.choice == 1,'CueTime'};
nogoIdx = trials{trials.choice == 0, 'CueTime'};

%% (LJ) Plot combined photometry PSTHs

if session.withPhotometry
    timeRange = [-1,3]; lick_binSize = 0.1;
    
    % 2. Plot traces
    initializeFig(1,1);
    tiledlayout(2,2);
    
    % 2.1 Plot photometry traces
    nexttile;
    [~,~] = plotTraces(waterIdx,timeRange,rollingGreenLP,blueWhiteRed(1,:),params);
    [~,~] = plotTraces(goIdx,timeRange,rollingGreenLP,[0.2, 0.2, 0.2],params);
    [~,~] = plotTraces(nogoIdx,timeRange,rollingGreenLP,[.75 .75 .75],params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score');
    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
        ['Lick trials (n=',num2str(length(goIdx)),')'],...
        ['No-lick trials (n=',num2str(length(nogoIdx)),')']};
    legend(taskLegend,'Location','northeast');

    % 2.2 Plot lick traces
    nexttile;
    plotLicks(waterIdx,timeRange,lick_binSize,blueWhiteRed(1,:),[],rightLick,params);
    plotLicks(goIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
    plotLicks(nogoIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('Licks/s');
    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
        ['Lick trials (n=',num2str(length(goIdx)),')'],...
        ['No-lick trials (n=',num2str(length(nogoIdx)),')']};
    legend(taskLegend,'Location','northeast');

    % 2.3 Plot photometry traces for trial subtype
    nexttile;
    [~,~] = plotTraces(hitIdx,timeRange,rollingGreenLP,blueWhiteRed(1,:),params);
    [~,~] = plotTraces(missIdx,timeRange,rollingGreenLP,blueWhiteRed(100,:),params);
    [~,~] = plotTraces(faIdx,timeRange,rollingGreenLP,blueWhiteRed(400,:),params);
    [~,~] = plotTraces(crIdx,timeRange,rollingGreenLP,blueWhiteRed(500,:),params);
    plotEvent('Cue',0.025,'r');
    xlabel('Time (s)'); ylabel('z-score');
    taskLegend = {['Hit trials (n=',num2str(length(hitIdx)),')'],...
        ['Miss trials (n=',num2str(length(missIdx)),')'],...
        ['False alarm trials (n=',num2str(length(faIdx)),')'],...
        ['Correct rejection trials (n=',num2str(length(crIdx)),')']};
    legend(taskLegend,'Location','northeast');

    % 2.2 Plot lick traces
    nexttile;
    plotLicks(hitIdx,timeRange,lick_binSize,blueWhiteRed(1,:),[],rightLick,params);
    plotLicks(missIdx,timeRange,lick_binSize,blueWhiteRed(100,:),[],rightLick,params);
    plotLicks(faIdx,timeRange,lick_binSize,blueWhiteRed(400,:),[],rightLick,params);
    plotLicks(crIdx,timeRange,lick_binSize,blueWhiteRed(500,:),[],rightLick,params);
    plotEvent('Cue',0.025,'r');
    xlabel('Time (s)'); ylabel('Licks/s');
    taskLegend = {['Hit trials (n=',num2str(length(hitIdx)),')'],...
        ['Miss trials (n=',num2str(length(missIdx)),')'],...
        ['False alarm trials (n=',num2str(length(faIdx)),')'],...
        ['Correct rejection trials (n=',num2str(length(crIdx)),')']};
    legend(taskLegend,'Location','northeast');
    
    
    saveas(gcf,strcat(sessionpath,filesep,'psth_lj_combined_',session.name,'.png'));
end

%% (NI) Plot combined PSTH

if session.ni_photometryON
    timeRange = [-1,3]; lick_binSize = 0.1;
    
    % 2. Plot traces
    initializeFig(1,1);
    tiledlayout(2,2);
    
    % 2.1 Plot photometry traces
    nexttile;
    [~,~] = plotTraces(waterIdx,timeRange,photometryNI,blueWhiteRed(1,:),params,photometrySystem='ni');
    [~,~] = plotTraces(goIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni');
    [~,~] = plotTraces(nogoIdx,timeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score');
    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
        ['Lick trials (n=',num2str(length(goIdx)),')'],...
        ['No-lick trials (n=',num2str(length(nogoIdx)),')']};
    legend(taskLegend,'Location','northeast');

    % 2.2 Plot lick traces
    nexttile;
    plotLicks(waterIdx,timeRange,lick_binSize,blueWhiteRed(1,:),[],rightLick,params);
    plotLicks(goIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
    plotLicks(nogoIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('Licks/s');
    taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
        ['Lick trials (n=',num2str(length(goIdx)),')'],...
        ['No-lick trials (n=',num2str(length(nogoIdx)),')']};
    legend(taskLegend,'Location','northeast');

    % 2.3 Plot photometry traces for trial subtype
    nexttile;
    [~,~] = plotTraces(hitIdx,timeRange,photometryNI,blueWhiteRed(1,:),params,photometrySystem='ni');
    [~,~] = plotTraces(missIdx,timeRange,photometryNI,blueWhiteRed(100,:),params,photometrySystem='ni');
    [~,~] = plotTraces(faIdx,timeRange,photometryNI,blueWhiteRed(400,:),params,photometrySystem='ni');
    [~,~] = plotTraces(crIdx,timeRange,photometryNI,blueWhiteRed(500,:),params,photometrySystem='ni');
    plotEvent('Cue',0.025,'r');
    xlabel('Time (s)'); ylabel('z-score');
    taskLegend = {['Hit trials (n=',num2str(length(hitIdx)),')'],...
        ['Miss trials (n=',num2str(length(missIdx)),')'],...
        ['False alarm trials (n=',num2str(length(faIdx)),')'],...
        ['Correct rejection trials (n=',num2str(length(crIdx)),')']};
    legend(taskLegend,'Location','northeast');

    % 2.2 Plot lick traces
    nexttile;
    plotLicks(hitIdx,timeRange,lick_binSize,blueWhiteRed(1,:),[],rightLick,params);
    plotLicks(missIdx,timeRange,lick_binSize,blueWhiteRed(100,:),[],rightLick,params);
    plotLicks(faIdx,timeRange,lick_binSize,blueWhiteRed(400,:),[],rightLick,params);
    plotLicks(crIdx,timeRange,lick_binSize,blueWhiteRed(500,:),[],rightLick,params);
    plotEvent('Cue',0.025,'r');
    xlabel('Time (s)'); ylabel('Licks/s');
    taskLegend = {['Hit trials (n=',num2str(length(hitIdx)),')'],...
        ['Miss trials (n=',num2str(length(missIdx)),')'],...
        ['False alarm trials (n=',num2str(length(faIdx)),')'],...
        ['Correct rejection trials (n=',num2str(length(crIdx)),')']};
    legend(taskLegend,'Location','northeast');
    
    saveas(gcf,strcat(sessionpath,filesep,'psth_ni_combined_',session.name,'.png'));
end

%% Plot lick raster plot

timeRange = [-1,5]; markerSize = 20;

initializeFig(1,.67);

% Get event time and number by trial type
hitTrialIdx = trials{strcmp('H',trials.outcome),["TrialNumber","CueTime","OutcomeTime"]};
missTrialIdx = trials{strcmp('M',trials.outcome),["TrialNumber","CueTime","OutcomeTime"]};
faTrialIdx = trials{strcmp('FA',trials.outcome),["TrialNumber","CueTime","OutcomeTime"]};
crTrialIdx = trials{strcmp('CR',trials.outcome),["TrialNumber","CueTime","OutcomeTime"]};
hitTrialIdx(:,3) = hitTrialIdx(:,3)./params.sync.behaviorFs;
missTrialIdx(:,3) = missTrialIdx(:,3)./params.sync.behaviorFs;
faTrialIdx(:,3) = faTrialIdx(:,3)./params.sync.behaviorFs;
crTrialIdx(:,3) = crTrialIdx(:,3)./params.sync.behaviorFs;

% getLicks by trial type
[~,~,hitLicks] = getLicks(timeRange,hitTrialIdx(:,2),lick_binSize,[],rightLick,...
                            params.sync.behaviorFs,params.sync.timeNI);
[~,~,missLicks] = getLicks(timeRange,missTrialIdx(:,2),lick_binSize,[],rightLick,...
                            params.sync.behaviorFs,params.sync.timeNI);
[~,~,faLicks] = getLicks(timeRange,faTrialIdx(:,2),lick_binSize,[],rightLick,...
                            params.sync.behaviorFs,params.sync.timeNI);
[~,~,crLicks] = getLicks(timeRange,crTrialIdx(:,2),lick_binSize,[],rightLick,...
                            params.sync.behaviorFs,params.sync.timeNI);


% Plot overall raster plot (color coded by trial type)
for i = 1:size(hitLicks,1)
    scatter(hitLicks{i},hitTrialIdx(i,1),markerSize,'filled','MarkerFaceColor',blueWhiteRed(1,:)); hold on
    scatter(hitTrialIdx(i,3),hitTrialIdx(i,1),markerSize+10,blueWhiteRed(1,:),'pentagram','filled'); hold on
end
for i = 1:size(missLicks,1)
    scatter(missLicks{i},missTrialIdx(i,1),markerSize,'filled','MarkerFaceColor',blueWhiteRed(100,:)); hold on
    scatter(missTrialIdx(i,3),missTrialIdx(i,1),markerSize+10,blueWhiteRed(1,:),'pentagram','filled'); hold on
end
for i = 1:size(faLicks,1)
    scatter(faLicks{i},faTrialIdx(i,1),markerSize,'filled','MarkerFaceColor',blueWhiteRed(400,:)); hold on
    %scatter(faTrialIdx(i,3),faTrialIdx(i,1),markerSize+10,blueWhiteRed(1,:),'pentagram','filled'); hold on
end
for i = 1:size(crLicks,1)
    scatter(crLicks{i},crTrialIdx(i,1),markerSize,'filled','MarkerFaceColor',blueWhiteRed(500,:)); hold on
    %scatter(crTrialIdx(i,3),crTrialIdx(i,1),markerSize+10,blueWhiteRed(1,:),'pentagram','filled'); hold on
end
xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Trial'); ylim([0,size(trials,1)]);
plotEvent("Cue",0.025,'r');


saveas(gcf,strcat(sessionpath,filesep,'lick_summary_',session.name,'.png'));


end