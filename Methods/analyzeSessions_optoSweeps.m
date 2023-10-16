function analyzeSessions_optoSweeps(sessionpath,task,stimPattern,options)

arguments
    sessionpath string
    task string
    stimPattern cell
    options.plot logical = true % Plot session figures
    options.redo logical = true % Recalculate trial table and all preprocessing
    options.stay logical = true % do not exit the method so I can do further plotting
end

%% Notes
% Shun_analyzeBehavior_optoSweeps
% Shun Li, 2023/07/25
% Rewrite from analyzeSessions_optoPair to plot traces for sweeps of opto
% stimulation

% Assume trial structure as follows:
% 1. random/pair opto trials
% 2. opto sweeps at different frequencies or power (same repetition)
% PS: some reward happens in between each optosweep pattern

% Variables to set:
% 3. Opto sweep pattern
% 4. Repetitions per opto sweep pattern

%% Load data

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
             
% 1. Select session via uigetdir
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; 
if ispc; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
elseif isunix; session.projectPath = strcat('/',fullfile(dirsplit{2:end-1}));
end
clear dirsplit

% Opto stim params
sweepFreqList = str2num(stimPattern{1}); % str2double doesn't work
repeatPerPattern = str2double(stimPattern{2});  

disp(strcat('********** ',sessionName,'**********'));
load(strcat(sessionpath,filesep,'sync_',sessionName,'.mat'));
if ~isfield(session,'name'); session.name = sessionName; end
disp(['Session ',sessionName,' loaded']);

% Set other trial related variables
timeRange = [-0.5,3]; minLicks = 1; reactionTimeSamp = 2 * params.sync.behaviorFs;
toneDuration = 0.5;

lick_binSize = 0.2; blink_thresh = 2.5; % in turns of z score

%% Preprocess outcome and opto data

% Reward/punishment params
rewardUnit = 0.012; % 8ms opening to dispense 1ul water
rewardList = [0 2 5]; % in ul
punishList = [0 0.1];
toneList = [0 0.5 1]; % in sec

disp('Ongoing: preprocess outcome and opto data');

if ~exist('rightSolenoid_rounded','var') || options.redo
    rightSolenoid = rightSolenoid ./ rewardUnit;
    
    % Round reward and tone
    leftTone_rounded = roundToTarget(leftTone, toneList); disp('Finished rounding: leftTone');
    rightSolenoid_rounded = roundToTarget(rightSolenoid, rewardList); disp('Finished rounding: rightSolenoid');
    airpuff_rounded = roundToTarget(airpuff,punishList); disp('Finished rounding: airpuff');
    
    % Save rounded data
    save(strcat(sessionpath,filesep,'sync_',session.name),...
        'leftTone_rounded','rightSolenoid_rounded','airpuff_rounded','-append');
    disp('Finished: rounding cue/outcome data');
end

if ~exist('firstPulse','var') || options.redo
    if ~isempty(find(redLaser, 1))
        allPulses = find(redLaser);
        intervalThreshold = 10000;
        temp_interval = [100000,diff(allPulses)];
        firstPulse = allPulses(temp_interval > intervalThreshold)';

        % If repeatPerPattern is none (i.e. only one freq), determine based on number of
        % laser pulses
        if (repeatPerPattern==0 || isnan(repeatPerPattern)) && length(sweepFreqList) == 1
            repeatPerPattern = length(firstPulse);
        end

        % Return error if total pulse mismatches
        if length(sweepFreqList) * repeatPerPattern > length(firstPulse)
            disp(['Supposed total pulse: ', num2str(length(sweepFreqList)), 'frequencies x ',...
                num2str(repeatPerPattern), ' repeatPerPattern = ', ...
                num2str(length(sweepFreqList) * repeatPerPattern)]);
            disp(['Actual total pulse: ', num2str(length(firstPulse))]);
            error('Error; mismatch between actual total pulse and supposed total pulse!');
        end
    
        % Add extra column to first pulse indicate pulse pattern of that pulse
        firstPulse(end-length(sweepFreqList)*repeatPerPattern+1:end,2) = repelem(sweepFreqList,repeatPerPattern);

        % Save first pulse data
        save(strcat(sessionpath,filesep,'sync_',session.name),"firstPulse",'-append');
        disp('Finished: first pulse data');
    else 
        firstPulse = [];
    end
end

%% Stop if dont want to plot figures

if ~options.plot; return; end

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


%% Event params
    
% Select event idx
waterIdx = find(rightSolenoid);  
airpuffIdx = find(airpuff_rounded);
toneIdx = find(leftTone);
stimIdx = firstPulse;

% Create task legend
taskLegend = {['Water (n=',num2str(length(waterIdx)),')'],...
        ['Airpuff (n=',num2str(length(airpuffIdx)),')'],...
        ['Tone (n=',num2str(length(toneIdx)),')'],...
        ['Avg stim (n=',num2str(length(stimIdx)),')']};

stimLegend = cell(length(sweepFreqList)+1,1);
for i= 1:length(sweepFreqList)
    stimLegend{i} = [num2str(sweepFreqList(i)),'Hz stim (n=',num2str(repeatPerPattern),')'];
    if i==length(sweepFreqList)
        stimLegend{i+1} = ['Airpuff (n=',num2str(length(airpuffIdx)),')'];
    end
end
% stimLegend = {[num2str(sweepFreqList(1)),'Hz stim (n=',num2str(repeatPerPattern),')'],...
%     [num2str(sweepFreqList(2)),'Hz stim (n=',num2str(repeatPerPattern),')'],...
%     [num2str(sweepFreqList(3)),'Hz stim (n=',num2str(repeatPerPattern),')'],...
%     [num2str(sweepFreqList(4)),'Hz stim (n=',num2str(repeatPerPattern),')'],...
%     [num2str(sweepFreqList(5)),'Hz stim (n=',num2str(repeatPerPattern),')'],...
%     ['Airpuff (n=',num2str(length(airpuffIdx)),')']};


%% (LJ) Plot combined photometry PSTHs

if session.withPhotometry
    timeRange = [-1,3]; lick_binSize = 0.2;
    
    % 2. Plot traces
    initializeFig(0.5,0.5);

    % 2.1 Plot photometry traces
    subplot(2,1,1)
    [~,~] = plotTraces(waterIdx,timeRange,rollingGreenLP,bluePurpleRed(1,:),params);
    [~,~] = plotTraces(airpuffIdx,timeRange,rollingGreenLP,[0.2, 0.2, 0.2],params);
    [~,~] = plotTraces(toneIdx,timeRange,rollingGreenLP,bluePurpleRed(350,:),params);
    [~,~] = plotTraces(stimIdx,timeRange,rollingGreenLP,bluePurpleRed(end,:),params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score');
    legend(taskLegend,'Location','northeast');

    % 2.2 Plot different stim frequencies with airpuff
    subplot(2,1,2)
    palette = floor(linspace(1,500,length(sweepFreqList)));
    for i = 1:length(sweepFreqList)
        patternIdx = firstPulse(firstPulse(:,2)==sweepFreqList(i),1);
        [~,~] = plotTraces(patternIdx,timeRange,rollingGreenLP,bluePurpleRed(palette(i),:),params);
    end
    [~,~] = plotTraces(airpuffIdx,timeRange,rollingGreenLP,[0.2, 0.2, 0.2],params);
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score');
    legend({stimLegend},'Location','northeast');
    
    % % 2.2 Plot lick traces
    % subplot(2,1,2)
    % plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
    % plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
    % plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
    % plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
    % plotEvent('',0,'r');
    % xlabel('Time (s)'); ylabel('Licks/s');  %legend('Shutter','Water','Stim'); 
    % legend(taskLegend,'Location','best');

    saveas(gcf,strcat(sessionpath,filesep,'psth_lj_combined_',session.name,'.png'));
end

%% (LJ) Plot single stimulus PSTH

if session.withPhotometry
    eventIdxes = {stimIdx,waterIdx,toneIdx,airpuffIdx};
    labels = {'Stim','Water','Tone','Airpuff'};
    eventDurations = [0.5,0,0.5,0.01];
    groupSizes = [10,30,10,30];
    longTimeRange = [-5,10];
    shortTimeRange = [-0.5,3]; 
   
    for event = 1:length(eventIdxes)
        eventIdx = eventIdxes{event};
        if isempty(eventIdx); continue; end
        label = labels{event}; eventDuration = eventDurations(event); groupSize = groupSizes(event); % num of trials to plot in one line
        
        initializeFig(0.5,1);
        if strcmp(label,'Stim')
            subplot(4,1,1)
            palette = floor(linspace(1,500,length(sweepFreqList)));
            for i = 1:length(sweepFreqList)
                patternIdx = firstPulse(firstPulse(:,2)==sweepFreqList(i),1);
                [traces,t] = plotTraces(patternIdx,longTimeRange,rollingGreenLP,bluePurpleRed(palette(i),:),params);
            end
            plotEvent('',0,'r');
            xlabel('Time (s)'); ylabel('z-score');
            legend(stimLegend,'Location','northeast');
        else
            subplot(4,1,1)
            [traces,t] = plotTraces(eventIdx,longTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
            % [~,~] = plotTraces(randShutterIdx,longTimeRange,rollingGreenLP,[.75 .75 .75],params);
            plotEvent(label,eventDuration,bluePurpleRed(end,:));
            xlabel('Time (s)'); ylabel('z-score');
            legend({[label,' (n=',num2str(length(eventIdx)),')']},...
                'Location','northeast'); 
        end
        
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
        
        if strcmp(label,'Stim')
            subplot(4,1,2)
            palette = floor(linspace(1,500,length(sweepFreqList)));
            for i = 1:length(sweepFreqList)
                patternIdx = firstPulse(firstPulse(:,2)==sweepFreqList(i),1);
                [traces,t] = plotTraces(patternIdx,shortTimeRange,rollingGreenLP,bluePurpleRed(palette(i),:),params);
            end
            plotEvent('',0,'r');
            xlabel('Time (s)'); ylabel('z-score');
            legend(stimLegend,'Location','northeast');
        else
            subplot(4,1,2)
            [traces,t] = plotTraces(eventIdx,shortTimeRange,rollingGreenLP,bluePurpleRed(1,:),params);
            % [~,~] = plotTraces(randShutterIdx,shortTimeRange,rollingGreenLP,[.75 .75 .75],params);
            plotEvent(label,eventDuration,bluePurpleRed(end,:)); 
            xlabel('Time (s)'); ylabel('z-score'); % legend('Shutter',label);
            legend({[label,' (n=',num2str(length(eventIdx)),')']},...
                'Location','northeast'); 
        end
        
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
        
        saveas(gcf,strcat(sessionpath,filesep,'psth_lj_',label,'_',session.name,'.png'));
    end
end

%% (NI) Plot combined PSTH

if session.ni_photometryON
    timeRange = [-1,3]; lick_binSize = 0.2;
    
    % 2. Plot traces
    initializeFig(0.5,0.5);

    % 2.1 Plot photometry traces
    subplot(2,1,1)
    [~,~] = plotTraces(waterIdx,timeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
    [~,~] = plotTraces(airpuffIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni');
    [~,~] = plotTraces(toneIdx,timeRange,photometryNI,bluePurpleRed(350,:),params,photometrySystem='ni');
    [~,~] = plotTraces(stimIdx,timeRange,photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni');
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score'); 
    legend(taskLegend,'Location','best'); 

    % 2.2 Plot different stim frequencies with airpuff
    subplot(2,1,2)
    palette = floor(linspace(1,500,length(sweepFreqList)));
    for i = 1:length(sweepFreqList)
        patternIdx = firstPulse(firstPulse(:,2)==sweepFreqList(i),1);
        [~,~] = plotTraces(patternIdx,timeRange,photometryNI,bluePurpleRed(palette(i),:),params,photometrySystem='ni');
    end
    [~,~] = plotTraces(airpuffIdx,timeRange,photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni');
    plotEvent('',0,'r');
    xlabel('Time (s)'); ylabel('z-score');
    legend(stimLegend,'Location','northeast');
    
    % % 2.2 Plot lick traces
    % subplot(2,1,2)
    % plotLicks(waterIdx,timeRange,lick_binSize,bluePurpleRed(1,:),[],rightLick,params);
    % plotLicks(airpuffIdx,timeRange,lick_binSize,[0.2, 0.2, 0.2],[],rightLick,params);
    % plotLicks(toneIdx,timeRange,lick_binSize,bluePurpleRed(350,:),[],rightLick,params);
    % plotLicks(stimIdx,timeRange,lick_binSize,bluePurpleRed(end,:),[],rightLick,params);
    % plotLicks(randShutterIdx,timeRange,lick_binSize,[.75 .75 .75],[],rightLick,params);
    % plotEvent('',0,'r');
    % xlabel('Time (s)'); ylabel('Licks/s'); 
    % legend(taskLegend,'Location','best');

    saveas(gcf,strcat(sessionpath,filesep,'psth_ni_combined_',session.name,'.png'));
end

%% (NI) Plot single stimulus PSTH

if session.ni_photometryON

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

        if strcmp(label,'Stim')
            subplot(4,1,1)
            palette = floor(linspace(1,500,length(sweepFreqList)));
            for i = 1:length(sweepFreqList)
                patternIdx = firstPulse(firstPulse(:,2)==sweepFreqList(i),1);
                [traces,t] = plotTraces(patternIdx,longTimeRange,photometryNI,bluePurpleRed(palette(i),:),params,photometrySystem='ni');
            end
            plotEvent('',0,'r');
            xlabel('Time (s)'); ylabel('z-score');
            legend(stimLegend,'Location','northeast');
        else
            subplot(4,1,1)
            [traces,t] = plotTraces(eventIdx,longTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
            plotEvent(label,eventDuration,bluePurpleRed(end,:));
            xlabel('Time (s)'); ylabel('z-score');
            legend({[label,' (n=',num2str(length(eventIdx)),')']},...
                'Location','northeast'); 

        end
        
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
        
        if strcmp(label,'Stim')
            subplot(4,1,2)
            palette = floor(linspace(1,500,length(sweepFreqList)));
            for i = 1:length(sweepFreqList)
                patternIdx = firstPulse(firstPulse(:,2)==sweepFreqList(i),1);
                [traces,t] = plotTraces(patternIdx,shortTimeRange,photometryNI,bluePurpleRed(palette(i),:),params,photometrySystem='ni');
            end
            plotEvent('',0,'r');
            xlabel('Time (s)'); ylabel('z-score');
            legend(stimLegend,'Location','northeast');
        else
            subplot(4,1,2)
            [traces,t] = plotTraces(eventIdx,shortTimeRange,photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni');
            % [~,~] = plotTraces(randShutterIdx,shortTimeRange,photometryNI,[.75 .75 .75],params,photometrySystem='ni');
            plotEvent(label,eventDuration,bluePurpleRed(end,:));
            xlabel('Time (s)'); ylabel('z-score');
            legend({[label,' (n=',num2str(length(eventIdx)),')']},...
                'Location','northeast'); 
        end
        
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
        
        saveas(gcf,strcat(sessionpath,filesep,'psth_ni_',label,'_',session.name,'.png'));
    end
end

return
%% Plot lick raster plot

timeRange = [-1,5]; markerSize = 20;

if contains(task,'Optopair')
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
    
    saveas(gcf,strcat(sessionpath,filesep,'lick_summary_',session.name,'.png'));
    
end

return

end