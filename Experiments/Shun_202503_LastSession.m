% Shun_analyzeExperiments_template
% 2023/12/04

% Template for multi-session analysis of an experiment

% To start, copy this matlab file and replace template with specific
% experiments. In theory, there will be a specific analyzeExperiments file
% for each individual experiments, as the specific needs for analysis
% varies between different experiments.

%% Analysis pipeline
% The pipeline in general is the following:

% 1. Select whether to load a previously animals struct (described below)
% or select individual session to combine.

% 2. After selecting ALL SESSIONS from an experiments, the pipeline will
% automatically concatenate analysis.mat for each recording sessions.
% Rename properties as needed in order to facilitate further analysis.

% 3. Run getAnimalStruct.m function to recreate animals struct from
% summary. This combines all sessions from the same animals together while
% cutoffs between individual sessions are also recorded.

% 4. Save animals struct if needed. Note: saving summary struct will take
% extremely long (>5hrs) so while saving animals struct is much shorter 
% (~2min). animals struct should contain information that satisfies MOST 
% plotting requirements so saving summary struct is not needed.

% 5. Data analysis and plotting. This part is designed to vary across
% experiments. Thus, following codes are just for demonstration of 
% essential functions.

%% Essential functions

% getAnimalStruct(summary)
% combine sessions of the same animal, from the same task,
% of the same event, recorded from the same signal (eg NAc, LHb, cam,
% Lick) together. As described above, animals struct will be the MOST
% IMPORTANT struct that stores information about the experiments.

% combineTraces(animals,options)
% combine traces and their relevant statstics of selected animals, 
% selected tasks, selected trialRange, selected totalTrialRange, 
% selected signals, and selected events together. 
% This is the MOST IMPORTANT and USED function in this script. 
% Important features are listed as follows:
    % 1. the function returns a structure with fields. data fields stores
    % the data (photometry, cam, lick rate traces) of selected sessions.
    % 2. Field stats stores stageAvg/Max/Min of each traces at selected stage
    % time (often determined when creating analysis.mat but can modify later).
    % 3. Field options contains following important variables:
        % options.empty: true if no session is found that fits the input criteria. 
        % Should skip during plotting or further analysis 
        % options.animalStartIdx: Records index (in field data) of the
        % first trace for each animals. Used in plotGroupTraces
        % options.sessionStartIdx: Records index (in field data) of the
        % first trace for each session.
    % 4. totalTrialRange and trialRange
        % totalTrialRange selects the ACTUAL trial number within each
        % session while trialRange selects the samples across selected
        % sessions. For example: I have 3 session where I inhibit CaMKII
        % activity for the first 60 trials of each session. Within these
        % first 60 trials, 30% of them are stim-only trials. If I want to
        % only plot the 50-100th stim-only trials with CaMKII inhibition 
        % across all sessions, I will set totalTrialRange=[1,60] and
        % trialRange=[50,100]. Detailed description and automatic handling
        % of edge cases is documented within the method.
    % 5. The function can take both animals and summary struct as inputs.

% plotGroupTraces(combined.data,combined.timestamp,options)
% While plotGroupTraces is also used in analyzeSessions.m; here, we can
% plot traces across all animal easily (see code below). Key options are as
% follows:
    % 1. groupSize and nGroups
        % You need to provide either groupSize or nGroups for the function to
        % run. If you provide both, plotGroupTraces will plot to the maximum
        % number of groups based on groupSize. Thus, for a input with 50
        % trials and groupSize = 10, the function will automatically plot 5
        % lines even when nGroups=10
    % 2. options.animalStartIdx
        % Use this to reorganize input data so that its plotted based on
        % animals. eg when I want to plot Trial 1-10, 11-20 for EACH animal
        % across all sessions
    % 3. options.remaining
        % There inevitably will be some traces that does not fully form a group
        % (eg 5 traces remaining for a groupSize of 50 traces). These traces,
        % if plotted separately, can induce lines with great variations and
        % error bars. To address this, one can either set remaining='include'
        % to include these traces to prev group; set to 'exclude' to not plot
        % these traces, or 'separate' if you really want to plot these traces
        % separately

%% Intro of sample data set

% The sample data set is recorded by Shun Li in 2023. It contains 4
% animals, with 1 animals with off-target expression ('SL137'). dLight
% signals in NAc, pupil/Eye area, and lick are simultaneously recorded for
% all sessions.

% There are 5 major phases:
    % 1. Random: water, airpuff, EP stim, and tone (75dB) are delivered
    % randomly
    % 2. Reward1/2: where EP stim and tone are paired with water
    % 3. Punish1/2: where EP stim and tone are paired with airpuff
    % 4. Timeline: Random (2 sessions) -> Reward1 (3 sessions) -> Punish1
    % (3 sessions) -> Reward2 (3 sessions) -> 1 week rest -> Punish2 (3 sessions but 3 animals)

%% Setup

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[~,~,~,~,~,~,bluePurpleRed] = loadColors;

% Define result directory
resultspath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'))';
    groupSessions = true;
    % Update resultspath
    projectName = 'Summary-LastSessions'; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
    fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results'))';
    groupSessions = false;
    % Update resultspath
    projectName = 'Summary-LastSessions'; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Load selected files
    for file = 1:length(fileList)
        dirsplit = strsplit(fileList{file},filesep);
        disp(['Ongoing: loading ',dirsplit{end}]);
        load(fileList{file});
        disp(['Finished: loaded ',dirsplit{end}]);
    end

elseif strcmpi(answer,'Load sample data')
    groupSessions = false;
    % Update resultspath
    dirsplit = strsplit(fileList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(osPathSwith('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Tutorials/Sample data/Results'),filesep,projectName);
    % Load selected files
    for file = 1:length(fileList)
        dirsplit = strsplit(fileList{file},filesep);
        disp(['Ongoing: loading ',dirsplit{end}]);
        load(fileList{file});
        disp(['Finished: loaded ',dirsplit{end}]);
    end
end

%% Optional: Create summary struct (only need to do this for initial loading)

if groupSessions    
    summary = concatAnalysis(sessionList);
    trialTables = loadTrialTables(sessionList);
end

% Check summary format (should all be chars NOT strings)
stringColumnsLabels = {'animal','date','session','task','event','name','system'};
for i = 1:length(stringColumnsLabels)
    for row = 1:length(summary)
        if isstring(summary(row).(stringColumnsLabels{1}))
            summary(row).(stringColumnsLabels{i}) = convertStringsToChars(summary(row).(stringColumnsLabels{i}));
        end
    end
end
disp('Finished: summary struct and trialtables loaded');

%% Optional: Make changes to summary for further analysis (first reward & punish sessions)

% Change some names if needed
for i = 1:length(summary)
    cur_task = summary(i).task;
    cur_event = summary(i).event;
    cur_name = summary(i).name;
    cur_animal = summary(i).animal;

    if strcmp('random',cur_task)
        summary(i).task = 'Random';
    elseif contains('reward pairing',cur_task,IgnoreCase=true)
        summary(i).task = 'Reward';
    elseif contains('punish pairing',cur_task,IgnoreCase=true)
        summary(i).task = 'Punish';
    end

    if contains('Stim only',cur_event,IgnoreCase=true)
        summary(i).event = 'Stim';
    elseif contains('Tone only',cur_event,IgnoreCase=true)
        summary(i).event = 'Tone';
    elseif contains('Pair',cur_event,IgnoreCase=true)
        summary(i).event = 'Pair';
    end

    if strcmpi('NAc-green',cur_name)
        summary(i).name = 'dLight';
    % elseif strcmpi('PMT',cur_name)
    %     disp(i);
    %     summary(i).name = 'GCaMP8m';
    end
end

% Remove some rows if needed
% rowIdx = cellfun(@(x) contains(x,'Blue stim',IgnoreCase=true), {summary.event});
% summary(rowIdx) = [];

%% Create animals struct

if isempty(dir(fullfile(resultspath,'animals*.mat'))) || groupSessions
    animals = getAnimalsStruct(summary);
end

%% Save animals struct

prompt = 'Enter database notes (animals_20230326_notes.mat):';
dlgtitle = 'Save animals struct'; fieldsize = [1 45]; definput = {''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
today = char(datetime('today','Format','yyyyMMdd'));
filename = strcat('animals_',today,'_',answer{1});

% Save animals.mat
if ~isempty(answer)
    disp(['Ongoing: saving animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
    save(strcat(resultspath,filesep,filename),'animals','trialTables','sessionList','-v7.3');
    disp(['Finished: saved animals.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
end


%% Test: Plot traces from summary/animals struct

event = 'Stim';
animal = 'All';
task = 'Reward';
signal = 'dLight';
session = 'all';

initializeFig(0.5,0.5);
combined = combineTraces(animals,timeRange=[-0.5,3],...
                            eventRange=event,...
                            animalRange=animal,...
                            taskRange=task,...
                            totalTrialRange='All',...
                            trialRange='All',...
                            signalRange=signal,...
                            sessionRange=session);

% plotTraces(combined.data{1},combined.timestamp,color=bluePurpleRed(500,:));
legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                        groupby='sessions',startIdx=combined.options.startIdx,remaining='include');

plotEvent('Stim',0.5,color=bluePurpleRed(500,:))
xlabel('Time (s)'); ylabel('z-score');


%% Calculate DA trend for each animal

eventRange = 'Stim';
animalRange = unique({animals.animal});
DAtrend = struct([]);
nboot = 1000;
trialConditions = 'trials.performing';
maxTrials = max(arrayfun(@(x) length(x.CueMax_slope), DAtrend));

for a = 1:length(animalRange)
    cur_animal = animalRange{a}; disp(['Ongoing: analyzing ',cur_animal]);
    cur_task = unique({animals(strcmp({animals.animal}, cur_animal)).task});
    DA_Max = getGroupedTrialStats(animals,'stageMax',...
                                eventRange=eventRange,...
                                animalRange=cur_animal,...
                                taskRange=cur_task,...
                                signalRange='dLight',...
                                trialConditions=trialConditions);
    DA_Min = getGroupedTrialStats(animals,'stageMin',...
                                eventRange=eventRange,...
                                animalRange=cur_animal,...
                                taskRange=cur_task,...
                                signalRange='dLight',...
                                trialConditions=trialConditions);
    DA_Avg = getGroupedTrialStats(animals,'stageAvg',...
                                eventRange=eventRange,...
                                animalRange=cur_animal,...
                                taskRange=cur_task,...
                                signalRange='dLight',...
                                trialConditions=trialConditions);
    % initializeFig(.7,.7); tiledlayout('flow');
    % DA_stageMax_result = plotGroupedTrialStats(DA_stageMax,'DA Max',groupSize=1,color=bluePurpleRed(500,:));

    % Save to DAtrend
    DAtrend(a).animal = cur_animal;
    DAtrend(a).task   = cur_task;
    DAtrend(a).CueMax = DA_Max.stats{1}{1}(:,2);
    DAtrend(a).CueMin = DA_Min.stats{1}{1}(:,2);
    DAtrend(a).CueAvg = DA_Avg.stats{1}{1}(:,2);
    DAtrend(a).CueAmp = getAmplitude(DA_Max.stats{1}{1}(:,2),DA_Min.stats{1}{1}(:,2));

    % Save smoothed DAtrend
    DAtrend(a).CueMax_smoothed = smoothdata(DAtrend(a).CueMax,'movmean',3);
    DAtrend(a).CueMin_smoothed = smoothdata(DAtrend(a).CueMin,'movmean',3);
    DAtrend(a).CueAvg_smoothed = smoothdata(DAtrend(a).CueAvg,'movmean',3);
    DAtrend(a).CueAmp_smoothed = smoothdata(DAtrend(a).CueAmp,'movmean',3);

    % Calculate DA trending stats
    fields = {'CueMax', 'CueMin', 'CueAvg', 'CueAmp',...
              'CueMax_smoothed','CueMin_smoothed','CueAvg_smoothed','CueAmp_smoothed'};
    nTrials = length(DAtrend(a).CueMax); % assuming all have the same length
    
    for f = 1:numel(fields)
        data = DAtrend(a).(fields{f});
        slopes = nan(nTrials, 1);
        pvals  = nan(nTrials, 1);
        
        % Loop for trending statistics: last n trials (n = 1:nTrials)
        for n = 3:nTrials
            Y = data(end-n+1:end);  % last n trials
            X = (1:n)';             % corresponding trial indices for regression
            fit_obs = polyfit(X, Y, 1); obsSlope = fit_obs(1);
            fit_boot = bootstrp(nboot, @(i) polyfit(X, Y(i), 1), 1:length(Y));
            slopes(n) = obsSlope;
            pvals(n)  = sum(abs(fit_boot(:,1)) >= abs(obsSlope))/nboot;   % p-value for the slope
        end
        
        % Store the trending statistics into DAtrend
        DAtrend(a).([fields{f} '_slope']) = slopes;
        DAtrend(a).([fields{f} '_pval'])  = pvals;
        
        % Compute difference between the last 10 trials and the last 30-40 trials
        if nTrials >= 40
            DAtrend(a).([fields{f} '_diff']) = mean(data(end-40+1:end-30)) - mean(data(end-10+1:end));
        else
            DAtrend(a).([fields{f} '_diff']) = NaN;  % Not enough trials
        end
    end
end

%% Save DAtrend struct
prompt = 'Enter database notes (DAtrend_20230326_notes.mat):';
dlgtitle = 'Save DAtrend struct'; fieldsize = [1 45]; definput = {''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
today = char(datetime('today','Format','yyyyMMdd'));
filename = strcat('DAtrend_',today,'_',answer{1});

disp(['Ongoing: saving DAtrend.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
save(strcat(resultspath,filesep,filename),'DAtrend','sessionList','-v7.3');
disp(['Finished: saved DAtrend.mat (',char(datetime('now','Format','HH:mm:ss')),')']);


%% Plot DA trend for each animal

animalRange = unique({animals.animal});
shortWindow = 20;%16:21;
longWindow  = 40;%38:43;

shortWindowColor = [128, 179, 255]./255;
longWindowColor = [143, 88, 219]./255;

for a = 1:length(animalRange)
    
    % Calculate p-value during plasticity window
    cur_animal = animalRange{a};

    % Plot DA trend summary for this animal
    close all; initializeFig(.5,1); tiledlayout(4,3);

    nexttile([2 1]);
    combined = combineTraces(animals,timeRange=[-0.5,1],...
                                eventRange='Stim',...
                                animalRange=cur_animal,...
                                signalRange='dLight');
    legendList = plotGroupTraces(combined.data{1},combined.timestamp,bluePurpleRed,...
                            groupby='trials',groupSize=10,remaining='include',...
                            LineWidth=5,plotPatch=true);
    plotEvent('Stim',0.5,color=bluePurpleRed(500,:));
    xlabel('Time (s)'); ylabel('dLight z-score'); legend(legendList);

    nexttile; 
    data = DAtrend(a).CueAmp; smoothed = DAtrend(a).CueAmp_smoothed;
    scatter(1:length(data),data,80,[0.3010 0.7450 0.9330],'filled'); hold on
    plot(smoothed,color=[0 0.4470 0.7410],LineWidth=5);
    plotEvent('Long window',range(longWindow),color=longWindowColor,x=numel(data)-longWindow(end));
    plotEvent('',longWindow(1),color=longWindowColor,x=numel(data)-longWindow(1));
    plotEvent('Short window',range(shortWindow),color=shortWindowColor,x=numel(data)-shortWindow(end));
    plotEvent('',shortWindow(1),color=shortWindowColor,x=numel(data)-shortWindow(1));
    xlabel('Trials'); ylabel('DA Amp during cue');
    xlim([1,numel(data)]);
    title('DA Amplitude'); box off;

    nexttile;
    plotScatterBar(1,mean(DAtrend(a).CueAmp_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(2,mean(DAtrend(a).CueAmp_smoothed_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(3,mean(DAtrend(a).CueAmp_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    plotScatterBar(4,mean(DAtrend(a).CueAmp_smoothed_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    % plotScatterBar(5,DAtrend(a).CueAmp_diff,style='bar'); hold on;
    text(1-0.15,mean(DAtrend(a).CueAmp_slope(shortWindow)),num2str(mean(DAtrend(a).CueAmp_pval(shortWindow)),'%.3f'));
    text(2-0.15,mean(DAtrend(a).CueAmp_smoothed_slope(shortWindow)),num2str(mean(DAtrend(a).CueAmp_smoothed_pval(shortWindow)),'%.3f'));
    text(3-0.15,mean(DAtrend(a).CueAmp_slope(longWindow)),num2str(mean(DAtrend(a).CueAmp_pval(longWindow)),'%.3f'));
    text(4-0.15,mean(DAtrend(a).CueAmp_smoothed_slope(longWindow)),num2str(mean(DAtrend(a).CueAmp_smoothed_pval(longWindow)),'%.3f'));
    xticks([1 2 3 4]);
    xticklabels({'Amp (short)','Amp smoothed (short)',...
                 'Amp (long)','Amp smoothed (long)'});
    ylabel('Slope');

    nexttile;
    data = DAtrend(a).CueAvg; smoothed = DAtrend(a).CueAvg_smoothed;
    scatter(1:length(data),data,80,[0.3010 0.7450 0.9330],'filled'); hold on
    plot(smoothed,color=[0 0.4470 0.7410],LineWidth=5);
    plotEvent('Long window',range(longWindow),color=longWindowColor,x=numel(data)-longWindow(end));
    plotEvent('',longWindow(1),color=longWindowColor,x=numel(data)-longWindow(1));
    plotEvent('Short window',range(shortWindow),color=shortWindowColor,x=numel(data)-shortWindow(end));
    plotEvent('',shortWindow(1),color=shortWindowColor,x=numel(data)-shortWindow(1));
    xlabel('Trials'); ylabel('DA Avg during cue');
    xlim([1,numel(data)]);
    title('DA Average'); box off;

    nexttile;
    plotScatterBar(1,mean(DAtrend(a).CueAvg_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(2,mean(DAtrend(a).CueAvg_smoothed_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(3,mean(DAtrend(a).CueAvg_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    plotScatterBar(4,mean(DAtrend(a).CueAvg_smoothed_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    text(1-0.15,mean(DAtrend(a).CueAvg_slope(shortWindow)),num2str(mean(DAtrend(a).CueAvg_pval(shortWindow)),'%.3f'));
    text(2-0.15,mean(DAtrend(a).CueAvg_smoothed_slope(shortWindow)),num2str(mean(DAtrend(a).CueAvg_smoothed_pval(shortWindow)),'%.3f'));
    text(3-0.15,mean(DAtrend(a).CueAvg_slope(longWindow)),num2str(mean(DAtrend(a).CueAvg_pval(longWindow)),'%.3f'));
    text(4-0.15,mean(DAtrend(a).CueAvg_smoothed_slope(longWindow)),num2str(mean(DAtrend(a).CueAvg_smoothed_pval(longWindow)),'%.3f'));
    xticks([1 2 3 4]);
    xticklabels({'Avg (short)','Avg smoothed (short)',...
                 'Avg (long)','Avg smoothed (long)'});
    ylabel('Slope');

    nexttile([2 1]);
    plotHeatmap(combined.data{1},combined.timestamp);
    plotEvent('Stim',0.5,color=bluePurpleRed(500,:));
    xlabel('Time (s)'); ylabel('dLight z-score');
    

    nexttile;
    data = DAtrend(a).CueMax; smoothed = DAtrend(a).CueMax_smoothed;
    scatter(1:length(data),data,80,[0.3010 0.7450 0.9330],'filled'); hold on
    plot(smoothed,color=[0 0.4470 0.7410],LineWidth=5);
    plotEvent('Long window',range(longWindow),color=longWindowColor,x=numel(data)-longWindow(end));
    plotEvent('',longWindow(1),color=longWindowColor,x=numel(data)-longWindow(1));
    plotEvent('Short window',range(shortWindow),color=shortWindowColor,x=numel(data)-shortWindow(end));
    plotEvent('',shortWindow(1),color=shortWindowColor,x=numel(data)-shortWindow(1));
    xlabel('Trials'); ylabel('DA Max during cue');
    xlim([1,numel(data)]);
    title('DA Max'); box off;

    nexttile;
    plotScatterBar(1,mean(DAtrend(a).CueMax_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(2,mean(DAtrend(a).CueMax_smoothed_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(3,mean(DAtrend(a).CueMax_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    plotScatterBar(4,mean(DAtrend(a).CueMax_smoothed_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    text(1-0.15,mean(DAtrend(a).CueMax_slope(shortWindow)),num2str(mean(DAtrend(a).CueMax_pval(shortWindow)),'%.3f'));
    text(2-0.15,mean(DAtrend(a).CueMax_smoothed_slope(shortWindow)),num2str(mean(DAtrend(a).CueMax_smoothed_pval(shortWindow)),'%.3f'));
    text(3-0.15,mean(DAtrend(a).CueMax_slope(longWindow)),num2str(mean(DAtrend(a).CueMax_pval(longWindow)),'%.3f'));
    text(4-0.15,mean(DAtrend(a).CueMax_smoothed_slope(longWindow)),num2str(mean(DAtrend(a).CueMax_smoothed_pval(longWindow)),'%.3f'));
    xticks([1 2 3 4]);
    xticklabels({'Max (short)','Max smoothed (short)',...
                 'Max (long)','Max smoothed (long)'});
    ylabel('Slope');

    nexttile([1 1]);
    data = DAtrend(a).CueMin; smoothed = DAtrend(a).CueMin_smoothed;
    scatter(1:length(data),data,80,[0.3010 0.7450 0.9330],'filled'); hold on
    plot(smoothed,color=[0 0.4470 0.7410],LineWidth=5);
    plotEvent('Long window',range(longWindow),color=longWindowColor,x=numel(data)-longWindow(end));
    plotEvent('',longWindow(1),color=longWindowColor,x=numel(data)-longWindow(1));
    plotEvent('Short window',range(shortWindow),color=shortWindowColor,x=numel(data)-shortWindow(end));
    plotEvent('',shortWindow(1),color=shortWindowColor,x=numel(data)-shortWindow(1));
    xlabel('Trials'); ylabel('DA Min during cue');
    xlim([1,numel(data)]);
    title('DA Min'); box off;

    nexttile;
    plotScatterBar(1,mean(DAtrend(a).CueMin_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(2,mean(DAtrend(a).CueMin_smoothed_slope(shortWindow)),color=shortWindowColor,style='bar'); hold on;
    plotScatterBar(3,mean(DAtrend(a).CueMin_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    plotScatterBar(4,mean(DAtrend(a).CueMin_smoothed_slope(longWindow)),color=longWindowColor,style='bar'); hold on;
    text(1-0.15,mean(DAtrend(a).CueMin_slope(shortWindow)),num2str(mean(DAtrend(a).CueMin_pval(shortWindow)),'%.3f'));
    text(2-0.15,mean(DAtrend(a).CueMin_smoothed_slope(shortWindow)),num2str(mean(DAtrend(a).CueMin_smoothed_pval(shortWindow)),'%.3f'));
    text(3-0.15,mean(DAtrend(a).CueMin_slope(longWindow)),num2str(mean(DAtrend(a).CueMin_pval(longWindow)),'%.3f'));
    text(4-0.15,mean(DAtrend(a).CueMin_smoothed_slope(longWindow)),num2str(mean(DAtrend(a).CueMin_smoothed_pval(longWindow)),'%.3f'));
    xticks([1 2 3 4]);
    xticklabels({'Min (short)','Min smoothed (short)',...
                 'Min (long)','Min smoothed (long)'});
    ylabel('Slope');

    saveFigures(gcf,cur_animal,strcat(resultspath,filesep,'DAtrend'),...
            saveFIG=false,savePDF=false,savePNG=true);
end


%%

animal = 'SL068';

initializeFig(.5,1); tiledlayout(3,1);

nexttile;
statsTypes = 'nAnticipatoryLicks';
ylabels = 'Anticipatory licks';
stats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange='Stim',...
                            animalRange=animal,...
                            signalRange='dLight');
results = plotGroupedTrialStats(stats,ylabels,groupSize=1,...
                            plot=true,plotNextTile=false);


nexttile;
statsTypes = 'nLicks';
ylabels = 'Total licks';
stats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange='Stim',...
                            animalRange=animal,...
                            signalRange='dLight');
results = plotGroupedTrialStats(stats,ylabels,groupSize=1,...
                            plot=true,plotNextTile=false);


nexttile;
statsTypes = 'nOutcomeLicks';
ylabels = 'Outcome licks';
stats = getGroupedTrialStats(animals,statsTypes,...
                            eventRange='Stim',...
                            animalRange=animal,...
                            signalRange='dLight');
results = plotGroupedTrialStats(stats,ylabels,groupSize=1,...
                            plot=true,plotNextTile=false);
