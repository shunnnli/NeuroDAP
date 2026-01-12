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
resultspath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results');

% Building summary struct from selected sessions
answer = questdlg('Group sessions or load combined data?','Select load sources',...
                  'Group single sessions','Load combined data','Load sample data','Load combined data');

if strcmpi(answer,'Group single sessions')
    sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));
    groupSessions = true;
    % Update resultspath
    dirsplit = strsplit(sessionList{1},filesep); projectName = dirsplit{end-1}; 
    resultspath = strcat(resultspath,filesep,projectName);
    % Create resultspath if necessary
    if isempty(dir(resultspath)); mkdir(resultspath); end

elseif strcmpi(answer,'Load combined data')
    fileList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Results'));
    groupSessions = false;
    % Update resultspath
    dirsplit = strsplit(fileList{1},filesep); projectName = dirsplit{end-1}; 
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

columnLabels = {'animal','date','session','system','name','type','pwm','duration','data','timestamp'};

if groupSessions    
    summary = concatAnalysis(sessionList, columnLabels=columnLabels);
end

disp('Finished: summary struct loaded');

%% Optional: Make changes to summary for further analysis (first reward & punish sessions)

% Add rows
[summary.task] = deal('random');    % char
[summary.finalFs] = deal(50);    % char
[summary.timeRange] = deal([-1,5]);    % char
[summary.options] = deal(struct());

% Add events row
events = arrayfun(@(s) sprintf('%s-%g-%g', s.type, s.pwm, s.duration), ...
                  summary, 'UniformOutput', false);
[summary.event] = events{:};   % adds the new field to every element

% Change some names if needed
for i = 1:length(summary)
    cur_session = summary(i).session;

    if contains(cur_session,'DAT-1')
        summary(i).animal = 'DAT1';
    elseif contains(cur_session,'DAT-2')
        summary(i).animal = 'DAT2';
    elseif contains(cur_session,'DAT-3')
        summary(i).animal = 'DAT3';
    elseif contains(cur_session,'DAT-4')
        summary(i).animal = 'DAT4';
    elseif contains(cur_session,'DAT-5')
        summary(i).animal = 'DAT5';
    end

    if contains(cur_session,'combine')
        summary(i).task = 'combine';
    end
end

% Remove some rows if needed
eventIdx = cellfun(@(x) contains(x,'dLight',IgnoreCase=true), {summary.name});
summary(eventIdx) = [];

eventIdx = cellfun(@(x) contains(x,'GCaMP8m',IgnoreCase=true), {summary.name});
summary(eventIdx) = [];


%% Save summary struct

prompt = 'Enter database notes (summary_20230326_notes.mat):';
dlgtitle = 'Save summary struct'; fieldsize = [1 45]; definput = {''};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);
today = char(datetime('today','Format','yyyyMMdd'));
filename = strcat('summary_',today,'_',answer{1});

% Save animals.mat
if ~isempty(answer)
    disp(['Ongoing: saving summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
    save(strcat(resultspath,filesep,filename),'summary','sessionList','-v7.3');
    disp(['Finished: saved summary.mat (',char(datetime('now','Format','HH:mm:ss')),')']);
end

%% Plot excitation vs water

waterColor = [.34, .76, .87];

animal = 'all';

initializeFig(0.5,0.5);
% combined = combineTraces(summary,timeRange=[-0.5,3],...
%                             eventRange='water',...
%                             animalRange=animal,...
%                             taskRange='random',...
%                             combineStats=false);
% plotTraces(combined.data{1},combined.timestamp,color=waterColor);

durationList = [0.1, 0.5, 1];
stimColor = 'red';
ev  = string({summary.event});
dur = [summary.duration];
alphas = linspace(0.2, 1, length(durationList));

for d = length(durationList):-1:1
    mask = startsWith(ev,stimColor) & abs(dur - durationList(d)) < 1e-9;
    eventRange = unique(ev(mask));
    combined = combineTraces(summary,timeRange=[-0.5,3],...
                            eventRange=eventRange,...
                            animalRange=animal,...
                            taskRange='random',...
                            combineStats=false);
    durationColor = addOpacity(bluePurpleRed(500,:),alphas(d));
    plotTraces(combined.data{1},combined.timestamp,color=durationColor);
end

plotEvent('Stim',0,color=bluePurpleRed(500,:));
xlabel('Time (s)'); ylabel('z-score');

%% Plot inhibition vs airpuff

airpuffColor = [.7, .7, .7];

animal = 'all';

initializeFig(0.5,0.5);
% combined = combineTraces(summary,timeRange=[-0.5,3],...
%                             eventRange='airpuff',...
%                             animalRange=animal,...
%                             taskRange='random',...
%                             combineStats=false);
% plotTraces(combined.data{1},combined.timestamp,color=airpuffColor);

durationList = [0.1, 0.5, 1];
stimColor = 'blue';
ev  = string({summary.event});
dur = [summary.duration];
alphas = linspace(0.2, 1, length(durationList));

for d = length(durationList):-1:1
    mask = startsWith(ev,stimColor) & abs(dur - durationList(d)) < 1e-9;
    eventRange = unique(ev(mask));
    combined = combineTraces(summary,timeRange=[-0.5,3],...
                            eventRange=eventRange,...
                            animalRange=animal,...
                            taskRange='random',...
                            combineStats=false);
    durationColor = addOpacity(bluePurpleRed(1,:),alphas(d));
    plotTraces(combined.data{1},combined.timestamp,color=durationColor);
end

plotEvent('Stim',0,color=bluePurpleRed(500,:));
xlabel('Time (s)'); ylabel('z-score');


%% Plot water collision

animal = 'all';
waterColor = bluePurpleRed(1,:);

initializeFig(0.5,0.5);
combined = combineTraces(summary,timeRange=[-0.5,3],...
                            eventRange='water (ctrl)',...
                            animalRange=animal,...
                            taskRange='combine',...
                            combineStats=false);
plotTraces(combined.data{1},combined.timestamp,color=waterColor);

combined = combineTraces(summary,timeRange=[-0.5,3],...
                            eventRange='water (laser)',...
                            animalRange=animal,...
                            taskRange='combine',...
                            combineStats=false);
plotTraces(combined.data{1},combined.timestamp,color=addOpacity(waterColor,0.5));

plotEvent('Water',0,color=bluePurpleRed(1,:));
xlabel('Time (s)'); ylabel('z-score');
legend({'Water (ctrl)','Water (laser)'});


%% Show sample excitation and inhibition

sessionName = '20250918-W-BiPOLES-DAT-2_g0';
rootPath    = '/Volumes/MICROSCOPE/Shun/Project clamping/Recordings/202509-BiPOLES/';
dataPath    = fullfile(rootPath, sessionName, "data_" + sessionName + ".mat");
load(dataPath);

% Low pass DA
Fs = 10000;      % sampling rate (Hz)
cutoff = 200;        % cutoff (Hz)
[b,a] = butter(4, cutoff/(Fs/2), 'low');           % 4th-order low-pass
dopamine_filtered = filtfilt(b, a, double(photometry_raw));

% Down sample
targetFs = 10000; 
if targetFs == Fs
    dopamine_finalFs = dopamine_filtered;
    red_laser = redClamp;
    blue_laser = blueClamp;
else
    downsample_factor = max(1, floor(Fs/targetFs));
    red_laser = redClamp(1:downsample_factor:end);
    blue_laser = blueClamp(1:downsample_factor:end);
    dopamine_finalFs = dopamine_filtered(1:downsample_factor:end);
end

% Z score DA
windowSizeSec = 180;
window_size = round(windowSizeSec * targetFs); 
DA = rollingZ(dopamine_finalFs,window_size);

% Threshold laser signals (0/1)
red_laser  = double(red_laser  > 0.5);
blue_laser = double(blue_laser > 0.5);

%% Plot sample window
data_window_sec = [180 220];                          % seconds
i1 = floor(data_window_sec(1)*targetFs) + 1;          % inclusive start (MATLAB)
i2 = floor(data_window_sec(2)*targetFs);              % exclusive end in Python -> here end
idx = i1:i2;
t = linspace(0,length(idx)/targetFs,length(idx));

close all;
initializeFig(1,.5); tiledlayout(1,1,'Padding','loose');
nexttile; box off;

yyaxis left
plot(t, red_laser(idx), color=addOpacity(bluePurpleRed(500,:),0.3), ...
    LineStyle='-', LineWidth=0.1, DisplayName='Red Laser'); hold on
plot(t, blue_laser(idx), color=addOpacity(bluePurpleRed(1,:),0.3), ...
    LineStyle='-', LineWidth=0.1, DisplayName='Blue Laser');
xlim([0,range(t)]); ylabel('laser')

yyaxis right
plot(t, DA(idx), color=[.2 .8 .3], LineWidth=3, DisplayName='Dopamine')
xlim([0,range(t)]); xlabel('time (s)'); 
ylabel('dopamine');


%% Show sample collision test with water

sessionName = '20250918-W-BiPOLES-DAT-2-combine_g0';
rootPath    = '/Volumes/MICROSCOPE/Shun/Project clamping/Recordings/202509-BiPOLES/';
dataPath    = fullfile(rootPath, sessionName, "data_" + sessionName + ".mat");
load(dataPath);

% Low pass DA
Fs = 10000;      % sampling rate (Hz)
cutoff = 200;        % cutoff (Hz)
[b,a] = butter(4, cutoff/(Fs/2), 'low');           % 4th-order low-pass
dopamine_filtered = filtfilt(b, a, double(photometry_raw));

% Down sample
targetFs = 10000; 
if targetFs == Fs
    dopamine_finalFs = dopamine_filtered;
    red_laser = redClamp;
    blue_laser = blueClamp;
else
    downsample_factor = max(1, floor(Fs/targetFs));
    red_laser = redClamp(1:downsample_factor:end);
    blue_laser = blueClamp(1:downsample_factor:end);
    dopamine_finalFs = dopamine_filtered(1:downsample_factor:end);
end

% Z score DA
windowSizeSec = 180;
window_size = round(windowSizeSec * targetFs); 
DA = rollingZ(dopamine_finalFs,window_size);

% Threshold laser signals (0/1)
red_laser  = double(red_laser  > 0.5);
blue_laser = double(blue_laser > 0.5);
lick       = double(rightLick > 0.5);
water      = double(rightSolenoid > 0.5);
lick  = find(diff([0, lick])  == 1);
water = find(diff([0, water]) == 1);

%% Plot sample window
data_window_sec = [270 310];                          % seconds
i1 = floor(data_window_sec(1)*targetFs) + 1;          % inclusive start (MATLAB)
i2 = floor(data_window_sec(2)*targetFs);              % exclusive end in Python -> here end
idx = i1:i2;
t = linspace(0,length(idx)/targetFs,length(idx));

close all;
initializeFig(1,.5); tiledlayout(1,1,'Padding','compact');
nexttile; box off;

yyaxis left
% Plot licking
lickON = lick(lick>=idx(1) & lick<=idx(end)) - idx(1) + 1;           % onsets within window
lickON = lickON / targetFs;
scatter(lickON, ones(size(lickON)), 50, bluePurpleRed(500,:), ...
        'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.8, ...
        'DisplayName','Lick onset'); hold on
% Plot water (not working for now)
waterON = water(water>=idx(1) & water<=idx(end)) - idx(1) + 1;   
waterON = waterON / targetFs;
scatter(waterON, ones(size(waterON)), 50, bluePurpleRed(1,:), ...
        'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.8, ...
        'DisplayName','Lick onset'); hold on
% Plot blue laser
plot(t, blue_laser(idx), color=addOpacity(bluePurpleRed(1,:),0.3), ...
    LineStyle='-', LineWidth=0.1, DisplayName='Blue Laser');
xlim([0,range(t)]); ylabel('Laser');

yyaxis right
plot(t, DA(idx), color=[.2 .8 .3], LineWidth=3, DisplayName='Dopamine')
xlim([0,range(t)]); 
xlabel('time (s)'); ylabel('Dopamine');



