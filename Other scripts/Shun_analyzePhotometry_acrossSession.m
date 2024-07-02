%% Load multiple sessionList

clear; close all;
% addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
[~,~,~,blueGreenYellow,blueWhiteRed,blueDarkPurpleRed,bluePurpleRed] = loadColors;

sessionList = uipickfiles('FilterSpec','\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');

resultspath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\project valence\Recordings\Results\';
% resultspath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\RCL\Results\';

% For pairing, sessions needs arrange in days then in animal number

%% Read files for opto pair sessions (have trial table)
paramsList = cell(size(sessionList));

photometryLJs = cell(size(sessionList));
photometryNIs = cell(size(sessionList));
trialTables = cell(size(sessionList));
randomShutter = cell(size(sessionList));
rightLickON = cell(size(sessionList));
timeNIs = cell(size(sessionList));
timePhotometrys = cell(size(sessionList));

% Loading all data into cell
for i = 1:length(sessionList)
    dirsplit = strsplit(sessionList{i},filesep); 
    sessionName = dirsplit{end};
    clear dirsplit

    load(strcat(sessionList{i},'\','sync_',sessionName,'.mat'),...
        'trials','rollingGreenLP','blueLaser','params','rightLick',...
        'timeNI','timePhotometry',...
        'rightSolenoid','redLaser','airpuff');
 
    paramsList{i} = params;

    photometryLJs{i} = rollingGreenLP;
    photometryNIs{i} = photometryNI;
    trialTables{i} = trials;
    randomShutter{i} = find(blueLaser);
    rightLickON{i} = rightLick;

    timeNIs{i} = timeNI;
    timePhotometrys{i} = timePhotometry;

    disp(['Session ',sessionName,' loaded']);
end

%% Find dLight traces around stim during stim only trials

timeRange = [-0.5,3];
binSize = params.finalTimeStep; groupSize = 20; % num of trials to plot in one line

eventIdx = cell(size(sessionList)); eventInLJ = cell(size(sessionList));
baselineInLJ = cell(size(sessionList));

traces = cell(size(sessionList)); baseline = cell(size(sessionList));

for i = 1:length(sessionList)
    % Find time points
    eventIdx{i} = trialTables{i}{trialTables{i}.isTone == 0 & trialTables{i}.isStim == 1,"CueTime"};
    eventInLJ{i} = zeros(size(eventIdx{i}));
    eventInLJ{i} = findCorrespondingTime(eventIdx{i},timeNIs{i},timePhotometrys{i});
    baselineInLJ{i} = findCorrespondingTime(randomShutter{i}',timeNIs{i},timePhotometrys{i});

    % Find traces
    [traces{i},t] = getTraces(eventInLJ{i}/params.finalFs,photometryLJs{i},timeRange,binSize);
    [baseline{i},~] = getTraces(baselineInLJ{i}/params.finalFs,photometryLJs{i},timeRange,binSize);
    disp(i);
end

%% Combine dataset 

% Opt 2: organized by mouse
% Cons: strictly defining trials means not taking account thirst state
% sl043Idx = [1,5,9,13]; sl043_traces = traces(1,sl043Idx); sl043_allTraces = vertcat(sl043_traces{:});
% sl044Idx = sl043Idx + 1; sl044_traces = traces(1,sl044Idx); sl044_allTraces = vertcat(sl044_traces{:});
% sl045Idx = sl043Idx + 2; sl045_traces = traces(1,sl045Idx); sl045_allTraces = vertcat(sl045_traces{:});
% sl046Idx = sl043Idx + 3; sl046_traces = traces(1,sl046Idx); sl046_allTraces = vertcat(sl046_traces{:});
% 
% trialRange = 1:30;
% trials1to30 = [sl043_allTraces(trialRange,:); sl044_allTraces(trialRange,:);...
%                sl045_allTraces(trialRange,:); sl046_allTraces(trialRange,:)];
% trialRange = 50:80;
% trials50to80 = [sl043_allTraces(trialRange,:); sl044_allTraces(trialRange,:);...
%                sl045_allTraces(trialRange,:); sl046_allTraces(trialRange,:)];

%% dLight traces during stim only trials across sessions (SL043) (reward only: first three traces)

initializeFig(.5,.5);
label = 'Stim'; eventDuration = 0.5;

% tiledlayout(1,2);

% For plotting first three traces
% nexttile;
plotSEM(t,baseline{1}(600:1100,:),[.75 .75 .75]);
plotSEM(t,traces{1}(1:30,:),blueWhiteRed(150,:));
plotSEM(t,traces{1}(50:80,:),blueWhiteRed(100,:));
plotSEM(t,traces{2}(1:30,:),blueWhiteRed(1,:));
ylim([-1.5,2]);
plotEvent(label,eventDuration,'r');
xlabel('Time (s)'); ylabel('z-score'); 
% legend({'Random shutter (n=500)',...
%     'D1, EP stim only #1-30'},...
%     'Location','northeast'); 
legend({'Random shutter (n=500)',...
    'D1, EP stim only #1-30','D1, EP stim only #50-80',...
    'D2, EP stim only #1-30','D3, EP stim only #1-30'},...
    'Location','northeast'); 

% nexttile;
% plotLicks(eventIdx{1}(1:30,:),timeRange,0.2,blueWhiteRed(500,:),[],rightLickON{1},paramsList{1});
% plotLicks(eventIdx{2}(1:30,:),timeRange,0.2,blueWhiteRed(100,:),[],rightLickON{1},paramsList{1});
% plotLicks(eventIdx{2}(1:30,:),timeRange,0.2,blueWhiteRed(1,:),[],rightLickON{2},paramsList{1});

% saveas(gcf,strcat(sessionList{1}.projectPath,'\psth_stimAcrossSession_3traces.png'));
% saveas(gcf,strcat('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\project valence\Recordings\Photometry-flexibleLearning-NAc-dLight\psth_stimAcrossSession_1traces.png'));


%% dLight traces during stim only trials across sessions (SL043) (reward only: all traces)

initializeFig(.5,.5);
label = 'Stim'; eventDuration = 0.5;

% SL043Idx = [1,5,9,13]; % traces_043 = traces(1,SL043Idx);

plotSEM(t,baseline{1}(600:1100,:),[.75 .75 .75]);
plotSEM(t,traces{1}(1:30,:),blueWhiteRed(210,:));
plotSEM(t,traces{1}(50:80,:),blueWhiteRed(210,:));
plotSEM(t,traces{2}(1:30,:),blueWhiteRed(210,:));
plotSEM(t,traces{3}(1:30,:),blueWhiteRed(90,:));
% plotSEM(t,traces{8}(1:30,:),blueWhiteRed(100,:));
plotSEM(t,traces{9}(1:30,:),blueWhiteRed(1,:));
ylim([-1.5,2]);
plotEvent(label,eventDuration,'r');
xlabel('Time (s)'); ylabel('z-score'); 

legend({'Random shutter (n=500)',...
    'S1, EP stim only #1-30','S1, EP stim only #50-80',...
    'S2, EP stim only #1-30','S3, EP stim only #1-30','S9, EP stim only #1-30'},...
    'Location','northeast'); 

scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
saveFigures(gcf,'SL043-allReward',resultspath);

% saveas(gcf,strcat(sessionList{1}.projectPath,'\psth_stimAcrossSession_3traces.png'));
% saveas(gcf,strcat(['\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\project valence\Recordings\Photometry-flexibleLearning-NAc-dLight\Results\' ...
%     'psth_stimAcrossSession_rewardOnly_5traces_background.png']));

%% SL043 second reward pairing

initializeFig(.5,.5);
label = 'Stim'; eventDuration = 0.5;

% SL043Idx = [1,5,9,13]; % traces_043 = traces(1,SL043Idx);

plotSEM(t,baseline{1}(600:1100,:),[.75 .75 .75]);
plotSEM(t,traces{13}(1:30,:),blueWhiteRed(500,:));
plotSEM(t,traces{14}(1:30,:),blueWhiteRed(1,:));
plotSEM(t,traces{15}(1:30,:),blueWhiteRed(100,:));
plotSEM(t,traces{16}(1:30,:),blueWhiteRed(150,:));
plotSEM(t,traces{17}(1:30,:),blueWhiteRed(200,:));
% plotSEM(t,traces{17}(1:30,:),blueWhiteRed(1,:));
ylim([-1.5,2]);
plotEvent(label,eventDuration,'r');
xlabel('Time (s)'); ylabel('z-score'); 

legend({'Random shutter (n=500)',...
    'S1, EP stim only #1-30','S2, EP stim only #1-30',...
    'S3, EP stim only #1-30','S4, EP stim only #1-30','S5, EP stim only #1-30'},...
    'Location','northeast'); 

scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
saveFigures(gcf,'SL043-secondPairing',resultspath);

%% dLight traces during stim only trials across sessions (SL043) (reward and punishment: first three traces)

initializeFig(.5,.5);
label = 'Stim'; eventDuration = 0.5;

% SL043Idx = [1,5,9,13]; % traces_043 = traces(1,SL043Idx);

plotSEM(t,baseline{1}(600:1100,:),[.75 .75 .75]);
plotSEM(t,traces{1}(1:30,:),blueWhiteRed(150,:));
plotSEM(t,traces{1}(50:80,:),blueWhiteRed(100,:));
plotSEM(t,traces{2}(1:30,:),blueWhiteRed(1,:));
plotSEM(t,traces{10}(1:30,:),blueWhiteRed(325,:));
plotSEM(t,traces{10}(80:110,:),blueWhiteRed(400,:));
plotSEM(t,traces{11}(1:30,:),blueWhiteRed(500,:));
% plotSEM(t,traces{12}(1:30,:),blueWhiteRed(1,:));
ylim([-2,2]);
plotEvent(label,eventDuration,'r');
xlabel('Time (s)'); ylabel('z-score'); 

legend({'Random shutter (n=500)',...
    'Reward D1, EP stim only #1-30','Reward D1, EP stim only #50-80',...
    'Reward D2, EP stim only #1-30',...
    'Punishment D1, EP stim only #1-30','Punishment D1, EP stim only #80-110',...
    'Punishment D2, EP stim only #1-30'},...
    'Location','northeast'); 

scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');

saveFigures(gcf,'SL043-optoPair',resultspath);

%% dLight traces during stim only trials across sessions across sessions (all mice)

% Opt 1: organize by day/session
% d1_traces = traces(1,1:4); d2_traces = traces(1,5:8);
% d3_traces = traces(1,9:12); d4_traces = traces(1,13:16);
% 
% d1_first_traces = cellfun(@(x) x(1:30,:),d1_traces,'UniformOutput',false);
% % d1_last_traces = cellfun(@(x) x(end-40:end-11,:),d1_traces,'UniformOutput',false);
% d2_first_traces = cellfun(@(x) x(1:30,:),d2_traces,'UniformOutput',false);
% d3_first_traces = cellfun(@(x) x(1:30,:),d3_traces,'UniformOutput',false);
% d4_first_traces = cellfun(@(x) x(1:30,:),d4_traces,'UniformOutput',false);
% 
% d1_allTraces = vertcat(d1_first_traces{:}); d2_allTraces = vertcat(d2_first_traces{:});
% d3_allTraces = vertcat(d3_first_traces{:}); d4_allTraces = vertcat(d4_first_traces{:});

d1_allTraces = combineTraces(traces,sessionRange=1:4,trialRange=1:30);
d2_allTraces = combineTraces(traces,sessionRange=5:8,trialRange=1:30);
d3_allTraces = combineTraces(traces,sessionRange=9:12,trialRange=1:30);
d4_allTraces = combineTraces(traces,sessionRange=13:16,trialRange=1:30);
d5_allTraces = combineTraces(traces,sessionRange=17:20,trialRange=1:30);
d6_allTraces = combineTraces(traces,sessionRange=21:24,trialRange=1:30);

% Plot traces
initializeFig(0.5,0.5);
label = 'Stim'; eventDuration = 0.5;

plotSEM(t,baseline{1}(600:1100,:),[.75 .75 .75]);
plotSEM(t,d1_allTraces,blueWhiteRed(150,:));
plotSEM(t,d2_allTraces,blueWhiteRed(100,:));
plotSEM(t,d3_allTraces,blueWhiteRed(1,:));

plotSEM(t,d4_allTraces,blueWhiteRed(325,:));
plotSEM(t,d5_allTraces,blueWhiteRed(400,:));
plotSEM(t,d6_allTraces,blueWhiteRed(500,:));
ylim([-1.5,2]);
plotEvent(label,eventDuration,'r');
xlabel('Time (s)'); ylabel('z-score'); 

legend({'Random shutter (n=500 trials)',...
    'Reward D1, EP stim only #1-30 (4 mice)','Reward D2, EP stim only #1-30 (4 mice)',...
    'Reward D3, EP stim only #1-30 (4 mice)','Airpuff D1, EP stim only #1-30 (4 mice)',...
    'Airpuff D2, EP stim only #1-30 (4 mice)','Airpuff D3, EP stim only #1-30 (4 mice)'},...
    'Location','northeast'); 

% saveas(gcf,strcat('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\project valence\Recordings\Photometry-flexibleLearning-NAc-dLight\psth_stimAcrossSession_3traces.png'));
saveFigures(gcf,'all-optoPair',resultspath,resolution=3000);

%% Compared area under curve (sum during time window) of baseline vs stim period across sessions

preStimRange = find(t<0); 
duringStimRange = preStimRange(end)+1:preStimRange(end)+length(preStimRange);
% d1_preStimSum = cellfun(@(x) sum(x(:,preStimRange),2),d1_first_traces,'UniformOutput',false);
% d2_preStimSum = cellfun(@(x) sum(x(:,preStimRange),2),d2_first_traces,'UniformOutput',false);
% d3_preStimSum = cellfun(@(x) sum(x(:,preStimRange),2),d3_first_traces,'UniformOutput',false);
% d1_duringStimSum = cellfun(@(x) sum(x(:,duringStimRange),2),d1_first_traces,'UniformOutput',false);
% d2_duringStimSum = cellfun(@(x) sum(x(:,duringStimRange),2),d2_first_traces,'UniformOutput',false);
% d3_duringStimSum = cellfun(@(x) sum(x(:,duringStimRange),2),d3_first_traces,'UniformOutput',false);
% 
% d1_preStimSum = vertcat(d1_preStimSum{:});
% d2_preStimSum = vertcat(d2_preStimSum{:});
% d3_preStimSum = vertcat(d3_preStimSum{:});
% d1_duringStimSum = vertcat(d1_duringStimSum{:});
% d2_duringStimSum = vertcat(d2_duringStimSum{:});
% d3_duringStimSum = vertcat(d3_duringStimSum{:});

d1_preStimSum = sum(d1_allTraces(:,preStimRange),2);
d2_preStimSum = sum(d2_allTraces(:,preStimRange),2);
d3_preStimSum = sum(d3_allTraces(:,preStimRange),2);
d4_preStimSum = sum(d4_allTraces(:,preStimRange),2);
d5_preStimSum = sum(d5_allTraces(:,preStimRange),2);
d6_preStimSum = sum(d6_allTraces(:,preStimRange),2);

d1_duringStimSum = sum(d1_allTraces(:,duringStimRange),2);
d2_duringStimSum = sum(d2_allTraces(:,duringStimRange),2);
d3_duringStimSum = sum(d3_allTraces(:,duringStimRange),2);
d4_duringStimSum = sum(d4_allTraces(:,duringStimRange),2);
d5_duringStimSum = sum(d5_allTraces(:,duringStimRange),2);
d6_duringStimSum = sum(d6_allTraces(:,duringStimRange),2);

% T-test
[~,d1_test_p,~] = kstest2(d1_preStimSum,d1_duringStimSum);
[~,d2_test_p,~] = kstest2(d2_preStimSum,d2_duringStimSum);
[~,d3_test_p,~] = kstest2(d3_preStimSum,d3_duringStimSum);
[~,d4_test_p,~] = kstest2(d4_preStimSum,d4_duringStimSum);
[~,d5_test_p,~] = kstest2(d5_preStimSum,d5_duringStimSum);
[~,d6_test_p,~] = kstest2(d6_preStimSum,d6_duringStimSum);

%% Plot bar plot (AUC prestim vs AUC during stim)
initializeFig(0.5,0.5);

x_baseline = [2*ones(120,1), 6*ones(120,1), 10*ones(120,1),...
              14*ones(120,1), 18*ones(120,1), 22*ones(120,1)];
y_baseline = [d1_preStimSum,d2_preStimSum,d3_preStimSum,...
              d4_preStimSum,d5_preStimSum,d6_preStimSum];
x_d1stim = 3*ones(120,1); x_d2stim = 7*ones(120,1); x_d3stim = 11*ones(120,1);
x_d4stim = 15*ones(120,1); x_d5stim = 19*ones(120,1); x_d6stim = 23*ones(120,1);

swarmchart(x_baseline,y_baseline,[],[.75 .75 .75],'filled'); hold on
swarmchart(x_d1stim,d1_duringStimSum,[],blueWhiteRed(150,:),'filled'); hold on
swarmchart(x_d2stim,d2_duringStimSum,[],blueWhiteRed(100,:),'filled'); hold on
swarmchart(x_d3stim,d3_duringStimSum,[],blueWhiteRed(1,:),'filled'); hold on
swarmchart(x_d4stim,d4_duringStimSum,[],blueWhiteRed(325,:),'filled'); hold on
swarmchart(x_d5stim,d5_duringStimSum,[],blueWhiteRed(400,:),'filled'); hold on
swarmchart(x_d6stim,d6_duringStimSum,[],blueWhiteRed(500,:),'filled'); hold on

barplot = gca;
barplot.XTick = [2.5 6.5 10.5 14.5 18.5 22.5];
barplot.XTickLabel = ["Reward D1", "Reward D2", "Reward D3",...
                      "Airpuff D1", "Airpuff D2", "Airpuff D3"];

ylabel("sum(z-score)");

%% Difference between AUC plotted across days (one animal per trace)

sl043_d1_diff = d1_duringStimSum(1:30) - d1_preStimSum(1:30);
sl044_d1_diff = d1_duringStimSum(31:60) - d1_preStimSum(31:60);
sl045_d1_diff = d1_duringStimSum(61:90) - d1_preStimSum(61:90);
sl046_d1_diff = d1_duringStimSum(91:120) - d1_preStimSum(91:120);

sl043_d2_diff = d2_duringStimSum(1:30) - d2_preStimSum(1:30);
sl044_d2_diff = d2_duringStimSum(31:60) - d2_preStimSum(31:60);
sl045_d2_diff = d2_duringStimSum(61:90) - d2_preStimSum(61:90);
sl046_d2_diff = d2_duringStimSum(91:120) - d2_preStimSum(91:120);

sl043_d3_diff = d3_duringStimSum(1:30) - d3_preStimSum(1:30);
sl044_d3_diff = d3_duringStimSum(31:60) - d3_preStimSum(31:60);
sl045_d3_diff = d3_duringStimSum(61:90) - d3_preStimSum(61:90);
sl046_d3_diff = d3_duringStimSum(91:120) - d3_preStimSum(91:120);

sl043_d4_diff = d4_duringStimSum(1:30) - d4_preStimSum(1:30);
sl044_d4_diff = d4_duringStimSum(31:60) - d4_preStimSum(31:60);
sl045_d4_diff = d4_duringStimSum(61:90) - d4_preStimSum(61:90);
sl046_d4_diff = d4_duringStimSum(91:120) - d4_preStimSum(91:120);

sl043_d5_diff = d5_duringStimSum(1:30) - d5_preStimSum(1:30);
sl044_d5_diff = d5_duringStimSum(31:60) - d5_preStimSum(31:60);
sl045_d5_diff = d5_duringStimSum(61:90) - d5_preStimSum(61:90);
sl046_d5_diff = d5_duringStimSum(91:120) - d5_preStimSum(91:120);

sl043_d6_diff = d6_duringStimSum(1:30) - d6_preStimSum(1:30);
sl044_d6_diff = d6_duringStimSum(31:60) - d6_preStimSum(31:60);
sl045_d6_diff = d6_duringStimSum(61:90) - d6_preStimSum(61:90);
sl046_d6_diff = d6_duringStimSum(91:120) - d6_preStimSum(91:120);

%% Plot line plot

initializeFig(0.5,0.5);
sl043_diff = [sl043_d1_diff,sl043_d2_diff,sl043_d3_diff,...
                sl043_d4_diff,sl043_d5_diff,sl043_d6_diff];
sl044_diff = [sl044_d1_diff,sl044_d2_diff,sl044_d3_diff,...
                sl044_d4_diff,sl044_d5_diff,sl044_d6_diff];
sl045_diff = [sl045_d1_diff,sl045_d2_diff,sl045_d3_diff,...
                sl045_d4_diff,sl045_d5_diff,sl045_d6_diff];
sl046_diff = [sl046_d1_diff,sl046_d2_diff,sl046_d3_diff,...
                sl046_d4_diff,sl046_d5_diff,sl046_d6_diff];

plotSEM(1:6,sl043_diff,blueWhiteRed(1,:));
plotSEM(1:6,sl044_diff,blueWhiteRed(100,:));
plotSEM(1:6,sl045_diff,blueWhiteRed(400,:));
plotSEM(1:6,sl046_diff,blueWhiteRed(500,:));

gca.XTick = [1 2 3 4 5 6];
xlabel('Session'); ylabel('\Delta AUC')

%% %%%%%%%%%%%%%%%%%%%%%%%%% Pre pairing %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read files for opto pair sessions (no trial table)
paramsList = cell(size(sessionList));

photometryLJs = cell(size(sessionList));
randomShutter = cell(size(sessionList));
timeNIs = cell(size(sessionList));
timePhotometrys = cell(size(sessionList));

rightSolenoids = cell(size(sessionList));
redLasers = cell(size(sessionList));
airpuffs = cell(size(sessionList));

trialTables = cell(size(sessionList));
stimOnlys = cell(size(sessionList));

% Loading all data into cell
for i = 1:length(sessionList)
    dirsplit = strsplit(sessionList{i},filesep); 
    sessionName = dirsplit{end};
    clear dirsplit

    trials = [];
    
    load(strcat(sessionList{i},'\','sync_',sessionName,'.mat'),...
        'rollingGreenLP','blueLaser','params','trials',...
        'timeNI','timePhotometry',...
        'rightSolenoid','redLaser','airpuff');

    paramsList{i} = params;

    photometryLJs{i} = rollingGreenLP;
    randomShutter{i} = find(blueLaser);

    timeNIs{i} = timeNI;
    timePhotometrys{i} = timePhotometry;

    rightSolenoids{i} = find(rightSolenoid);
    redLasers{i} = find(redLaser);
    airpuffs{i} = find(airpuff);

    if isempty(trials); stimOnlys{i} = find(redLaser);
    else
        trialTables{i} = trials;
        stimOnlys{i} = trials{trials.isReward == 0 & trials.isStim == 1,"CueTime"};
    end

    disp(['Session ',sessionName,' loaded']);
end
%% Find dLight traces around water/laser/airpuff

timeRange = [-0.5,3];
binSize = params.finalTimeStep;

% waterTraces = cell(size(sessionList));
% airpuffTraces = cell(size(sessionList));
stimOnlyTraces = cell(size(sessionList));
baselineInLJ = cell(size(sessionList)); baseline = cell(size(sessionList));

for i = 1:length(sessionList)
    eventIdx = cell(size(sessionList)); eventInLJ = cell(size(sessionList));
    % Find time points
    %eventIdx{i} = trialTables{i}{trialTables{i}.isReward == 0 & trialTables{i}.isStim == 1,"CueTime"};
    %eventIdx{i} = airpuffs{i};
    eventIdx{i} = stimOnlys{i};

    eventInLJ{i} = findCorrespondingTime(eventIdx{i},timeNIs{i},timePhotometrys{i});
    % baselineInLJ{i} = findCorrespondingTime(randomShutter{i}',timeNIs{i},timePhotometrys{i});

    % Find traces
    [stimOnlyTraces{i},t] = getTraces(eventInLJ{i}/params.finalFs,photometryLJs{i},timeRange,binSize);
    %[baseline{i},~] = getTraces(baselineInLJ{i}/params.finalFs,photometryLJs{i},timeRange,binSize);

    disp(i);
end


%% Combine event traces into one

waterCombined = combineTraces(waterTraces);
stimCombined = combineTraces(stimOnlyTraces,sessionRange=[1 2 3 4 5 7 10 38:43]);
airpuffCombined = combineTraces(airpuffTraces,sessionRange=[11 12 38:43]);
baselineCombined = combineTraces(baseline);

%% dLight response of water vs airpuff vs baseline

initializeFig(0.5,0.5);

plotSEM(t,waterCombined,blueWhiteRed(1,:));
plotSEM(t,stimCombined,blueWhiteRed(500,:));
plotSEM(t,airpuffCombined,[.2 .2 .2]);
% plotEvent('',0,'r');
xlabel('Time (s)'); ylabel('z-score'); 
legend({['Water (n=',num2str(size(waterCombined,1)),', 6 mice)'],...
    ['Stim (n=',num2str(size(stimCombined,1)),', 4 mice)'],...
    ['Airpuff (n=',num2str(size(airpuffCombined,1)),', 3 mice)']},...
    'Location','best'); 

scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');

% saveFigures(gcf,'psth_reward&airpuff&stim',resultspath);


%% dLight response of tone across sessions across mice

%% %%%%%%%%%%%%%%%%%%%%%% EP stim + LHb iGluSnFR %%%%%%%%%%%%%%%%%%%%%%%
%% Read files
paramsList = cell(size(sessionList));

lj_photometry = cell(size(sessionList));
trialTables = cell(size(sessionList));
randomShutter = cell(size(sessionList));
rightLickON = cell(size(sessionList));
timeNIs = cell(size(sessionList));
timePhotometrys = cell(size(sessionList));
ni_photometry = cell(size(sessionList));

% Loading all data into cell
for i = 1:length(sessionList)
    dirsplit = strsplit(sessionList{i},filesep); 
    sessionName = dirsplit{end};
    clear dirsplit

    load(strcat(sessionList{i},'\','sync_',sessionName,'.mat'),...
        'trials','rollingGreenLP','blueLaser','params','rightLick',...
        'timeNI','timePhotometry','photometryNI',...
        'rightSolenoid','redLaser','airpuff');

    paramsList{i} = params;

    lj_photometry{i} = rollingGreenLP;
    trialTables{i} = trials;
    randomShutter{i} = find(blueLaser);
    rightLickON{i} = rightLick;
    ni_photometry{i} = photometryNI;

    timeNIs{i} = timeNI;
    timePhotometrys{i} = timePhotometry;

    disp(['Session ',sessionName,' loaded']);
end

%% Extract stim only traces
timeRange = [-0.5,3];
binSize = params.finalTimeStep;

eventIdx = cell(size(sessionList)); eventInLJ = cell(size(sessionList));

LJtraces = cell(size(sessionList)); 
NItraces = cell(size(sessionList)); 

for i = 1:length(sessionList)
    % Find time points
    eventIdx{i} = trialTables{i}{trialTables{i}.isTone == 0 & trialTables{i}.isStim == 1,"CueTime"};
    eventInLJ{i} = zeros(size(eventIdx{i}));
    eventInLJ{i} = findCorrespondingTime(eventIdx{i},timeNIs{i},timePhotometrys{i});

    % Find traces
    [LJtraces{i},t_lj] = getTraces(eventInLJ{i}/params.finalFs,lj_photometry{i},timeRange,binSize);
    [NItraces{i},t_ni] = getTraces(eventIdx{i}/params.sync.behaviorFs,ni_photometry{i},timeRange,1/50);
    disp(i);
end

%% Combine traces

allNI = combineTraces(NItraces,merge=true);

initializeFig(0.8,0.5);

tiledlayout(1,2); nexttile;
plotSEM(t_ni,allNI(1:30,:),blueWhiteRed(1,:));
plotSEM(t_ni,allNI(31:60,:),blueWhiteRed(100,:));
plotSEM(t_ni,allNI(71:100,:),blueWhiteRed(300,:));
plotSEM(t_ni,allNI(111:140,:),blueWhiteRed(400,:));
plotSEM(t_ni,allNI(148:178,:),blueWhiteRed(500,:));
plotEvent('Stim',0.5,'r');
xlabel('Time (s)'); ylabel('z-score'); 
legend({'Reward D1, EP stim only #1-30','Reward D2, EP stim only #1-30',...
    'Reward D3, EP stim only #1-30','Reward D4, EP stim only #1-30',...
    'Reward D5, EP stim only #1-30'},...
    'Location','northeast'); 

nexttile;
plotLicks(eventIdx{1}(1:30,:),timeRange,0.2,blueWhiteRed(1,:),[],rightLickON{1},paramsList{1});
plotLicks(eventIdx{2}(1:30,:),timeRange,0.2,blueWhiteRed(100,:),[],rightLickON{2},paramsList{1});
plotLicks(eventIdx{3}(1:30,:),timeRange,0.2,blueWhiteRed(300,:),[],rightLickON{3},paramsList{1});
plotLicks(eventIdx{4}(1:30,:),timeRange,0.2,blueWhiteRed(400,:),[],rightLickON{4},paramsList{1});
plotLicks(eventIdx{5}(1:30,:),timeRange,0.2,blueWhiteRed(500,:),[],rightLickON{5},paramsList{1});
plotEvent('Stim',0.5,'r');
xlabel('Time (s)'); ylabel('Licks/s'); 

saveFigures(gcf,'stimOnly-optoPair',resultspath,resolution=3000);

%% First two day subtraces

initializeFig(0.8,0.5);

tiledlayout(1,2);
nexttile;
plotSEM(t_ni,allNI(1:10,:),blueWhiteRed(1,:));
plotSEM(t_ni,allNI(11:20,:),blueWhiteRed(100,:));
plotSEM(t_ni,allNI(21:30,:),blueWhiteRed(200,:));

plotSEM(t_ni,allNI(31:40,:),blueWhiteRed(300,:));
plotSEM(t_ni,allNI(41:50,:),blueWhiteRed(400,:));
plotSEM(t_ni,allNI(51:60,:),blueWhiteRed(500,:));
plotEvent('Stim',0.5,'r');
xlabel('Time (s)'); ylabel('z-score'); 
legend({'Reward D1, EP stim only #1-10','Reward D1, EP stim only #11-20',...
    'Reward D1, EP stim only #21-30','Reward D2, EP stim only #1-10',...
    'Reward D2, EP stim only #11-20','Reward D2, EP stim only #21-30'},...
    'Location','northeast'); 

nexttile;
plotLicks(eventIdx{1}(1:10,:),timeRange,0.2,blueWhiteRed(1,:),[],rightLickON{1},paramsList{1});
plotLicks(eventIdx{1}(11:20,:),timeRange,0.2,blueWhiteRed(100,:),[],rightLickON{1},paramsList{1});
plotLicks(eventIdx{1}(21:30,:),timeRange,0.2,blueWhiteRed(200,:),[],rightLickON{1},paramsList{1});

plotLicks(eventIdx{2}(1:10,:),timeRange,0.2,blueWhiteRed(300,:),[],rightLickON{2},paramsList{1});
plotLicks(eventIdx{2}(11:20,:),timeRange,0.2,blueWhiteRed(400,:),[],rightLickON{2},paramsList{1});
plotLicks(eventIdx{2}(21:30,:),timeRange,0.2,blueWhiteRed(500,:),[],rightLickON{2},paramsList{1});
plotEvent('Stim',0.5,'r');
xlabel('Time (s)'); ylabel('Licks/s'); 
% legend({'Reward D1, EP stim only #1-10','Reward D1, EP stim only #11-20',...
%     'Reward D1, EP stim only #21-30','Reward D2, EP stim only #1-10',...
%     'Reward D2, EP stim only #11-20','Reward D2, EP stim only #21-30'},...
%     'Location','northeast'); 

saveFigures(gcf,'stimOnly-firstTwoSession',resultspath,resolution=3000);

