% Shun_DAClamping_analysis

% 2024/10/10

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

%% (Optional) Load sessions

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'))';
errorSessionIdx = [];

% Select anlaysis params
[analysisParams,canceled] = inputAnalysisParams(sessionList,...
                                reloadAll=false,...
                                recordLJ='[0 0 0]',...
                                plotPhotometry=false,...
                                plotBehavior=true,...
                                withPhotometryNI=true);
if canceled; return; end
for s = 1:length(sessionList)
    analysisParams(s).rollingWindowTime = str2double(analysisParams(s).rollingWindowTime);
    analysisParams(s).recordLJ = eval(analysisParams(s).recordLJ);
end 
% Select session params
[sessionParams,canceled] = inputSessionParams(sessionList,...
                                paradigm=2,...
                                redStim=true,pavlovian=false,...
                                reactionTime=2,...
                                redPulseFreq=50,redPulseDuration=5,redStimDuration=500,...
                                bluePulseFreq=30,bluePulseDuration=10,blueStimDuration=500,...
                                includeOtherStim=true);
if canceled; return; end
taskList = cell(size(sessionList));
taskOptions = {'random','reward pairing','punish pairing'};
redStimPatternList = cell(size(sessionList));
blueStimPatternList = cell(size(sessionList));
for s = 1:length(sessionList)
    taskList{s} = taskOptions{sessionParams(s).Paradigm};
    redStimPatternList{s} = {sessionParams(s).RedPulseFreq,sessionParams(s).RedPulseDuration,sessionParams(s).RedStimDuration};
    blueStimPatternList{s} = {sessionParams(s).BluePulseFreq,sessionParams(s).BluePulseDuration,sessionParams(s).BlueStimDuration};
    sessionParams(s).ReactionTime = str2double(sessionParams(s).ReactionTime);
    sessionParams(s).minLicks = str2double(sessionParams(s).minLicks);
end


% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList errorSessionIdx analysisParams sessionParams taskList redStimPatternList blueStimPatternList withPhotometryNI plotPhotometry reloadAll
    
    idx = s; % if loop through all
    % idx = errorSessionIdx(s); % if loop through error sessions

    dirsplit = strsplit(sessionList{idx},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{idx},reloadAll=analysisParams(idx).reloadAll,...
            invertStim=sessionParams(idx).OptoInverted,...
            withPhotometryNI=analysisParams(idx).withPhotometryNI,photometryNI_mod=false,...
            recordLJ=analysisParams(idx).recordLJ,...
            rollingWindowTime=analysisParams(idx).rollingWindowTime,...
            followOriginal=false,getConsecutive=false,...
            saveDigitalNI=true,saveDigitalNIChannelIdx=8);
        analyzeSessions_optoPair(sessionList{idx},...
            task=taskList{idx},...
            redStimPattern=redStimPatternList{idx},...
            blueStimPattern=blueStimPatternList{idx},...
            redStim=sessionParams(idx).redStim,...
            redo=true,round=false,performing=false,...
            analyzeTraces=true,...
            plotPhotometry=analysisParams(idx).plotPhotometry,...
            plotBehavior=analysisParams(idx).plotBehavior,...
            pavlovian=sessionParams(idx).Pavlovian,...
            reactionTime=sessionParams(idx).ReactionTime,...
            combineOmission=true,...
            shortTimeRange=[-1,5],longTimeRange=[-5,10],...
            includeOtherStim=sessionParams(idx).IncludeOtherStim);
    catch ME
        errorSessionIdx = [errorSessionIdx;idx];
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end
close all

%% Settings
originalFs = 10000;
Fs_5ms = 1 / 0.005;
detrendWindow = 180; % in sec

% Plot options
baselineColor = [.7 .7 .7];
pidColor = [.232 .76 .58];
blueColor = [130 189 252]./255;

%% Load baseline sessions

% Select sessions via uipickfiles
sessionPath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));
sessionpath = sessionPath{1};

% Load sessions
% Get animal name and session datec
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; clear dirsplit
load(strcat(sessionpath,filesep,'data_',sessionName,'.mat'));

baseline_raw = voltage2arduino(photometry_raw);
baseline_lp = lowpass(baseline_raw',500,originalFs);
% baseline_zscore = rollingZ(baseline_lp,detrendWindow*originalFs);

baseline_blueLaser = blueLaser;
baseline_airpuffIdx = find(airpuff);
baseline_rewardIdx = find(rightSolenoid);
baseline_lickON = find(rightLick);
baseline_tone = find(leftTone);

% Find water lick (first lick in response to water)
baseline_waterLickIdx = nan(size(baseline_rewardIdx));
for i = 1:length(baseline_rewardIdx)
    nextLick = baseline_lickON(find(baseline_lickON>=baseline_rewardIdx(i),1));
    if ~isempty(nextLick); baseline_waterLickIdx(i) = nextLick; end
end
baseline_waterLickIdx = rmmissing(baseline_waterLickIdx);

% Downsample photometry signal
if ~exist('baseline_5ms','var')
    downsampled = downsampleSignal(baseline_lp,targetFs=Fs_5ms,originalFs=originalFs,rollingZ=false);
    baseline_5ms = downsampled.dsData;
    save(strcat(sessionpath,filesep,'data_',sessionName,'.mat'),'baseline_5ms','-append');
    disp('Finished: downsampled and saved baseline raw photometry');
end

% Get duration of blueLaser
if ~exist('baseline_blueLaserON','var')
    baseline_blueLaserON = ~savedDigitalNI .* getConsecutive(~savedDigitalNI)./originalFs;
    save(strcat(sessionpath,filesep,'data_',sessionName,'.mat'),'baseline_blueLaserON','-append');
    disp('Finished: saved baseline blueLaserON');
end

%% Load PID sessions

% Select sessions via uipickfiles
sessionPath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));
sessionpath = sessionPath{1};

% Load sessions
% Get animal name and session date
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; clear dirsplit
load(strcat(sessionpath,filesep,'data_',sessionName,'.mat'));

pid_raw = voltage2arduino(photometry_raw);
pid_lp = lowpass(pid_raw',500,originalFs);
% pid_zscore = rollingZ(pid_lp,detrendWindow*originalFs);

pid_blueLaser = blueLaser;
pid_airpuffIdx = find(airpuff);
pid_rewardIdx = find(rightSolenoid);
pid_lickON = find(rightLick);
pid_tone = find(leftTone);

% Find water lick (first lick in response to water)
pid_waterLickIdx = nan(size(pid_rewardIdx));
for i = 1:length(pid_rewardIdx)
    nextLick = pid_lickON(find(pid_lickON>=pid_rewardIdx(i),1));
    if ~isempty(nextLick); pid_waterLickIdx(i) = nextLick; end
end
pid_waterLickIdx = rmmissing(pid_waterLickIdx);

% Downsample photometry signal
if ~exist('pid_5ms','var')
    downsampled = downsampleSignal(pid_lp,targetFs=Fs_5ms,originalFs=originalFs,rollingZ=false);
    pid_5ms = downsampled.dsData;
    save(strcat(sessionpath,filesep,'data_',sessionName,'.mat'),'pid_5ms','-append');
    disp('Finished: downsampled and saved PID raw photometry');
end

% Get duration of blueLaser
if ~exist('pid_blueLaserON','var')
    pid_blueLaserON = ~savedDigitalNI .* getConsecutive(~savedDigitalNI)./originalFs;
    save(strcat(sessionpath,filesep,'data_',sessionName,'.mat'),'pid_blueLaserON','-append');
    disp('Finished: saved pid blueLaserON');
end

%% Plot DA session summary trace

baseline_photometry = baseline_5ms;
pid_photometry = pid_5ms;
Fs = Fs_5ms;
fftWindow = 10;

initializeFig(1,1); tiledlayout(2,4);

nexttile(1,[1 2]);
% minLength = min(length(pid_photometry),length(baseline_photometry));
% t = (1:minLength)/Fs; plotWindow = 1:minLength;
timeToPlot = [0,0]; % in sec
% Plot photometry trace
if timeToPlot(2) <= 0
    plotWindow = 1 : min([length(baseline_photometry),length(pid_photometry)]);
else
    plotWindow = timeToPlot(1)*Fs : timeToPlot(2)*Fs;
end
t = (timeToPlot(1)*Fs : timeToPlot(2)*Fs) / Fs;
% yyaxis left;
plot(t,baseline_photometry(plotWindow),Color=baselineColor,LineWidth=3,LineStyle='-'); hold on;
plot(t,pid_photometry(plotWindow),Color=pidColor,LineWidth=3,LineStyle='-'); hold on;
% yyaxis right;
% % Plot blue laser
% plotWindow = timeToPlot(1)*originalFs : timeToPlot(2)*originalFs;
% t = (timeToPlot(1)*originalFs : timeToPlot(2)*originalFs) / originalFs;
% ylabel('Intensity');
% plot(t,pid_blueLaserON(plotWindow)*1000,Color=blueColor,LineWidth=3);
xlabel('Time (s)'); box off;
legend('Baseline','PID');


nexttile(4);
[baseline_fft,baseline_fftpower] = plotFFT(baseline_photometry,Fs=Fs,color=baselineColor); 
[pid_fft,pid_fftpower] = plotFFT(pid_photometry,Fs=Fs,color=pidColor); 
legend('Baseline','PID');
title('FFT');

nexttile(8);
yyaxis left; histogram(baseline_photometry,200,EdgeColor=baselineColor,FaceColor=baselineColor); hold on
yyaxis right; histogram(pid_photometry,200,EdgeColor=pidColor,FaceColor=pidColor); hold on
xlabel('Intensity'); legend({'Baseline','PID'}); box off;
title('dLight');
subtitle({strcat("Basline std: ",num2str(std(baseline_lp)),", PID std: ",num2str(std(pid_lp)))});


nexttile(5);
% yyaxis left
[baseline_rewardTraces,~] = plotTraces(baseline_waterLickIdx,[-1,5],baseline_lp,...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=baselineColor);
[pid_rewardTraces,~] = plotTraces(pid_waterLickIdx,[-1,5],pid_lp,...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=pidColor);
% yyaxis right
% [~,~] = plotTraces(pid_waterLickIdx,[-1,5],pid_blueLaserON*1000,...
%                      signalFs=originalFs,sameSystem=true,...
%                      plotIndividual=false,color=blueColor,...
%                      xlabel='Time (s)',ylabel='Stim duration (ms)');
legend({['Baseline (n=',num2str(size(baseline_rewardTraces,1)),')'],...
        ['PID (n=',num2str(size(pid_rewardTraces,1)),')']});
xlabel('Time (s)'); ylabel('Intensity');
plotEvent('Water',0);


nexttile(6);
% yyaxis left
[baseline_rewardTraces_delta,~] = plotTraces(baseline_waterLickIdx,[-1,5],baseline_lp,...
                         baseline=[-1,0],...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=baselineColor);
[pid_rewardTraces_delta,~] = plotTraces(pid_waterLickIdx,[-1,5],pid_lp,...
                         baseline=[-1,0],...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=pidColor);
% yyaxis right
% [~,~] = plotTraces(pid_waterLickIdx,[-1,5],pid_blueLaserON*1000,...
%                      signalFs=originalFs,sameSystem=true,...
%                      plotIndividual=false,color=blueColor,...
%                      xlabel='Time (s)',ylabel='Stim duration (ms)');
legend({['Baseline (n=',num2str(size(baseline_rewardTraces_delta,1)),')'],...
        ['PID (n=',num2str(size(pid_rewardTraces_delta,1)),')']});
xlabel('Time (s)'); ylabel('\Delta Intensity');
plotEvent('Water',0);


nexttile(3,[2 1]);
plotHeatmap(pid_rewardTraces,timestamp);
xlabel('Time (s)'); ylabel('Trials');
title('PID: water');

% Save
% saveFigures(gcf,'Summary',sessionpath); close all;
return

%% Plot DA aligned trace

baseline_photometry = baseline_5ms;
pid_photometry = pid_5ms;
Fs = Fs_5ms;
fftWindow = 10;

initializeFig(1,1); tiledlayout(2,4);

nexttile(1,[1 2]);
yyaxis left
[baseline_toneTraces,~] = plotTraces(baseline_tone(1:100),[-1,5],baseline_raw,...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=baselineColor,...
                         xlabel='Time (s)',ylabel='Intensity');
yyaxis right
[pid_toneTraces,~] = plotTraces(pid_tone,[-1,5],pid_raw,...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=pidColor,...
                         xlabel='Time (s)',ylabel='Intensity');
yyaxis right
[~,~] = plotTraces(pid_tone,[-1,5],pid_blueLaserON*1000,...
                     signalFs=originalFs,sameSystem=true,...
                     plotIndividual=false,color=blueColor,...
                     xlabel='Time (s)',ylabel='Stim duration (ms)');
legend({['Baseline (n=',num2str(size(baseline_toneTraces,1)),')'],...
        ['PID (n=',num2str(size(pid_toneTraces,1)),')'],...
        ['Blue laser (n=',num2str(size(pid_toneTraces,1)),')']});
plotEvent('Tone',0);


nexttile(4);
[baseline_fft,baseline_fftpower] = plotFFT(baseline_photometry,Fs=Fs,color=baselineColor); 
[pid_fft,pid_fftpower] = plotFFT(pid_photometry,Fs=Fs,color=pidColor); 
legend('Baseline','PID');
title('FFT');

nexttile(8);
yyaxis left; histogram(baseline_photometry,200,EdgeColor=baselineColor,FaceColor=baselineColor); hold on
yyaxis right; histogram(pid_photometry,200,EdgeColor=pidColor,FaceColor=pidColor); hold on
xlabel('Intensity'); legend({'Baseline','PID'}); box off;
title('dLight');
subtitle({strcat("Basline std: ",num2str(std(baseline_raw)),", PID std: ",num2str(std(pid_raw)))});


nexttile(5,[1 2]);
yyaxis left
[baseline_rewardTraces,~] = plotTraces(baseline_waterLickIdx(1:100),[-1,5],baseline_raw,...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=baselineColor,...
                         xlabel='Time (s)',ylabel='Intensity');
[pid_rewardTraces,timestamp] = plotTraces(pid_waterLickIdx,[-1,5],pid_raw,...
                         signalFs=originalFs,sameSystem=true,...
                         plotIndividual=false,color=pidColor,...
                         xlabel='Time (s)',ylabel='Intensity');
yyaxis right
[~,~] = plotTraces(pid_waterLickIdx,[-1,5],pid_blueLaserON*1000,...
                     signalFs=originalFs,sameSystem=true,...
                     plotIndividual=false,color=blueColor,...
                     xlabel='Time (s)',ylabel='Stim duration (ms)');
legend({['Baseline (n=',num2str(size(baseline_rewardTraces,1)),')'],...
        ['PID (n=',num2str(size(pid_rewardTraces,1)),')'],...
        ['Blue laser (n=',num2str(size(pid_rewardTraces,1)),')']});
plotEvent('Water',0);


nexttile(3,[2 1]);
plotHeatmap(pid_toneTraces,timestamp);
xlabel('Time (s)'); ylabel('Trials');
title('PID: water');

% Save
% saveFigures(gcf,'Summary',sessionpath); close all;

%% Align to reward and airpuff

initializeFig(1,0.5); tiledlayout(1,1,Padding='loose',TileSpacing='tight'); 

nexttile;
timeToPlot = [170,190]; % in sec
% Plot photometry trace
plotWindow = timeToPlot(1)*Fs : timeToPlot(2)*Fs;
t = (timeToPlot(1)*Fs : timeToPlot(2)*Fs) / Fs;
yyaxis left; ylabel('Intensity'); 
plot(t,pid_photometry(plotWindow),Color=pidColor,LineWidth=3,LineStyle='-'); hold on;
% Plot blue laser
plotWindow = timeToPlot(1)*originalFs : timeToPlot(2)*originalFs;
t = (timeToPlot(1)*originalFs : timeToPlot(2)*originalFs) / originalFs;
yyaxis right; ylabel('Time (ms)');
plot(t,pid_blueLaserON(plotWindow)*1000,Color=blueColor,LineWidth=3);
% Plot settings
xlabel('Time (s)'); box off;
legend('PID','Blue laser');

%% Load !fake! PID sessions

% Select sessions via uipickfiles
sessionPath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));
sessionpath = sessionPath{1};
clear pid_5ms pid_blueLaserON

% Load sessions
% Get animal name and session date
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; clear dirsplit
load(strcat(sessionpath,filesep,'data_',sessionName,'.mat'));

pid_raw = photometry_raw;
pid_blueLaser = blueLaser;

% Downsample photometry signal
if ~exist('pid_5ms','var')
    downsampled = downsampleSignal(pid_raw,targetFs=Fs_5ms,originalFs=originalFs,rollingZ=false);
    pid_5ms = downsampled.dsData;
    save(strcat(sessionpath,filesep,'data_',sessionName,'.mat'),'pid_5ms','-append');
    disp('Finished: downsampled and saved PID raw photometry');
end

% Get duration of blueLaser
if ~exist('pid_blueLaserON','var')
    pid_blueLaserON = ~savedDigitalNI .* getConsecutive(~savedDigitalNI)./originalFs;
    save(strcat(sessionpath,filesep,'data_',sessionName,'.mat'),'pid_blueLaserON','-append');
    disp('Finished: saved pid blueLaserON');
end

% Plot session
Fs = Fs_5ms; pid_photometry = pid_5ms;
initializeFig(1,0.5); tiledlayout(1,1,Padding='loose',TileSpacing='tight');

nexttile;
title(sessionName(1:end-3));
timeToPlot = [50,80]; % in sec
% Plot photometry trace
plotWindow = timeToPlot(1)*Fs : timeToPlot(2)*Fs;
t = (timeToPlot(1)*Fs : timeToPlot(2)*Fs) / Fs;
yyaxis left; ylabel('Intensity'); 
plot(t,pid_photometry(plotWindow),Color=pidColor,LineWidth=3,LineStyle='-'); hold on;
% Plot blue laser
plotWindow = timeToPlot(1)*originalFs : timeToPlot(2)*originalFs;
t = (timeToPlot(1)*originalFs : timeToPlot(2)*originalFs) / originalFs;
yyaxis right; ylabel('Time (ms)');
plot(t,pid_blueLaserON(plotWindow)*1000,Color=blueColor,LineWidth=3);
% Plot settings
xlabel('Time (s)'); box off;
legend('PID','Blue laser');

saveFigures(gcf,'Summary',sessionpath);

% nexttile;
% yyaxis left;
% [pid_rewardTraces,timestamp] = plotTraces(pid_blueLaserON,[-1,5],pid_raw,...
%                          signalFs=originalFs,sameSystem=true,...
%                          plotIndividual=false,color=pidColor,...
%                          xlabel='Time (s)',ylabel='Intensity');
% yyaxis right
% [~,~] = plotTraces(pid_blueLaser,[-1,5],pid_blueLaserON*1000,...
%                      signalFs=originalFs,sameSystem=true,...
%                      plotIndividual=false,color=blueColor,...
%                      xlabel='Time (s)',ylabel='Stim duration (ms)');


%% Generate dataset for reward lick decoder analysis

% Sample: n x m
% n = number of windows (1 sec/window)
% m = number of samples within each window

% Target: n x 1
% whether reward lick happens in that window

% Settings
originalFs = 10000; windowLength = 1; % sec
nSamplesPerWindow = windowLength * originalFs;

% Step 0: select files
sessionPath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));

% Loop over each files
for i = 1:length(sessionPath)
    sessionpath = sessionPath{i};

    % Step 1: load photometry & reward lick
    dirsplit = strsplit(sessionpath,filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    load(strcat(sessionpath,filesep,'data_',sessionName,'.mat'));
    lickON = find(rightLick);
    rewardIdx = find(rightSolenoid);
    waterLickIdx = nan(size(rewardIdx));
    for i = 1:length(rewardIdx)
        nextLick = lickON(find(lickON>=rewardIdx(i),1));
        if ~isempty(nextLick); waterLickIdx(i) = nextLick; end
    end
    waterLickIdx = rmmissing(waterLickIdx);
    waterLick = zeros(size(photometry_raw)); waterLick(waterLickIdx) = 1;
    
    % Step 2: truncate photometry into windows
    maxSample = nSamplesPerWindow * floor(length(photometry_raw)/nSamplesPerWindow);
    photometry_reshape = reshape(photometry_raw(1:maxSample),[],nSamplesPerWindow);
    
    % Step 3: truncate waterLickIdx into windows
    waterLick_reshape = reshape(waterLick(1:maxSample),[],nSamplesPerWindow);
    waterLick_reshape = any(waterLick_reshape,2);
    
    % Step 4: truncate reward into windows
    reward_reshape = reshape(rightSolenoid(1:maxSample),[],nSamplesPerWindow);
    reward_reshape = any(reward_reshape,2);
    
    % Step 5: truncate lick int windows
    lick_reshape = reshape(rightLick(1:maxSample),[],nSamplesPerWindow);
    lick_reshape = any(lick_reshape,2);
    
    % Save data for this session
    save(strcat(sessionpath,filesep,'decoder_',sessionName,'.mat'),...
        'photometry_reshape','waterLick_reshape',...
        'lick_reshape','reward_reshape','-v7.3');
    disp(['Decoder dataset saved: ',sessionName]);
end