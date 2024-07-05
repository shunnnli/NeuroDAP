%% Shun_loadSiceData
% Modified from Shun_analyzeSlice

% 09/13/23
% Separated from Shun_analyzeSlice, the idea is to plot individual and
% average trace from each epoch without referencing Excel data

% 09/14/23
% Package loading part and anlaysis part into separate function

%% Define data path
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
saveDataPath = 'default'; % strcat(parentPath,filesep,'20231221_ally');

% Set comman params
[sessionParams,canceled] = inputSessionParams_singleSlice(expPath,...
                                paradigm=2,redStim=true,...
                                reload=false,calculateQC=false,...
                                timeRange='[-20,100]',nArtifactSamples='10');
taskOptions = {'random','reward pairing','punish pairing'};
task = taskOptions{sessionParams.Paradigm};
timeRange = eval(sessionParams.timeRange);
nArtifactSamples = str2double(sessionParams.nArtifactSamples);

if ~isscalar(expPath); error('Multiple sessions were selected, for multi-session analysis see analyzeSlice pipeline!'); end
expPath = expPath{1};

%% (Optional) Just load epochs.mat

[epochs_old,cells] = loadSlices(expPath,reload=sessionParams.reload);

% Reprocess post QC epochs.mat
[epochs,cells] = loadSlices(epochs_old,reload=sessionParams.reload,...
                animal=sessionParams.Animal,task=task,...
                timeRange=timeRange,...
                filterSignal=false,filterSweeps=true,...
                calculateQC=sessionParams.calculateQC,...
                nArtifactSamples=nArtifactSamples,...
                saveDataPath=saveDataPath,...
                save=true);
return

%% Load epoch for single session

dirsplit = split(expPath,filesep); expName = dirsplit{end};

[epochs,cells] = loadSlices(expPath,reload=sessionParams.reload,...
                animal=sessionParams.Animal,task=task,...
                timeRange=timeRange,...
                filterSignal=false,filterSweeps=true,...
                calculateQC=sessionParams.calculateQC,...
                nArtifactSamples=nArtifactSamples,...
                saveDataPath=saveDataPath);

analyzeSlice_OptoPair(expPath,...
            timeRange=timeRange,...
            nArtifactSamples=nArtifactSamples);

save(strcat(expPath,filesep,'PreQC',filesep,'epochs_',expName),'epochs','-v7.3');
disp(strcat("Saved: ",expName," in PreQC folder"));
close all

% Message for editing quality checks
f = msgbox("Edit the epochs table and run following block after finished",...
    "Quality check","help");
return

%% Quality checks: save modified epoch file

% Message for editing quality checks
answer = questdlg('Confirm and save quality check results?', ...
    'Quality check confirmation','Yes','Not yet','Not yet');
switch answer
    case 'Yes'; qc = true;
    case 'Not yet'; qc = false;
end

if qc
    dirsplit = split(expPath,filesep); expName = dirsplit{end};
    save(strcat(expPath,filesep,'epochs_',expName),'epochs','-v7.3');
    save(strcat(expPath,filesep,'PostQC',filesep,'epochs_',expName),'epochs','-v7.3');
    disp(strcat("Saved: ",expName," in PostQC folder"));

    analyzeSlice_OptoPair(expPath,...
                timeRange=timeRange,nArtifactSamples=nArtifactSamples,...
                resultPath='PostQC',plotAll=false);
    close all
end

return

%% Useful code to plot raw sweeps

close all
initializeFig(0.67, 0.5); tiledlayout(1,3);
row = 13;

% Find event window
timeRange = [-20,100];
timeRangeStartSample = 10000 + 10000*timeRange(1)/1000;
timeRangeEndSample = 10000 + 10000*timeRange(2)/1000;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = (plotWindow-1*10000) ./ (10000/1000);
analysisWindow = (10000+nArtifactSamples)-timeRangeStartSample : length(plotWindow);

% Plot all traces
nexttile;
included = ones(size(epochs{row,'Raw sweeps'}{1},1),1);
traces = epochs{row,'Raw sweeps'}{1}(included==1,plotWindow);
if ~isempty(traces)
    plotSEM(timeRangeInms,traces,[0.343, 0.75, 0.232],...
            meanOnly=true,plotIndividual=true);
    xlabel('Time (ms)');
    ylabel('Current (pA)');
    yMin = min(traces(:,analysisWindow),[],"all");
    yMax = max(traces(:,analysisWindow),[],"all");
    yPad = abs(yMax-yMin)*0.1;
    ylim([yMin-yPad,yMax+yPad]);
end
title(strcat('Epochs #',num2str(epochs{row,'Epoch'}),' (all traces)'));

% Plot included trace only
nexttile;
included = epochs{row,'Included'}{1};
traces = epochs{row,'Raw sweeps'}{1}(included==1,plotWindow);
if ~isempty(traces)
    plotSEM(timeRangeInms,traces,[0.343, 0.75, 0.232],...
            meanOnly=true,plotIndividual=true);
    xlabel('Time (ms)');
    ylabel('Current (pA)');
    yMin = min(traces(:,analysisWindow),[],"all");
    yMax = max(traces(:,analysisWindow),[],"all");
    yPad = abs(yMax-yMin)*0.1;
    ylim([yMin-yPad,yMax+yPad]);
end
title(strcat('Epochs #',num2str(epochs{row,'Epoch'}),' (included traces)'));

% Plot excluded trace only
nexttile;
included = epochs{row,'Included'}{1};
traces = epochs{row,'Raw sweeps'}{1}(included~=1,plotWindow);
if ~isempty(traces)
    plotSEM(timeRangeInms,traces,[0.343, 0.75, 0.232],...
            meanOnly=true,plotIndividual=true);
    xlabel('Time (ms)');
    ylabel('Current (pA)');
    yMin = min(traces(:,analysisWindow),[],"all");
    yMax = max(traces(:,analysisWindow),[],"all");
    yPad = abs(yMax-yMin)*0.1;
    ylim([yMin-yPad,yMax+yPad]);
end
title(strcat('Epochs #',num2str(epochs{row,'Epoch'}),' (excluded traces)'));

%% Useful code to move sweeps between epochs

originalRow = 7;
newRow = 10;
sweepIdx = 1:7;

% Move sweep names
sweepNames = epochs{originalRow,"Sweep names"}{1}(sweepIdx);
disp(["Moving sweeps: ", sweepNames]);
epochs{newRow,"Sweep names"}{1} = [epochs{newRow,"Sweep names"}{1}, sweepNames];

% Edit included to new row, set true in new row, false in old row
epochs{newRow,"Included"}{1} = [epochs{newRow,"Included"}{1}; ones(length(sweepIdx),1)];
epochs{originalRow,"Included"}{1}(sweepIdx) = zeros(length(sweepIdx),1);

% Move raw sweeps
moveSweeps = epochs{originalRow,"Raw sweeps"}{1}(sweepIdx,:);
epochs{newRow,"Raw sweeps"}{1} = [epochs{newRow,"Raw sweeps"}{1}; moveSweeps];
% Move processed sweeps
moveSweeps = epochs{originalRow,"Processed sweeps"}{1}(sweepIdx,:);
epochs{newRow,"Processed sweeps"}{1} = [epochs{newRow,"Processed sweeps"}{1}; moveSweeps];

% Move Vhold sweeps
moveSweeps = epochs{originalRow,"Vhold sweep trace"}{1}(sweepIdx,:);
epochs{newRow,"Vhold sweep trace"}{1} = [epochs{newRow,"Vhold sweep trace"}{1}; moveSweeps];
% Recalculate Vhold epoch trace
epochs{newRow,"Vhold epoch trace"}{1} = mean(epochs{newRow,"Vhold sweep trace"}{1});
% Recalculate Vhold epoch mean
epochs{newRow,"Vhold epoch mean"} = mean(epochs{newRow,"Vhold sweep trace"}{1},'all');

% Move Peaks
moveSweeps = epochs{originalRow,"Peaks"}{1}(sweepIdx);
epochs{newRow,"Peaks"}{1} = [epochs{newRow,"Peaks"}{1}; moveSweeps];
% Move AUCs
moveSweeps = epochs{originalRow,"AUCs"}{1}(sweepIdx);
epochs{newRow,"AUCs"}{1} = [epochs{newRow,"AUCs"}{1}; moveSweeps];
% Move Rin
moveSweeps = epochs{originalRow,"Rin"}{1}(sweepIdx);
epochs{newRow,"Rin"}{1} = [epochs{newRow,"Rin"}{1}; moveSweeps];
% Move Rs
moveSweeps = epochs{originalRow,"Rs"}{1}(sweepIdx);
epochs{newRow,"Rs"}{1} = [epochs{newRow,"Rs"}{1}; moveSweeps];
% Move Cm
moveSweeps = epochs{originalRow,"Cm"}{1}(sweepIdx);
epochs{newRow,"Cm"}{1} = [epochs{newRow,"Cm"}{1}; moveSweeps];

disp('Moving finished');

%% Useful code to undo moving

% originalRow = 7;
% newRow = 10;
% sweepIdx = 1:7;
% 
% % Move sweep names
% sweepNames = epochs{originalRow,"Sweep names"}{1}(sweepIdx);
% disp(["Undo moving sweeps: ", sweepNames]);
% epochs{newRow,"Sweep names"}{1}(end-length(sweepIdx)) = [];
% 
% % Remove included to new row, set true in old row
% epochs{newRow,"Included"}{1}(end-length(sweepIdx)) = [];
% epochs{originalRow,"Included"}{1}(sweepIdx) = ones(length(sweepIdx),1);
% 
% % Move raw sweeps
% epochs{newRow,"Raw sweeps"}{1}(end-length(sweepIdx)) = [];
% % Move processed sweeps
% epochs{newRow,"Processed sweeps"}{1}(end-length(sweepIdx)) = [];
% 
% % Move Vhold sweeps
% epochs{newRow,"Vhold sweep trace"}{1}(end-length(sweepIdx)) = [];
% % Recalculate Vhold epoch trace
% epochs{newRow,"Vhold epoch trace"}{1} = mean(epochs{newRow,"Vhold sweep trace"}{1});
% % Recalculate Vhold epoch mean
% epochs{newRow,"Vhold epoch mean"} = mean(epochs{newRow,"Vhold sweep trace"}{1},'all');
% 
% % Move Peaks
% epochs{newRow,"Peaks"}{1}(end-length(sweepIdx)) = [];
% % Move AUCs
% epochs{newRow,"AUCs"}{1}(end-length(sweepIdx)) = [];
% % Move Rin
% epochs{newRow,"Rin"}{1}(end-length(sweepIdx)) = [];
% % Move Rs
% epochs{newRow,"Rs"}{1}(end-length(sweepIdx)) = [];
% % Move Cm
% epochs{newRow,"Cm"}{1}(end-length(sweepIdx)) = [];
% 
% disp('Undo moving finished');

%% Run Quality Control
disp('Running Quality Control Check');
QCs = zeros(length(raw_concatenated_traces), length(raw_concatenated_traces(1).data));
for QCCheck = 1:length(raw_concatenated_traces)
    QCs(QCCheck,:) = raw_concatenated_traces(QCCheck).data;
end

tau = [];
rs = [];
rin = [];
cm = [];
baseline = [];

for l = 1:size(QCs,1)
    dt = 0.1; % ms
    T = 3000; % ms
    tvec = dt:dt:T; % time vector
    NS = length(tvec); % number samples

    % You'd provide this I'm just showing you how to make it
    data = QCs(l,:);
    amplitudeRC = -5; % mV
    
    if isequal(recorder, 'KM') == 1
        rcStart = 250;
        rcEnd = 550;
    elseif isequal(recorder, 'WW') == 1 
        rcStart = 2800;
        rcEnd = 2900;
    end 

    rcDelay = .05; % ms to wait after start of rc check to do exponential fit and a few other things

    startBaselineSample = find(tvec>=rcStart-(rcEnd-rcStart)*0.1,1,'first'); % calculate baseline with same amount of time
    rcStartSample = find(tvec>=rcStart,1,'first'); 
    startEndlineSample = find(tvec>=rcStart + (rcEnd-rcStart)*0.9,1,'first') - 1; % endline start sample
    rcEndSample = find(tvec>=rcEnd,1,'first');
    numDelaySamples = round(rcDelay / dt);

    % Compute peak and get peak location
    if amplitudeRC > 0
        [pk, pkIdx] = max(data(rcStartSample:rcEndSample));
    else
        [pk, pkIdx] = min(data(rcStartSample:rcEndSample));
    end
    pkIdx = pkIdx + rcStartSample - 1; 

    baseline(l) = mean(data(startBaselineSample:rcStartSample-1));
    endline = mean(data(startEndlineSample:rcEndSample));


    % Dan- play with this. Higher percentage biases fit to the fast component
    % which is more accurate in dendritic cells
    percentageDecayCutoff = 15; 
    decayCutoffValue = (pk-endline)*percentageDecayCutoff/100 + endline;

    if amplitudeRC > 0
        decayEnd = find(data(pkIdx+numDelaySamples:rcEndSample)<decayCutoffValue,1,'first'); % Get's first below 20%
    else
        decayEnd = find(data(pkIdx+numDelaySamples:rcEndSample)>decayCutoffValue,1,'first'); % Get's first below 20%
    end

    decayTime = tvec(pkIdx+numDelaySamples:pkIdx+numDelaySamples+decayEnd-1);
    
    if length(decayTime)>1
        decayTime = decayTime - decayTime(1); %start at 0
        decayData = data(pkIdx+numDelaySamples:pkIdx+numDelaySamples+decayEnd-1) - endline; 
        if amplitudeRC < 0, decayData = -decayData; end

        % Fit Options
        fitType = fittype('peak*exp(-x/tau)');
        fitOptions = fitoptions(fitType);
        fitOptions.StartPoint = [pk-endline decayEnd*dt/3];
        fitOptions.Lower = [0 0];
        fitOptions.Upper = [abs(pk)*10 T];

        % Do Fit
    
        
        decayFit = fit(decayTime(:),decayData(:),fitType,fitOptions);
        tau(l) = decayFit.tau; % in whatever units you made tvec (should be ms)
    else
        tau(l) = NaN;
    end
    
    relativeStart = (pkIdx+numDelaySamples - rcStartSample)*dt;

    % Extrapolate back for exponential estimate of peak
    pkEstimateExponential = decayFit.peak * exp(relativeStart/decayFit.tau) * (-1*(amplitudeRC<0)+1*(amplitudeRC>0)); 

    choiceOfPeak = pk; % Or use "pkEstimateExponential"
    % - using the true pk, rather than pkEstimateExponential, will probably be
    % more accurate if you're patching dendritic cells
    rs(l) = 1000 * amplitudeRC / (choiceOfPeak-endline);
    rin(l) = (1000 * amplitudeRC/(endline-baseline(l))) - rs(l);
    cm(l) = 1000 * tau(l) / rs(l);
end


%Plot everything
figure
for m = 1:size(QCs,1)
    subplot(3,2,1);
    plot(QCs(m, (rcStart*10)-100:(rcStart*10)+200));
    ylabel('QC Pulse');
    xlim([0 300]);
    ylim([(min(QCs(m, (rcStart*10)-100:(rcStart*10)+200))-300) (max(QCs(m, (rcStart*10)-100:(rcStart*10)+200))+300)]);
    hold on
end
subplot(3,2,2);
plot(baseline, 'ko');
ylabel('Baseline');
xlim([0 size(QCs,1)])
subplot(3,2,3);
plot(rs,'ko');
xlim([0 size(QCs,1)])
ylabel('Rs (MOhm)')
subplot(3,2,4);
plot(rin,'bo');
xlim([0 size(QCs,1)])
ylabel('Rin (MOhm)')
subplot(3,2,5);
plot(tau,'ko');
xlim([0 size(QCs,1)])
ylabel('{\tau}_m (ms)')
subplot(3,2,6);
plot(cm,'ko');
xlim([0 size(QCs,1)])
ylabel('Cm (pF)')
xlabel('Sweep #')

%Make it look pretty
sgtitle(strcat(mouseID, ' cell ', celll, ', epoch ', epochh));
set(gcf, 'Color', 'w');

%Adjust and save
print(fullfile(savePath1, strcat(mouseID, '_cell', celll, '_epoch', epochh, '_QualityControl')), '-dpdf', '-fillpage', '-r1000');
print(fullfile(savePath, strcat(mouseID, 'QualityControl')), '-dpsc', '-fillpage', '-append', '-r1000'); 

%Close everything
close all
clear tau tauCap tauParams Cm Rm Raccess extrapolPkRC fittedVal f fittingTrace2 fittingTrace x x2 steadyRC baselineRC RCtrace locOnset
clear AverageTraces AllTracesTable