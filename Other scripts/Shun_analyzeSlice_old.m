% Shun_analyzeSlice
% Shun Li, 12/23/2022

%% Load data
clear; close all;
% addpath(genpath('D:\Shun\Analysis\Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
[~,colors,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed] = loadColors;

resultspath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\project valence\Recordings\Photometry-flexibleLearning-NAc-dLight\Results\';

% epochList = uipickfiles('FilterSpec','\\research.files.med.harvard.edu\neurobio\MICROSCOPE\wengang\Exp_withShun');
sessionpath = uigetdir('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\wengang\Exp_withShun');
epochList = struct2cell(dir(fullfile(sessionpath,['AD0_e*','p1avg.mat'])))';

%% Build a cell array for all cell, load raw data of each sweep

sweepList = cell(length(epochList),1);

for i = 1:length(epochList)
    % Load epoch file to find individual sweep.mat
    load(fullfile(sessionpath,epochList{i,1}));
    namesplit = strsplit(epochList{i,1},{'e','p1avg'}); epoch = str2num(namesplit{2});
    sweepNums = evalin('base',['AD0_e',num2str(epoch),'p1avg.UserData.Components']);

    % Load data of individual sweeps
    for j = 1:length(sweepNums)
        load(fullfile(sessionpath,strcat(sweepNums{j},'.mat'))); 
        disp(['Loading ',sweepNums{j},'.mat for epoch ',num2str(epoch)]);
        sweeps(j,:) = eval([sweepNums{j},'.data']);
        sweepList{epoch} = sweeps;
    end
    clearvars AD0_* sweeps
end

%% Perform mean subtraction (baseline as 0-1000ms)

processed = cell(size(sweepList));

for i = 1:length(sweepList)
    epochTraces = sweepList{i};

    if isempty(epochTraces); continue; end

    % Pre-processing
    % Calculate baseline current amplitude
    base = mean(epochTraces(:,1:1999),2); baseM = repmat(base,1,size(epochTraces,2));
    epoch_subtracted = epochTraces-baseM;  
    
    % Low pass filter at 2kHz
    Fs = 10000; % Sampling frequency  
    LP = lowpass(epoch_subtracted',2000,Fs);
    % Notch filter
    % d = designfilt('bandstopiir','FilterOrder',2, ...
    %                'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    %                'DesignMethod','butter','SampleRate',Fs);
    % Notch = filtfilt(d,LP);   
    
    % Smooth data using sgolay filter
    yT = sgolayfilt(LP,5,27); % polynomial order of 5 and framelength of 27
    y = yT';
    
    % Median filter using 0.5ms window
    y = movmedian(y,6,2);

    % figure;
    % plot(epoch_subtracted(2,:)); hold on; plot(y(2,:));

    % Subtract mean again (in SeulAh's code)
    base2 = mean(y(:,1:1999),2); baseM2 = repmat(base2,1,size(y,2));
    processed{i} = y - baseM2;
    
end

%%  Eliminate error epochs based on spreadsheet

deleteThese = [1,2,7,10,11,12]; % SL043
% deleteThese = []; % SL044
% deleteThese = [1:13 16:43]; % SL045
% deleteThese = [1 2 3 9 12 15 16 22:24]; % SL046

for i = 1:length(deleteThese)
    processed{deleteThese(i)} = [];
end

% SL046: remove some sweep in epoch 17
% temp = processed{17};
% processed{17} = temp(end-5:end,:);
% temp = processed{10};
% processed{10} = temp(end-5:end,:);

%% test plot
timeRange = (950:1200)*10;
trace = processed{10};

figure;
plotSEM(timeRange,trace(:,timeRange),blueWhiteRed(1,:));

%% Extract EPSC or IPSC

% SL043
mouseID = 'SL043';
cellEPSC = [1 3 5 8 10 13 15 17 19 21]; 
cellIPSC = [2 4 6 9 11 14 16 18 20 22];

% SL044
% mouseID = 'SL044';
% cellEPSC = 1:2:length(processed); 
% cellIPSC = cellEPSC + 1;

% SL045
% mouseID = 'SL045';
% cellEPSC = 14; cellIPSC = 15;

% SL046
% mouseID = 'SL046';
% cellEPSC = [1 4 6 8 11 14 18 20];
% cellIPSC = [2 5 7 10 13 17 19 21];

% Initialize peak arrays
peakEPSC = nan(size(cellEPSC)); peakIPSC = nan(size(cellIPSC));

initializeFig(1,1);
tiledlayout('flow');
for i = 1:length(cellEPSC)
    EPSC = processed{cellEPSC(i)}; IPSC = processed{cellIPSC(i)};

    if isempty(EPSC)||isempty(IPSC); continue; end

    % Build noise model

    % Extract peak
    pTimeRange = (1000:1050)*10;
    peakEPSC(i) = min(mean(EPSC(:,pTimeRange)));
    peakIPSC(i) = max(mean(IPSC(:,pTimeRange)));
    
    % Plot cells
    timeRange = (950:1200)*10;
    nexttile
    plotSEM(timeRange,EPSC(:,timeRange),blueWhiteRed(1,:));
    plotSEM(timeRange,IPSC(:,timeRange),blueWhiteRed(500,:));
    xlabel("Time (ms)"); ylabel("Current (pA)");
    title(['Cell',num2str(i)]);

end

nexttile
scatter(abs(peakEPSC),peakIPSC,'filled'); refline(1,0);
xlabel("Peak EPSC (pA)"); ylabel("Peak IPSC (pA)");

saveFigures(gcf,strcat(mouseID,"-EPSCvsIPSC"),resultspath);
save(fullfile(resultspath,strcat(mouseID,"-EPSCvsIPSC"),mouseID),...
    'processed','sweepList','deleteThese',...
    'cellEPSC','cellIPSC','peakEPSC','peakIPSC','-v7.3');


%% ************** Multiple animals ******************

mouseIDs = {'SL043','SL044','SL045','SL046'};
peakEPSCs = cell(size(mouseIDs)); peakIPSCs = cell(size(mouseIDs));

colorList = [1, 500, 400, 100];
initializeFig(0.5,0.5);
for i = 1:length(mouseIDs)
    disp(mouseIDs{i});
    load(fullfile(resultspath,strcat(mouseIDs{i},'-EPSCvsIPSC'),strcat(mouseIDs{i},'.mat')));

    peakEPSCs{i} = peakEPSC; peakIPSCs{i} = peakIPSC;
    scatter(abs(peakEPSC),peakIPSC,[],blueWhiteRed(colorList(i),:),'filled'); hold on
end

hline = refline(1,0); hline.LineStyle= '--';
xlabel("Peak EPSC (pA)"); ylabel("Peak IPSC (pA)");
legend(mouseIDs);

saveFigures(gcf,"all-EPSCvsIPSC",resultspath);
save(fullfile(resultspath,"all-EPSCvsIPSC","params"),...
    'peakEPSCs','peakIPSCs','-v7.3');
