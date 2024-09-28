%% Shun_loadSiceData
% Modified from Shun_analyzeSlice

% 09/13/23
% Separated from Shun_analyzeSlice, the idea is to plot individual and
% average trace from each epoch without referencing Excel data

% 09/14/23
% Package loading part and anlaysis part into separate function

%% Define data path
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
saveDataPath = 'default'; % strcat(parentPath,filesep,'20231221_ally');

% Set comman params
[sessionParams,canceled] = inputSessionParams_singleSlice(expPath,...
                                paradigm=1,redStim=true,...
                                reload=false,calculateQC=true,...
                                timeRange='[-10,50]',nArtifactSamples='10');
taskOptions = {'random','reward pairing','punish pairing'};

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
today = char(datetime('today','Format','yyyyMMdd')); 

%% (Optional) Just load epochs.mat

combined_epochs = []; combined_cells = []; close all;

for i = 1:length(expPath)
    timeRange = eval(sessionParams(i).timeRange);
    nArtifactSamples = str2double(sessionParams(i).nArtifactSamples);

    % (Optional) Define inclusion criteria if necessary
    QCThreshold.include = {};
    QCThreshold.Rs = 30;
    QCThreshold.Verror = 10;
    QCThreshold.Ibaseline = -300;
    QCThreshold.Ibaseline_std = 20;

    % Reprocess post QC epochs.mat
    [epochs] = loadSlices(expPath{i},reload=true,...
                                timeRange=timeRange,...
                                filterSignal=false,filterSweeps=true,...
                                calculateQC=sessionParams(i).calculateQC,...
                                nArtifactSamples=nArtifactSamples,...
                                saveDataPath=saveDataPath,...
                                save=true,...
                                QCThreshold=QCThreshold,...
                                getCellTable=false,...
                                animal=sessionParams(i).Animal,...
                                task=taskOptions{sessionParams(i).Paradigm});

    % Load & combine epochs.mat
    % [epochs,cells] = loadSlices(expPath{i},reload=sessionParams(i).reload,reloadCell=true);
    % combined_epochs = [combined_epochs; epochs];
    % combined_cells = [combined_cells; cells];
end

return

%% Save current epochs to session folder

notes = 'QC';
sessionPath = epochs{1,'Session'};

% Save current epochs to the newest results folder
resultsFolders = sortrows(struct2cell(dir(fullfile(sessionPath,"Epochs-*")))',[1 3]);
resultFolder = resultsFolders{end,1};

dirsplit = split(sessionPath,filesep); expName = dirsplit{end};
dirsplit = split(resultFolder,'-'); folderDate = dirsplit{end};
savePath = fullfile(sessionPath,resultFolder);

save(strcat(sessionPath,filesep,'epochs_',folderDate,'_',notes),'epochs','-v7.3');
disp(strcat("Saved: ",expName," in session folder"));

save(strcat(savePath,filesep,'epochs_',folderDate,'_',notes),'epochs','-v7.3');
disp(strcat("Saved: ",expName," in results folder"));

%% Plot epoch summary

for rowIdx = 1:size(epochs,1)
    close all;
    plotEpochSummary(epochs,rowIdx,save=true);
end
close all;

%% Useful code to plot raw sweeps

close all
initializeFig(0.67, 0.5); tiledlayout(1,3);
row = 6;

plotWholeTrace = false;

if plotWholeTrace
    plotWindow = 1:30000;
    timeRangeInms = (plotWindow-1*10000) ./ (10000/1000); 
    analysisWindow = 1:30000;
else
    % Find event window
    timeRange = [-10,50];
    timeRangeStartSample = 10000 + 10000*timeRange(1)/1000;
    timeRangeEndSample = 10000 + 10000*timeRange(2)/1000;
    plotWindow = timeRangeStartSample : timeRangeEndSample;
    timeRangeInms = (plotWindow-1*10000) ./ (10000/1000);
    analysisWindow = (10000+10)-timeRangeStartSample : length(plotWindow);
end

% Plot all traces
nexttile;
included = ones(size(epochs{row,'Raw sweeps'}{1},1),1);
traces = epochs{row,'Raw sweeps'}{1}(:,plotWindow);
if ~isempty(traces)
    plotSEM(timeRangeInms,traces,[0.343, 0.75, 0.232],...
            plotPatch=false,plotIndividual=true,plotMean=false);
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
            plotPatch=false,plotIndividual=true);
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
            plotPatch=false,plotIndividual=true);
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

%% Misc: plot histogram for Rs and voltage error

close all;
initializeFig(1,1); tiledlayout('flow');

nexttile;
allRs = cell2mat(cellfun(@(x) x.Rs, epochs.QC,UniformOutput=false));
histogram(allRs,100);
title('Rs (MOhm)');

nexttile;
allVerror = cell2mat(cellfun(@(x) abs(x.Verror), epochs.QC,UniformOutput=false));
histogram(allVerror,100);
title('|Verror| (mV)');

nexttile;
allIbaseline = cell2mat(cellfun(@(x) x.Ibaseline, epochs.QC,UniformOutput=false));
histogram(allIbaseline,100);
title('Ibaseline (pA)');

nexttile;
allIbaseline_std = cell2mat(cellfun(@(x) x.Ibaseline_std, epochs.QC,UniformOutput=false));
histogram(allIbaseline_std,100);
title('Ibaseline STD');