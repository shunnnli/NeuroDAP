% Shun_loadSliceData_DMD

% Load and analyze DMD slice data (single session)

%% Load sessions
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
saveDataPath = 'default';

% Set comman params
[sessionParams,canceled] = inputSessionParams_singleSlice(expPath,...
                                paradigm=1,redStim=true,...
                                reload=false,calculateQC=true,...
                                timeRange='[-10,50]',nArtifactSamples='0');
taskOptions = {'random','reward pairing','punish pairing'};
task = taskOptions{sessionParams.Paradigm};
timeRange = eval(sessionParams(1).timeRange);
nArtifactSamples = str2double(sessionParams(1).nArtifactSamples);

if ~isscalar(expPath); error('Multiple sessions were selected, for multi-session analysis see analyzeSlice pipeline!'); end
expPath = expPath{1};

%% Process data at epoch level

% This is for epoch-level analysis. It assumes all sweeps within the epoch 
% have the same length (not true for randomSearch)
% Therefore, for randomSearch epochs, raw sweeps and processed sweeps are
% meaningless (should be all zeros). These epochs will be address
% separately in later (see below).

% sessionParams.calculateQC = true;
epochs = loadSlices(expPath,reload=sessionParams.reload,...
                    reloadCell=false,...
                    animal=sessionParams.Animal,task=task,...
                    timeRange=timeRange,...
                    vholdChannel='AD1',...
                    filterSignal=false,filterSweeps=true,...
                    calculateQC=sessionParams.calculateQC,...
                    nArtifactSamples=nArtifactSamples,...
                    getCellTable=false,...
                    saveDataPath=saveDataPath,...
                    plot=false);

%% Process data for random search epoch
% This will take the most amount of time and is the most important step

% For each search depth, it separately save a file that stores all the
% params and also responses, QC, and statistics of that search depth
close all;

% loadSlicesDMD(epochs,reload=false,reloadCells=true,reloadCellAnalysis=true);
loadSlicesDMD(epochs,reload=true);

%% Plot DMD results
% Plot results for all searches in this session
close all;
cells = analyzeSlice_DMD(expPath,resultsPathDate='newest',...
                        timeRange=timeRange,...
                        redStim=true,...
                        plotSearch=true,plotPairs=true,...
                        savePNG=false,savePDF=true,saveFIG=false);
return

%% Plot search summary
% Just plot results for a specific search in a session

cellNum = 4; searchIdx = 2; pairIdx = 1; close all;

disp(['Ongoing: plotting searches for cell',num2str(cellNum)]);
curCell = cells(cells.Cell == cellNum,:);
searchPerCell = length(curCell.Epochs{1});

traceColor = [127 182 227]./255;
whiteBlue = getColormap([238 240 241],[48 154 209],500,midCol=[150 198 227]);

if isnan(searchIdx) || ~isnumeric(searchIdx)
    for searchIdx = 1:searchPerCell
        analyzeDMDSearch(curCell,searchIdx,...
                     redStim=true,...
                     savePNG=false,savePDF=true,saveFIG=false);
    end
else
    analyzeDMDSearch(curCell, searchIdx, ...
                     redStim=true,...
                     savePNG=false, savePDF=true, saveFIG=false);
end

% analyzeDMDSearchPair(curCell,pairIdx,...
%              redStim=true,...
%              savePNG=false,savePDF=true,saveFIG=false);