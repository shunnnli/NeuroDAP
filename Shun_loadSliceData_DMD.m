% Shun_loadSliceData_DMD

% Load and analyze DMD slice data (single session)

%% Load sessions
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
saveDataPath = 'default';

% Set comman params
[sessionParams,canceled] = inputSessionParams_singleSlice(expPath,...
                                paradigm=1,redStim=true,...
                                reload=false,calculateQC=true,...
                                timeRange='[-10,50]',nArtifactSamples='0');
taskOptions = {'random','reward pairing','punish pairing'};
task = taskOptions{sessionParams.Paradigm};
timeRange = eval(sessionParams.timeRange);
nArtifactSamples = str2double(sessionParams.nArtifactSamples);

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
                    filterSignal=false,filterSweeps=true,...
                    calculateQC=sessionParams.calculateQC,...
                    nArtifactSamples=nArtifactSamples,...
                    getCellTable=false,...
                    saveDataPath=saveDataPath,...
                    plot=false);

%% Process data for random search epoch
% This will take the most amount of time and is the most important step
close all;

loadSlicesDMD(epochs,reload=false,reloadCells=true,reloadCellAnalysis=true);
% loadSlicesDMD(epochs,reload=true);

%% Plot DMD results

close all;
cells = analyzeSlice_DMD(expPath,resultsPathDate='newest',...
                        timeRange=timeRange,...
                        redStim=true,...
                        plotSearch=true,plotPairs=true,...
                        savePNG=false,savePDF=true,saveFIG=false);
return

%% Plot search summary

c = 3; searchIdx = 4; pairIdx = 1; close all;
disp(['Ongoing: plotting searches for cell',num2str(c)]);
curCell = cells(cells.Cell == c,:);
searchPerCell = length(curCell.Vhold{1});

traceColor = [127 182 227]./255;
whiteBlue = getColormap([t238 240 241],[48 154 209],500,midCol=[150 198 227]);

for searchIdx = 1:searchPerCell
    analyzeDMDSearch(curCell,searchIdx,...
                 redStim=false,...
                 color=traceColor,...
                 colormap=whiteBlue,...
                 savePNG=false,savePDF=true,saveFIG=false);
end

% analyzeDMDSearchPair(curCell,pairIdx,...
%              redStim=true,...
%              savePNG=false,savePDF=true,saveFIG=false);

%% Test

c = 11; searchIdx = 3; close all;
disp(['Ongoing: plotting searches for cell',num2str(c)]);
curCell = cells(cells.Cell == c,:);

curCell = getHotspots(curCell,);