% Shun_loadSliceData_DMD

% Load and analyze DMD slice data (single session)

%% Load sessions
clear; close all;
addpath(genpath(osPathSwitch('\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Shun\Analysis\NeuroDAP')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/wengang/Exp_withShun/');
parentPath = osPathSwitch('\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Paolo');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
saveDataPath = 'default';

% Set comman params
[sessionParams,canceled] = inputSessionParams_singleSlice(expPath,...
                                paradigm=2,redStim=true,...
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

% For each search depth, it separately save a file that stores all the
% params and also responses, QC, and statistics of that search depth
close all;

loadSlicesDMD(epochs,reload=false,reloadCells=true,reloadCellAnalysis=true);
% loadSlicesDMD(epochs,reload=true);

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

cellIdx = 3; searchIdx = 1; pairIdx = 1; close all;
disp(['Ongoing: plotting searches for cell',num2str(cellIdx)]);
curCell = cells(cells.Cell == cellIdx,:); 
searchPerCell = length(curCell.Vhold{1});

traceColor = [127 182 227]./255;
whiteBlue = getColormap([238 240 241],[48 154 209],500,midCol=[150 198 227]);

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