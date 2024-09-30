%% Shun_analyzeCellEI

% 2024/09/25

%% Define data path

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
parentPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
resultPath = '/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/Combined';

[~,~,~,~,~,~,bluePurpleRed] = loadColors;

%% (Optional) Get cell table

% combined_cells = getCellTable(combined_epochs,save=true,...
%                          timeRange=[-10,50]);
% 
% save(fullfile(resultPath,'combined_cells.mat'),"combined_cells","-v7.3");
% save(fullfile(resultPath,'combined_epochs.mat'),"combined_epochs","-v7.3");

%% Extract EI information from combined_epochs

nCells = size(combined_cells,1);

EPSC_peaks = cellfun(@(x) x.summary.EPSC.peakAvg, combined_cells.Stats);
IPSC_peaks = cellfun(@(x) x.summary.IPSC.peakAvg, combined_cells.Stats);

EPSC_aucs = cellfun(@(x) x.summary.EPSC.aucAvg, combined_cells.Stats);
IPSC_aucs = cellfun(@(x) x.summary.IPSC.aucAvg, combined_cells.Stats);

% EI statistics
% 1. EPSC + IPSC
EIsum_peaks = IPSC_peaks + EPSC_peaks;
EIsum_aucs = IPSC_aucs + EPSC_aucs;

% 2. EPSC/IPSC
EIratio_peaks = abs(EPSC_peaks) ./ abs(IPSC_peaks);
EIratio_aucs = abs(EPSC_aucs) ./ abs(IPSC_aucs);

% 3. EI index (-1: all excitatory & 1: all inhibitory)
EIindex_peaks = (abs(IPSC_peaks)-abs(EPSC_peaks)) ./ (abs(IPSC_peaks)+abs(EPSC_peaks));
EIindex_aucs = (abs(IPSC_aucs)-abs(EPSC_aucs)) ./ (abs(IPSC_aucs)+abs(EPSC_aucs));


% Remove from analysis if both currents are less than 5pA
combined_cells.Analyze = combined_cells.Learned;
removeIdx = find(abs(EPSC_peaks) < 0 & abs(IPSC_peaks) < 0);
combined_cells.Analyze(removeIdx) = 0;

% Set task range
randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardIdx = find(strcmpi('Reward pairing',combined_cells.Task) & combined_cells.Analyze);
punishIdx = find(strcmpi('Punish pairing',combined_cells.Task) & combined_cells.Analyze);

%% Plot multiple conditions

close all;

% Select data to plot
groupIdx = {randomIdx,rewardIdx,punishIdx}; 
plotGroup = [1,1,1];
groupColors = {[.7 .7 .7],bluePurpleRed(1,:),bluePurpleRed(end,:)};

fig = initializeFig(1,1); master = tiledlayout(2,4);
master.TileSpacing = 'compact'; master.Padding = 'compact';
dotSize = 200;

% Plot peak EPSC vs IPSC
nexttile(master,1);
for group = 1:length(plotGroup)
    curGroup = length(plotGroup) - group + 1;
    if plotGroup(curGroup)
        idx = groupIdx{curGroup};
        groupColor = groupColors{curGroup};
        scatter(abs(EPSC_peaks(idx)),abs(IPSC_peaks(idx)),dotSize,groupColor,"filled"); hold on; 
    end
end
diagonal = refline(1,0); diagonal.Color=[.9 .9 .9]; diagonal.LineWidth=4; diagonal.LineStyle='--';
set(gca, 'Children', flipud(get(gca,'Children'))); % Flip drawing order
xlabel('EPSC amplitude (pA)');
ylabel('IPSC amplitude (pA)');

% Plot AUC EPSC vs IPSC
nexttile(master,1+4);
for group = 1:length(plotGroup)
    curGroup = length(plotGroup) - group + 1;
    if plotGroup(curGroup)
        idx = groupIdx{curGroup};
        groupColor = groupColors{curGroup};
        scatter(abs(EPSC_aucs(idx)),abs(IPSC_aucs(idx)),dotSize,groupColor,"filled"); hold on; 
    end
end
diagonal = refline(1,0); diagonal.Color=[.9 .9 .9]; diagonal.LineWidth=4; diagonal.LineStyle='--';
set(gca, 'Children', flipud(get(gca,'Children'))); % Flip drawing order
xlabel('EPSC charge (pC)');
ylabel('IPSC charge (pC)');

% Plot amplitude EPSC+IPSC index
nexttile(master,2);
plotDistribution(EIsum_peaks,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                masterlayout=master,tile=2,xlabel='EPSC+IPSC (pA)');

% Plot charge EPSC+IPSC index
nexttile(master,6);
plotDistribution(EIsum_aucs,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                masterlayout=master,tile=6,xlabel='EPSC+IPSC (pC)');

% Plot amplitude EPSC/IPSC index
nexttile(master,3);
plotDistribution(EIratio_peaks,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                masterlayout=master,tile=3,xlabel='EPSC/IPSC (amplitude)');

% Plot charge EPSC/IPSC index
nexttile(master,7);
plotDistribution(EIratio_aucs,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                masterlayout=master,tile=7,xlabel='EPSC/IPSC (charge)');

% Plot EI amplitude index
nexttile(master,4);
plotDistribution(EIindex_peaks,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                masterlayout=master,tile=4,xlabel='EI amplitude index');

% Plot EI charge index
nexttile(master,8);
plotDistribution(EIindex_aucs,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                masterlayout=master,tile=8,xlabel='EI charge index');

saveFigures(fig,'Summary',resultPath,...
            saveFIG=true,savePDF=true,savePNG=true);