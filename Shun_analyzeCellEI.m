%% Shun_analyzeCellEI

% 2024/09/25

%% Define data path

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
parentPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/');
expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');
resultPath = '/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/OptoPair-EP/';

%% (Optional) Get cell table

combined_cells = getCellTable(combined_epochs,save=true,...
                         timeRange=[-10,50]);

save(fullfile(resultPath,'combined_cells.mat'),"combined_cells","-v7.3");
save(fullfile(resultPath,'combined_epochs.mat'),"combined_epochs","-v7.3");

%% Extract EI information from combined_epochs

nCells = size(combined_cells,1);

EPSC_peaks = cellfun(@(x) x.summary.EPSC.peakAvg, combined_cells.Stats);
IPSC_peaks = cellfun(@(x) x.summary.IPSC.peakAvg, combined_cells.Stats);

EPSC_aucs = cellfun(@(x) x.summary.EPSC.aucAvg, combined_cells.Stats);
IPSC_aucs = cellfun(@(x) x.summary.IPSC.aucAvg, combined_cells.Stats);

EIindex_peaks = (abs(IPSC_peaks)-abs(EPSC_peaks)) ./ (abs(IPSC_peaks)+abs(EPSC_peaks));
EIindex_aucs = (abs(IPSC_aucs)-abs(EPSC_aucs)) ./ (abs(IPSC_aucs)+abs(EPSC_aucs));

%% Set task range



%% Plot EI scatter plot

fig = initializeFig(0.5,1); master = tiledlayout(2,4); 
master.TileSpacing = 'compact'; master.Padding = 'compact';

% Plot params
dotSize = 200;
color = [12, 173, 74]./255;

% Bootstrap data set using random index
bootIdx = round((nCells-1).*rand(5000,1) + 1);

nexttile(master,1,[1 2]);
scatter(abs(EPSC_peaks),abs(IPSC_peaks),dotSize,color,"filled");
diagonal = refline(1,0); diagonal.Color=[.7 .7 .7]; diagonal.LineWidth=4; diagonal.LineStyle='--';
set(gca, 'Children', flipud(get(gca,'Children'))); % Flip drawing order
xlabel('EPSC amplitude (pA)');
ylabel('IPSC amplitude (pA)');

nexttile(master,5,[1 2]);
scatter(abs(EPSC_aucs),abs(IPSC_aucs),dotSize,color,"filled");
diagonal = refline(1,0); diagonal.Color=[.7 .7 .7]; diagonal.LineWidth=4; diagonal.LineStyle='--';
set(gca, 'Children', flipud(get(gca,'Children'))); % Flip drawing order
xlabel('EPSC charge (pC)');
ylabel('IPSC charge (pC)');

nexttile(master,3,[1 2]); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 3; children.Layout.TileSpan = [1 2];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;
nexttile(children,1);
plotScatterBar(EIindex_peaks,1,color=color,dotSize=100,orientation='horizontal');
h = gca; h.YAxis.Visible = 'off';
nexttile(children,2,[3 1]);
histogram(EIindex_peaks(bootIdx),FaceColor=color,EdgeColor=color); 
xlabel('EI amplitude index'); ylabel('Count'); box off;
xlabel('EI peak index');

nexttile(master,7,[1 2]);
children = tiledlayout(master,4,1);
children.Layout.Tile = 7; children.Layout.TileSpan = [1 2];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;
nexttile(children,1);
plotScatterBar(EIindex_aucs,1,color=color,dotSize=100,orientation='horizontal');
h = gca; h.YAxis.Visible = 'off';
nexttile(children,2,[3 1]);
histogram(EIindex_aucs(bootIdx),FaceColor=color,EdgeColor=color); 
xlabel('EI charge index'); ylabel('Count'); box off;
xlabel('EI charge index');

saveFigures(fig,'Summary-baseline',resultPath,...
            saveFIG=true,savePDF=true,savePNG=true);