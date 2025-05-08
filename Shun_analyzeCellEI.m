%% Shun_analyzeCellEI

% 2024/09/25

%% Define data path

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

rootPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/Combined');
[~,~,~,~,blueWhiteRed,~,bluePurpleRed] = loadColors;
today = char(datetime('today','Format','yyyyMMdd'));

% Select cell and DA trend data
expPath = uipickfiles('FilterSpec',rootPath,'Prompt','Select combined mat files');
DAtrend_path = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results/Summary-LastSessions'),...
                           'Prompt','Select DAtrend.mat')';

% Load combined_cells and combined_epochs
for i = 1:length(length(expPath))
    dirsplit = split(expPath{i},filesep); 
    filename = dirsplit{end};
    disp(['Loading: ', filename]);
    load(fullfile(expPath{i}));
    disp(['Finished: ',filename, ' loaded']);
end
resultPath = fileparts(expPath{1});
disp(strcat('resultPath: ',resultPath));

% Load DAtrend struct
load(DAtrend_path{1});
disp('Finished: DAtrend loaded');

%% (Optional) change included
% If all the sweeps included is 0, change all to 1

combined_cells = changeIncluded(combined_cells);

%% (Optional) Add epochs to combined_epochs

combined_epochs = [combined_epochs;epochs];
disp('Finished: epochs added to combined_epochs.mat');

%% (Optional) Extract response trace

close all;
animalRange = 'SL317';
[~,~] = getResponseTraces(combined_cells,animalRange=animalRange,plot=true);

%% (Optional) Get cell table

combined_cells = getCellTable(combined_epochs,save=true,...
                         timeRange=[-10,50]);

combined_cells.Learned = ones(size(combined_cells,1),1);
notLearnedAnimals = {'SL044','SL232','SL212','SL254','SL271'}; 
notLearnedIdx = find(ismember(combined_cells.Animal,notLearnedAnimals));
combined_cells.Analyze = logical(combined_cells.Learned);
combined_cells.Analyze(notLearnedIdx) = 0;

disp('Saving combined_cells and combined_epochs...');
rootPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/Combined');
resultPath = strcat(rootPath,filesep,today);
if ~isfolder(resultPath); mkdir(resultPath); end
save(fullfile(resultPath,strcat('combined_cells_',today,'.mat')),"combined_cells","-v7.3");
save(fullfile(resultPath,strcat('combined_epochs_',today,'.mat')),"combined_epochs","-v7.3");
nCells = size(combined_cells,1);
disp(strcat("Finished: saving ",num2str(nCells),' cells to combined_cells'));

%% Fig 3: select sample cells for each condition

randomSampleAnimal = {'SL255','SL285','SL333'}; randomSampleCell = [7,2,10];
rewardSampleAnimal = {'SL344','SL343','SL354'}; rewardSampleCell = [5,1,4];
punishSampleAnimal = {'SL207','SL207','SL213'}; punishSampleCell = [3,6,10];

randomCellIdx = cellfun(@(anm,cnum) ...
                find(strcmp(combined_cells.Animal, anm) & combined_cells.Cell == cnum), ...
                randomSampleAnimal, num2cell(randomSampleCell));
rewardCellIdx = cellfun(@(anm,cnum) ...
                find(strcmp(combined_cells.Animal, anm) & combined_cells.Cell == cnum), ...
                rewardSampleAnimal, num2cell(rewardSampleCell));
punishCellIdx = cellfun(@(anm,cnum) ...
                find(strcmp(combined_cells.Animal, anm) & combined_cells.Cell == cnum), ...
                punishSampleAnimal, num2cell(punishSampleCell));

randomColor = [.7 .7 .7];
rewardColor = [0 158 115]./255;
punishColor = [135 104 247]./255;

close all; initializeFig(0.5,1); t = tiledlayout('flow');
t.TileSpacing = 'none';
t.Padding = 'none';
[~,~] = getResponseTraces(combined_cells,cellList=randomCellIdx,color=randomColor,...
                          plot=true,ylim=[-1550,1000],plotIndividual=true,newFig=false);
[~,~] = getResponseTraces(combined_cells,cellList=rewardCellIdx,color=rewardColor,...
                          plot=true,ylim=[-1500,1550],plotIndividual=true,newFig=false);
[~,~] = getResponseTraces(combined_cells,cellList=punishCellIdx,color=punishColor,...
                          plot=true,ylim=[-1100,500],plotIndividual=true,newFig=false);

% saveFigures(gcf,'Fig3-SampleCells',resultPath,savePDF=true,savePNG=true);

%% Define reversal animals (fig 3)

rewardAnimals = {'SL206','SL343','SL344','SL345','SL353','SL354','SL355'};
punishAnimals = {'SL207','SL213','SL253','SL342','SL044'};

randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardIdx = find(ismember(combined_cells.Animal, rewardAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
punishIdx = find(ismember(combined_cells.Animal, punishAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Reward pairing',combined_cells.Task) & ~combined_cells.Analyze)];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Punish pairing',combined_cells.Task) & ~combined_cells.Analyze)];
baselineIdx = [randomIdx; rewardCtrlIdx; punishCtrlIdx];


%% Plot cell EI based on tasks
% Make sure to run the block above first!

figureName = 'CellEI-Reversal';

% Select data to plot
groupIdx = {randomIdx,rewardIdx,punishIdx,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [1,1,1,0,0];
groupNames = {'Random','Reward','Punish','Reward Ctrl','Punish Ctrl'};

% Set color
opacity = 0.5;
rewardColor = [0 158 115]./255; rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = [135 104 247]./255; punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};

% Plot figure
close all;
plotCellEI(combined_cells,groupIdx,nboot=10000,centerEI=false,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);

%% (Analysis) Calculate animal avg EI charge & amplitude index

% EPSC_peaks = cellfun(@(x) x.summary.EPSC.peakAvg, combined_cells.Stats);
% IPSC_peaks = cellfun(@(x) x.summary.IPSC.peakAvg, combined_cells.Stats);
% EPSC_aucs = cellfun(@(x) x.summary.EPSC.aucAvg, combined_cells.Stats);
% IPSC_aucs = cellfun(@(x) x.summary.IPSC.aucAvg, combined_cells.Stats);

EPSC_stats = getCellStats(combined_cells,'exci',type='min',average=true);
IPSC_stats = getCellStats(combined_cells,'inhi',type='max',average=true);
EPSC_peaks = vertcat(EPSC_stats.response);
IPSC_peaks = vertcat(IPSC_stats.response);

EPSC_stats = getCellStats(combined_cells,'exci',type='auc',average=true);
IPSC_stats = getCellStats(combined_cells,'inhi',type='auc',average=true);
EPSC_aucs = vertcat(EPSC_stats.response);
IPSC_aucs = vertcat(IPSC_stats.response);

% EI statistics
% 1. EPSC + IPSC
EIsum_peaks = IPSC_peaks + EPSC_peaks;
EIsum_aucs = IPSC_aucs + EPSC_aucs;

% 2. EPSC/IPSC
% EIratio_peaks = abs(EPSC_peaks) ./ abs(IPSC_peaks);
% EIratio_aucs = abs(EPSC_aucs) ./ abs(IPSC_aucs);

% 3. EI index (-1: all excitatory & 1: all inhibitory)
EIindex_peaks = (abs(IPSC_peaks)-abs(EPSC_peaks)) ./ (abs(IPSC_peaks)+abs(EPSC_peaks));
EIindex_aucs = (abs(IPSC_aucs)-abs(EPSC_aucs)) ./ (abs(IPSC_aucs)+abs(EPSC_aucs));

% 4. z-scored EI index
bootstrapping = @(B,nBoot) reshape(B(randi(numel(B), numel(B), nBoot)), [], 1);
centering = @(X,B) (X - mean(B)) ./ std(B);
% Step 1: calculate baseline distribution
randomIdx = strcmpi('Random',combined_cells.Task) & combined_cells.Analyze;
baselineEPSC_peaks = bootstrapping(EPSC_peaks(randomIdx), 5000); %EPSC_peaks(randomIdx);
baselineIPSC_peaks = bootstrapping(IPSC_peaks(randomIdx), 5000);%IPSC_peaks(randomIdx);
baselineEPSC_aucs = bootstrapping(EPSC_aucs(randomIdx), 5000);%EPSC_aucs(randomIdx);
baselineIPSC_aucs = bootstrapping(IPSC_aucs(randomIdx), 5000);%IPSC_aucs(randomIdx);
% Step 2: calculate zscore of EPSC and IPSC stats based on baseline
EPSC_peaks_centered = centering(EPSC_peaks,baselineEPSC_peaks);
IPSC_peaks_centered = centering(IPSC_peaks,baselineIPSC_peaks);
EPSC_aucs_centered = centering(EPSC_aucs,baselineEPSC_aucs);
IPSC_aucs_centered = centering(IPSC_aucs,baselineIPSC_aucs);
% Step 3: calculate EI index
EIindex_peaks_centered = (abs(IPSC_peaks_centered)-abs(EPSC_peaks_centered)) ./ (abs(IPSC_peaks_centered)+abs(EPSC_peaks_centered));
EIindex_aucs_centered = (abs(IPSC_aucs_centered)-abs(EPSC_aucs_centered)) ./ (abs(IPSC_aucs_centered)+abs(EPSC_aucs_centered));


% Calculate animal average
animalList = unique({DAtrend.animal})';
animalEIindex_peaks = cellfun(@(a) mean(EIindex_peaks(strcmp(combined_cells.Animal, a))), animalList);
animalEIindex_aucs = cellfun(@(a) mean(EIindex_aucs(strcmp(combined_cells.Animal, a))), animalList);

% Calculate cell count per animal
cellsWithDAIdx = ismember(combined_cells.Animal, animalList);
cellsWithDA = combined_cells(cellsWithDAIdx,:);
nCells = groupcounts(cellsWithDA.Animal);

%% Fig 4: compare last 30 trial, first 30 trial, value

trialWindow = 20;
endingSlope = zeros(length(DAtrend),1); 
startingSlope = zeros(length(DAtrend),1);
sessionValue = zeros(length(DAtrend),1);

for a = 1:length(DAtrend)
    % Get last n trial slope
    startTrial = max(1,DAtrend(a).nTrials-trialWindow);
    endingSlope(a) = DAtrend(a).amp.slopeMap.smoothed.map(startTrial,end);

    % Get first n trial slope
    endingTrial = min(DAtrend(a).nTrials,trialWindow);
    startingSlope(a) = DAtrend(a).amp.slopeMap.smoothed.map(1,endingTrial);

    % Get session average DA value
    sessionValue(a) = mean(DAtrend(a).amp.raw,'omitnan');
end


% Plot figure
close all;
initializeFig(0.8,1); 
master = tiledlayout('flow');
endingColor = [0, 153, 153] / 255; %[128, 179, 255]./255;
startingColor = [102, 204, 204] / 255; %[152, 201, 163]./255;
valueColor = [.7 .7 .7];
sweepColor = [184 140 1]/255;
LineWidth = 5;
nboot = 1000;

% Plot first trial slope vs animal EI
X = animalEIindex_peaks; Y = startingSlope; color = startingColor;
[model,p_value,fit_boot] = fitScatter(X,Y,type='weighted',weights=nCells); slope = model(1);
[r_spearman_starting, p_spearman_starting] = corr(X(:), Y(:), 'Type', 'Spearman');

nexttile(master,1); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 1; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;

nexttile(children,1,[1 1]);
histogram(fit_boot(:,1),200,EdgeColor=addOpacity(color,0.5),FaceColor=addOpacity(color,0.5)); hold on; 
xline(slope,'-',['p=',num2str(p_value)],Color=color,LineWidth=3);
box off; c = gca; c.YAxis(1).Visible = 'off';

nexttile(children,2,[3 1]);
scatter(X,Y,200,color,"filled"); hold on;
xFit = linspace(min(X), max(X), 100);
yFit = polyval(model, xFit);
plot(xFit, yFit, color=color, LineWidth=LineWidth);
xlabel('Animal EI index');
ylabel('Distant DA slope');



% Plot last trial slope vs animal EI
X = animalEIindex_peaks; Y = endingSlope; color = endingColor;
[model,p_value,fit_boot] = fitScatter(X,Y,type='weighted',weights=nCells); slope = model(1);
[r_spearman_ending, p_spearman_ending] = corr(X(:), Y(:), 'Type', 'Spearman');

nexttile(master,2); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 2; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;

nexttile(children,1,[1 1]);
histogram(fit_boot(:,1),200,EdgeColor=addOpacity(color,0.5),FaceColor=addOpacity(color,0.5)); hold on; 
xline(slope,'-',['p=',num2str(p_value)],Color=color,LineWidth=3);
box off; c = gca; c.YAxis(1).Visible = 'off';

nexttile(children,2,[3 1]);
scatter(X,Y,200,color,"filled"); hold on;
xFit = linspace(min(X), max(X), 100);
yFit = polyval(model, xFit);
plot(xFit, yFit, color=color, LineWidth=LineWidth);
xlabel('Animal EI index');
ylabel('Recent DA slope');



% Plot average trial value vs animal EI
X = animalEIindex_peaks; Y = sessionValue; color = valueColor;
[model,p_value,fit_boot] = fitScatter(X,Y,type='weighted',weights=nCells); slope = model(1);
[r_spearman_value, p_spearman_value] = corr(X(:), Y(:), 'Type', 'Spearman');

nexttile(master,3); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 3; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;

nexttile(children,1,[1 1]);
histogram(fit_boot(:,1),200,EdgeColor=addOpacity(color,0.5),FaceColor=addOpacity(color,0.5)); hold on; 
xline(slope,'-',['p=',num2str(p_value)],Color=color,LineWidth=3);
box off; c = gca; c.YAxis(1).Visible = 'off';

nexttile(children,2,[3 1]);
scatter(X,Y,200,color,"filled"); hold on;
xFit = linspace(min(X), max(X), 100);
yFit = polyval(model, xFit);
plot(xFit, yFit, color=color, LineWidth=LineWidth);
xlabel('Animal EI index');
ylabel('Average DA amplitude');


nexttile(master,4);
plotScatterBar(1,r_spearman_starting,style='bar',color=startingColor,LineWidth=LineWidth);
text(2,r_spearman_starting,['p = ',num2str(p_spearman_starting)]);
plotScatterBar(2,r_spearman_ending,style='bar',color=endingColor,LineWidth=LineWidth);
text(1,r_spearman_ending,['p = ',num2str(p_spearman_ending)]);
plotScatterBar(3,r_spearman_value,style='bar',color=valueColor,LineWidth=LineWidth);
text(3,r_spearman_value,['p = ',num2str(p_spearman_value)]);
xticks([1 2 3]);
xticklabels({'Distant DA slope','Recent DA slope','Average DA amplitude'});
ylabel('Spearman correlation coefficient');


% See trial window effects
trialWindowLength = 10:80;
ending_slopes = zeros(length(trialWindowLength),1);
ending_pvals = zeros(length(trialWindowLength),1);
ending_spearman = zeros(length(trialWindowLength),1);
ending_spearman_p = zeros(length(trialWindowLength),1);
starting_slopes = zeros(length(trialWindowLength),1);
starting_pvals = zeros(length(trialWindowLength),1);
starting_spearman = zeros(length(trialWindowLength),1);
starting_spearman_p = zeros(length(trialWindowLength),1);
ending_vals = zeros(length(trialWindowLength),1);
ending_vals_pvals = zeros(length(trialWindowLength),1);
ending_spearman_val = zeros(length(trialWindowLength),1);
ending_spearman_val_p = zeros(length(trialWindowLength),1);
starting_vals = zeros(length(trialWindowLength),1);
starting_vals_pvals = zeros(length(trialWindowLength),1);
starting_spearman_val = zeros(length(trialWindowLength),1);
starting_spearman_val_p = zeros(length(trialWindowLength),1);

for i = 1:length(trialWindowLength)
    % Get slopes
    trialWindow = trialWindowLength(i);
    endingSlope = zeros(length(DAtrend),1); 
    startingSlope = zeros(length(DAtrend),1);
    sessionValue = zeros(length(DAtrend),1);
    for a = 1:length(DAtrend)
        % Get last n trial slope
        startTrial = max(1,DAtrend(a).nTrials-trialWindow);
        endingSlope(a) = DAtrend(a).amp.slopeMap.raw.map(startTrial,end);
        % Get first n trial slope
        endingTrial = min(DAtrend(a).nTrials,trialWindow);
        startingSlope(a) = DAtrend(a).amp.slopeMap.raw.map(1,endingTrial);

        % Get last n trial slope
        startTrial = max(1,DAtrend(a).nTrials-trialWindow);
        endingVal(a) = mean(DAtrend(a).amp.raw(startTrial:end));
        % Get first n trial slope
        endingTrial = min(DAtrend(a).nTrials,trialWindow);
        startingVal(a) = mean(DAtrend(a).amp.raw(1:endingTrial));
    end

    % Calculate last n trial fitted slopes
    Y = endingSlope;
    [model,ending_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    ending_slopes(i) = model(1);
    [ending_spearman(i),ending_spearman_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');

    % Calculate first n trial fitted slopes
    Y = startingSlope;
    [model,starting_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    starting_slopes(i) = model(1);
    [starting_spearman(i),starting_spearman_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');

    % Calculate last n trial fitted slopes
    Y = endingVal;
    [model,ending_vals_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    ending_vals(i) = model(1);
    [ending_spearman_val(i),ending_spearman_val_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');

    % Calculate first n trial fitted slopes
    Y = startingVal;
    [model,starting_vals_pvals(i),~] = fitScatter(X,Y,weights=nCells);
    starting_vals(i) = model(1);
    [starting_spearman_val(i),starting_spearman_val_p(i)] = corr(X(:), Y(:), 'Type', 'Spearman');
end

% Plot slope vs window
nexttile(master,5);
plot(trialWindowLength,ending_slopes,LineWidth=LineWidth,color=endingColor); hold on;
plot(trialWindowLength,starting_slopes,LineWidth=LineWidth,color=startingColor); hold on;
% plot(trialWindowLength,ending_vals,LineWidth=LineWidth,color=[.5 .5 .5]); hold on;
% plot(trialWindowLength,starting_vals,LineWidth=LineWidth,color=valueColor); hold on;
xlabel('trial window length');
ylabel('fitted slope');
box off

nexttile(master,6);
plot(trialWindowLength,ending_pvals,LineWidth=LineWidth,color=endingColor); hold on;
plot(trialWindowLength,starting_pvals,LineWidth=LineWidth,color=startingColor); hold on;
plot(trialWindowLength,ending_vals_pvals,LineWidth=LineWidth,color=[.5 .5 .5]); hold on;
plot(trialWindowLength,starting_vals_pvals,LineWidth=LineWidth,color=valueColor); hold on;
yline(0.05,'--',color=[0.7 0.7 0.7],LineWidth=LineWidth);
xlabel('trial window length');
ylabel('p value');
box off

% Plot correlation vs window
nexttile(master,7);
plot(trialWindowLength,ending_spearman,LineWidth=LineWidth,color=endingColor); hold on;
plot(trialWindowLength,starting_spearman,LineWidth=LineWidth,color=startingColor); hold on;
plot(trialWindowLength,ending_spearman_val,LineWidth=LineWidth,color=[.5 .5 .5]); hold on;
plot(trialWindowLength,starting_spearman_val,LineWidth=LineWidth,color=valueColor); hold on;
xlabel('trial window length');
ylabel('Spearman correlation coefficient');
box off

nexttile(master,8);
plot(trialWindowLength,ending_spearman_p,LineWidth=LineWidth,color=endingColor); hold on;
plot(trialWindowLength,starting_spearman_p,LineWidth=LineWidth,color=startingColor); hold on;
plot(trialWindowLength,ending_spearman_val_p,LineWidth=LineWidth,color=[.5 .5 .5]); hold on;
plot(trialWindowLength,starting_spearman_val_p,LineWidth=LineWidth,color=valueColor); hold on;
yline(0.05,'--',color=[0.7 0.7 0.7],LineWidth=LineWidth);
xlabel('trial window length');
ylabel('p value');
box off

% Plot moving 20 trial window
sweepWindow = 20;
nWindows  = 80;                                % # of windows you want
nAnimals  = numel(DAtrend);
movSlopes = nan(nAnimals, nWindows);
movAvgs   = nan(nAnimals, nWindows);
% precompute your slope‐kernel once
xc = (1:sweepWindow)' - mean(1:sweepWindow);
b  = flipud(xc)/sum(xc.^2);
for a = 1:nAnimals
    raw = DAtrend(a).amp.smoothed(:); % make sure it’s a column
    % truncate or pad to exactly nWindow samples
    data = raw(1 : min(end,nWindows));
    data(end+1:nWindows,1) = nan;
    
    movSlopes(a,:) = filter(b, 1, data);
    movAvgs(a,:)   = movmean(data, sweepWindow);
end

sweeping_slopes = zeros(length(nWindows),1);
sweeping_pvals = zeros(length(nWindows),1);
sweeping_spearman = zeros(length(nWindows),1);
sweeping_spearman_p = zeros(length(nWindows),1);
sweeping_vals = zeros(length(nWindows),1);
sweeping_vals_p = zeros(length(nWindows),1);
sweeping_vals_spearman = zeros(length(nWindows),1);
sweeping_vals_spearman_p = zeros(length(nWindows),1);

for w = 1:nWindows
    % Calculate fitted slopes with animal EI
    Y = movSlopes(:,w);
    [model,sweeping_pvals(w),~] = fitScatter(X,Y,weights=nCells);
    sweeping_slopes(w) = model(1);
    [X_clean,Y_clean] = removeNaNs(X,Y);
    [sweeping_spearman(w),sweeping_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
    % Calculate average val with animal EI
    Y = movAvgs(:,w);
    [model,sweeping_vals_p(w),~] = fitScatter(X,Y,weights=nCells);
    sweeping_vals(w) = model(1);
    [X_clean,Y_clean] = removeNaNs(X,Y);
    [sweeping_vals_spearman(w),sweeping_vals_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
end

nexttile(master,9);
plot(1:nWindows,sweeping_slopes,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals,LineWidth=LineWidth,color=valueColor); hold on;
xlabel('trials');
ylabel('fitted slopes');
box off

nexttile(master,10);
plot(1:nWindows,sweeping_pvals,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals_p,LineWidth=LineWidth,color=valueColor); hold on;
yline(0.05,'--',color=[0.7 0.7 0.7],LineWidth=LineWidth);
xlabel('trials');
ylabel('p value');
box off

nexttile(master,11);
plot(1:nWindows,sweeping_spearman,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals_spearman,LineWidth=LineWidth,color=valueColor); hold on;
xlabel('trials');
ylabel('correlation coefficient');
box off

nexttile(master,12);
plot(1:nWindows,sweeping_spearman_p,LineWidth=LineWidth,color=sweepColor); hold on;
plot(1:nWindows,sweeping_vals_spearman_p,LineWidth=LineWidth,color=valueColor); hold on;
yline(0.05,'--',color=[0.7 0.7 0.7],LineWidth=LineWidth);
xlabel('trials');
ylabel('p value');
box off

% saveFigures(gcf,'Fig4-all',resultPath,saveFIG=true,savePDF=true,savePNG=true);

%% (Fig 4) Calculate slopes of last n trials

trialWindow = 15:40;

recentSlope_max = zeros(length(DAtrend),1);
recentSlope_min = zeros(length(DAtrend),1);
recentSlope_avg = zeros(length(DAtrend),1);
recentSlope_amp = zeros(length(DAtrend),1);

distantSlope_max = zeros(length(DAtrend),1);
distantSlope_min = zeros(length(DAtrend),1);
distantSlope_avg = zeros(length(DAtrend),1);
distantSlope_amp = zeros(length(DAtrend),1);

sessionValue_max = zeros(length(DAtrend),1);
sessionValue_min = zeros(length(DAtrend),1);
sessionValue_avg = zeros(length(DAtrend),1);
sessionValue_amp = zeros(length(DAtrend),1);

% Calculate average slopes in trialWindow
for a = 1:length(DAtrend)
    % Get last n trial slope
    startTrial = max(1,DAtrend(a).nTrials-trialWindow);
    recentSlope_max(a) = mean(DAtrend(a).max.slopeMap.smoothed.map(startTrial,end));
    recentSlope_min(a) = mean(DAtrend(a).min.slopeMap.smoothed.map(startTrial,end));
    recentSlope_avg(a) = mean(DAtrend(a).avg.slopeMap.smoothed.map(startTrial,end));
    recentSlope_amp(a) = mean(DAtrend(a).amp.slopeMap.smoothed.map(startTrial,end));

    % Get first n trial slope
    endingTrial = min(DAtrend(a).nTrials,trialWindow);
    distantSlope_max(a) = mean(DAtrend(a).max.slopeMap.smoothed.map(1,endingTrial));
    distantSlope_min(a) = mean(DAtrend(a).min.slopeMap.smoothed.map(1,endingTrial));
    distantSlope_avg(a) = mean(DAtrend(a).avg.slopeMap.smoothed.map(1,endingTrial));
    distantSlope_amp(a) = mean(DAtrend(a).amp.slopeMap.smoothed.map(1,endingTrial));

    sessionValue_max(a) = mean(DAtrend(a).max.raw,'omitnan');
    sessionValue_min(a) = mean(DAtrend(a).min.raw,'omitnan');
    sessionValue_avg(a) = mean(DAtrend(a).avg.raw,'omitnan');
    sessionValue_amp(a) = mean(DAtrend(a).amp.raw,'omitnan');
end

% Define the eight slope vectors in a cell array and label them
slopeVectors = {recentSlope_max, recentSlope_min, recentSlope_avg, recentSlope_amp};
slopeNames = {'CueMax Slope (Short)', 'CueMin Slope (Short)', 'CueAvg Slope (Short)', 'CueAmp Slope (Short)'};

% Set color
upColor = [0, 158, 115]/255; %[122, 201, 67]/255;
downColor = [135, 104, 247]/255; %[84, 137, 45]/255; 
stableColor = [198, 156, 109]/255;

% Set group names
% numGroups = 2;
% groupNames = {'Down DA','Up DA'};
% groupColors = {downColor,upColor};

numGroups = 3;
groupNames = {'Down DA','Stable DA','Up DA'};
groupColors = {downColor,stableColor,upColor};


%% (Analysis) kmeans clustering to define up/down/stable animals

% Preallocate cell arrays to save grouping results for each slope vector
downAnimals_kmeans = cell(numel(slopeVectors),1);
stableAnimals_kmeans = cell(numel(slopeVectors),1);
upAnimals_kmeans = cell(numel(slopeVectors),1);

initializeFig(0.5,1); tiledlayout(2,2);

for i = 1:numel(slopeVectors)
    nexttile;
    data = slopeVectors{i};
    
    % Perform k-means clustering with replicates for stability
    [groupIdx, groupCenters] = kmeans(data, numGroups, Replicates=100);
    
    % Sort group centers so that:
    %   Group 1: obvious negative (lowest center)
    %   Group 2: ambiguous (middle center)
    %   Group 3: obvious positive (highest center)
    [~, sortOrder] = sort(groupCenters);
    mapping = zeros(numGroups,1);
    mapping(sortOrder) = 1:numGroups;
    groupLabels = mapping(groupIdx);
    
    % Save grouping results
    if numGroups == 2
        downAnimals_kmeans{i} = find(groupLabels == 1);
        upAnimals_kmeans{i} = find(groupLabels == 2);
    elseif numGroups == 3
        downAnimals_kmeans{i} = find(groupLabels == 1);
        stableAnimals_kmeans{i} = find(groupLabels == 2);
        upAnimals_kmeans{i} = find(groupLabels == 3);
    end
    
    hold on;
    for j = 1:numGroups
        idx = groupLabels == j;
        scatter(find(idx), data(idx), 200, groupColors{j}, 'filled', DisplayName=groupNames{j});
    end
    hold off;
    
    title(slopeNames{i});
    xlabel('Animal');
    ylabel('Slope');
    legend('Location','southeast');
end
sgtitle('K-means classification');
saveFigures(gcf,'K-means',resultPath,saveFIG=true,savePDF=true,savePNG=true);

%% (Analysis) Plot cell EI clustered by kmeans

% Set up getCellIndices function
getCellIndices = @(indices) cell2mat(arrayfun(@(idx) find(strcmpi(combined_cells.Animal, animalList{idx})), indices, 'UniformOutput', false));

for i = 4:4%numel(slopeNames)

    % Get the animal indices from quant
    animalIdxDown = downAnimals_kmeans{i}; 
    cellIdxDown   = getCellIndices(animalIdxDown);
    animalIdxUp = upAnimals_kmeans{i};
    cellIdxUp     = getCellIndices(animalIdxUp);
    if numGroups == 3
        animalIdxStable = stableAnimals_kmeans{i}; 
        cellIdxStable = getCellIndices(animalIdxStable);
    end
    
    if numGroups == 2
        groupIdx = {cellIdxDown,cellIdxUp};
        plotGroup = [1, 1];
    elseif numGroups == 3
        groupIdx = {cellIdxStable,cellIdxDown,cellIdxUp};
        plotGroup = [1, 1, 1];
        groupNames = {'Stable DA','Down DA','Up DA'};
        groupColors = {stableColor,downColor,upColor};
    end
    
    % Create a figure name based on the current criterion (e.g., "CellEI-CueMaxSlopeShort")
    % Removing spaces and parentheses for a clean filename.
    figureName = ['CellEI-' regexprep(slopeNames{i}, '[ ()]', ''),'-kmeans'];
    close all;
    plotCellEI(combined_cells,groupIdx,centerEI=true,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);
end

%% (Fig 4 supp) plot cell EI removing reversal cells

% Find reversal cells
rewardAnimals = {'SL206','SL343','SL344','SL345','SL354','SL355'};
punishAnimals = {'SL207','SL213','SL253','SL342','SL044'};

% Exclude reversal cells
rewardIdx = find(strcmpi('Reward pairing',combined_cells.Task) & combined_cells.Analyze);
punishIdx = find(strcmpi('Punish pairing',combined_cells.Task) & combined_cells.Analyze);
rewardReversalIdx = find(ismember(combined_cells.Animal, rewardAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
punishReversalIdx = find(ismember(combined_cells.Animal, punishAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
allReversalIdx = [rewardReversalIdx;punishReversalIdx];
rewardIdx_excluded = setdiff(rewardIdx,allReversalIdx);
punishIdx_excluded = setdiff(punishIdx,allReversalIdx);
cellIdxStable_excluded = setdiff(cellIdxStable,allReversalIdx);
cellIdxUp_excluded = setdiff(cellIdxUp,allReversalIdx);
cellIdxDown_excluded = setdiff(cellIdxDown,allReversalIdx);


% Other categories
randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Reward pairing',combined_cells.Task) & ~combined_cells.Analyze)];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Punish pairing',combined_cells.Task) & ~combined_cells.Analyze)];
baselineIdx = [randomIdx; rewardCtrlIdx; punishCtrlIdx];


% Set color
opacity = 0.5;
rewardColor = [0 158 115]./255; rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = [135 104 247]./255; punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};


% Plot cell EI based on task
groupIdx = {randomIdx,rewardIdx_excluded,punishIdx_excluded,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [0,1,1,0,0];
groupNames = {'Random','Reward','Punish','Reward Ctrl','Punish Ctrl'};
% Plot figure
close all;
figureName = 'CellEI-Task-ReversalExcluded';
plotCellEI(combined_cells,groupIdx,nboot=10000,centerEI=false,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=true,figureName=figureName,resultPath=resultPath,print=true);


% Plot cell EI based on DA slope
groupIdx = {cellIdxStable_excluded,cellIdxDown_excluded,cellIdxUp_excluded};
plotGroup = [1, 1, 1];
groupNames = {'Stable DA','Down DA','Up DA'};
groupColors = {stableColor,downColor,upColor};
close all;
figureName = 'CellEI-DAslope-ReversalExcluded';
plotCellEI(combined_cells,groupIdx,centerEI=true,...
       plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
       save=true,figureName=figureName,resultPath=resultPath,print=true);

%% (Fig 4) Plot sample DA stat vs trials and animal EPSC IPSC trace

% In paper: SL066, SL068, SL316
rowIdx = [8 9 6];
dotColor = [243 186 0]/255;
earlyColor = [102 204 204]/255;
endingColor = [0 153 153]/255;

close all; initializeFig(0.4,0.7);
tiledlayout(length(rowIdx),3);

for i = 1:length(rowIdx)
    r = rowIdx(i);
    trace = DAtrend(r).amp.raw;
    trace_smoothed = DAtrend(r).amp.smoothed;
    
    nexttile(i*3-2,[1 2]);
    scatter(1:length(trace),trace,150,dotColor,'filled',MarkerFaceAlpha=0.8); hold on
    plot(trace_smoothed,color=[dotColor,0.3],LineWidth=4);
    xlim([1,length(trace)]);
    xticks([1, 20, length(trace)-19, length(trace)]);
    xticklabels({'1','20','n-19','n'}); yticks([]);
    xlabel('trials'); ylabel('DA amp');
    
    % Add shading
    yl = ylim; % get current y-limits
    % shade first 20 points
    rectangle('Position',[1,yl(1),19,yl(2)-yl(1)], ...
        'FaceColor', earlyColor, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    % shade last 20 points
    rectangle('Position',[(length(trace)-19),yl(1),20,yl(2)-yl(1)], ...
        'FaceColor', endingColor,  'FaceAlpha', 0.1, 'EdgeColor', 'none');
    % draw vertical boundary lines
    xline(1, 'LineWidth',3, 'Color', earlyColor);
    xline(20, 'LineWidth',3, 'Color', earlyColor);
    xline(length(trace)-19, 'LineWidth',3, 'Color', endingColor);
    xline(length(trace), 'LineWidth',3, 'Color', endingColor);
    % push the rectangles to the bottom so they don't cover your lines/dots
    uistack(findobj(gca,'Type','rectangle'),'bottom');
    
    % Fit and plot best fit line for the first 20 dots
    x1 = 2:19;%1:20;
    y1 = trace_smoothed(x1);
    p1 = polyfit(x1, y1, 1);
    yfit1 = polyval(p1, x1);
    plot(x1, yfit1, color=[earlyColor,0.8], LineWidth=8);
    
    % Fit and plot best fit line for the last 20 dots
    x2 = (length(trace)-18):length(trace)-1;%(length(trace)-19):length(trace);
    y2 = trace_smoothed(x2);
    p2 = polyfit(x2, y2, 1);
    yfit2 = polyval(p2, x2);
    plot(x2, yfit2, color=[endingColor,0.8], LineWidth=8);

    % Plot EPSC and IPSC of each cell
    nexttile(i*3,[1 1]);
    [~,~] = getResponseTraces(combined_cells,animalRange=animalList(r),...
                              plot=true,overlay=true,newFig=false,...
                              plotLegend=false,plotNormalized=false,plotIndividual=false);
end
saveFigures(gcf,'Fig4-schematic',resultPath,saveFIG=true,savePDF=true,savePNG=true);

%% (Fig 1 supp) For TRN patching

% Plot example cell
close all;
initializeFig(0.5,1); tiledlayout('flow');
[~,~] = getResponseTraces(combined_cells,cellList=18,...
                          plot=true,overlay=false,newFig=false,ylim=[-20,20],...
                          plotLegend=false,plotNormalized=false,plotIndividual=false);

% Plot currents
% Step 1: get all currents
epscRespAmp    = []; epscRespCharge = [];
epscBaseAmp    = []; epscBaseCharge = [];
ipscRespAmp    = []; ipscRespCharge = [];
ipscBaseAmp    = []; ipscBaseCharge = [];

for i = 1:height(combined_cells)
    % pull out the Vhold vector for this row:
    vhold = combined_cells.Vhold{i};    % e.g. [-70; 0] or [0; -70], etc.
    S     = combined_cells.Stats{i};    % struct with .response and .baseline
    
    for j = 1:numel(vhold)
        if vhold(j) == -70
            % EPSC
            epscRespAmp   (end+1,1) = mean(S.response.min{j});
            epscRespCharge(end+1,1) = mean(S.response.auc{j});    % or .auc
            epscBaseAmp   (end+1,1) = mean(S.baseline.min{j});
            epscBaseCharge(end+1,1) = mean(S.baseline.auc{j});
            
        elseif vhold(j) == 0
            % IPSC
            ipscRespAmp   (end+1,1) = mean(S.response.max{j});
            ipscRespCharge(end+1,1) = mean(S.response.auc{j});
            ipscBaseAmp   (end+1,1) = mean(S.baseline.max{j});
            ipscBaseCharge(end+1,1) = mean(S.baseline.auc{j});
        end
    end
end

exciColor = [255 157 33]/255; %[122, 201, 67]/255;
inhiColor = [71 144 253]/255; %[84, 137, 45]/255; 
ctrlColor = [.7 .7 .7];

nexttile;
plotScatterBar(1,epscRespAmp,color=exciColor);
plotScatterBar(2,epscBaseAmp,color=ctrlColor);
xticks([1,2]); xticklabels({'EPSC amp','Baseline amp'});
ylim([-50,50]); ylabel('Current (pA)');

nexttile;
plotScatterBar(1,ipscRespAmp,color=inhiColor);
plotScatterBar(2,ipscBaseAmp,color=ctrlColor);
xticks([1,2]); xticklabels({'IPSC amp','Baseline amp'});
ylim([-50,50]); ylabel('Current (pA)');

nexttile;
plotScatterBar(1,epscRespCharge,color=exciColor);
plotScatterBar(2,epscBaseCharge,color=ctrlColor);
xticks([1,2]); xticklabels({'EPSC charge','Baseline charge'});
ylim([-0.5,0.5]); ylabel('Current (pA)');

nexttile;
plotScatterBar(1,ipscRespCharge,color=inhiColor);
plotScatterBar(2,ipscBaseCharge,color=ctrlColor);
xticks([1,2]); xticklabels({'EPSC charge','Baseline charge'});
ylim([-0.5,0.5]); ylabel('Current (pA)');

%% (Fig 3 Supp) Day of patching summary

rewardColor = [0 158 115]./255;
punishColor = [135 104 247]./255;

% days = {'1st pairing, D1','1st pairing, D2','1st pairing, D3+',...
%         '2nd pairing, D1','2nd pairing, D2','2nd pairing, D3+',...
%         '3rd pairing, D1','3rd pairing, D2','3rd pairing, D3+'}; 
% y = [0 3; 1 0; 5 2;...
%      7 1; 0 0; 8 5;...
%      1 1; 0 0; 0 0];

days = {'Pairing D1','Pairing D2-3','Pairing D3+'}; 
y = [8 5; 7 5; 7 2];

b = bar(y, EdgeColor='none');
b(1).FaceColor = rewardColor;
b(2).FaceColor = punishColor;
xticklabels(days);
box off

%% (Fig 3 Supp, optional) QC summary

EPSC_qc = getCellStats(combined_cells,'exci',type='qc',average=true);
IPSC_qc = getCellStats(combined_cells,'inhi',type='qc',average=true);

% Plot Rs, Cm, Rm for each cell
EPSC_Rs = cellfun(@(x) x.Rs, {EPSC_qc.QC})'; 
EPSC_Rs_header = cellfun(@(x) x.Rs_headerString, {EPSC_qc.QC})'; 
EPSC_Rs(90) = EPSC_Rs_header(90);
EPSC_Rs_final = min(EPSC_Rs,EPSC_Rs_header);
IPSC_Rs = cellfun(@(x) x.Rs, {IPSC_qc.QC})';
IPSC_Rs_header = cellfun(@(x) x.Rs_headerString, {IPSC_qc.QC})';
IPSC_Rs_final = min(IPSC_Rs,IPSC_Rs_header);
Rs = min(EPSC_Rs_final,IPSC_Rs_final); Rs(60) = EPSC_Rs(60);

EPSC_Cm = cellfun(@(x) x.Cm, {EPSC_qc.QC})';
EPSC_Cm_header = cellfun(@(x) x.Cm_headerString, {EPSC_qc.QC})'; 
EPSC_Cm_final = min(EPSC_Cm,EPSC_Cm_header);
IPSC_Cm = cellfun(@(x) x.Cm, {IPSC_qc.QC})';
IPSC_Cm_header = cellfun(@(x) x.Cm_headerString, {IPSC_qc.QC})';
IPSC_Cm_final = min(IPSC_Cm,IPSC_Cm_header);
Cm = min(EPSC_Cm_final,IPSC_Cm_final);

EPSC_Rm = cellfun(@(x) x.Rm, {EPSC_qc.QC})';
EPSC_Rm_header = cellfun(@(x) x.Rm_headerString, {EPSC_qc.QC})'; 
EPSC_Rm_final = min(EPSC_Rm,EPSC_Rm_header);
IPSC_Rm = cellfun(@(x) x.Rm, {IPSC_qc.QC})';
IPSC_Rm_header = cellfun(@(x) x.Rm_headerString, {IPSC_qc.QC})';
IPSC_Rm_final = min(IPSC_Rm,IPSC_Rm_header);
Rm = min(EPSC_Rm_final,IPSC_Rm_final);

% Plot QC verus response
EPSCcolor = [255 157 33]/255; %[122, 201, 67]/255;
IPSCcolor = [71 144 253]/255; %[84, 137, 45]/255; 
ctrlColor = [.7 .7 .7];
dotSize = 200;

initializeFig(0.4,1); tiledlayout('flow');

Rs = Rs;
Rm = EPSC_Rm;
Cm = EPSC_Cm_final;

nexttile;
lmE   = fitlm(Rs, abs(EPSC_peaks));
slopeE = lmE.Coefficients.Estimate(2);
pvalE  = lmE.Coefficients.pValue(2);
lmI   = fitlm(Rs, abs(IPSC_peaks));
slopeI = lmI.Coefficients.Estimate(2);
pvalI  = lmI.Coefficients.pValue(2);
scatter(Rs,abs(EPSC_peaks),dotSize,EPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
scatter(Rs,abs(IPSC_peaks),dotSize,IPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
xFit = linspace(min(Rs), max(Rs), 100)';
yFitE = slopeE*xFit + lmE.Coefficients.Estimate(1);
yFitI = slopeI*xFit + lmI.Coefficients.Estimate(1);
plot(xFit, yFitE, 'Color', EPSCcolor, 'LineWidth', 5);
plot(xFit, yFitI, 'Color', IPSCcolor, 'LineWidth', 5);
xlabel('Rs'), ylabel('Amplitude');
title('Rs vs amplitude');

nexttile;
lmE   = fitlm(Rs, abs(EPSC_aucs));
slopeE = lmE.Coefficients.Estimate(2);
pvalE  = lmE.Coefficients.pValue(2);
lmI   = fitlm(Rs, abs(IPSC_aucs));
slopeI = lmI.Coefficients.Estimate(2);
pvalI  = lmI.Coefficients.pValue(2);
scatter(Rs,abs(EPSC_aucs),dotSize,EPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
scatter(Rs,abs(IPSC_aucs),dotSize,IPSCcolor,"filled",MarkerFaceAlpha=0.5); hold on;
xFit = linspace(min(Rs), max(Rs), 100)';
yFitE = slopeE*xFit + lmE.Coefficients.Estimate(1);
yFitI = slopeI*xFit + lmI.Coefficients.Estimate(1);
plot(xFit, yFitE, 'Color', EPSCcolor, 'LineWidth', 5);
plot(xFit, yFitI, 'Color', IPSCcolor, 'LineWidth', 5);
xlabel('Rs'), ylabel('Charge');
title('Rs vs charge');

% nexttile;
% scatter(Rm,abs(EPSC_peaks),dotSize,EPSCcolor,"filled"); hold on; lsline;
% % diagonal = lsline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% scatter(Rm,abs(IPSC_peaks),dotSize,IPSCcolor,"filled"); hold on; lsline;
% % diagonal = lsline; diagonal.Color=IPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% xlabel('Rm'), ylabel('Amplitude');
% title('Rm vs amplitude');
% 
% nexttile;
% scatter(Rm,abs(EPSC_aucs),dotSize,EPSCcolor,"filled"); hold on; lsline;
% % diagonal = refline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% scatter(Rm,abs(IPSC_aucs),dotSize,IPSCcolor,"filled"); hold on; lsline;
% % diagonal = refline; diagonal.Color=[.8 .8 .8]; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% xlabel('Rm'), ylabel('Charge');
% title('Rm vs charge');
% 
% nexttile;
% scatter(Cm,abs(EPSC_peaks),dotSize,EPSCcolor,"filled"); hold on; lsline;
% % diagonal = lsline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% scatter(Cm,abs(IPSC_peaks),dotSize,IPSCcolor,"filled"); hold on; lsline;
% % diagonal = lsline; diagonal.Color=IPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% xlabel('Cm'), ylabel('Amplitude');
% title('Cm vs amplitude');
% 
% nexttile;
% scatter(Cm,abs(EPSC_aucs),dotSize,EPSCcolor,"filled"); hold on; lsline;
% % diagonal = refline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% scatter(Cm,abs(IPSC_aucs),dotSize,IPSCcolor,"filled"); hold on; lsline;
% % diagonal = refline; diagonal.Color=[.8 .8 .8]; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
% xlabel('Cm'), ylabel('Charge');
% title('Cm vs charge');














%% (Analysis) Calculate relationship between DA slope and patch data
% Calculate the slope of pairwise trials vs animal avg EI index

nTrials = 70; % look at first 30 trials + last 30 trials
metric = 'slope';
skipWholeSession = true;

% Get DA vs EI map
DAvsEIaucs_slopeMap  = getDAvsEImap(DAtrend,animalEIindex_aucs,mapType='slope',...
                                    nTrials=nTrials,metric=metric,...
                                    skipWholeSession=skipWholeSession,...
                                    weights=nCells);
DAvsEIaucs_diffMap   = getDAvsEImap(DAtrend,animalEIindex_aucs,mapType='diff',...
                                    nTrials=nTrials,metric=metric,...
                                    skipWholeSession=skipWholeSession, ...
                                    weights=nCells);
DAvsEIaucs_avgMap    = getDAvsEImap(DAtrend,animalEIindex_aucs,mapType='avg',...
                                    nTrials=nTrials,metric=metric,...
                                    skipWholeSession=skipWholeSession, ...
                                    weights=nCells);

DAvsEIpeaks_slopeMap = getDAvsEImap(DAtrend,animalEIindex_peaks,mapType='slope',...
                                    nTrials=nTrials,metric=metric,...
                                    skipWholeSession=skipWholeSession, ...
                                    weights=nCells);
DAvsEIpeaks_diffMap  = getDAvsEImap(DAtrend,animalEIindex_peaks,mapType='diff',...
                                    nTrials=nTrials,metric=metric,...
                                    skipWholeSession=skipWholeSession, ...
                                    weights=nCells);
DAvsEIpeaks_avgMap   = getDAvsEImap(DAtrend,animalEIindex_peaks,mapType='avg',...
                                    nTrials=nTrials,metric=metric,...
                                    skipWholeSession=skipWholeSession, ...
                                    weights=nCells);

% save(strcat(resultPath,filesep,'DAvsEImap'),...
%     'DAvsEIaucs_slopeMap','DAvsEIaucs_diffMap','DAvsEIaucs_avgMap',...
%     'DAvsEIpeaks_slopeMap','DAvsEIpeaks_diffMap','DAvsEIpeaks_avgMap',...
%     '-v7.3');

%% test map

% close all;
initializeFig(1,1); masterLayout = tiledlayout(1,1);

nexttile(masterLayout,1); 
heatmapLayout = tiledlayout(masterLayout,1,2); 
heatmapLayout.Layout.Tile = 1;
plotDAvsEImap(DAvsEIaucs_slopeMap,statType='amp',dataType='smoothed',nTrials=50,...
              layout=heatmapLayout,startingTrialAxis='x')

%% (Analysis) Plot DA vs EI heatmap

mapsToPlot = {DAvsEIaucs_slopeMap; DAvsEIaucs_diffMap;  DAvsEIaucs_avgMap;...
             DAvsEIpeaks_slopeMap; DAvsEIpeaks_diffMap; DAvsEIpeaks_avgMap};
EItypeList = {'EI charge index'; 'EI charge index'; 'EI charge index';...
              'EI amplitude index'; 'EI amplitude index'; 'EI amplitude index'};
mapTypeList = {'slope';'diff';'avg';'slope';'diff';'avg'};
figureNameList = {'DAvsEIaucs_slopeMap'; 'DAvsEIaucs_diffMap'; 'DAvsEIaucs_avgMap';...
                  'DAvsEIpeaks_slopeMap'; 'DAvsEIpeaks_diffMap'; 'DAvsEIpeaks_avgMap'};
fields = {'max','min','avg','amp'};
nTrials = 50;

for i = 1:length(mapsToPlot)
    
    DAvsEImap = mapsToPlot{i}; figureName = figureNameList{i};
    EItype = EItypeList{i}; mapType = mapTypeList{i};
    
    initializeFig(1,1); masterLayout = tiledlayout(2,4);
    
    for f = 1:length(fields)
        cur_field = fields{f};

        nexttile(masterLayout,f); 
        heatmapLayout = tiledlayout(masterLayout,1,2); 
        heatmapLayout.Layout.Tile = f;
        plotDAvsEImap(DAvsEImap,statType=cur_field,dataType='raw',nTrials=nTrials,...
                      layout=heatmapLayout,...
                      title=['DA ',mapType,' (raw ',cur_field,') vs ', EItype]);
    
        nexttile(masterLayout,f+length(fields)); 
        heatmapLayout = tiledlayout(masterLayout,1,2); 
        heatmapLayout.Layout.Tile = f+length(fields);
        plotDAvsEImap(DAvsEImap,statType=cur_field,dataType='smoothed',nTrials=nTrials,...
                      layout=heatmapLayout,...
                      title=['DA ',mapType,' (smoothed ',cur_field,') vs ', EItype]);
    end
    
    saveFigures(gcf,figureName,resultPath,...
                saveFIG=true,savePDF=true,savePNG=true);
    close all;
end


%% Generate DAtrend_manim for plotting
DAtrend_manim = struct([]);

for a = 1:30
    DAtrend_manim(a).animal = DAtrend(a).animal;
    DAtrend_manim(a).nTrials = DAtrend(a).nTrials;
    DAtrend_manim(a).slopeMap_raw = DAtrend(a).amp.slopeMap.raw.map;
    DAtrend_manim(a).slopeMap_smoothed = DAtrend(a).amp.slopeMap.smoothed.map;
end

DAvsEImap_manim = DAvsEIpeaks_slopeMap.amp;

% Save DAtrend_manim, DAvsEImap_manim, animalEIinex_peaks
%%













%% (Old) Animals that undergo reward/punish and have significant DA trend
% increase/decrease

% Get DA trend of last 30 trials for each animal
lastTrialWindow = 30;
animalRange = {DAtrend.animal};
DAslope = zeros(length(DAtrend),1); 
DAslope_pval = zeros(length(DAtrend),1);
rewardCandidates = [];
punishCandidates = [];

for a = 1:length(DAtrend)
    startTrial = max(1,DAtrend(a).nTrials-lastTrialWindow);
    DAslope(a) = DAtrend(a).amp.slopeMap.smoothed.map(startTrial,end);
    DAslope_pval(a) = DAtrend(a).amp.slopeMap.smoothed.pval(startTrial,end);

    if DAslope(a) > 0 %&& contains(DAtrend(a).task,'Reward')
        rewardCandidates = [rewardCandidates;a];
    elseif DAslope(a) <= 0 %&& contains(DAtrend(a).task,'Punish')
        punishCandidates = [punishCandidates;a];
    end
end

% rewardAnimals = {'SL046','SL063','SL068','SL206','SL231','SL316','SL343','SL344','SL345','SL354','SL355'};
% punishAnimals = {'SL207','SL213','SL320','SL323','SL342','SL352'};

% rewardAnimals = animalRange(rewardCandidates);
% punishAnimals = animalRange(punishCandidates);


%% (Old) Define categories

% Set DA trend range
% Animals classified by DA amplitude trend in the last session
% Stable includes no amplitude change or net change is close to zero
upAnimals = {'SL043','SL063','SL068',...
             'SL206','SL208','SL229','SL231',...
             'SL316','SL317',...
             'SL321','SL322'};
stableAnimals = {'SL044','SL060','SL207','SL208',...
                 'SL212','SL213','SL229','SL231',...
                 'SL253','SL254','SL271','SL317',...
                 'SL320','SL321','SL323'};
downAnimals = {'SL046','SL062','SL064','SL066',...
               'SL208','SL232','SL253','SL254',...
               'SL320'};
overlapAnimals = {'SL208','SL229','SL231','SL253','SL254','SL317','SL320'};

upDAIdx = find(ismember(combined_cells.Animal, upAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
stableDAIdx = find(ismember(combined_cells.Animal, stableAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
downDAIdx = find(ismember(combined_cells.Animal, downAnimals) & ...
               contains(combined_cells.Task, 'pairing'));

% Set task range
% Remove not learned animals
% combined_cells.Learned = ones(size(combined_cells,1),1);
notLearnedAnimals = {'SL044','SL232','SL212','SL254','SL271'}; 
notLearnedIdx = find(ismember(combined_cells.Animal,notLearnedAnimals));
combined_cells.Analyze = logical(combined_cells.Learned);
combined_cells.Analyze(notLearnedIdx) = 0;
removeIdx = find(abs(EPSC_peaks) <= 0 & abs(IPSC_peaks) <= 0);
combined_cells.Analyze(removeIdx) = 0;

randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardIdx = find(strcmpi('Reward pairing',combined_cells.Task) & combined_cells.Analyze);
punishIdx = find(strcmpi('Punish pairing',combined_cells.Task) & combined_cells.Analyze);
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Reward pairing',combined_cells.Task) & ~combined_cells.Analyze)];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Punish pairing',combined_cells.Task) & ~combined_cells.Analyze)];


%% Plot response traces

% Extract respons trace
[EPSC_traces,IPSC_traces] = getResponseTraces(combined_cells,animalRange='All');

% Plot params
initializeFig(.7,1); tiledlayout(2,3);

plotIndividual = false;
nexttile;
plotSEM(timeRangeInms,EPSC_traces(randomIdx,:),bluePurpleRed(end,:),plotIndividual=plotIndividual,label='EPSC');
plotSEM(timeRangeInms,IPSC_traces(randomIdx,:),bluePurpleRed(1,:),plotIndividual=plotIndividual,label='IPSC');
xlabel('Time (ms)'); ylabel('Current (pA)'); legend; ylim([-400,500]);
title('Baseline');

nexttile;
plotSEM(timeRangeInms,EPSC_traces(rewardIdx,:),bluePurpleRed(end,:),plotIndividual=plotIndividual,label='EPSC');
plotSEM(timeRangeInms,IPSC_traces(rewardIdx,:),bluePurpleRed(1,:),plotIndividual=plotIndividual,label='IPSC');
xlabel('Time (ms)'); ylabel('Current (pA)'); legend; ylim([-400,500]);
title('Reward pairing');

nexttile;
plotSEM(timeRangeInms,EPSC_traces(punishIdx,:),bluePurpleRed(end,:),plotIndividual=plotIndividual,label='EPSC');
plotSEM(timeRangeInms,IPSC_traces(punishIdx,:),bluePurpleRed(1,:),plotIndividual=plotIndividual,label='IPSC');
xlabel('Time (ms)'); ylabel('Current (pA)'); legend; ylim([-400,500]);
title('Punish pairing');

plotIndividual = true;
nexttile;
plotSEM(timeRangeInms,EPSC_traces(randomIdx,:),bluePurpleRed(end,:),plotIndividual=plotIndividual,label='EPSC');
plotSEM(timeRangeInms,IPSC_traces(randomIdx,:),bluePurpleRed(1,:),plotIndividual=plotIndividual,label='IPSC');
xlabel('Time (ms)'); ylabel('Current (pA)'); legend;
title('Baseline');

nexttile;
plotSEM(timeRangeInms,EPSC_traces(rewardIdx,:),bluePurpleRed(end,:),plotIndividual=plotIndividual,label='EPSC');
plotSEM(timeRangeInms,IPSC_traces(rewardIdx,:),bluePurpleRed(1,:),plotIndividual=plotIndividual,label='IPSC');
xlabel('Time (ms)'); ylabel('Current (pA)'); legend;
title('Reward pairing');

nexttile;
plotSEM(timeRangeInms,EPSC_traces(punishIdx,:),bluePurpleRed(end,:),plotIndividual=plotIndividual,label='EPSC');
plotSEM(timeRangeInms,IPSC_traces(punishIdx,:),bluePurpleRed(1,:),plotIndividual=plotIndividual,label='IPSC');
xlabel('Time (ms)'); ylabel('Current (pA)'); legend;
title('Punish pairing');

saveFigures(gcf,'Trace',resultPath,...
            saveFIG=true,savePDF=true,savePNG=true);

%% (Old)(Analysis) Quantile based method to define up/down/stable animals

% % Preallocate cell arrays to save grouping results for each slope vector
% downAnimals_quant = cell(numel(slopeVectors),1);
% stableAnimals_quant = cell(numel(slopeVectors),1);
% upAnimals_quant = cell(numel(slopeVectors),1);
% 
% initializeFig(0.5,1); tiledlayout(2,2);
% 
% for i = 1:numel(slopeVectors)
%     nexttile;
%     data = slopeVectors{i};
% 
%     % Use quantiles to classify data into groups
%     % Save grouping results
%     if numGroups == 2
%         edges = quantile(data, [0, 1/numGroups 1]);
%         groupLabels = discretize(data, edges);
%         downAnimals_quant{i} = find(groupLabels == 1);
%         upAnimals_quant{i} = find(groupLabels == 2);
%     elseif numGroups == 3
%         edges = quantile(data, [0, 1/numGroups, 2/numGroups, 1]);
%         groupLabels = discretize(data, edges);
%         downAnimals_quant{i} = find(groupLabels == 1);
%         stableAnimals_quant{i} = find(groupLabels == 2);
%         upAnimals_quant{i} = find(groupLabels == 3);
%     end
% 
%     hold on;
%     for j = 1:numGroups
%         idx = groupLabels == j;
%         scatter(find(idx), data(idx), 200, groupColors{j}, 'filled', DisplayName=groupNames{j});
%     end
%     hold off;
% 
%     title(slopeNames{i});
%     xlabel('Animal');
%     ylabel('Slope');
%     legend('Location','southeast');
% end
% sgtitle('Quantiles classification');
% saveFigures(gcf,'Quantile',resultPath,...
%                 saveFIG=true,savePDF=true,savePNG=true);

%% (Analysis) Plot cell EI clustered by quantile

% % Set up getCellIndices function
% getCellIndices = @(indices) cell2mat(arrayfun(@(idx) find(strcmpi(combined_cells.Animal, animalList{idx})), indices, 'UniformOutput', false));
% 
% for i = 1:numel(slopeNames)
% 
%     % Get the animal indices  from quant
%     animalIdxDown = downAnimals_quant{i}; 
%     cellIdxDown   = getCellIndices(animalIdxDown);
%     animalIdxUp = upAnimals_quant{i};
%     cellIdxUp     = getCellIndices(animalIdxUp);
%     if numGroups == 3
%         animalIdxStable = stableAnimals_quant{i}; 
%         cellIdxStable = getCellIndices(animalIdxStable);
%     end
% 
%     if numGroups == 2
%         groupIdx = {cellIdxDown,cellIdxUp};
%         plotGroup = [1, 1];
%     elseif numGroups == 3
%         groupIdx = {cellIdxStable,cellIdxDown,cellIdxUp};
%         plotGroup = [1, 1, 1];
%         groupNames = {'Stable DA','Down DA','Up DA'};
%         groupColors = {stableColor,downColor,upColor};
%     end
% 
%     DAslopePerCell = getDAperCell(combined_cells,recentSlope_amp,animalList);
% 
%     % Create a figure name based on the current criterion (e.g., "CellEI-CueMaxSlopeShort")
%     % Removing spaces and parentheses for a clean filename.
%     figureName = ['CellEI-' regexprep(slopeNames{i}, '[ ()]', ''),'-quantile'];
%     close all;
%     plotCellEI(combined_cells,groupIdx,...
%            plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
%            save=true,figureName=figureName,resultPath=resultPath,print=true);
% end

%% (Fig 4) Calculate DA based on animal EI

exciColor = [255 157 33]/255; %[122, 201, 67]/255;
inhiColor = [71 144 253]/255; %[84, 137, 45]/255; 
neutralColor = [.5 .5 .5];

numGroups = 3;
groupNames = {'Exci','Neutral','Inhi'};
groupColors = {exciColor,neutralColor,inhiColor};
separationMethod = 'quant';

% data = EIindex_peaks(cellsWithDAIdx);
data = animalEIindex_peaks;

% Separate animal EI into three groups
close all; initializeFig(0.5,1); tiledlayout('flow');

% Option 1: use kmeans
[groupIdx, groupCenters] = kmeans(data, numGroups, Replicates=100);
[~, sortOrder] = sort(groupCenters);
mapping = zeros(numGroups,1);
mapping(sortOrder) = 1:numGroups;
groupLabels = mapping(groupIdx);
% Save grouping results
if numGroups == 2
    exciIdx_kmeans = find(groupLabels == 1);
    inhiIdx_kmeans = find(groupLabels == 2);
elseif numGroups == 3
    exciIdx_kmeans = find(groupLabels == 1);
    neutralIdx_kmeans = find(groupLabels == 2);
    inhiIdx_kmeans = find(groupLabels == 3);
end
nexttile; hold on;
for j = 1:numGroups
    idx = groupLabels == j;
    scatter(find(idx), data(idx), 200, groupColors{j}, 'filled', DisplayName=groupNames{j});
end
xlabel('Cells'); ylabel('EI index');
legend('Location','southeast'); title('K-means classification');

% Option 2: use quantile
if numGroups == 2
    edges = quantile(data, [0, 1/numGroups 1]);
    groupLabels = discretize(data, edges);
    exciIdx_quant = find(groupLabels == 1);
    inhiIdx_quant = find(groupLabels == 2);
elseif numGroups == 3
    edges = quantile(data, [0, 1/numGroups, 2/numGroups, 1]);
    groupLabels = discretize(data, edges);
    exciIdx_quant = find(groupLabels == 1);
    neutralIdx_quant = find(groupLabels == 2);
    inhiIdx_quant = find(groupLabels == 3);
end
nexttile; hold on;
for j = 1:numGroups
    idx = groupLabels == j;
    scatter(find(idx), data(idx), 200, groupColors{j}, 'filled', DisplayName=groupNames{j});
end
xlabel('Cells'); ylabel('EI index');
legend('Location','southeast'); title('Quantile classification');

% Find the mean DA slope, and mean DA value
recentSlopePerCell = getDAperCell(combined_cells,recentSlope_amp,animalList);
distantSlopePerCell = getDAperCell(combined_cells,distantSlope_amp,animalList);
sessionValuePerCell = getDAperCell(combined_cells,sessionValue_amp,animalList);

if strcmpi(separationMethod,'kmeans')
    avgRecentSlope_exci = recentSlopePerCell(exciIdx_kmeans);
    avgRecentSlope_inhi = recentSlopePerCell(inhiIdx_kmeans);
    avgRecentSlope_neutral = recentSlopePerCell(neutralIdx_kmeans);
    avgDistantSlope_exci = distantSlopePerCell(exciIdx_kmeans);
    avgDistantSlope_inhi = distantSlopePerCell(inhiIdx_kmeans);
    avgDistantSlope_neutral = distantSlopePerCell(neutralIdx_kmeans);
    avgSessionValue_exci = sessionValuePerCell(exciIdx_kmeans);
    avgSessionValue_inhi = sessionValuePerCell(inhiIdx_kmeans);
    avgSessionValue_neutral = sessionValuePerCell(neutralIdx_kmeans);
else
    avgRecentSlope_exci = recentSlopePerCell(exciIdx_quant);
    avgRecentSlope_inhi = recentSlopePerCell(inhiIdx_quant);
    avgRecentSlope_neutral = recentSlopePerCell(neutralIdx_quant);
    avgDistantSlope_exci = distantSlopePerCell(exciIdx_quant);
    avgDistantSlope_inhi = distantSlopePerCell(inhiIdx_quant);
    avgDistantSlope_neutral = distantSlopePerCell(neutralIdx_quant);
    avgSessionValue_exci = sessionValuePerCell(exciIdx_quant);
    avgSessionValue_inhi = sessionValuePerCell(inhiIdx_quant);
    avgSessionValue_neutral = sessionValuePerCell(neutralIdx_quant);
end

% (Fig 4) Plot bar plot
nexttile;
plotScatterBar(1,avgRecentSlope_exci,color=exciColor);
plotScatterBar(2,avgRecentSlope_neutral,color=neutralColor);
plotScatterBar(3,avgRecentSlope_inhi,color=inhiColor);
plotStats(avgRecentSlope_exci,avgRecentSlope_neutral,[1 2],testType='kstest');
plotStats(avgRecentSlope_neutral,avgRecentSlope_inhi,[2 3],testType='kstest');
plotStats(avgRecentSlope_exci,avgRecentSlope_inhi,[1 3],testType='kstest');
xticks([1 2 3]); xticklabels(groupNames); ylabel('Recent DA slope');

nexttile
plotScatterBar(1,avgDistantSlope_exci,color=exciColor);
plotScatterBar(2,avgDistantSlope_neutral,color=neutralColor);
plotScatterBar(3,avgDistantSlope_inhi,color=inhiColor);
plotStats(avgDistantSlope_exci,avgDistantSlope_neutral,[1 2],testType='kstest');
plotStats(avgDistantSlope_neutral,avgDistantSlope_inhi,[2 3],testType='kstest');
plotStats(avgDistantSlope_exci,avgDistantSlope_inhi,[1 3],testType='kstest');
xticks([1 2 3]); xticklabels(groupNames); ylabel('Distant DA slope');