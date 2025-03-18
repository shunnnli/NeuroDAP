%% Shun_analyzeCellEI

% 2024/09/25

%% Define data path

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

rootPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/Combined');
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
today = char(datetime('today','Format','yyyyMMdd'));

% Select cell and DA trend data
expPath = uipickfiles('FilterSpec',rootPath,'Prompt','Select combined folder');
DAtrend_path = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Results/Summary-LastSessions'),...
                           'Prompt','Select DAtrend.mat')';

% Load combined_cells and combined_epochs
files = dir(fullfile(expPath{1}, 'combined_*.mat'));
for k = 1:length(files)
    filename = files(k).name;
    disp(['Loading: ', filename]);
    load(fullfile(expPath{1}, filename));
    disp(['Finished: ',filename, ' loaded']);
end
resultPath = expPath{1};
disp(strcat('resultPath: ',resultPath));

% Load DAtrend struct
load(DAtrend_path{1});
disp('Finished: DAtrend loaded');

%% (Optional) Add epochs to combined_epochs

combined_epochs = [combined_epochs;epochs];

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

%% Define categories

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

% Animals with prolonged reward
longRewardAnimals = {'SL043','SL208',...
                     'SL060','SL062','SL063','SL064','SL066','SL068',...
                     'SL208','SL229','SL231','SL232','SL316'};

% Animals with prolonged punish
longPunishAnimals = {'SL044','SL045',...
                     'SL254','SL271','SL317',...
                     'SL320','SL321','SL322','SL323'};

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


% Set DA trend range
upDAIdx = find(ismember(combined_cells.Animal, upAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
stableDAIdx = find(ismember(combined_cells.Animal, stableAnimals) & ...
               contains(combined_cells.Task, 'pairing'));
downDAIdx = find(ismember(combined_cells.Animal, downAnimals) & ...
               contains(combined_cells.Task, 'pairing'));

%% (Optional) Extract response trace

animalRange = 'SL043';
[~,~] = getResponseTraces(combined_cells,animalRange=animalRange,plot=true);

%% Plot cell EI based on tasks

figureName = 'CellEI-Task';

% Select data to plot
groupIdx = {randomIdx,rewardIdx,punishIdx,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [0,1,1,0,0];
groupNames = {'Random','Reward','Punish','Reward Ctrl','Punish Ctrl'};

% Set color
rewardColor = bluePurpleRed(1,:); rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = bluePurpleRed(end,:); punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};

% Plot figure
close all;
plotCellEI(combined_cells,groupIdx,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=true,figureName=figureName,resultPath=resultPath,print=true);

%% Plot cell EI based on DA amplitude trend

figureName = 'CellEI-DAtrend';

% Select data to plot
groupIdx = {stableDAIdx,upDAIdx,downDAIdx,randomIdx,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [1,1,1,0,0,0];
groupNames = {'Stable DA','Up DA','Down DA','Random','Reward Ctrl','Punish Ctrl'};

% Set color
upColor = bluePurpleRed(1,:); rewardCtrlColor = 1 - opacity*(1-upColor);
downColor = bluePurpleRed(end,:); punishCtrlColor = 1 - opacity*(1-downColor);
groupColors = {[.2 .2 .2],upColor,downColor,...
               [.7 .7 .7],rewardCtrlColor,punishCtrlColor};

% Plot figure
close all;
plotCellEI(combined_cells,groupIdx,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=true,figureName=figureName,resultPath=resultPath,print=true);

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

%% Plot DA vs EI charge index

color = [27 227 101]./255;

initializeFig(.6,.5); tiledlayout(1,2);

nexttile;
DA = combined_cells.Learned([rewardIdx;punishIdx]);
EI = EIindex_aucs([rewardIdx;punishIdx]);
scatter(DA,EI,dotSize,color,'filled'); hold on;
p = polyfit(DA, EI, 1);
yfit = polyval(p, DA);
plot(DA, yfit, LineWidth=3);

nexttile;
animals = combined_cells.Animal([rewardIdx;punishIdx]);
DA_animal = splitapply(@mean, DA, findgroups(animals));
EI_animal = splitapply(@mean, EI, findgroups(animals));
scatter(DA_animal,EI_animal,dotSize,color,'filled'); hold on;
p = polyfit(DA, EI, 1);
yfit = polyval(p, DA);
plot(DA, yfit, LineWidth=3);


%% Extract cell QC data

qcIdx = cellfun(@(x) find(x==-70,1), combined_cells.Vhold);
QCs = struct2table(cellfun(@(q, idx) q{idx}, combined_cells.QC, num2cell(qcIdx)));

included = cellfun(@(q, idx) q{idx}, combined_cells.Included, num2cell(qcIdx),UniformOutput=false);
included = cellfun(@(x) x + (~any(x))*ones(size(x)), included, UniformOutput=false);

Rs = cellfun(@(x,idx) mean(x(find(idx))), QCs.Rs, included);
Rm = cellfun(@(x,idx) mean(x(find(idx))), QCs.Rm, included);
Cm = cellfun(@(x,idx) mean(x(find(idx))), QCs.Cm, included);

%% Plot QC vs response

EPSCcolor = bluePurpleRed(end,:);
IPSCcolor = bluePurpleRed(1,:);

initializeFig(0.4,1); tiledlayout('flow');

nexttile;
scatter(Rs,abs(EPSC_peaks),dotSize,EPSCcolor,"filled"); hold on; lsline;
% diagonal = lsline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
scatter(Rs,abs(IPSC_peaks),dotSize,IPSCcolor,"filled"); hold on; lsline;
% diagonal = lsline; diagonal.Color=IPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
xlabel('Rs'), ylabel('Amplitude');
title('Rs vs IPSC amplitude');

nexttile;
scatter(Rs,abs(EPSC_aucs),dotSize,EPSCcolor,"filled"); hold on; lsline;
% diagonal = refline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
scatter(Rs,abs(IPSC_aucs),dotSize,IPSCcolor,"filled"); hold on; lsline;
% diagonal = refline; diagonal.Color=[.8 .8 .8]; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
xlabel('Rs'), ylabel('Charge');
title('Rs vs IPSC charge');

nexttile;
scatter(Rm,abs(EPSC_peaks),dotSize,EPSCcolor,"filled"); hold on; lsline;
% diagonal = lsline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
scatter(Rm,abs(IPSC_peaks),dotSize,IPSCcolor,"filled"); hold on; lsline;
% diagonal = lsline; diagonal.Color=IPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
xlabel('Rm'), ylabel('Amplitude');
title('Rm vs IPSC amplitude');

nexttile;
scatter(Rm,abs(EPSC_aucs),dotSize,EPSCcolor,"filled"); hold on; lsline;
% diagonal = refline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
scatter(Rm,abs(IPSC_aucs),dotSize,IPSCcolor,"filled"); hold on; lsline;
% diagonal = refline; diagonal.Color=[.8 .8 .8]; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
xlabel('Rm'), ylabel('Charge');
title('Rm vs IPSC charge');

nexttile;
scatter(Cm,abs(EPSC_peaks),dotSize,EPSCcolor,"filled"); hold on; lsline;
% diagonal = lsline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
scatter(Cm,abs(IPSC_peaks),dotSize,IPSCcolor,"filled"); hold on; lsline;
% diagonal = lsline; diagonal.Color=IPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
xlabel('Cm'), ylabel('Amplitude');
title('Cm vs IPSC amplitude');

nexttile;
scatter(Cm,abs(EPSC_aucs),dotSize,EPSCcolor,"filled"); hold on; lsline;
% diagonal = refline; diagonal.Color=EPSCcolor; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
scatter(Cm,abs(IPSC_aucs),dotSize,IPSCcolor,"filled"); hold on; lsline;
% diagonal = refline; diagonal.Color=[.8 .8 .8]; diagonal.LineWidth=4; diagonal.LineStyle='--'; hold on;
xlabel('Cm'), ylabel('Charge');
title('Cm vs IPSC charge');


%% Plot summay trace (for TRN-LHb)

%{

load('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/TRN-LHb/combined_epochs_20241003.mat');
resultPath = '/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/TRN-LHb';
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
today = char(datetime('today','Format','yyyyMMdd')); 

%}

%% (Analysis) Calculate animal avg EI charge & amplitude index

EPSC_peaks = cellfun(@(x) x.summary.EPSC.peakAvg, combined_cells.Stats);
IPSC_peaks = cellfun(@(x) x.summary.IPSC.peakAvg, combined_cells.Stats);

EPSC_aucs = cellfun(@(x) x.summary.EPSC.aucAvg, combined_cells.Stats);
IPSC_aucs = cellfun(@(x) x.summary.IPSC.aucAvg, combined_cells.Stats);

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

% Calculate animal average
animalList = unique({DAtrend.animal})';
animalEIindex_peaks = cellfun(@(a) mean(EIindex_peaks(strcmp(combined_cells.Animal, a))), animalList);
animalEIindex_aucs = cellfun(@(a) mean(EIindex_aucs(strcmp(combined_cells.Animal, a))), animalList);

%% (Analysis) Calculate relationship between DA slope and patch data
% Calculate the slope of last n trials vs animal avg EI index

nboot = 1000;
maxTrials = max(arrayfun(@(x) length(x.CueMax_slope), DAtrend));
statsTypes = {'CueMax', 'CueMin', 'CueAvg', 'CueAmp',...
              'CueMax_smoothed','CueMin_smoothed','CueAvg_smoothed','CueAmp_smoothed'};

% Preallocate structures to store slopes and p-values for each cue type and each EI measure
for st = 1:length(statsTypes)
    slopes_aucs.(statsTypes{st})  = nan(numTrials,1);
    pvals_aucs.(statsTypes{st})   = nan(numTrials,1);
    slopes_peaks.(statsTypes{st}) = nan(numTrials,1);
    pvals_peaks.(statsTypes{st})  = nan(numTrials,1);
end

% Loop over n
for n = 1:maxTrials

    % For each cue type, extract the nth slope from each animal in DAtrend
    nthSlopes = getNthSlope(DAtrend,n);

    % Skip fitting if animal number is less than 5
    validIdx = ~isnan(nthSlopes.CueMax);
    if sum(validIdx) <= 5; continue; end
    disp(['Ongoing: calculating fits for last ',num2str(n),' trials']);
    
    % For each cue type, fit regressions for both EI measures and store the results
    for st = 1:length(statsTypes)
        % Regression with animalEIindex_aucs
        Y = nthSlopes.(statsTypes{st})(validIdx); 
        X = animalEIindex_aucs(validIdx);
        fit_obs = polyfit(X, Y, 1);
        fit_boot = bootstrp(nboot, @(i) polyfit(X, Y(i), 1), 1:length(Y));
        slopes_aucs.(statsTypes{st})(n) = fit_obs(1);
        pvals_aucs.(statsTypes{st})(n)  = sum(abs(fit_boot(:,1)) >= abs(fit_obs(1)))/nboot;   % p-value for the slope

        % stats = regstats(animalEIindex_aucs, nthSlopes.(statsTypes{st}), 'linear', {'tstat'});
        % slopes_aucs.(statsTypes{st})(n) = stats.tstat.beta(2);
        % pvals_aucs.(statsTypes{st})(n)  = stats.tstat.pval(2);
        
        % Regression with animalEIindex_peaks
        Y = nthSlopes.(statsTypes{st})(validIdx); 
        X = animalEIindex_peaks(validIdx);
        fit_obs = polyfit(X, Y, 1);
        fit_boot = bootstrp(nboot, @(i) polyfit(X, Y(i), 1), 1:length(Y));
        slopes_aucs.(statsTypes{st})(n) = fit_obs(1);
        pvals_aucs.(statsTypes{st})(n)  = sum(abs(fit_boot(:,1)) >= abs(fit_obs(1)))/nboot;   % p-value for the slope
    end
end
disp('Finished: slopes fitted for animal EI vs DA trend');


%% Plot summary scatters for each statstype with top slope with
% significant p values

nTrialsUsed = 5:50;
minTrials = min(arrayfun(@(x) length(x.CueMax_slope), DAtrend)); 

initializeFig(1,1); tiledlayout(2,8);
for st = 1:length(statsTypes)
    cue = statsTypes{st};
    
    % Find top slope with significant p values
    significantIdx = pvals_aucs.(cue) <= 0.05;
    [maxSlope,idx] = max(slopes_aucs.(cue)(significantIdx));
    significantTrials = find(significantIdx); maxIdx = significantTrials(idx);

    % Get slopes and calculate fit
    topSlopes = getNthSlope(DAtrend,maxIdx);
    validIdx = ~isnan(topSlopes.CueMax);
    Y = topSlopes.(cue)(validIdx); 
    X = animalEIindex_aucs(validIdx);
    fit = polyfit(X, Y, 1);

    % Plot scatter plot and fitted line
    nexttile;
    scatter(X,Y,100,'filled'); hold on;
    plot(X,polyval(fit,X),LineWidth=5);
    xlabel('Animal EI charge index');
    ylabel(['Slope for ',cue]);
    title(['Last ',num2str(maxIdx),' trials']);


    % Get slopes and calculate fit
    ctrlIdx = 10;
    topSlopes = getNthSlope(DAtrend,ctrlIdx);
    validIdx = ~isnan(topSlopes.CueMax);
    Y = topSlopes.(cue)(validIdx); 
    X = animalEIindex_aucs(validIdx);
    fit = polyfit(X, Y, 1);

    % Plot scatter plot and fitted line
    nexttile;
    scatter(X,Y,100,'filled'); hold on;
    plot(X,polyval(fit,X),LineWidth=5);
    xlabel('Animal EI charge index');
    ylabel(['Slope for ',cue]);
    title(['Last ',num2str(ctrlIdx),' trials']);
end

%% (Analysis) Plot the relationship between DA trend and animal EI

% Positive slope here means that increase in DA amplitude in the last n trials
% corresponds to more cells being inhibitory
nTrialsUsed = 5:77;
initializeFig(1,1); tiledlayout(2,4);
for st = 1:length(statsTypes)
    cue = statsTypes{st};
    nexttile;
    
    % Left y-axis: plot slopes
    yyaxis left;
    plot(nTrialsUsed, slopes_aucs.(cue)(nTrialsUsed), '-o', 'LineWidth', 1.5);
    ylabel('Slope');
    xlabel('n trials used'); xlim([nTrialsUsed(1),nTrialsUsed(end)]);
    title(sprintf('%s: Slope and p-value vs n', cue));
    
    % Right y-axis: plot p-values
    yyaxis right;
    plot(nTrialsUsed, pvals_aucs.(cue)(nTrialsUsed), '-s', 'LineWidth', 1.5);
    xlim([nTrialsUsed(1),nTrialsUsed(end)]);
    ylabel('p-value');
    ylim([0 1]); hold on;
    plot(nTrialsUsed, repmat(0.05, size(nTrialsUsed)), '--k', 'LineWidth', 1);
    hold off; box off;
end

%% (Analysis) Calculate slopes of last n trials

shortWindow = 14:24;
longWindow  = 36:40;

% Calculate average slopes in shortWindow
avgCueMaxSlope_short  = arrayfun(@(d) mean(d.CueMax_slope(shortWindow)),  DAtrend)';
avgCueMinSlope_short  = arrayfun(@(d) mean(d.CueMin_slope(shortWindow)),  DAtrend)';
avgCueAvgSlope_short  = arrayfun(@(d) mean(d.CueAvg_slope(shortWindow)),  DAtrend)';
avgCueAmpSlope_short = arrayfun(@(d) mean(d.CueAmp_slope(shortWindow)), DAtrend)';

% Calculate average slopes in longWindow
avgCueMaxSlope_long  = arrayfun(@(d) mean(d.CueMax_slope(longWindow)),  DAtrend)';
avgCueMinSlope_long  = arrayfun(@(d) mean(d.CueMin_slope(longWindow)),  DAtrend)';
avgCueAvgSlope_long  = arrayfun(@(d) mean(d.CueAvg_slope(longWindow)),  DAtrend)';
avgCueAmpSlope_long = arrayfun(@(d) mean(d.CueAmp_slope(longWindow)), DAtrend)';

% % Calculate average p value in lastTrialWindow
% avgCueMaxPval  = arrayfun(@(d) mean(d.CueMax_pval(lastTrialWindow)),  DAtrend)';
% avgCueMinPval  = arrayfun(@(d) mean(d.CueMin_pval(lastTrialWindow)),  DAtrend)';
% avgCueAvgPval  = arrayfun(@(d) mean(d.CueAvg_pval(lastTrialWindow)),  DAtrend)';
% avgCueAreaPval = arrayfun(@(d) mean(d.CueAmp_pval(lastTrialWindow)), DAtrend)';

% Define the eight slope vectors in a cell array and label them
slopeVectors = {avgCueMaxSlope_short, avgCueMinSlope_short, avgCueAvgSlope_short, avgCueAmpSlope_short, ...
                avgCueMaxSlope_long,  avgCueMinSlope_long,  avgCueAvgSlope_long,  avgCueAmpSlope_long};
slopeNames = {'CueMax Slope (Short)', 'CueMin Slope (Short)', 'CueAvg Slope (Short)', 'CueAmp Slope (Short)', ...
              'CueMax Slope (Long)',  'CueMin Slope (Long)',  'CueAvg Slope (Long)',  'CueAmp Slope (Long)'};


% Set group names
numGroups = 2;
groupNames = {'Down DA','Stable DA','Up DA'};

% Set color
upColor = bluePurpleRed(1,:); 
downColor = bluePurpleRed(end,:); 
groupColors = {downColor,[.2 .2 .2],upColor};

%% (Analysis) TODO: Histograms of slope EI

%% (Analysis) kmeans clustering to define up/down/stable animals

% Preallocate cell arrays to save grouping results for each slope vector
downAnimals_kmeans = cell(numel(slopeVectors),1);
stableAnimals_kmeans = cell(numel(slopeVectors),1);
upAnimals_kmeans = cell(numel(slopeVectors),1);

figure; tiledlayout(2,4);

for i = 1:numel(slopeVectors)
    nexttile;
    data = slopeVectors{i};
    
    % Perform k-means clustering with replicates for stability
    [groupIdx, groupCenters] = kmeans(data, numGroups, 'Replicates', 100);
    
    % Sort group centers so that:
    %   Group 1: obvious negative (lowest center)
    %   Group 2: ambiguous (middle center)
    %   Group 3: obvious positive (highest center)
    [~, sortOrder] = sort(groupCenters);
    mapping = zeros(numGroups,1);
    mapping(sortOrder) = 1:numGroups;
    groupLabels = mapping(groupIdx);
    
    % Save grouping results
    downAnimals_kmeans{i} = find(groupLabels == 1);
    stableAnimals_kmeans{i} = find(groupLabels == 2);
    upAnimals_kmeans{i} = find(groupLabels == 3);
    
    hold on;
    for j = 1:numGroups
        idx = groupLabels == j;
        scatter(find(idx), data(idx), 100, groupColors{j}, 'filled', DisplayName=groupNames{j});
    end
    hold off;
    
    title(slopeNames{i});
    xlabel('Animal');
    ylabel('Slope');
    legend('Location','southeast');
end
sgtitle('K-means classification');

%% (Analysis) Quantile based method to define up/down/stable animals

% Preallocate cell arrays to save grouping results for each slope vector
downAnimals_quant = cell(numel(slopeVectors),1);
stableAnimals_quant = cell(numel(slopeVectors),1);
upAnimals_quant = cell(numel(slopeVectors),1);

figure; tiledlayout(2,4);

for i = 1:numel(slopeVectors)
    nexttile;
    data = slopeVectors{i};
    
    % Use quantiles to classify data into three groups
    edges = quantile(data, [0, 1/numGroups, 2/numGroups, 1]);
    groupLabels = discretize(data, edges);
    
    % Save grouping results
    downAnimals_quant{i} = find(groupLabels == 1);
    stableAnimals_quant{i} = find(groupLabels == 2);
    upAnimals_quant{i} = find(groupLabels == 3);
    
    hold on;
    for j = 1:numGroups
        idx = groupLabels == j;
        scatter(find(idx), data(idx), 100, groupColors{j}, 'filled', DisplayName=groupNames{j});
    end
    hold off;
    
    title(slopeNames{i});
    xlabel('Animal');
    ylabel('Slope');
    legend('Location','southeast');
end
sgtitle('Quantiles classification');


%% (Analysis) Plot cell EI

% Set up getCellIndices function
getCellIndices = @(indices) cell2mat(arrayfun(@(idx) find(strcmpi(combined_cells.Animal, animalList{idx})), indices, 'UniformOutput', false));

for i = 4:4%numel(slopeNames)
    % Get the animal indices from k-means
    animalIdxDown = downAnimals_kmeans{i}; 
    animalIdxStable = stableAnimals_kmeans{i};
    animalIdxUp = upAnimals_kmeans{i};

    % Convert animal indices to cell indices
    cellIdxDown   = getCellIndices(animalIdxDown);
    cellIdxStable = getCellIndices(animalIdxStable);
    cellIdxUp     = getCellIndices(animalIdxUp);
    groupIdx = {cellIdxDown,cellIdxStable,cellIdxUp};
    plotGroup = [1, 1, 0];
    
    % Create a figure name based on the current criterion (e.g., "CellEI-CueMaxSlopeShort")
    % Removing spaces and parentheses for a clean filename.
    figureName = ['CellEI-' regexprep(slopeNames{i}, '[ ()]', ''),'-kmeans'];
    close all;
    plotCellEI(combined_cells,groupIdx,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);


    % Get the animal indices from k-means
    animalIdxDown = downAnimals_quant{i}; 
    animalIdxStable = stableAnimals_quant{i};
    animalIdxUp = upAnimals_quant{i};

    % Convert animal indices to cell indices
    cellIdxDown   = getCellIndices(animalIdxDown);
    cellIdxStable = getCellIndices(animalIdxStable);
    cellIdxUp     = getCellIndices(animalIdxUp);
    groupIdx = {cellIdxDown,cellIdxStable,cellIdxUp};
    plotGroup = [1, 1, 1];
    
    % Create a figure name based on the current criterion (e.g., "CellEI-CueMaxSlopeShort")
    % Removing spaces and parentheses for a clean filename.
    figureName = ['CellEI-' regexprep(slopeNames{i}, '[ ()]', ''),'-quantile'];
    close all;
    plotCellEI(combined_cells,groupIdx,...
           plotGroup=plotGroup,groupColors=groupColors,groupNames=groupNames,...
           save=false,figureName=figureName,resultPath=resultPath,print=true);
end

%% (Analysis) Check relationship between DA start & end diff and patch data