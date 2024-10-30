%% Shun_analyzeCellEI

% 2024/09/25

%% Define data path

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions for analysis
% parentPath = osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/');
% expPath = uipickfiles('FilterSpec',parentPath,'Prompt','Select experiment folders');

resultPath = '/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/Combined';
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
today = char(datetime('today','Format','yyyyMMdd')); 

%% Drag the desired combined.mat file

% load('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Patch/Combined/combined_cells_20241010.mat');

%% (Optional) Get cell table

combined_cells = getCellTable(combined_epochs,save=true,...
                         timeRange=[-10,50]);

combined_cells.Learned = ones(size(combined_cells,1),1);
notLearnedAnimals = {'SL044','SL232','SL212','SL254','SL271'}; 
notLearnedIdx = find(ismember(combined_cells.Animal,notLearnedAnimals));
combined_cells.Analyze = logical(combined_cells.Learned);
combined_cells.Analyze(notLearnedIdx) = 0;

disp('Saving combined_cells and combined_epochs...');
save(fullfile(resultPath,strcat('combined_cells_',today,'.mat')),"combined_cells","-v7.3");
save(fullfile(resultPath,strcat('combined_epochs_',today,'.mat')),"combined_epochs","-v7.3");
nCells = size(combined_cells,1);
disp(strcat("Finished: saving ",num2str(nCells),' cells to combined_cells'));

%% Extract EI information from combined_epochs

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


% Remove from analysis if both currents are less than 5pA
% combined_cells.Learned = ones(size(combined_cells,1),1);
notLearnedAnimals = {'SL044','SL232','SL212','SL254','SL271'}; 
notLearnedIdx = find(ismember(combined_cells.Animal,notLearnedAnimals));
combined_cells.Analyze = logical(combined_cells.Learned);
combined_cells.Analyze(notLearnedIdx) = 0;
removeIdx = find(abs(EPSC_peaks) <= 0 & abs(IPSC_peaks) <= 0);
combined_cells.Analyze(removeIdx) = 0;

% Set task range
randomIdx = find(strcmpi('Random',combined_cells.Task) & combined_cells.Analyze);
rewardIdx = find(strcmpi('Reward pairing',combined_cells.Task) & combined_cells.Analyze);
punishIdx = find(strcmpi('Punish pairing',combined_cells.Task) & combined_cells.Analyze);
rewardCtrlIdx = [find(strcmpi('Reward control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Reward pairing',combined_cells.Task) & ~combined_cells.Analyze)];
punishCtrlIdx = [find(strcmpi('Punish control',combined_cells.Task) & combined_cells.Analyze);...
                 find(strcmpi('Punish pairing',combined_cells.Task) & ~combined_cells.Analyze)];


% Pool data and bootstrap
% Pool all cells together
pooledIdx = [randomIdx;rewardIdx;punishIdx;rewardCtrlIdx;punishCtrlIdx];
bootIdx = mat2cell(1:length(pooledIdx),1,[length(randomIdx),length(rewardIdx),length(punishIdx),length(rewardCtrlIdx),length(punishCtrlIdx)]);

% Generate shuffled indices
nboot = 5000;
pooledIdx_shuffled = arrayfun(@(x) pooledIdx(randperm(length(pooledIdx))), 1:nboot, UniformOutput=false);
pooledIdx_shuffled = cell2mat(pooledIdx_shuffled);

% Find statistics for pooled values
bs_EIsum_aucs = arrayfun(@(x) EIsum_aucs(x), pooledIdx_shuffled);
bs_EIsum_peaks = arrayfun(@(x) EIsum_peaks(x), pooledIdx_shuffled);
bs_EIindex_aucs = arrayfun(@(x) EIindex_aucs(x), pooledIdx_shuffled);
bs_EIindex_peaks = arrayfun(@(x) EIindex_peaks(x), pooledIdx_shuffled);

%% Plot multiple conditions

close all;

% Select data to plot
groupIdx = {randomIdx,rewardIdx,punishIdx,rewardCtrlIdx,punishCtrlIdx}; 
plotGroup = [1,0,0,0,1];

if sum(plotGroup) == 1; plotBootstrap = true; 
else; plotBootstrap = false; end
% if sum(plotGroup) > 2; error('Only plot two groups at max!'); end

% Set plotting params
dotSize = 200; opacity = 0.3;

% Set colors
rewardColor = bluePurpleRed(1,:); rewardCtrlColor = 1 - opacity*(1-rewardColor);
punishColor = bluePurpleRed(end,:); punishCtrlColor = 1 - opacity*(1-punishColor);
groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                rewardCtrlColor,punishCtrlColor};

fig = initializeFig(1,1); master = tiledlayout(2,4);
master.TileSpacing = 'compact'; master.Padding = 'compact';

% Plot AUC EPSC vs IPSC
nexttile(master,1);
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

% Plot peak EPSC vs IPSC
nexttile(master,5);
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

% Plot charge EPSC+IPSC index
nexttile(master,2);
plotDistribution(EIsum_aucs,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                bootstrap=plotBootstrap,bootData=bs_EIsum_aucs,bootIdx=bootIdx,...
                masterlayout=master,tile=2,xlabel='EPSC+IPSC (pC)');

% Plot amplitude EPSC+IPSC index
nexttile(master,6);
plotDistribution(EIsum_peaks,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                bootstrap=plotBootstrap,bootData=bs_EIsum_peaks,bootIdx=bootIdx,...
                masterlayout=master,tile=6,xlabel='EPSC+IPSC (pA)');

% Plot EI charge index
nexttile(master,3);
plotDistribution(EIindex_aucs,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                bootstrap=plotBootstrap,bootData=bs_EIindex_aucs,bootIdx=bootIdx,...
                masterlayout=master,tile=3,xlabel='EI charge index',limit=[-1,1]);

% Plot EI amplitude index
nexttile(master,7);
plotDistribution(EIindex_peaks,groupIdx=groupIdx,plotGroup=plotGroup,color=groupColors,...
                bootstrap=plotBootstrap,bootData=bs_EIindex_peaks,bootIdx=bootIdx,...
                masterlayout=master,tile=7,xlabel='EI amplitude index',limit=[-1,1]);

% Plot EI charge index
nexttile(master,4); 
children = tiledlayout(master,4,1);
children.Layout.Tile = 4; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;
nexttile(children,1,[1 1]);
for group = 1:length(plotGroup)
    curGroup = length(plotGroup) - group + 1;
    if plotGroup(curGroup)
        idx = groupIdx{curGroup};
        fakeIdx = bootIdx{curGroup};
        groupColor = groupColors{curGroup};
        bootColor = 1 - opacity*(1-groupColor);

        % Calculate area diff: observed - bootstrap_sum
        diffArea_observed = getCDF_diffArea(EIindex_aucs(idx),bs_EIindex_aucs(fakeIdx,:),separate=false);
        % Calculate area diff: bootstrap_sum - individual bootstrap stim
        diffArea_bs = getCDF_diffArea(bs_EIindex_aucs(fakeIdx,:),bs_EIindex_aucs(fakeIdx,:));
        % Calculate p value
        p_cdf = min([sum(diffArea_bs <= diffArea_observed)/nboot,...
                    sum(diffArea_bs >= diffArea_observed)/nboot]);

        % Plot histogram and p value
        histogram(diffArea_bs,200,EdgeColor=bootColor,FaceColor=bootColor); hold on; 
        xline(diffArea_observed,'-',['p=',num2str(p_cdf)],Color=groupColor,LineWidth=3);
        box off; c = gca; c.YAxis(1).Visible = 'off';
    end
end
nexttile(children,2,[3,1]);
for group = 1:length(plotGroup)
    curGroup = length(plotGroup) - group + 1;
    if plotGroup(curGroup)
        idx = groupIdx{curGroup};
        groupColor = groupColors{curGroup};
        cdf = cdfplot(EIindex_aucs(idx)); hold on; 
        cdf.Color = groupColor; cdf.LineWidth = 3;

        if plotBootstrap
            idx = bootIdx{curGroup};
            bootColor = 1 - opacity*(1-groupColor);
            bs_sum = reshape(bs_EIindex_aucs(idx,:), [], 1);
            cdf = cdfplot(bs_sum); hold on; 
            cdf.Color = bootColor; cdf.LineWidth = 3;
        end
    end
end
xlabel('EI charge index'); ylabel('CDF'); box off; grid off;

% Plot EI amplitude index
nexttile(master,8);
children = tiledlayout(master,4,1);
children.Layout.Tile = 8; children.Layout.TileSpan = [1 1];
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;
nexttile(children,1,[1 1]);
for group = 1:length(plotGroup)
    curGroup = length(plotGroup) - group + 1;
    if plotGroup(curGroup)
        idx = groupIdx{curGroup};
        fakeIdx = bootIdx{curGroup};
        groupColor = groupColors{curGroup};
        bootColor = 1 - opacity*(1-groupColor);

        % Calculate area diff: observed - bootstrap_sum
        diffArea_observed = getCDF_diffArea(EIindex_peaks(idx),bs_EIindex_peaks(fakeIdx,:),separate=false);
        % Calculate area diff: bootstrap_sum - individual bootstrap stim
        diffArea_bs = getCDF_diffArea(bs_EIindex_peaks(fakeIdx,:),bs_EIindex_peaks(fakeIdx,:));
        % Calculate p value
        p_cdf = min([sum(diffArea_bs <= diffArea_observed)/nboot,...
                    sum(diffArea_bs >= diffArea_observed)/nboot]);

        % Plot histogram and p value
        histogram(diffArea_bs,200,EdgeColor=bootColor,FaceColor=bootColor); hold on; 
        xline(diffArea_observed,'-',['p=',num2str(p_cdf)],Color=groupColor,LineWidth=3);
        box off; c = gca; c.YAxis(1).Visible = 'off';
    end
end
nexttile(children,2,[3,1]);
for group = 1:length(plotGroup)
    curGroup = length(plotGroup) - group + 1;
    if plotGroup(curGroup)
        idx = groupIdx{curGroup};
        groupColor = groupColors{curGroup};
        cdf = cdfplot(EIindex_peaks(idx)); hold on; 
        cdf.Color = groupColor; cdf.LineWidth = 3;

        if plotBootstrap
            idx = bootIdx{curGroup};
            bootColor = 1 - opacity*(1-groupColor);
            bs_sum = reshape(bs_EIindex_peaks(idx,:), [], 1);
            cdf = cdfplot(bs_sum); hold on; 
            cdf.Color = bootColor; cdf.LineWidth = 3;
        end
    end
end
xlabel('EI amplitude index'); ylabel('CDF'); box off; grid off;

% saveFigures(fig,'Summary-reward-punish',resultPath,...
%             saveFIG=true,savePDF=true,savePNG=true);

%% Extract response trace (stupid for-loop way)

timeRange = [-10,50]; 
time2sample = 10000/1000; 
timeRangeSamples = timeRange * time2sample;
timeWindowLength = timeRangeSamples(2) - timeRangeSamples(1) + 1;
timeRangeInms = linspace(-10,50,timeWindowLength);
nCells = size(combined_cells,1);

EPSC_traces = []; IPSC_traces = [];

for c = 1:nCells
    curCell = combined_cells(c,:); 

    % Get rows
    EPSCrows = curCell.Vhold{1} == -70;
    IPSCrows = curCell.Vhold{1} >= 0;

    % Get included
    if sum(EPSCrows)>0; EPSCincluded = curCell.Included{1}{EPSCrows};
    else; EPSCincluded = []; end
    if sum(IPSCrows)>0; IPSCincluded = curCell.Included{1}{IPSCrows};
    else; IPSCincluded = []; end
    if sum(EPSCincluded) == 0; EPSCincluded = ones(length(EPSCincluded),1); end
    if sum(IPSCincluded) == 0; IPSCincluded = ones(length(IPSCincluded),1); end
    EPSCincluded = logical(EPSCincluded);
    IPSCincluded = logical(IPSCincluded);

    % Concatenate traces
    if sum(EPSCrows) ~= 0
        EPSCIdx = find(EPSCrows); curTrace = [];
        for epsc = EPSCIdx'
            stimOnset = curCell.Protocol{1}{epsc}.stimOnset(1);
            timeRangeWindow = stimOnset+timeRangeSamples(1) : stimOnset+timeRangeSamples(2);
            curTrace = [curTrace;curCell.('Processed sweeps'){1}{epsc}(:,timeRangeWindow)];
            curTrace(~EPSCincluded,:) = [];
        end
        EPSC_traces = [EPSC_traces;mean(curTrace,1)];
    else
        EPSC_traces = [EPSC_traces;nan(1,size(EPSC_traces,2))];
    end

    if sum(IPSCrows) ~= 0
        IPSCIdx = find(IPSCrows); curTrace = [];
        for ipsc = IPSCIdx'
            stimOnset = curCell.Protocol{1}{ipsc}.stimOnset(1);
            timeRangeWindow = stimOnset+timeRangeSamples(1) : stimOnset+timeRangeSamples(2);
            curTrace = [curTrace;curCell.('Processed sweeps'){1}{ipsc}(:,timeRangeWindow)];
            curTrace(~IPSCincluded,:) = [];
        end
        IPSC_traces = [IPSC_traces;mean(curTrace,1)];
    else
        IPSC_traces = [IPSC_traces;nan(1,size(IPSC_traces,2))];
    end
end


%% Plot response traces

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