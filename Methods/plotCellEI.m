function plotCellEI(combined_cells, groupIdx,options)

arguments
    combined_cells table
    groupIdx cell

    options.nboot {mustBeNumeric} = 5000

    options.plotGroup double
    options.groupColors cell
    options.groupNames cell

    options.save logical = true
    options.resultPath string
    options.figureName string
    options.print logical = true
end

%% Set up

% Select data to plot
if ~isfield(options,'plotGroup')
    plotGroup = ones(length(groupIdx));
else
    plotGroup = options.plotGroup;
end

if sum(plotGroup) == 1; plotBootstrap = true; 
else; plotBootstrap = false; end

% Set plotting params
dotSize = 200; opacity = 0.3;
[~,~,~,~,~,~,bluePurpleRed] = loadColors;

% Set colors
if ~isfield(options,'groupColors')
    rewardColor = bluePurpleRed(1,:); rewardCtrlColor = 1 - opacity*(1-rewardColor);
    punishColor = bluePurpleRed(end,:); punishCtrlColor = 1 - opacity*(1-punishColor);
    groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                    rewardCtrlColor,punishCtrlColor};
else
    if length(options.groupColors) ~= length(groupIdx)
        error('groupColors should have the same length as groupIdx!');
    end
    groupColors = options.groupColors;
end

% Set groupNames
if options.print
    if ~isfield(options,'groupNames')
        warning('Option print is true but groupNames are not provided!');
        options.print = false;
    elseif length(options.groupNames) ~= length(groupIdx)
        warning('groupNames provided should have the same length as groupIdx.');
        options.print = false;
    end
end

%% Extract EI information

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

%% Bootstrap data

% Pool all cells together
pooledIdx = vertcat(groupIdx{:}); 
pooledLengths = cellfun(@length, groupIdx);
bootIdx = mat2cell(1:length(pooledIdx),1,pooledLengths);

% Generate shuffled indices
pooledIdx_shuffled = arrayfun(@(x) pooledIdx(randperm(length(pooledIdx))), 1:options.nboot, UniformOutput=false);
pooledIdx_shuffled = cell2mat(pooledIdx_shuffled);

% Find statistics for pooled values
bs_EIsum_aucs = arrayfun(@(x) EIsum_aucs(x), pooledIdx_shuffled);
bs_EIsum_peaks = arrayfun(@(x) EIsum_peaks(x), pooledIdx_shuffled);
bs_EIindex_aucs = arrayfun(@(x) EIindex_aucs(x), pooledIdx_shuffled);
bs_EIindex_peaks = arrayfun(@(x) EIindex_peaks(x), pooledIdx_shuffled);

%% Plot figure

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
        p_cdf = min([sum(diffArea_bs <= diffArea_observed)/options.nboot,...
                    sum(diffArea_bs >= diffArea_observed)/options.nboot]);

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
        p_cdf = min([sum(diffArea_bs <= diffArea_observed)/options.nboot,...
                    sum(diffArea_bs >= diffArea_observed)/options.nboot]);

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

if options.save
    saveFigures(fig,options.figureName,options.resultPath,...
                saveFIG=true,savePDF=true,savePNG=true);
end

if options.print
    % Loop through each group flagged for plotting.
    for i = 1:length(groupIdx)
        if plotGroup(i)
            idx = groupIdx{i};
            totalCells = numel(idx);  % Total number of cells (rows)
            totalAnimals = numel(unique(combined_cells.Animal(idx)));  % Unique animals count
            fprintf('%s: n = %d, N = %d\n', options.groupNames{i}, totalCells, totalAnimals);
        end
    end
end

end