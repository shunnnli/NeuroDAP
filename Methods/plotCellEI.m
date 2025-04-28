function plotCellEI(combined_cells, groupIdx,options)

arguments
    combined_cells table
    groupIdx cell

    options.nboot {mustBeNumeric} = 5000

    options.plotGroup double
    options.groupColors cell
    options.groupNames cell

    options.colorScale double
    options.colormap

    options.statsType string = 'mean'
    options.scatterScale string = 'original'

    options.save logical = true
    options.resultPath string
    options.figureName string
    options.print logical = true

    options.centerEI logical = false
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

% Set colors
if ~isfield(options,'groupColors')
    rewardColor = [0 158 115]./255; rewardCtrlColor = 1 - opacity*(1-rewardColor);
    punishColor = [135 104 247]./255; punishCtrlColor = 1 - opacity*(1-punishColor);
    groupColors = {[.7 .7 .7],rewardColor,punishColor,...
                    rewardCtrlColor,punishCtrlColor};
else
    if length(options.groupColors) ~= length(groupIdx)
        error('groupColors should have the same length as groupIdx!');
    end
    groupColors = options.groupColors;
end
if isfield(options,'colorScale')
    if ~isfield(options,'colormap')
        options.colormap = getColormap([135 104 247],[0 158 115],500,'midCol',[255 255 255]);
    end
    [~,colors] = mapValueToColormap(options.colorScale,options.colormap);
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

% Calculate EI statistics
% 1. EPSC + IPSC
EIsum_peaks = IPSC_peaks + EPSC_peaks;
EIsum_aucs = IPSC_aucs + EPSC_aucs;

% 3. EI index (-1: all excitatory & 1: all inhibitory)
EIindex_peaks = (abs(IPSC_peaks)-abs(EPSC_peaks)) ./ (abs(IPSC_peaks)+abs(EPSC_peaks));
EIindex_aucs = (abs(IPSC_aucs)-abs(EPSC_aucs)) ./ (abs(IPSC_aucs)+abs(EPSC_aucs));

% 4. z-scored EI index
bootstrapping = @(B,nBoot) reshape(B(randi(numel(B), numel(B), nBoot)), [], 1);
centering = @(X,B) X/mean(B); %@(X,B) (X - mean(B)) ./ std(B);
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

if options.centerEI
    EPSC_peaks = EPSC_peaks_centered;
    IPSC_peaks = IPSC_peaks_centered;
    EPSC_aucs = EPSC_aucs_centered;
    IPSC_aucs = IPSC_aucs_centered;

    EIsum_peaks = IPSC_peaks_centered + EPSC_peaks_centered;
    EIsum_aucs = IPSC_aucs_centered + EPSC_aucs_centered;
    EIindex_peaks = EIindex_peaks_centered;
    EIindex_aucs = EIindex_aucs_centered;
end

%% Bootstrap data for p value calculation
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
        if ~isfield(options,'colorScale')
            scatter(abs(EPSC_aucs(idx)),abs(IPSC_aucs(idx)),dotSize,groupColor,"filled",MarkerFaceAlpha=0.8,MarkerEdgeAlpha=0.8); hold on; 
        else
            scatter(abs(EPSC_aucs(idx)),abs(IPSC_aucs(idx)),dotSize,colors(idx),"filled"); hold on; 
            colormap(options.colormap);
        end
    end
end
diagonal = refline(1,0); diagonal.Color=[.9 .9 .9]; diagonal.LineWidth=4; diagonal.LineStyle='--';
if strcmpi(options.scatterScale,'log')
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
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
        if ~isfield(options,'colorScale')
            scatter(abs(EPSC_peaks(idx)),abs(IPSC_peaks(idx)),dotSize,groupColor,"filled",MarkerFaceAlpha=0.8,MarkerEdgeAlpha=0.8); hold on; 
        else
            scatter(abs(EPSC_peaks(idx)),abs(IPSC_peaks(idx)),dotSize,colors(idx),"filled"); hold on; 
            colormap(options.colormap);
        end
    end
end
diagonal = refline(1,0); diagonal.Color=[.9 .9 .9]; diagonal.LineWidth=4; diagonal.LineStyle='--';
if strcmpi(options.scatterScale,'log')
    set(gca, 'XScale', 'log', 'YScale', 'log');
end
set(gca, 'Children', flipud(get(gca,'Children'))); % Flip drawing order
xlabel('EPSC amplitude (pA)');
ylabel('IPSC amplitude (pA)');

%% 
% cellIdx = 84;
% scatter(abs(EPSC_peaks(cellIdx)),abs(IPSC_peaks(cellIdx)),dotSize,'b',"filled")
%%

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

        if contains(options.statsType,'diff')
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
        elseif contains(options.statsType,'mean')
            meanEIauc_observed = mean(EIindex_aucs(idx));
            meanEIauc_bs = mean(bs_EIindex_aucs(fakeIdx,:));
            p_cdf = min([sum(meanEIauc_bs <= meanEIauc_observed)/options.nboot,...
                        sum(meanEIauc_bs >= meanEIauc_observed)/options.nboot]);
            % Plot histogram and p value
            histogram(meanEIauc_bs,200,EdgeColor=bootColor,FaceColor=bootColor); hold on; 
            xline(meanEIauc_observed,'-',['p=',num2str(p_cdf)],Color=groupColor,LineWidth=3);
            box off; c = gca; c.YAxis(1).Visible = 'off';
        end
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

        if contains(options.statsType,'diff')
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
        elseif contains(options.statsType,'mean')
            meanEIpeaks_observed = mean(EIindex_peaks(idx));
            meanEIpeaks_bs = mean(bs_EIindex_peaks(fakeIdx,:));
            p_cdf = min([sum(meanEIpeaks_bs <= meanEIpeaks_observed)/options.nboot,...
                        sum(meanEIpeaks_bs >= meanEIpeaks_observed)/options.nboot]);
            % Plot histogram and p value
            histogram(meanEIpeaks_bs,200,EdgeColor=bootColor,FaceColor=bootColor); hold on; 
            xline(meanEIpeaks_observed,'-',['p=',num2str(p_cdf)],Color=groupColor,LineWidth=3);
            box off; c = gca; c.YAxis(1).Visible = 'off';
        end

        
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
    disp('---------- plotCellEI ----------');
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