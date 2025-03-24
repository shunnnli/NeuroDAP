function plotDAvsEImap(DAvsEImap,options)

arguments
    DAvsEImap struct

    options.dataType string = 'raw'
    options.statType string = 'amp'
    options.nTrials double

    options.layout
    options.tileSize double = [1,1]

    options.colorlim double
    options.colorBarLabel string = 'Slope'
    options.colormap
    options.colorBarPosition double = [0, 0, 0.025, 0.5]

    % p value map options
    options.pval logical = true
    options.pvalColor double = [.3 .3 .3]
    options.pvalLineWidth double = 3

    options.title
end

%% Check input

% Define dataType
if contains(options.dataType,"raw")
    dataType = 'raw';
elseif contains(options.dataType,"smooth")
    dataType = 'smoothed';
end

% Define statType
if contains(options.statType,"max")
    statType = 'max';
elseif contains(options.statType,"min")
    statType = 'min';
elseif sum(contains(options.statType,["avg","mean","average"]))
    statType = 'avg';
elseif contains(options.statType,"amp")
    statType = 'amp';
end

% Define colormap
if ~isfield(options,'colormap')
    [~,~,~,~,blueWhiteRed,~,~] = loadColors;
    options.colormap = blueWhiteRed;
end

%% Load in map data

% nTrials
if ~isfield(options,'nTrials')
    if ~isfield(DAvsEImap,'options')
        options.nTrials = 30;
    else
        options.nTrials = DAvsEImap.options.nTrials;
    end
end
    

% See whether whole session correlations needs to be plotted (default is false)
skipWholeSession = DAvsEImap.options.skipWholeSession;
if ~isfield(options,'layout')
    error('Need to provide tiledlayout instance to the method!');
end

if skipWholeSession
    % Get map range
    earlyIdx = 1:options.nTrials;
    lateIdx = options.nTrials+1:options.nTrials*2;
    % Get tick marks
    earlyTicks = 1:options.nTrials;
    lateTicks = -options.nTrials:-1;
    trialTickLoc = floor(linspace(1,options.nTrials,3));
else
    % Get tick marks
    trialIdx = [1:options.nTrials,-options.nTrials:-1];
    trialTickLoc = floor(linspace(1,options.nTrials,3));
    options.colorBarPosition = [0.02, 0.01, 0.025, 0.5];
end

%% Extract heatmap

map = DAvsEImap.(statType).(dataType);
pvalMap = DAvsEImap.(statType).(strcat('pval_',dataType));

if skipWholeSession
    earlyMap = map(earlyIdx,earlyIdx);
    lateMap = map(lateIdx,lateIdx);
    earlyPvalMap = pvalMap(earlyIdx,earlyIdx);
    latePvalMap = pvalMap(lateIdx,lateIdx);
end

% Define colorlim
if ~isfield(options,'colorlim')
    if skipWholeSession
        climit = max([abs(earlyMap),abs(lateMap)]);
    else
        climit = max(abs(map),[],'all');
    end
    options.colorlim = [-climit,climit];
end

%% Plot heatmap
if skipWholeSession
    heatmapLayout = options.layout;
    heatmapLayout.Layout.TileSpan = options.tileSize;
    heatmapLayout.TileSpacing = 'tight'; heatmapLayout.Padding = 'tight'; axis off;

    % Plot early session triange
    nexttile(heatmapLayout,1);
    plotHeatmap(earlyMap,[],dataType='heatmap',...
                colormap=options.colormap,...
                colorlim=options.colorlim,...
                colorBarLabel=options.colorBarLabel);
    % Plot pval
    if options.pval
        significantMap = earlyPvalMap <= 0.05;
        if sum(significantMap,"all") > 0
            contour(significantMap, [0.5 0.5],...
                    Color=options.pvalColor,LineWidth=options.pvalLineWidth);
        end
    end
    % Plot axis label
    set(gca, 'XAxisLocation', 'top','YAxisLocation', 'right');
    xticks(trialTickLoc); yticks(trialTickLoc);
    xticklabels(earlyTicks(trialTickLoc)); yticklabels(earlyTicks(trialTickLoc));
    xlabel('Ending trial'); ylabel('Starting trial');


    % Plot late session triange
    nexttile(heatmapLayout,2);
    plotHeatmap(lateMap,[],dataType='heatmap',...
                colormap=options.colormap,...
                colorlim=options.colorlim,...
                plotColorBar=false);
    % Plot pval
    if options.pval
        significantMap = latePvalMap <= 0.05;
        if sum(significantMap,"all") > 0
            contour(significantMap, [0.5 0.5],...
                    Color=options.pvalColor,LineWidth=options.pvalLineWidth);
        end
    end
    % Plot axis label
    set(gca, 'XAxisLocation', 'top','YAxisLocation', 'right');
    xticks(trialTickLoc); yticks(trialTickLoc);
    xticklabels(lateTicks(trialTickLoc)); yticklabels(lateTicks(trialTickLoc));
    xlabel('Ending trial'); ylabel('Starting trial');
    
    if isfield(options,'title')
        sgtitle(heatmapLayout,options.title);
    end

else
    ax = nexttile(options.tileSize);
    plotHeatmap(map,[],dataType='heatmap',...
                colormap=options.colormap,...
                colorlim=options.colorlim,...
                colorBarLabel=options.colorBarLabel,...
                tile=ax,...
                colorBarPosition=options.colorBarPosition);

    % Plot pval
    if options.pval
        pvalMap = DAvsEImap.(statType).(strcat('pval_',dataType));
        significantMap = pvalMap <= 0.05;
        if sum(significantMap,"all") > 0
            contour(significantMap, [0.5 0.5],...
                    Color=options.pvalColor,LineWidth=options.pvalLineWidth);
        end
    end

    % Plot axis label
    set(gca, 'XAxisLocation', 'top','YAxisLocation', 'right');
    xticks(trialTickLoc); yticks(trialTickLoc);
    xticklabels(trialIdx(trialTickLoc)); yticklabels(trialIdx(trialTickLoc));
    xlabel('Ending trial'); ylabel('Starting trial');

    if isfield(options,'title')
        sgtitle(options.title);
    end
end

end