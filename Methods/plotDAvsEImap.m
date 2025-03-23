function plotDAvsEImap(DAvsEImap,options)

arguments
    DAvsEImap struct

    options.dataType string = 'raw'
    options.statType string = 'amp'
    options.nTrials double = 30

    options.tile
    options.colorlim double
    options.colorBarLabel string = 'Slope'
    options.colormap
    options.colorBarPosition double = [0.05, 0.02, 0.025, 0.5]

    % p value map options
    options.pval logical = true
    options.pvalColor double = [.3 .3 .3]
    options.pvalLineWidth double = 3
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

% Get tick marks
trialIdx = [1:options.nTrials,-options.nTrials:-1];
trialTickLoc = round(linspace(1,options.nTrials*2,5)); 

%% Plot heatmap

map = DAvsEImap.(statType).(dataType);

% Define colorlim
if ~isfield(options,'colorlim')
    % Get the max and min of colormap
    climit = max([abs(max(map,[],'all')),abs(min(map,[],'all'))]);
    options.colorlim = [-climit,climit];
end

plotHeatmap(map,[],dataType='heatmap',...
            colormap=options.colormap,...
            colorlim=options.colorlim,...
            colorBarLabel=options.colorBarLabel,...
            tile=options.tile,...
            colorBarPosition=options.colorBarPosition);

%% Plot pval
if options.pval
    pvalMap = DAvsEImap.(statType).(strcat('pval_',dataType));
    significantMap = pvalMap <= 0.05;
    if sum(significantMap,"all") > 0
        contour(significantMap, [0.5 0.5],...
                Color=options.pvalColor,LineWidth=options.pvalLineWidth);
    end
end

%% Plot axis label
set(gca, 'XAxisLocation', 'top','YAxisLocation', 'right');
xticks(trialTickLoc); yticks(trialTickLoc);
xticklabels(trialIdx(trialTickLoc)); yticklabels(trialIdx(trialTickLoc));
xlabel('Trials'); ylabel('Trials');

end