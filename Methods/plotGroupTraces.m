function legendList = plotGroupTraces(traces,timestamp,colormap,options)

%% Notes
% plotGroupTraces.m, Shun Li
% 2023/11/30

%% Function
arguments
    traces double
    timestamp double
    colormap

    options.nGroups double
    options.groupSize double

    options.animalStartIdx double

    % Plot options
    options.plot logical = true % whether or not to plot traces
    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 2
end

% Check input
if ~isfield(options,'nGroups') && ~isfield(options,'groupSize')
    error('plotGroupTraces: need to define either nGroups or groupSize');
elseif isfield(options,'nGroups') && ~isfield(options,'groupSize')
    options.groupSize = ceil(size(traces,1)/options.nGroups);
elseif ~isfield(options,'nGroups') && isfield(options,'groupSize')
    options.nGroups = ceil(size(traces,1)/options.groupSize);
elseif isfield(options,'nGroups') && isfield(options,'groupSize')
    maxNumGroups = ceil(size(traces,1)/options.groupSize);
    if options.nGroups > maxNumGroups; options.nGroups = maxNumGroups; end
end

% Reorganize traces if animalStartIdx is provided (interleave animals)
if isfield(options,'animalStartIdx')
    traces_animals = cell(1,length(options.animalStartIdx));
    for i = 1:length(options.animalStartIdx)
        startIdx = options.animalStartIdx(i);
        if i == length(options.animalStartIdx); endIdx = size(traces,1);
        else; endIdx = options.animalStartIdx(i+1)-1; end
        traces_animals{i} = traces(startIdx:endIdx,:);
    end
    traces_animals_dim = cellfun('size',traces_animals,1);
end

%% Plot data
legendList = cell(options.nGroups,1);
nColors = round(linspace(1,size(colormap,1),options.nGroups));

for i = 1:options.nGroups
    startTrial = (i-1)*options.groupSize+1; 
    if i == options.nGroups; endTrial = min(size(traces,1),options.groupSize*options.nGroups);
    else; endTrial = i*options.groupSize; end

    if isfield(options,'animalStartIdx')
        plotData = cell2mat(cellfun(@(x,dim) x(startTrial:min(endTrial,dim),:),traces_animals,num2cell(traces_animals_dim),'UniformOutput',false)');
        plotSEM(timestamp,plotData,colormap(nColors(i),:),LineStyle=options.LineStyle,LineWidth=options.LineWidth);
        legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial),' (n=',num2str(size(plotData,1)),')'];
    else
        plotData = traces(startTrial:endTrial,:);
        plotSEM(timestamp,plotData,colormap(nColors(i),:),LineStyle=options.LineStyle,LineWidth=options.LineWidth);
        legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial),' (n=',num2str(size(plotData,1)),')'];
    end
    
    
end

end