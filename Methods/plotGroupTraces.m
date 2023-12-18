function legendList = plotGroupTraces(traces,timestamp,colormap,options)

%% Notes
% plotGroupTraces.m, Shun Li
% 2023/11/30

%% Function
arguments
    traces double
    timestamp double
    colormap

    options.groupby string = 'trials'
    options.nGroups double
    options.groupSize double

    options.startIdx struct % from combineTraces()
    options.animalStartIdx double
    options.sessionStartIdx double

    % Plot options
    options.plot logical = true % whether or not to plot traces
    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 2

    options.remaining string = 'include' % can also be 'exclude' or 'separate'
        % include: include the remaining traces to the last group
        % exclude: do not plot the remaining traces
        % separate: plot the remaining traces separately
end

%% Parse inputs

if isfield(options,'startIdx') 
    options.animalStartIdx = options.startIdx.animal;
    options.sessionStartIdx = options.startIdx.session;

    if iscell(options.animalStartIdx); options.animalStartIdx = options.animalStartIdx{1}; end
    if iscell(options.sessionStartIdx); options.sessionStartIdx = options.sessionStartIdx{1}; end
end

% Reorganize traces if animalStartIdx is provided (interleave animals)
if isfield(options,'animalStartIdx') && length(options.animalStartIdx)>1 && sum(strcmpi(options.groupby,["trial","trials"]))
    traces_animals = cell(1,length(options.animalStartIdx));
    for i = 1:length(options.animalStartIdx)
        startIdx = options.animalStartIdx(i);
        if i == length(options.animalStartIdx); endIdx = size(traces,1);
        else; endIdx = options.animalStartIdx(i+1)-1; end
        traces_animals{i} = traces(startIdx:endIdx,:);
    end
    traces_animals_dim = cellfun('size',traces_animals,1);
    nTraces = min(traces_animals_dim);
else
    nTraces = size(traces,1);
end

%% Calculate start and end trials

% Check input and update groupSize and nGroups
if sum(strcmpi(options.groupby,["trial","trials"]))
    if ~isfield(options,'nGroups') && ~isfield(options,'groupSize')
        error('plotGroupTraces: need to define either nGroups or groupSize');
    elseif isfield(options,'nGroups') && ~isfield(options,'groupSize')
        options.groupSize = ceil(nTraces/options.nGroups);
    elseif ~isfield(options,'nGroups') && isfield(options,'groupSize')
        options.nGroups = ceil(nTraces/options.groupSize);
    elseif isfield(options,'nGroups') && isfield(options,'groupSize')
        maxNumGroups = ceil(nTraces/options.groupSize);
        if options.nGroups > maxNumGroups; options.nGroups = maxNumGroups; end
    end

    % Calculate start and end trials
    groupIdx = 1:options.nGroups;
    startTrials = (groupIdx-1)*options.groupSize+1;
    endTrials = groupIdx * options.groupSize;
    if endTrials(end) > nTraces; endTrials(end) = nTraces; end

    % Update based on remaining number of traces that can not be grouped
    if strcmpi(options.remaining,'include') && nTraces > options.groupSize
        nTraces_lastGroup = endTrials(end) - startTrials(end);
        if nTraces_lastGroup > 0 && nTraces_lastGroup < options.groupSize/2
            endTrials(end-1) = nTraces;
            endTrials(end) = [];
            options.nGroups = options.nGroups - 1;
        end
    elseif strcmpi(options.remaining,'exclude') && nTraces > options.groupSize
        nTraces_lastGroup = endTrials(end) - startTrials(end);
        if nTraces_lastGroup > 0 && nTraces_lastGroup < options.groupSize/2
            endTrials(end) = [];
            options.nGroups = options.nGroups - 1;
        end
    end

elseif sum(strcmpi(options.groupby,["animal","animals"]))
    if ~isfield(options,'animalStartIdx') && ~isfield(options.startIdx)
        error('startIdx required if groupby=animal');
    end
    options.nGroups = length(options.animalStartIdx);

    % Calculate start and end trials
    startTrials = options.animalStartIdx;
    endTrials = [options.animalStartIdx(2:options.nGroups); nTraces];
    if endTrials(end) > nTraces; endTrials(end) = nTraces; end

elseif sum(strcmpi(options.groupby,["session","sessions"]))
    if ~isfield(options,'sessionStartIdx') && ~isfield(options.startIdx)
        error('startIdx required if groupby=session');
    end
    % Reorganize to find start and end trials
    options.nGroups = round(length(options.sessionStartIdx)/length(options.animalStartIdx)); % nSessions per animal
    trialStart_sessions = nan(length(options.nGroups),length(options.animalStartIdx));
    trialEnd_sessions = nan(length(options.nGroups),length(options.animalStartIdx));

    sessionInAnimal = 0; animalIdx = 0;
    for i = 1:length(options.sessionStartIdx)
        if ismember(options.sessionStartIdx(i),options.animalStartIdx)
            animalIdx = animalIdx + 1;
            sessionInAnimal = 0;
        end
        sessionInAnimal = sessionInAnimal + 1;
        trialStart_sessions(sessionInAnimal,animalIdx) = options.sessionStartIdx(i);
        if i == length(options.sessionStartIdx); trialEnd_sessions(sessionInAnimal,animalIdx) = nTraces;
        else; trialEnd_sessions(sessionInAnimal,animalIdx) = options.sessionStartIdx(i+1)-1; end
    end
else
    error('Input to options.groupby is wrong!');
end


%% Plot data
legendList = cell(options.nGroups,1);
nColors = round(linspace(1,size(colormap,1),options.nGroups));

for i = 1:options.nGroups
    if isfield(options,'animalStartIdx') && length(options.animalStartIdx)>1 && sum(strcmpi(options.groupby,["trial","trials"]))
        startTrial = startTrials(i); endTrial = endTrials(i);
        plotData = cell2mat(cellfun(@(x,dim) x(startTrial:min(endTrial,dim),:),traces_animals,num2cell(traces_animals_dim),'UniformOutput',false)');
        plotSEM(timestamp,plotData,colormap(nColors(i),:),LineStyle=options.LineStyle,LineWidth=options.LineWidth);
        legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial),' (n=',num2str(size(plotData,1)),')'];
    elseif sum(strcmpi(options.groupby,["session","sessions"]))
        trialWindow = [];
        for animal = 1:size(trialStart_sessions,2)
            trialWindow = [trialWindow,trialStart_sessions(i,animal):trialEnd_sessions(i,animal)];
        end
        plotData = traces(trialWindow,:);
        plotSEM(timestamp,plotData,colormap(nColors(i),:),LineStyle=options.LineStyle,LineWidth=options.LineWidth);
        legendList{i} = ['Session ', num2str(i),' (n=',num2str(size(plotData,1)),')'];
    else
        startTrial = startTrials(i); endTrial = endTrials(i);
        plotData = traces(startTrial:endTrial,:);
        plotSEM(timestamp,plotData,colormap(nColors(i),:),LineStyle=options.LineStyle,LineWidth=options.LineWidth);

        if sum(strcmpi(options.groupby,["trial","trials"]))
            legendList{i} = ['Trial ', num2str(startTrial),'-',num2str(endTrial),' (n=',num2str(size(plotData,1)),')'];
        elseif sum(strcmpi(options.groupby,["animal","animals"]))
            legendList{i} = ['Animal ', num2str(i),' (n=',num2str(size(plotData,1)),')'];
        end
    end
end

end