function results = plotGroupedTrialStats(stats,ylabels,options)

arguments
    stats
    ylabels

    options.stage double = 2

    options.groupSize double = 10
    options.color
    options.eventRange   
    options.inTrialTable logical
    options.xlimIdx double = 0
    options.xlim double
    options.ylim double
    options.dotSize double = 200

    % plotTraces options
    options.plot logical = true
    options.plotPatch logical = true
    options.plotIndividual logical = false
    options.individualColor = 'gray'
    options.individualAlpha double = 0
    options.plotNextTile logical = true
end

%% Parse inputs
if ~iscell(ylabels); ylabels = repelem(string(ylabels),length(stats),1); end

if isstruct(stats)
    options.inTrialTable = stats.options.inTrialTable;
    options.eventRange = stats.options.eventRange;
    stats = stats.stats;    
else
    if ~isfield(options,'inTrialTable')
        error('Need to provide inTrialTable information!');
    end
    if ~isfield(options,'eventRange')
        error('Need to provide eventRange information!');
    end
end

% Set color
if ~isfield(options,'color')
    options.color = {[.23 .34 .45]};
end
if ~iscell(options.color)
    options.color = {options.color};
end

% Set y labels
if isempty(ylabels)
    warning('y label is not provided!');
    ylabels = {'Stats'};
end

%% Initialize statistics
traces = cell(length(stats),1);
avgStats = cell(length(stats),1);

%% Loop through tasks
for task = 1:length(stats)
    % Plot grouped traces
    if options.plot 
        if options.plotNextTile; nexttile; end
    end
    stats_combined = stats{task};

    % Initialize data matrix
    xlimList = nan(length(options.eventRange),1);
    groupedTraces = cell(length(options.eventRange),1);
    groupedStats = cell(length(options.eventRange),1);

    for event = 1:length(options.eventRange)
        % Find the animal with the least number of trials
        event_stats = stats_combined(:,event);
        event_stats = event_stats(~cellfun('isempty',event_stats));
        minTrials = min(cellfun(@(x) x(1),cellfun(@size,event_stats,'UniformOutput',false)));
        nGroups = floor(minTrials/options.groupSize);
        nCommonTrials = options.groupSize * nGroups;
        
        % Initialize grouped array
        x = linspace(options.groupSize,nCommonTrials,nGroups);
        grouped = nan(length(event_stats),nGroups);
        animalStats = nan(length(event_stats),2);

        % Loop through animals
        for animal = 1:length(event_stats)
            % Calculate grouped average per animal
            if options.inTrialTable; data = event_stats{animal}(1:nCommonTrials);
            else; data = event_stats{animal}(1:nCommonTrials,options.stage); end
            grouped(animal,:) = mean(reshape(data,options.groupSize,[]),1);

            % Extract slopes and intercepts
            animalStats(animal,:) = polyfit(x,grouped(animal,:),1);
        end

        % Save event specific results
        groupedTraces{event} = grouped;
        if options.plot; xlimList(event) = nCommonTrials; end
        groupedStats{event} = animalStats;

        % Plot average data
        if options.plot
            if iscell(options.individualColor); individualColor = options.individualColor{event};
            else; individualColor = options.individualColor; end
            scatter(x,mean(grouped,1),options.dotSize,options.color{event},'filled',HandleVisibility='off'); hold on
            plotTraces(grouped,x,color=options.color{event},extract=false,...
                plotPatch=options.plotPatch,...
                plotIndividual=options.plotIndividual,...
                individualColor=individualColor);
        end
    end

    if options.plot
        % Adjust ylim
        if isfield(options,'ylim')
            if size(options.ylim,1)==1; ylim(options.ylim); 
            else; ylim(options.ylim(task,:)); end
        end
    
        % Adjust xlim
        if isfield(options,'xlim') && options.xlimIdx == 0
            if size(options.xlim,1)==1; xlimit = options.xlim; 
            else; xlimit = options.xlim(task,:); end
            xlim(xlimit);
        elseif options.xlimIdx ~= 0 && ~isfield(options,'xlim')
            xlim([0,xlimList(options.xlimIdx)]);
        elseif ~isfield(options,'xlim') && options.xlimIdx == 0
            options.xlimIdx = 1;
            xlim([0,xlimList(options.xlimIdx)]);
        else
            if size(options.xlim,1)==1; xlimit = options.xlim; 
            else; xlimit = options.xlim(task,:); end
            xlim([0,min(xlimList(options.xlimIdx),xlimit(2))]);
        end
        xlabel('Trials'); ylabel(ylabels{task});
        legend(options.eventRange);
    end

    % Save results
    traces{task} = groupedTraces;
    avgStats{task} = groupedStats;

%% Save results
results.traces = traces;
results.stats = avgStats;

end

end