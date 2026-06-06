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
    options.ylim double
    options.dotSize double = 200

    % plotTraces options
    options.meanOnly logical = false
    options.plotIndividual logical = false
    options.individualColor = [0.8, 0.8, 0.8]
end

%% Parse inputs
if ~iscell(ylabels); ylabels = repelem(string(ylabels),length(options.taskRange),1); end

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

%% Initialize statistics
traces = cell(length(stats),1);
avgStats = cell(length(stats),1);

%% Loop through tasks
for task = 1:length(stats)
    % Plot grouped traces
    nexttile; 
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
        xlimList(event) = nCommonTrials;
        groupedStats{event} = animalStats;

        % Plot average data
        if iscell(options.individualColor); individualColor = options.individualColor{event};
        else; individualColor = options.individualColor; end
        scatter(x,mean(grouped,1),options.dotSize,options.color{event},'filled',HandleVisibility='off'); hold on
        plotTraces(grouped,x,color=options.color{event},extract=false,...
            meanOnly=options.meanOnly,...
            plotIndividual=options.plotIndividual,...
            individualColor=individualColor);
    end

    % Adjust ylim
    if isfield(options,'ylim'); ylim(options.ylim); end

    % Adjust xlim
    if options.xlimIdx ~= 0
        xlim([0,xlimList(options.xlimIdx)]);
    end 
    xlabel('Trials'); ylabel(ylabels{task});
    legend(options.eventRange);

    % Save results
    traces{task} = groupedTraces;
    avgStats{task} = groupedStats;

%% Save results
results.traces = traces;
results.stats = avgStats;

end

end