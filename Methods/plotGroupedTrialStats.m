function plotGroupedTrialStats(stats,ylabels,options)

arguments
    stats
    ylabels

    options.stage double = 2

    options.groupSize double = 10
    options.color
    options.eventRange   
    options.inTrialTable logical

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

%% Plot
for task = 1:length(stats)
    % Plot
    nexttile; 
    stats_combined = stats{task};

    for event = 1:length(options.eventRange)
        % Find the animal with the least number of trials
        event_stats = stats_combined(:,event);
        event_stats = event_stats(~cellfun('isempty',event_stats));
        minTrials = min(cellfun(@(x) x(1),cellfun(@size,event_stats,'UniformOutput',false)));
        nGroups = floor(minTrials/options.groupSize);
        nCommonTrials = options.groupSize * nGroups;
        x = linspace(options.groupSize,nCommonTrials,nGroups);

        % Initialize grouped array
        grouped = nan(length(event_stats),nGroups);

        % Calculate grouped average per animal
        for animal = 1:length(event_stats)
            if options.inTrialTable; data = event_stats{animal}(1:nCommonTrials);
            else; data = event_stats{animal}(1:nCommonTrials,options.stage); end
            grouped(animal,:) = mean(reshape(data,options.groupSize,[]),1);
        end

        % Plot average data
        scatter(x,mean(grouped,1),200,options.color{event},'filled',HandleVisibility='off'); hold on
        plotTraces(grouped,x,color=options.color{event},extract=false,...
            meanOnly=options.meanOnly,...
            plotIndividual=options.plotIndividual,...
            individualColor=options.individualColor);
    end

    xlabel('Trials'); ylabel(ylabels{task});
    legend(options.eventRange);
end

end