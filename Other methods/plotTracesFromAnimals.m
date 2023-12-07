function [traces,timestamp] = plotTracesFromAnimals(animals,animalList,region,eventName,options)

% plot traces from struct animal
arguments
    animals struct
    animalList cell
    region string
    eventName string
    options.task string = ''
    
    % traceRange refers to the whole matrix of the animal 
    % IF option.task = ''
    % Otherwise, it refers to trials counting from the start of the same
    % type (ie 1:30 trials for reward->punish task)
    options.traceRange double = []
    options.plot logical = true % plot traces or not
    options.color
    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 2
    options.baselineWindow double = 0 % x secs before event time
end


% Get traces name
name_trace = region + '_' + eventName;
name_time = region + '_timestamp';
% Initialize matrix
traces = [];

% Get animal traces
for i = 1:length(animalList)
    % Find row number for that animal
    row = contains({animals.mouse},animalList{i});
    
    % Find timestemp
    timestamp = animals(row).traces.(name_time);
    
    % Extract traces
    sessionTraces = animals(row).traces.(name_trace);
    
    % Find trial range
    if isempty(options.task)
        if isempty(options.traceRange)
            traceRange = 1:size(sessionTraces,1);
        else
            if options.traceRange(end) > size(sessionTraces,1); traceRange = 1:size(sessionTraces,1);
            else; traceRange = options.traceRange; end
        end
    else
        taskRange = animals(row).taskRange.(options.task).(eventName);
        if isempty(options.traceRange) || options.traceRange(end) >= length(taskRange)
            if options.traceRange(1) > length(taskRange); traceRange = [];
            else; traceRange = options.traceRange(1):taskRange(end); end
        else
            traceRange = taskRange(options.traceRange);
        end
    end
    
    traces = [traces;sessionTraces(traceRange,:)];
end

traces = rmmissing(traces);
% Plot traces
if options.plot
    plotSEM(timestamp,traces,options.color,LineStyle=options.LineStyle,LineWidth=options.LineWidth); hold on
end

end