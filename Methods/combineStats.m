function combined = combineStats(summary,options)

arguments
    summary struct
    options.eventRange string = 'All'
    options.animalRange string = 'All'
    options.taskRange string = 'All'
    options.sessionRange string = 'All'
    options.trialRange = 'All'
    options.signalRange string = 'All'

    options.statsType string = 'stageAvg'
end

% Select rows based on animalRange
if ~strcmpi(options.animalRange,'All')
    animalIdx = find(cellfun(@(x) contains(x,options.animalRange), {summary.animal}));
    finalIdx = animalIdx;
else
    finalIdx = 1:size(summary,2);
end

% Select rows based on options.taskRange
if ~strcmpi(options.taskRange,'All')
    taskIdx = cellfun(@(x) strcmpi(x,options.taskRange), {summary(finalIdx).task});
    finalIdx = finalIdx(taskIdx);
end

% Select rows based on sessionRange
if ~strcmpi(options.sessionRange,'All')
    sessionIdx = cellfun(@(x) contains(x,options.sessionRange), {summary(finalIdx).session});
    finalIdx = finalIdx(sessionIdx);
end

% Select rows based on events
if ~strcmpi(options.eventRange,'All')
    eventIdx = cellfun(@(x) contains(x,options.eventRange,IgnoreCase=true), {summary(finalIdx).event}); % change to event
    finalIdx = finalIdx(eventIdx);
end

% Select rows based on signal
if strcmpi(options.signalRange,'All')
    options.signalRange = unique({summary(finalIdx).name});
else
    if sum(cellfun(@(x) contains(x,options.signalRange), {summary(finalIdx).name})) == 0
        error('analysis struct does not contain input range');
    end
end

% Concat .data in each selected rows into a array for plotting
data = cell(length(options.signalRange),1);
for signal = 1:length(options.signalRange)

    % Select rows for the current signal
    signalRows = cellfun(@(x) contains(x,options.signalRange{signal}), {summary(finalIdx).name});
    signalIdx = finalIdx(signalRows);

    % Loop through to combine data
    for i = 1:length(signalIdx)
        row = summary(signalIdx(i));

        % Find trialRange
        if strcmpi(options.trialRange,'All')
            rowData = row.(statsType).data;
            trialRange = 1:size(rowData,1); 
        else
            rowData = row.(statsType).data;
            trialRange = max(1,options.trialRange(1)):min(options.trialRange(end),size(rowData,1)); 
        end
    
        % Combined .data
        data{signal} = [data{signal}; rowData(trialRange,:)]; 
        options.finalFs = uniqueFs;
        options.system = row.system;
    end
end

combined.data = data;
combined.timestamp = t;
combined.options = options;

end