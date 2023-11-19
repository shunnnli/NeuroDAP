function combined = combineTraces(summary,options)

arguments
    summary struct
    options.timeRange = 'All'
    options.eventRange string = 'All'
    options.animalRange string = 'All'
    options.taskRange string = 'All'
    options.sessionRange string = 'All'
    options.trialRange = 'All'
    options.signalRange string = 'All'
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
    taskIdx = cellfun(@(x) contains(x,options.taskRange), {summary(finalIdx).task});
    finalIdx = finalIdx(taskIdx);
end

% Select rows based on sessionRange
if ~strcmpi(options.sessionRange,'All')
    sessionIdx = cellfun(@(x) contains(x,options.sessionRange), {summary(finalIdx).session});
    finalIdx = finalIdx(sessionIdx);
end

% Select rows based on events
if ~strcmpi(options.eventRange,'All')
    eventIdx = cellfun(@(x) contains(x,options.eventRange), {summary(finalIdx).event}); % change to event
    finalIdx = finalIdx(eventIdx);
end

% Select rows based on signal
if strcmpi(options.signalRange,'All')
    options.signalRange = unique({summary(finalIdx).name});
else
    if sum(cellfun(@(x) contains(x,options.signalRange), {summary(finalIdx).name})) == 0
        error('     combineTraces: analysis struct does not contain input signalRange');
    end
end

% Concat .data in each selected rows into a array for plotting
data = cell(length(options.signalRange),1);
for signal = 1:length(options.signalRange)

    % Select rows for the current signal
    signalRows = cellfun(@(x) contains(x,options.signalRange{signal}), {summary(finalIdx).name});
    signalIdx = finalIdx(signalRows);

    % Check whether Fs is the same
    uniqueFs = unique(cell2mat({summary(signalIdx).finalFs}));
    if length(uniqueFs) > 1
        error('     combineTraces: selected Fs are not consistent!');

        % can add automatic resampling later
    end

    % Loop through to combine data
    for i = 1:length(signalIdx)
        row = summary(signalIdx(i));
        
        % Check timeRange is valid
        if row.timeRange(1) > options.timeRange(1) || row.timeRange(2) < options.timeRange(2)
            error('     combineTraces: options.timeRange exceeds analysis.timeRange for this session');
        end

        % Find eventSample: option 1
        eventSample = abs(row.timeRange(1))*row.finalFs + 1;
        firstSample = eventSample + options.timeRange(1)*row.finalFs;
        lastSample = eventSample + options.timeRange(2)*row.finalFs;
        t = row.timestamp(firstSample:lastSample);

        % Find trialRange
        if strcmpi(options.trialRange,"All"); options.trialRange = 1:size(row.data,1);
        else; options.trialRange = 1:min(options.trialRange(end),size(row.data,1)); end
    
        % Combined .data
        data{signal} = [data{signal}; row.data(options.trialRange,firstSample:lastSample)]; 
    end
end

combined.data = data;
combined.timestamp = t;
combined.options = options;

end