function results = getTracesStats(traces,timeWindows,options)

arguments
    traces double % traces to analyze (row is trial, column is time)
    timeWindows double % time windows to analyze stat on
    options.returnStruct logical = true
    options.timeRange double
    options.timestamp double
end

% getTracesStat() gets a trace and return relevant stats
average = nan(size(traces,1),5); variance = nan(size(traces,1),5);
maximum = nan(size(traces,1),5); minimum = nan(size(traces,1),5);
maxloc = nan(size(traces,1),5); minloc = nan(size(traces,1),5);

% Calculate time window index within one animal
if isempty(options.timestamp)
    if ~isempty(options.timeRange)
        % if without time stamp, calculate directly
        timestamp = linspace(options.timeRange(1),options.timeRange(2),size(traces,2));
    else 
        error('timestamp or timeRange not found!');
    end
else
    timestamp = options.timestamp;
end

% Loop through all windows
for w = 1:size(timeWindows,1)
    
    % Calculate timeWindow range
    cur_window = find(timestamp >= timeWindows(w,1) & timestamp < timeWindows(w,2));
    
    % Calculate window stat
    average(:,w) = mean(traces(:,cur_window),2); 
    variance(:,w) = var(traces(:,cur_window),0,2);
    [maximum(:,w),maxloc(:,w)] = max(traces(:,cur_window),[],2); 
    [minimum(:,w),minloc(:,w)] = min(traces(:,cur_window),[],2);
end

if options.returnStruct
    results.avg = average;
    results.var = variance;
    results.max = maximum;
    results.min = minimum;
    results.maxloc = maxloc;
    results.minloc = minloc;
else
    results = {average;variance;maximum;minimum;maxloc;minloc};
end

end