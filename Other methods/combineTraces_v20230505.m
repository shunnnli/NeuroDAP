function combined = combineTraces(traces,options)

arguments
    traces cell
    options.sessionRange double = 1:size(traces,2)
    options.trialRange double = nan
    options.timeRange double = 1:size(traces{1},2)
    options.merge logical = true % whether or not to merge into one big matrix
end

% params
sessionRange = options.sessionRange;
trialRange = options.trialRange;
timeRange = options.timeRange;

% 1. select traces based on session
session_traces = traces(1,sessionRange);

% 2. select traces based on trial range
if isnan(options.trialRange) 
    trial_traces = cellfun(@(x) x(:,timeRange),session_traces,'UniformOutput',false);
else
    trial_traces = cellfun(@(x) x(trialRange,timeRange),session_traces,'UniformOutput',false);
end

% 3. merge to one matrix if needed
if options.merge
    combined = vertcat(trial_traces{:});
else
    combined = trial_traces;
end

end