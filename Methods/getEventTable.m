function eventTable = getEventTable(events,params,options)

arguments
    events cell
    params struct
    options.typeName cell = {"trialStart","airpuff","water","lick","tone","redStim"} 
end

% Load events as separte vectors
allTrials = events{1};
typeName = options.typeName;

% Find total events
totalEvents = 0;
for i = 2:length(events)
    totalEvents = totalEvents + length(events{i});
end

% Initialize event table
varTypes = {'string','double','logical','double'};
varNames = {'type','trial','inTrial','time'};
eventTable = table('Size',[totalEvents length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

% Fill in event table
timeNI = params.sync.timeNI;
cur_row = 0;
for i = 0:length(allTrials)
    % Find NI sample range
    if i == 0; niRange = [0,allTrials(1)];
    elseif i == length(allTrials)
        niRange = [allTrials(i),length(timeNI)+1];
    else; niRange = [allTrials(i),allTrials(i+1)]; 
    end

    % Find events within NI sample range
    eventsInRange = cell(length(events),1);
    for j = 2:length(events)
        eventsInRange{j} = events{j}(events{j} >= niRange(1) & events{j} < niRange(2));
    end

    % Fill in these events
    for j = 2:length(eventsInRange)
        if isempty(eventsInRange{j}); continue; end
        rowIdx = cur_row+1 : cur_row+length(eventsInRange{j});
        eventTable{rowIdx,1} = repmat(typeName{j},length(eventsInRange{j}),1);
        eventTable{rowIdx,2} = i;
        eventTable{rowIdx,3} = true;
        eventTable{rowIdx,4} = timeNI(eventsInRange{j})';
        cur_row = cur_row + length(eventsInRange{j});
    end
end

% Sort based on time
eventTable = sortrows(eventTable,"time");

end