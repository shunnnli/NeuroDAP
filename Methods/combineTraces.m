function combined = combineTraces(database,options)

arguments
    database struct
    
    options.timeRange double = [-15,15]
    options.eventRange string = 'All'
    options.animalRange string = 'All'
    options.taskRange string = 'All'
    options.sessionRange string = 'All'
    options.signalRange string = 'All'
    options.totalTrialRange = 'All'
    options.trialRange = 'All' % index from totalTrialRange

    % trialRange idx explainer:
    % if options.totalTrialRange = [10,50]; found 10 trials [13,15,29,32,..,36,40,60]
    % Normal usecase:
    % if options.trialRange = [1,3]; get trials [13,15,29]
    % if options.trialRange = [-3,0]; get trials [36,40,60]
    % if options.trialRange = [-3,-1]; get trials [36,40]
    % Edge scenarios:
    % if options.trialRange = [2,20]; get trials [15,29,...60]
    % if options.trialRange = [20,23]; get 0 trials
    % if options.trialRange = [20,53]; get 0 trials
    % if options.trialRange = [0/-1,3]; get trials [13,15,29,32] (ie get [1,3-0/3-(-1)])

    options.statsType string = 'All'
    options.withStageArea logical = true
    options.empty logical = false

    options.trialConditions string
    options.historyConditions string
    options.trialTables table % total trial table generated from loadTrialTables()
end

%% Check input
if options.timeRange(1) >= options.timeRange(2)
    error('Input timeRange error: timeRange(1) >= timeRange(2)');
end

if ~strcmpi(options.trialRange,'All')
    if options.trialRange(1) >= options.trialRange(2)
        error('Input trialRange error: trialRange(1) >= trialRange(2)');
    elseif length(options.trialRange) > 2
        options.trialRange(2) = options.trialRange(end);
        options.trialRange = options.trialRange(1:2);
        disp(['combineTraces: trialRange should be [start,end], changed to trialRange=[',num2str(options.trialRange),']']);
    end
end

if ~strcmpi(options.totalTrialRange,'All')
    if options.totalTrialRange(1) >= options.totalTrialRange(2)
        error('Input totalTrialRange error: totalTrialRange(1) >= totalTrialRange(2)');
    elseif length(options.totalTrialRange) > 2
        options.totalTrialRange(2) = options.totalTrialRange(end);
        options.totalTrialRange = options.totalTrialRange(1:2);
        disp(['combineTraces: totalTrialRange should be [start,end], changed to totalTrialRange=[',num2str(options.totalTrialRange),']']);
    end
end

if isfield(options,'historyConditions') && ~isfield(options,'trialTables')
    error('You need to provide trialTables when using historyConditions!');
end

%% Select rows based on input range

% Select rows based on animalRange
if ~strcmpi(options.animalRange,'All')
    animalIdx = find(cellfun(@(x) contains(x,options.animalRange,IgnoreCase=true), {database.animal}));
    finalIdx = animalIdx;
else
    options.animalRange = unique({database.animal});
    finalIdx = 1:size(database,2);
end

% Select rows based on options.taskRange
if ~strcmpi(options.taskRange,'All')
    taskIdx = cellfun(@(x) contains(x,options.taskRange,IgnoreCase=true), {database(finalIdx).task});
    finalIdx = finalIdx(taskIdx);
else
    options.taskRange = unique({database.task});
end

% Select rows based on sessionRange if database have field session (ie summary struct)
if isfield(database,'session')
    if  ~any((strcmpi(options.sessionRange,'All')))
        sessionIdx = cellfun(@(x) contains(x,options.sessionRange,IgnoreCase=true), {database(finalIdx).session});
        finalIdx = finalIdx(sessionIdx);
        % will save options.sessionList later
    else
        options.sessionList = unique({database.session});
    end
else
    options.sessionList = {};
end


% Select rows based on events
if ~strcmpi(options.eventRange,'All')
    eventIdx = cellfun(@(x) contains(x,options.eventRange,IgnoreCase=true), {database(finalIdx).event});
    % eventIdx = cellfun(@(x) strcmpi(x,options.eventRange), {database(finalIdx).event});
    finalIdx = finalIdx(eventIdx);
else
    options.eventRange = unique({database.event});
end

% Select rows based on signal
if strcmpi(options.signalRange,'All')
    options.signalRange = unique({database(finalIdx).name});
else
    if sum(cellfun(@(x) contains(x,options.signalRange,IgnoreCase=true), {database(finalIdx).name})) == 0
        warning('Did NOT find rows that fits the input range, skipped');
        options.empty = true;
        combined.options = options;
        disp(options);
        return
    end
end

% Update options.statsType
if strcmpi(options.statsType,'All')
    % if ~options.withStageArea; options.statsType = {'stageAvg','stageMax','stageMin'}; 
    % else; options.statsType = {'stageAvg','stageMax','stageMin','stageArea'}; end
    database_fieldnames = fieldnames(database);
    options.statsType = database_fieldnames(contains(database_fieldnames,'stage'));
end
if ~iscell(options.statsType); options.statsType = {options.statsType}; end

%% Initialize data/stats/trialNum arrays

data = cell(length(options.signalRange),1);
for type = 1:length(options.statsType)
    stats.(options.statsType{type}) = cell(length(options.signalRange),1);
end
trialNumData = cell(length(options.signalRange),1);
trialTableData = cell(length(options.signalRange),1);
options.signalRows = [];

% To record session/animal changes
options.startIdx.animal = cell(length(options.signalRange),1); prev_animal = '';
options.startIdx.session = cell(length(options.signalRange),1);
options.startIdx.stage = cell(length(options.signalRange),1);

% To record stageTime
options.stageTime = {};

%% Concat .data in each selected rows into a array for plotting
for signal = 1:length(options.signalRange)

    % Select rows for the current signal
    signalRows = cellfun(@(x) contains(x,options.signalRange{signal},IgnoreCase=true), {database(finalIdx).name});
    signalIdx = finalIdx(signalRows);
    options.signalRows = [options.signalRows; signalIdx];

    % Update sessionList (for generating animal struct specifically)
    if isfield(database,'session')
        options.sessionList = unique({database(signalIdx).session});
    end

    % Check whether Fs is the same
    uniqueFs = unique(cell2mat({database(signalIdx).finalFs}));
    if length(uniqueFs) > 1
        error('Selected Fs are not consistent!');
        % can add automatic resampling later
    end

    % Loop through to combine data
    for i = 1:length(signalIdx)
        row = database(signalIdx(i));
        if (iscell(row.data) || isstruct(row.data)) && strcmpi(row.system,'Lick')
            rowData = row.data.lickRate;
        else
            rowData = row.data; 
        end

        % If enter a new animal, update animalStartIdx
        if ~strcmpi(prev_animal,row.animal)
            options.startIdx.animal{signal} = [options.startIdx.animal{signal}; size(data{signal},1)+1];
            prev_animal = row.animal;
        end

        % To record stage time
        if any(cellfun(@(c) isequal(c, row.options.stageTime), options.stageTime))
            options.stageTime{end+1} = row.options.stageTime;
            options.startIdx.stage{signal} = [options.startIdx.stage{signal},size(data{signal},1)+1]; % Not tested
        end

        % Check timeRange is valid
        if row.timeRange(1) > options.timeRange(1) || row.timeRange(2) < options.timeRange(2)
            error('options.timeRange exceeds analysis.timeRange for this session');
        end

        % Find eventSample: option 1
        eventSample = abs(row.timeRange(1))*row.finalFs + 1;
        firstSample = eventSample + options.timeRange(1)*row.finalFs;
        lastSample = eventSample + options.timeRange(2)*row.finalFs;
        t = row.timestamp(firstSample:lastSample);

        % Find trialRange
        if strcmpi(options.trialRange,'All') & strcmpi(options.totalTrialRange,'All')
            trialRange = 1:size(rowData,1); 
        elseif strcmpi(options.trialRange,'All') & ~strcmpi(options.totalTrialRange,'All')
            trialRange = find(row.trialInfo.trialNumber >= options.totalTrialRange(1) & row.trialInfo.trialNumber <= options.totalTrialRange(2));
        elseif ~strcmpi(options.trialRange,'All') & strcmpi(options.totalTrialRange,'All')
            trialRange = max(1,options.trialRange(1)):min(options.trialRange(2),size(rowData,1));
        elseif ~strcmpi(options.trialRange,'All') & ~strcmpi(options.totalTrialRange,'All')
            totalTrialRange = find(row.trialInfo.trialNumber >= options.totalTrialRange(1) & row.trialInfo.trialNumber <= options.totalTrialRange(2));
            if options.trialRange(1) <= 0
                if options.trialRange(2) <= 0
                    trialRange = totalTrialRange(end+options.trialRange(1)+1:end+options.trialRange(2));
                else
                    warning('trialRange format error: input something like [-4,1], changed to [1,6]');
                    options.trialRange = options.trialRange + abs(options.trialRange(1)) + 1;
                    trialRange = totalTrialRange(options.trialRange(1):min(options.trialRange(2),length(totalTrialRange)));
                end
            else
                if options.trialRange(1) >= length(totalTrialRange)
                    if options.trialRange(2)-options.trialRange(1) <= length(totalTrialRange)
                        warning('trialRange format error: input something like [end,end+3], changed to [end-3,0]');
                        options.trialRange = [options.trialRange(1)-options.trialRange(2),0];
                        trialRange = totalTrialRange(end+options.trialRange(1)+1:end+options.trialRange(2));
                    else
                        warning('trialRange format error: input something like [end,end+10000], changed to All');
                        trialRange = totalTrialRange;
                    end
                    continue
                elseif options.trialRange(2) > length(totalTrialRange)
                    warning('trialRange format error: input something like [2,end+10], changed to [2,end]');
                    trialRange = totalTrialRange(options.trialRange(1):length(totalTrialRange));
                else
                    trialRange = totalTrialRange(options.trialRange(1):options.trialRange(2));
                end
            end
        end

        % Check sessionRange if database do not have field sessions (ie animal struct)
        if ~isfield(database,'session') && ~any((strcmpi(options.sessionRange,'All')))
            if ~isfield(row.options,'sessionList')
                disp('No options.sessionList for this row found, set to empty!'); 
                row.options.sessionList = {};
            end
            sessionIdx = cellfun(@(x) contains(x,options.sessionRange,IgnoreCase=true), row.options.sessionList);
            options.sessionList = [options.sessionList,row.options.sessionList(sessionIdx)];
            sessionTrials = [];
            for s = 1:length(sessionIdx)
                if sessionIdx(s)
                    startTrial = row.options.startIdx.session{signal}(s);
                    if s == length(sessionIdx); endTrial = size(rowData,1);
                    else; endTrial = row.options.startIdx.session{signal}(s+1)-1; end
                    sessionTrials = [sessionTrials,startTrial:endTrial];
                end
            end
            trialRange = intersect(trialRange,sessionTrials);
            if isempty(trialRange); warning('Did not find trials that fits the sessionRange!'); end
        else
            if ~isfield(row.options,'sessionList')
                row.options.sessionList = options.sessionList;
            end
            options.sessionList = [options.sessionList,row.options.sessionList];
        end

        % Check trialConditions in trialTable
        if isfield(options,'trialConditions')
            trials = row.trialInfo.trialTable; % required!!
            trialRange = intersect(trialRange,find(eval(options.trialConditions)));
            if isempty(trialRange); warning('Did not find trials that fits trialConditions! Check expression!'); end
        end

        % Check historyConditions in trialTable
        % if isfield(options,'historyConditions')
        %     sessionTrials = options.trialTable.data(strcmpi(options.trialTable.name,row.session));
        % end
    
        % Combined .data
        data{signal} = [data{signal}; rowData(trialRange,firstSample:lastSample)];  
        options.finalFs = uniqueFs;
        options.system = row.system;

        % Combine stats
        for type = 1:length(options.statsType)
            statsData = row.(options.statsType{type}).data;
            stats.(options.statsType{type}){signal} = [stats.(options.statsType{type}){signal}; statsData(trialRange,:)];
        end

        % Combine getTrialNum
        if ~isempty(row.trialInfo.trialNumber) && ~isempty(row.trialInfo.trialTable)
            cur_trialNumber = row.trialInfo.trialNumber(trialRange);
            trialNumData{signal} = [trialNumData{signal}; reshape(cur_trialNumber,[length(cur_trialNumber),1])];
            trialTableData{signal} = [trialTableData{signal}; row.trialInfo.trialTable(trialRange,:)];
        else
            warning('combineTraces: empty trialNumber or trialTable found, skipped for now');
            disp(options);
        end
    end
    options.startIdx.session{signal} = find([-1; diff(trialNumData{signal})]<0); % Not working for animals
end

%% Save
combined.data = data;
combined.stats = stats;
combined.timestamp = t;

combined.trialNumber = trialNumData;
combined.trialTable = trialTableData;

combined.options = options;


end