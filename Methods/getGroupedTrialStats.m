function combinedStats = getGroupedTrialStats(animals,statsTypes,options)

arguments
    animals struct
    statsTypes
    
    options.concatSessions logical = true
    
    % combineTraces options
    options.eventRange string = 'All'
    options.animalRange string = 'All'
    options.taskRange string = 'All'
    options.sessionRange string = 'All'
    options.signalRange string = 'All'
    options.totalTrialRange = 'All'
    options.trialRange = 'All' % index from totalTrialRange
    % options.trialConditions

    options.inTrialTable
end

%% Parse combineTraces inputs
if strcmpi(options.animalRange,'All'); options.animalRange = unique({animals.animal}); end
if strcmpi(options.eventRange,'All'); options.eventRange = unique({animals.event}); end

if strcmpi(options.taskRange,'All'); options.taskRange = unique({animals.task});
elseif length(options.taskRange)==1; options.taskRange = {options.taskRange}; end

if strcmpi(options.sessionRange,'All'); options.sessionRange = repelem("All",length(options.taskRange),1); end
if strcmpi(options.totalTrialRange,'All'); options.totalTrialRange = repelem("All",length(options.taskRange),1); end
if strcmpi(options.trialRange,'All'); options.trialRange = repelem("All",length(options.taskRange),1); end

if ~iscell(statsTypes); statsTypes = repelem(string(statsTypes),length(options.taskRange),1); end

%% get stats
stats = cell(length(options.taskRange),1);

for task = 1:length(options.taskRange)
    stats_combined = cell(length(options.animalRange),length(options.eventRange));
    animalList_task = unique({animals(contains({animals.task},options.taskRange{task},IgnoreCase=true)).animal});
    animalList = intersect(animalList_task,options.animalRange);

    % Get statsType
    statsType = convertStringsToChars(statsTypes{task});
    if contains(statsType,"stage")
        options.inTrialTable = false;
        if ~strcmpi(statsType,{'stageAvg','stageMax','stageMin'})
            warning('Not a valid statsType input, check for typos! Changed to stageAvg by default!');
            statsType = 'stageAvg';
        end
    else
        options.inTrialTable = true;
    end
    
    for event = 1:length(options.eventRange)
        for animal = 1:length(animalList)
            combined = combineTraces(animals,...
                                    eventRange=options.eventRange{event},...
                                    animalRange=animalList{animal},...
                                    taskRange=options.taskRange{task},...
                                    signalRange=options.signalRange,...
                                    totalTrialRange=options.totalTrialRange,...
                                    trialRange=options.trialRange);
            if combined.options.empty; continue; end

            % Extract stats
            if options.inTrialTable
                if ~any(strcmpi(statsType,combined.trialTable{1}.Properties.VariableNames))
                    error('Not a valid statsType input, check for typos!');
                end
                statsData = combined.trialTable{1}.(statsType);
            else
                statsData = combined.stats.(statsType){1};
            end

            % Separate sessions if necessary
            if options.concatSessions
                stats_combined{animal,event} = statsData;
            else
                sessionStartIdx = combined.options.startIdx.session{1};
                sessionStats = cell(length(sessionStartIdx),1);
        
                for session = 1:length(sessionStartIdx)
                    if session == length(sessionStartIdx); lastTrial = length(combined.trialNumber{1}); 
                    else; lastTrial = sessionStartIdx(session+1)-1; end
                    sessionWindow = sessionStartIdx(session):lastTrial;
        
                    if options.inTrialTable
                        sessionStats{session} = [sessionStats{session}; statsData(sessionWindow)];
                    else
                        sessionStats{session} = [sessionStats{session}; statsData(sessionWindow,:)];
                    end
                end
                stats_combined{animal,event} = [stats_combined{animal,event}, sessionStats];
            end
        end 
    end
    stats{task} = stats_combined;
end

combinedStats.stats = stats;
combinedStats.options = options;

end