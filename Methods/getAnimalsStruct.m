function animals = getAnimalsStruct(summary)

arguments
    summary struct
end

animals = struct([]);
animalList = unique({summary.animal});
disp('Finished: animal.mat not found, created a new one');

for a = 1:length(animalList)
    animalIdx = find(cellfun(@(x) contains(x,animalList{a}), {summary.animal}));
    taskList = unique({summary(animalIdx).task});
    
    for task = 1:length(taskList)
        taskIdx = cellfun(@(x) strcmpi(x,taskList{task}), {summary(animalIdx).task});
        taskRows = animalIdx(taskIdx);
        eventList = unique({summary(taskRows).event});

        for event = 1:length(eventList)
            eventIdx = cellfun(@(x) contains(x,eventList{event},IgnoreCase=true), {summary(taskRows).event});
            eventRows = taskRows(eventIdx);
            signalList = unique({summary(eventRows).name});

            for signal = 1:length(signalList)
                disp(['Ongoing: ',animalList{a},' -> ',taskList{task},' -> ',eventList{event},...
                    ' -> ',signalList{signal}]);

                combined = combineTraces(summary,animalRange=animalList{a},...
                                eventRange=eventList{event},...
                                taskRange=taskList{task},...
                                signalRange=signalList{signal},...
                                statsType='All');

                row = size(animals,2) + 1;
                animals(row).animal = animalList{a};
                animals(row).task = taskList{task};
                animals(row).event = eventList{event};
                animals(row).name = signalList{signal};
                animals(row).system = combined.options.system;
                animals(row).data = combined.data{1};
                animals(row).stageAvg.data = combined.stats.stageAvg{1};
                animals(row).stageMax.data = combined.stats.stageMax{1};
                animals(row).stageMin.data = combined.stats.stageMin{1};
                animals(row).timestamp = combined.timestamp;
                animals(row).timeRange = combined.options.timeRange;
                animals(row).finalFs = combined.options.finalFs;
                animals(row).trialInfo.trialNumber = combined.trialNumber{1};
                animals(row).trialInfo.trialTable = combined.trialTable{1};
                animals(row).options = combined.options;
            end
        end
    end
end

end