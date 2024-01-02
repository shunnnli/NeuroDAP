function analysis = analyzeTraces(timeSeries,lick,analysisEvents,analysisLabels,params,options)

arguments
    timeSeries struct
    lick double
    analysisEvents cell
    analysisLabels cell
    params 
    
    options.task string = "NaN"
    options.timeRange double = [-15, 15]
    options.stageTime double = [-2,0;0,2]
    options.lick_binSize double = 0.1
    options.nboot double = 10000
    
    options.trialTable table
    options.trialNumber cell

end

disp('Ongoing: analyze traces and saved in anlaysis struct');
options.analysisDate = char(datetime('today','Format','yyyyMMdd'));

% Remove empty indices
if isfield(options,'trialNumber')
    options.trialNumber = options.trialNumber(~cellfun('isempty',analysisEvents));
end
analysisLabels = analysisLabels(~cellfun('isempty',analysisEvents));
analysisEvents = analysisEvents(~cellfun('isempty',analysisEvents));

%% Define analysis params
analysis = struct([]);
timeRange = options.timeRange;
stageTime = options.stageTime;

% Define tasks and other properties
if isfield(params.session,'task')
    options.task = params.session.task;
end

% Make sure trialTable and trialNumber is provided
if ~isfield(options,'trialNumber') || ~isfield(options,'trialTable')
    warning('analyzeTraces: Strongly recommend to include trialTable and trialNumber!')
    disp('Ongoing: starting a 30s timer, stop and modify the code if needed!');
    pause(30);
    disp('Finished: 30s timer passed, code resume as is');
end

% Change all trialNumber to nx1 vector
for i = 1:length(options.trialNumber)
    options.trialNumber{i} = reshape(options.trialNumber{i},[length(options.trialNumber{i}),1]);
end

%% Loop through all events
for i = 1:length(analysisEvents)
    if isempty(analysisEvents{i})
        disp(['     Ongoing: ',analysisLabels{i},' is empty, skipped!']);
        continue
    end
    %% 1. Loop through all timeSeries
    for signal = 1:size(timeSeries,2)
        disp(['     Ongoing: analyze ',analysisLabels{i},' in ',timeSeries(signal).name,' signal']);

        % Load signal of interest
        data = timeSeries(signal).data;
        finalFs = round(timeSeries(signal).finalFs);
        system = timeSeries(signal).system;

        % Save overall traces
        [trace,t] = plotTraces(analysisEvents{i},timeRange,data,params,...
                        signalFs=finalFs,signalSystem=system,plot=false);
        if isempty(trace)
            disp(['     Ongoing: ',analysisLabels{i},' is empty, skipped!']); 
            continue; 
        end

        % analyzeStages
        stats = analyzeStages(trace,stageTime,nboot=options.nboot,finalFs=finalFs);

        % Save traces and anlaysis data
        row = size(analysis,2) + 1;
        analysis(row).animal = params.session.animal;
        analysis(row).date = params.session.date;
        analysis(row).session = params.session.name;
        analysis(row).task = options.task;
        analysis(row).event = analysisLabels{i};
        analysis(row).name = timeSeries(signal).name;
        analysis(row).system = system;
        analysis(row).data = trace;
        analysis(row).timestamp = t;
        analysis(row).timeRange = timeRange;
        analysis(row).finalFs = finalFs;

        analysis(row).stageAvg = stats.stageAvg;
        analysis(row).stageMax = stats.stageMax;
        analysis(row).stageMin = stats.stageMin;
        
        analysis(row).trialInfo.trialNumber = options.trialNumber{i};
        if ~isempty(options.trialTable)
            analysis(row).trialInfo.trialTable = options.trialTable(options.trialNumber{i},:);
        end
        analysis(row).options = options;
    end
    
    %% 2. Get lick events
    disp(['     Ongoing: analyze ',analysisLabels{i},' in lick signal']);
    [lickRate,lickEvents,t_licks] = plotLicks(analysisEvents{i},timeRange,...
                    options.lick_binSize,[1 1 1],[],lick,params,plot=false);

    % analyzeStages
    stats_lick = analyzeStages(lickRate,stageTime,nboot=options.nboot,finalFs=1/options.lick_binSize);

    % Save lick traces
    row = size(analysis,2) + 1;
    analysis(row).animal = params.session.animal;
    analysis(row).date = params.session.date;
    analysis(row).session = params.session.name;
    analysis(row).task = options.task;
    analysis(row).event = analysisLabels{i};
    analysis(row).name = 'Lick';
    analysis(row).system = 'Lick';
    analysis(row).data.lickRate = lickRate;
    analysis(row).data.lickEvents = lickEvents;
    analysis(row).timestamp = t_licks;
    analysis(row).timeRange = timeRange;
    analysis(row).finalFs = 1/options.lick_binSize;

    analysis(row).stageAvg = stats_lick.stageAvg;
    analysis(row).stageMax = stats_lick.stageMax;
    analysis(row).stageMin = stats_lick.stageMin;
    
    analysis(row).trialInfo.trialNumber = options.trialNumber{i};
    if ~isempty(options.trialTable)
        analysis(row).trialInfo.trialTable = options.trialTable(options.trialNumber{i},:);
    end
    analysis(row).options = options;
end

% Save
save(strcat(params.session.path,filesep,'analysis_',params.session.name),'analysis','-append');
disp('Finished: analysis struct created and saved');

end