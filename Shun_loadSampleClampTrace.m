% Shun_loadSampleClampTrace
% Load a processed session and plot example DA and clampIndex traces.
% Edit timeOption and the settings in the corresponding section below.
% Session time is in seconds since the first NI sample; trial times use CueTime.
% Results remain in the workspace as sampleClampTrace. No session files are changed.

clear; close all; loadNeuroDAP;
timeOption = 2; % 1: enter session time in sec; 2: select trials from the trial table

%% Load sessions of choice

% Select session via uipickfiles
sessionPath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'), ...
                         'NumFiles',1,'Prompt','Select a processed clamp session');
if isempty(sessionPath) || isequal(sessionPath,0); return; end
sessionpath = sessionPath{1};
if ~isfolder(sessionpath); error('Select a session folder.'); end
dirsplit = strsplit(sessionpath,filesep);
sessionName = dirsplit{end}; clear dirsplit

% Select the output name if it differs from the session folder
if ~isfile(fullfile(sessionpath,['timeseries_',sessionName,'.mat']))
    fileList = dir(fullfile(sessionpath,'timeseries_*.mat'));
    if isempty(fileList)
        error('No timeseries_*.mat found. Run loadSessions with withClamp=true first.');
    end
    fileIdx = 1;
    if length(fileList) > 1
        [fileIdx,selected] = listdlg('PromptString','Select session output:', ...
            'SelectionMode','single','ListString',{fileList.name},'ListSize',[600,250]);
        if ~selected; return; end
    end
    [~,sessionName] = fileparts(fileList(fileIdx).name);
    sessionName = sessionName(length('timeseries_')+1:end);
end

% Load the same processed signals and sync params as analyzeSessions_clamp
disp(strcat('**********',sessionName,'**********'));
load(strcat(sessionpath,filesep,'timeseries_',sessionName,'.mat'),'timeSeries');
load(strcat(sessionpath,filesep,'sync_',sessionName,'.mat'),'params');

% Trial table is optional when specifying the time window in seconds
trials = table;
behavior = struct();
behaviorFile = fullfile(sessionpath,['behavior_',sessionName,'.mat']);
if isfile(behaviorFile)
    behavior = load(behaviorFile);
    if isfield(behavior,'trials'); trials = behavior.trials;
    elseif isfield(behavior,'trialTable'); trials = behavior.trialTable;
    end
end
timeNI = params.sync.timeNI;
trialRows = [];
disp(['Finished: Session ',sessionName,' loaded']);

% Check all recorded photometry channels, excluding clamp commands and camera
photometryIdx = find(cellfun(@(x) any(strcmpi(x,{'NI','LJ','labjack'})),{timeSeries.system}) & ...
    ~contains({timeSeries.name},'clamp',IgnoreCase=true) & ~cellfun(@isempty,{timeSeries.data}));
if isempty(photometryIdx); error('No photometry channels found in timeSeries.'); end
daIdx = photometryIdx(1); % editable channel for manual time/trial selection; automatic search checks all

%% Define time to plot (opt 1: just enter in sec)

if timeOption == 1
    sampleStart_sec = 112;
    sampleEnd_sec = 172;

    % Convert session seconds to the native NI reference sample
    validateattributes([sampleStart_sec,sampleEnd_sec],{'numeric'},{'real','finite','nonnegative'});
    if sampleEnd_sec <= sampleStart_sec || sampleEnd_sec > timeNI(end)-timeNI(1)
        error('Choose an increasing time window within the session.');
    end
    [~,referenceIdx] = min(abs(timeNI-(timeNI(1)+sampleStart_sec)));
    timeRange = [sampleStart_sec,sampleEnd_sec] - (timeNI(referenceIdx)-timeNI(1));
    trialRows = [];
elseif timeOption ~= 2
    error('timeOption must be 1 (seconds) or 2 (trial table).');
end

%% Define time to plot (opt 2: use trial table / random-session events)

if timeOption == 2
    findBestClampWindow = true; % false: select start/end trials manually

    if ~findBestClampWindow
        firstTrialIdx = 1; %#ok<UNRCH> % enabled for manual selection
        secondTrialIdx = 1;
        trialWindow_sec = [-5,10];
        requiredColumns = {'TrialNumber','CueTime','isReward','isClamp'};
        if ~istable(trials) || ~all(ismember(requiredColumns,trials.Properties.VariableNames))
            error('Manual selection requires the trial table saved by analyzeSessions_clamp.');
        end
        startCondition = trials.isReward == 1 & trials.isClamp == 0;
        endCondition = trials.isReward == 1 & trials.isClamp == 1;
        % startCondition = trials.TrialNumber == 17;
        % endCondition = trials.TrialNumber == 18;
        validateattributes([firstTrialIdx,secondTrialIdx],{'numeric'},{'finite','integer','positive'});
        validateattributes(trialWindow_sec,{'numeric'},{'size',[1,2],'real','finite'});
        validTrials = trials.TrialNumber > 0 & isfinite(trials.CueTime) & ...
            trials.CueTime >= 1 & trials.CueTime <= numel(timeNI);
        firstTrialRows = find(startCondition & validTrials);
        if length(firstTrialRows) < firstTrialIdx; error('Not enough trials matching startCondition.'); end
        firstTrialRow = firstTrialRows(firstTrialIdx);
        secondTrialRows = find(endCondition & validTrials & trials.CueTime >= trials.CueTime(firstTrialRow));
        if length(secondTrialRows) < secondTrialIdx; error('Not enough later trials matching endCondition.'); end
        secondTrialRow = secondTrialRows(secondTrialIdx);
        trialRows = unique([firstTrialRow,secondTrialRow],'stable');
        referenceIdx = trials.CueTime(firstTrialRow);
        firstCue_sec = timeNI(referenceIdx)-timeNI(1);
        secondCue_sec = timeNI(trials.CueTime(secondTrialRow))-timeNI(1);
        sampleStart_sec = firstCue_sec + trialWindow_sec(1);
        sampleEnd_sec = secondCue_sec + trialWindow_sec(2);
        timeRange = [sampleStart_sec,sampleEnd_sec]-firstCue_sec;
        disp(trials(trialRows,{'TrialNumber','CueTime','isReward','isClamp'}));
    end

    % Find a continuous example with flat clamped and strong unclamped responses
    if findBestClampWindow
        windowDurations_sec = 60; % default: a 60 s example; optionally use 30:5:60
        windowStep_sec = 1;
        baselineWindow_sec = [-2,0]; % mean baseline; excludes the endpoint
        responseWindow_sec = [0,2]; % relative to actual water / airpuff / tone delivery
        clampWeight = 2;

        validateattributes(windowDurations_sec,{'numeric'},{'vector','nonempty','finite','>=',30,'<=',60});
        validateattributes(windowStep_sec,{'numeric'},{'scalar','finite','positive'});
        validateattributes(baselineWindow_sec,{'numeric'},{'size',[1,2],'finite','real'});
        validateattributes(responseWindow_sec,{'numeric'},{'size',[1,2],'finite','real'});
        validateattributes(clampWeight,{'numeric'},{'scalar','finite','nonnegative'});
        if diff(baselineWindow_sec) <= 0 || baselineWindow_sec(2) > 0 || diff(responseWindow_sec) <= 0
            error('Use increasing windows, with the baseline ending at or before event delivery.');
        end

        % Reuse the event classification saved by analyzeSessions_clamp.
        % Its clampON_event marks clamp activity within 3 s after each delivery.
        eventFields = {'waterClampIdx','waterUnclampIdx','airpuffClampIdx', ...
                       'airpuffUnclampIdx','toneClampIdx','toneUnclampIdx'};
        if ~all(isfield(behavior,eventFields))
            error('Run analyzeSessions_clamp for this random session to save water/airpuff/tone clamp classifications.');
        end
        eventTypes = ["Water","Water","Airpuff","Airpuff","Tone","Tone"];
        eventSamples = []; eventType = strings(0,1); eventIsClamp = false(0,1);
        for i = 1:length(eventFields)
            idx = unique(double(behavior.(eventFields{i})(:)));
            idx = idx(isfinite(idx) & idx >= 1 & idx <= numel(timeNI) & idx == fix(idx));
            eventSamples = [eventSamples;idx]; %#ok<AGROW>
            eventType = [eventType;repmat(eventTypes(i),length(idx),1)]; %#ok<AGROW>
            eventIsClamp = [eventIsClamp;repmat(mod(i,2) == 1,length(idx),1)]; %#ok<AGROW>
        end
        eventTable = table(eventSamples,eventType,eventIsClamp, ...
            'VariableNames',{'Sample','Type','IsClamp'});
        eventTable = sortrows(eventTable,'Sample');
        if isempty(eventTable)
            disp('No qualifying window exists: no water, airpuff or tone events were found.');
            return
        end
        eventTable.Time_sec = reshape(timeNI(eventTable.Sample)-timeNI(1),[],1);

        candidates = zeros(0,16);
        eventScores = cell(size(timeSeries));
        for daIdx = reshape(photometryIdx,1,[])
            % RMS is direction-independent: peaks and dips both count as responses.
            signal = timeSeries(daIdx).data;
            finalFs = timeSeries(daIdx).finalFs;
            system = timeSeries(daIdx).system;
            extractRange = [min(baselineWindow_sec(1),responseWindow_sec(1)), ...
                            max(baselineWindow_sec(2),responseWindow_sec(2))];
            [eventTraces,eventTime] = plotTraces(eventTable.Sample,extractRange,signal,params, ...
                signalFs=finalFs,signalSystem=system,eventSystem='ni',plot=false,rmmissing=false);
            baselineIdx = eventTime >= baselineWindow_sec(1)-1e-9 & eventTime < baselineWindow_sec(2)-1e-9;
            responseIdx = eventTime >= responseWindow_sec(1)-1e-9 & eventTime <= responseWindow_sec(2)+1e-9;
            if ~any(baselineIdx) || sum(responseIdx) < 2; error('Not enough DA samples in the scoring windows.'); end
            response = eventTraces(:,responseIdx)-mean(eventTraces(:,baselineIdx),2);
            eventTable.RMS = sqrt(mean(response.^2,2));
            eventTable.Valid = all(isfinite(eventTraces),2) & isfinite(eventTable.RMS);

            % Score only events whose full baseline and response fit in the window.
            for duration_sec = reshape(windowDurations_sec,1,[])
                for start_sec = 0:windowStep_sec:(timeNI(end)-timeNI(1)-duration_sec)
                    end_sec = start_sec+duration_sec;
                    included = eventTable.Valid & eventTable.Time_sec+extractRange(1) >= start_sec & ...
                               eventTable.Time_sec+extractRange(2) <= end_sec;
                    clamped = included & eventTable.IsClamp;
                    unclamped = included & ~eventTable.IsClamp;
                    if ~any(clamped) || ~any(unclamped); continue; end
                    unclampRMS = mean(eventTable.RMS(unclamped));
                    clampRMS = max(eventTable.RMS(clamped));
                    nWater = sum(included & eventTable.Type == "Water");
                    nAirpuff = sum(included & eventTable.Type == "Airpuff");
                    nTone = sum(included & eventTable.Type == "Tone");
                    nWaterClamp = sum(clamped & eventTable.Type == "Water");
                    nWaterUnclamp = sum(unclamped & eventTable.Type == "Water");
                    nAirpuffClamp = sum(clamped & eventTable.Type == "Airpuff");
                    nAirpuffUnclamp = sum(unclamped & eventTable.Type == "Airpuff");
                    score = unclampRMS-clampWeight*clampRMS;
                    candidates(end+1,:) = [start_sec,end_sec,duration_sec,score,unclampRMS, ...
                        clampRMS,sum(clamped),sum(unclamped),nWater,nAirpuff,nTone,daIdx, ...
                        nWaterClamp,nWaterUnclamp,nAirpuffClamp,nAirpuffUnclamp]; %#ok<SAGROW>
                end
            end
            eventScores{daIdx} = eventTable;
        end
        windowCandidates = array2table(candidates,'VariableNames',{'Start_sec','End_sec','Duration_sec', ...
            'Score','MeanUnclampRMS','WorstClampRMS','NClamp','NUnclamp','NWater','NAirpuff','NTone','ChannelIndex', ...
            'NWaterClamp','NWaterUnclamp','NAirpuffClamp','NAirpuffUnclamp'});
        windowCandidates.ChannelName = reshape(string({timeSeries(windowCandidates.ChannelIndex).name}),[],1);
        if isempty(windowCandidates)
            disp('No qualifying window exists: no recorded photometry channel has a candidate of the requested duration with complete clamped and unclamped events.');
            return
        end
        % Trial composition comes first, regardless of the response score.
        combinationCounts = windowCandidates{:,{'NWaterClamp','NWaterUnclamp','NAirpuffClamp','NAirpuffUnclamp'}};
        windowCandidates.NCombinations = sum(combinationCounts > 0,2);
        windowCandidates.NMatchedTypes = double(all(combinationCounts(:,1:2) > 0,2)) + ...
            double(all(combinationCounts(:,3:4) > 0,2));
        windowCandidates.HasBothVersions = windowCandidates.NCombinations == 4;
        windowCandidates = sortrows(windowCandidates, ...
            {'NCombinations','NMatchedTypes','Score','Duration_sec','Start_sec'}, ...
            {'descend','descend','descend','ascend','ascend'});

        % Ensure the entire displayed DA / clamp interval is recorded and finite.
        selectedWindowIdx = [];
        for i = 1:height(windowCandidates)
            daIdx = windowCandidates.ChannelIndex(i);
            sampleStart_sec = windowCandidates.Start_sec(i);
            sampleEnd_sec = windowCandidates.End_sec(i);
            [~,referenceIdx] = min(abs(timeNI-(timeNI(1)+sampleStart_sec)));
            timeRange = [sampleStart_sec,sampleEnd_sec]-(timeNI(referenceIdx)-timeNI(1));
            try
                candidateTrace = getSampleClampTrace(timeSeries,params,referenceIdx,timeRange, ...
                    daIdx=daIdx,removeClampBaseline=false);
            catch ME
                if strcmp(ME.identifier,'getSampleClampTrace:WindowOutsideRecording'); continue; end
                rethrow(ME)
            end
            if any(~isfinite([candidateTrace.DA,candidateTrace.blueClamp,candidateTrace.redClamp])); continue; end
            selectedWindowIdx = i;
            break
        end
        if isempty(selectedWindowIdx)
            disp('No qualifying window exists with continuous, finite DA and clamp recordings.');
            return
        end
        selectedWindow = windowCandidates(selectedWindowIdx,:);
        eventTable = eventScores{daIdx};
        selectedEvents = eventTable(eventTable.Time_sec >= sampleStart_sec & eventTable.Time_sec <= sampleEnd_sec,:);
        selectedEvents.InScore = selectedEvents.Valid & selectedEvents.Time_sec+extractRange(1) >= sampleStart_sec & ...
            selectedEvents.Time_sec+extractRange(2) <= sampleEnd_sec;
        trialRows = [];
        disp('Selected random-session window:'); disp(selectedWindow);
        if ~selectedWindow.HasBothVersions
            fprintf('No valid window contains all four combinations. Selected the best available coverage: %d/4 combinations.\n',selectedWindow.NCombinations);
        end
        disp(selectedEvents);
    end
end

%% Plot the trace

% Same extraction helpers as session analysis; retain each signal's stored Fs
removeClampBaseline = timeOption == 2 && ~findBestClampWindow; % manual trials only
sampleClampTrace = getSampleClampTrace(timeSeries,params,referenceIdx,timeRange, ...
                        daIdx=daIdx,removeClampBaseline=removeClampBaseline);
sampleClampTrace.sessionpath = sessionpath;
sampleClampTrace.sessionName = sessionName;
sampleClampTrace.trialRow = trialRows;
sampleClampTrace.trialNumber = [];
if ~isempty(trialRows)
    sampleClampTrace.trial = trials(trialRows,:);
    sampleClampTrace.trialNumber = trials.TrialNumber(trialRows)';
end
if timeOption == 2 && findBestClampWindow
    sampleClampTrace.windowCandidates = windowCandidates;
    sampleClampTrace.photometryIdx = photometryIdx;
    sampleClampTrace.selectedChannelIdx = daIdx;
    sampleClampTrace.eventScores = eventScores;
    sampleClampTrace.selectedWindow = selectedWindow;
    sampleClampTrace.events = selectedEvents;
    sampleClampTrace.baselineWindow_sec = baselineWindow_sec;
    sampleClampTrace.responseWindow_sec = responseWindow_sec;
    sampleClampTrace.clampWeight = clampWeight;
    sampleClampTrace.rankingOrder = {'NCombinations','NMatchedTypes','Score'};
end

% Plot against session time so the interval between trials remains visible
DA = sampleClampTrace.DA;
clampIndex = sampleClampTrace.clampIndex;
timeDA = sampleClampTrace.sessionTimeDA_sec;
timeClamp = sampleClampTrace.sessionTimeClamp_sec;
clampColor = [.232 .76 .58];

initializeFig(0.65,0.65,ScreenSize=get(groot,'ScreenSize'),Visible='on'); tiledlayout(2,1);
axDA = nexttile;
plotTraces(DA,timeDA,color=clampColor,plotPatch=false);
title(sampleClampTrace.DAName,'Interpreter','none');
xlabel('Time (s)'); ylabel(sampleClampTrace.DAUnits);

axClamp = nexttile;
plotTraces(clampIndex,timeClamp,color=[.3 .3 .3],plotPatch=false);
if removeClampBaseline; title('Clamp index (red - blue, baseline subtracted)');
else; title('Clamp index (red - blue)'); end
xlabel('Time (s)'); ylabel('Clamp index (%)');

for i = 1:length(trialRows)
    cueTime_sec = timeNI(trials.CueTime(trialRows(i)))-timeNI(1);
    trialLabel = sprintf('Trial %g',trials.TrialNumber(trialRows(i)));
    xline(axDA,cueTime_sec,':',trialLabel,'HandleVisibility','off');
    xline(axClamp,cueTime_sec,':',trialLabel,'HandleVisibility','off');
end
if timeOption == 2 && findBestClampWindow
    for i = 1:height(selectedEvents)
        if selectedEvents.IsClamp(i); state = 'clamp'; eventColor = clampColor;
        else; state = 'unclamp'; eventColor = [.6 .3 .1]; end
        eventLabel = sprintf('%s (%s)',selectedEvents.Type(i),state);
        if ~selectedEvents.InScore(i); eventLabel = sprintf('%s *',eventLabel); end
        xline(axDA,selectedEvents.Time_sec(i),':',eventLabel,'Color',eventColor,'HandleVisibility','off');
        xline(axClamp,selectedEvents.Time_sec(i),':',eventLabel,'Color',eventColor,'HandleVisibility','off');
    end
    % * marks events displayed but excluded from scoring (partial or missing windows).
end
linkaxes([axDA,axClamp],'x'); xlim([sampleStart_sec,sampleEnd_sec]);
sgtitle(sprintf('%s | %.3f to %.3f s',sessionName,sampleStart_sec,sampleEnd_sec),'Interpreter','none');
sampleClampTrace.figure = gcf;
disp('Finished: example DA and clampIndex traces plotted');

return
