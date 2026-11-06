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

%% (Optional) Define time to plot (opt 1: just enter in sec)


sampleStart_sec = 100;
sampleEnd_sec = 165;

% Convert session seconds to the native NI reference sample
validateattributes([sampleStart_sec,sampleEnd_sec],{'numeric'},{'real','finite','nonnegative'});
if sampleEnd_sec <= sampleStart_sec || sampleEnd_sec > timeNI(end)-timeNI(1)
    error('Choose an increasing time window within the session.');
end
[~,referenceIdx] = min(abs(timeNI-(timeNI(1)+sampleStart_sec)));
timeRange = [sampleStart_sec,sampleEnd_sec] - (timeNI(referenceIdx)-timeNI(1));
trialRows = [];


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

    % Find a continuous example with flat clamped and strong unclamped responses,
    % favouring windows whose clamped water and airpuff responses are the flattest.
    if findBestClampWindow
        windowDurations_sec = 60; % default: a 60 s example; optionally use 30:5:60
        windowStep_sec = 1;
        baselineWindow_sec = [-2,0]; % mean baseline; excludes the endpoint
        responseWindow_sec = [0,2]; % relative to actual water / airpuff / tone delivery
        clampWeight = 2;
        waterClampWeight = 4; % extra penalty on the worst clamped water response
        airpuffClampWeight = 4; % extra penalty on the worst clamped airpuff response

        validateattributes(windowDurations_sec,{'numeric'},{'vector','nonempty','finite','>=',30,'<=',60});
        validateattributes(windowStep_sec,{'numeric'},{'scalar','finite','positive'});
        validateattributes(baselineWindow_sec,{'numeric'},{'size',[1,2],'finite','real'});
        validateattributes(responseWindow_sec,{'numeric'},{'size',[1,2],'finite','real'});
        validateattributes(clampWeight,{'numeric'},{'scalar','finite','nonnegative'});
        validateattributes(waterClampWeight,{'numeric'},{'scalar','finite','nonnegative'});
        validateattributes(airpuffClampWeight,{'numeric'},{'scalar','finite','nonnegative'});
        if diff(baselineWindow_sec) <= 0 || baselineWindow_sec(2) > 0 || diff(responseWindow_sec) <= 0
            error('Use increasing windows, with the baseline ending at or before event delivery.');
        end

        % Reuse the event classification saved by analyzeSessions_clamp.
        % Its clampON_event marks clamp activity within 3 s after each delivery.
        eventTable = getClampEventTable(behavior,timeNI);
        if isempty(eventTable)
            disp('No qualifying window exists: no water, airpuff or tone events were found.');
            return
        end

        candidates = zeros(0,18);
        eventScores = cell(size(timeSeries));
        for daIdx = reshape(photometryIdx,1,[])
            % RMS is direction-independent: peaks and dips both count as responses.
            [eventTable,extractRange] = scoreClampEventRMS(eventTable,timeSeries(daIdx),params, ...
                baselineWindow_sec=baselineWindow_sec,responseWindow_sec=responseWindow_sec);

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
                    % Score the clamped water and airpuff traces separately so a window
                    % that is flat for both outranks one that is flat only for tone or
                    % for a single event type. 0 when the window holds no clamped event
                    % of that type, which the composition ranking below deprioritises.
                    waterClamped = clamped & eventTable.Type == "Water";
                    if any(waterClamped); waterClampRMS = max(eventTable.RMS(waterClamped));
                    else; waterClampRMS = 0; end
                    airpuffClamped = clamped & eventTable.Type == "Airpuff";
                    if any(airpuffClamped); airpuffClampRMS = max(eventTable.RMS(airpuffClamped));
                    else; airpuffClampRMS = 0; end
                    nWater = sum(included & eventTable.Type == "Water");
                    nAirpuff = sum(included & eventTable.Type == "Airpuff");
                    nTone = sum(included & eventTable.Type == "Tone");
                    nWaterClamp = sum(clamped & eventTable.Type == "Water");
                    nWaterUnclamp = sum(unclamped & eventTable.Type == "Water");
                    nAirpuffClamp = sum(clamped & eventTable.Type == "Airpuff");
                    nAirpuffUnclamp = sum(unclamped & eventTable.Type == "Airpuff");
                    score = unclampRMS-clampWeight*clampRMS- ...
                        waterClampWeight*waterClampRMS-airpuffClampWeight*airpuffClampRMS;
                    candidates(end+1,:) = [start_sec,end_sec,duration_sec,score,unclampRMS, ...
                        clampRMS,waterClampRMS,airpuffClampRMS,sum(clamped),sum(unclamped),nWater,nAirpuff,nTone,daIdx, ...
                        nWaterClamp,nWaterUnclamp,nAirpuffClamp,nAirpuffUnclamp]; %#ok<SAGROW>
                end
            end
            eventScores{daIdx} = eventTable;
        end
        windowCandidates = array2table(candidates,'VariableNames',{'Start_sec','End_sec','Duration_sec', ...
            'Score','MeanUnclampRMS','WorstClampRMS','WorstWaterClampRMS','WorstAirpuffClampRMS', ...
            'NClamp','NUnclamp','NWater','NAirpuff','NTone','ChannelIndex', ...
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

        % Screening step: keep the best non-overlapping candidates whose DA and
        % clamp commands are recorded in full, then show DA alone for each one.
        nTopCandidates = 5;
        windowCandidates.Candidate = (1:height(windowCandidates))';
        [topCandidates,topCandidateTraces] = selectTopClampWindows(windowCandidates, ...
            timeSeries,params,timeNI,nTop=nTopCandidates);
        if isempty(topCandidates)
            disp('No qualifying window exists with continuous, finite DA and clamp recordings.');
            return
        end
        topCandidates.Label = string(compose('%d/4 combos | clamp RMS %.3f', ...
            topCandidates.NCombinations,topCandidates.WorstClampRMS));
        candidateClampEpochs = getClampEpochs(timeSeries);
        % Type / IsClamp / Time_sec are the same on every channel, so any copy labels
        labelEvents = eventTable(:,{'Sample','Type','IsClamp','Time_sec'});
        candidateFig = plotClampWindowGrid(topCandidates,topCandidateTraces, ...
            events=labelEvents,clampEpochs=candidateClampEpochs,channelIdx=photometryIdx, ...
            titleText=sprintf('%s | best trial-type candidates (DA only)',sessionName));
        disp('Top trial-type candidates:');
        disp(topCandidates(:,{'ChannelName','Rank','Start_sec','End_sec','NCombinations', ...
            'NMatchedTypes','Score','WorstClampRMS','WorstWaterClampRMS','WorstAirpuffClampRMS'}));

        % The best candidate overall becomes the window the detail cell plots
        [~,bestRow] = min(topCandidates.Candidate);
        selectedWindow = topCandidates(bestRow,:);
        daIdx = selectedWindow.ChannelIndex;
        sampleStart_sec = selectedWindow.Start_sec;
        sampleEnd_sec = selectedWindow.End_sec;
        referenceIdx = selectedWindow.ReferenceIdx;
        timeRange = selectedWindow.TimeRange;
        eventTable = eventScores{daIdx};
        selectedEvents = eventTable(eventTable.Time_sec >= sampleStart_sec & ...
            eventTable.Time_sec <= sampleEnd_sec,:);
        selectedEvents.InScore = selectedEvents.Valid & ...
            selectedEvents.Time_sec+extractRange(1) >= sampleStart_sec & ...
            selectedEvents.Time_sec+extractRange(2) <= sampleEnd_sec;
        trialRows = [];
        disp('Selected random-session window:'); disp(selectedWindow);
        if ~selectedWindow.HasBothVersions
            fprintf('No valid window contains all four combinations. Selected the best available coverage: %d/4 combinations.\n',selectedWindow.NCombinations);
        end
        disp(selectedEvents);

        fprintf(['Current window: %s %g-%g s. Run the detail cell for DA + clamp index, or copy ' ...
                 'another window into the opt 1 cell first.\n'], ...
            selectedWindow.ChannelName,sampleStart_sec,sampleEnd_sec);

        clampWindowSearch.sessionName = sessionName;
        clampWindowSearch.windowCandidates = windowCandidates;
        clampWindowSearch.topCandidates = topCandidates;
        clampWindowSearch.traces = topCandidateTraces;
        clampWindowSearch.eventScores = eventScores;
        clampWindowSearch.selectedWindow = selectedWindow;
        clampWindowSearch.events = selectedEvents;
        clampWindowSearch.photometryIdx = photometryIdx;
        clampWindowSearch.rankingOrder = {'NCombinations','NMatchedTypes','Score'};
        clampWindowSearch.settings = struct('windowDurations_sec',windowDurations_sec, ...
            'windowStep_sec',windowStep_sec,'baselineWindow_sec',baselineWindow_sec, ...
            'responseWindow_sec',responseWindow_sec,'clampWeight',clampWeight, ...
            'waterClampWeight',waterClampWeight,'airpuffClampWeight',airpuffClampWeight, ...
            'nTopCandidates',nTopCandidates);
        clampWindowSearch.figure = candidateFig;
    end
end

%% Plot the detail trace: DA and clamp index for the window set above

% The detail step of the screening flow: pick a window from a screening figure
% (or type one into the opt 1 cell), then run this cell to see the clamp index
% under DA. It uses whatever referenceIdx / timeRange is currently set, so it
% works after opt 1, after the trial-type search and after the duration search.

removeClampBaseline = false; % true subtracts each clamp command's [-0.2,0) s median
detailDaIdx = daIdx;         % photometry channel to show

sampleClampTrace = getSampleClampTrace(timeSeries,params,referenceIdx,timeRange, ...
                        daIdx=detailDaIdx,removeClampBaseline=removeClampBaseline);
sampleClampTrace.sessionpath = sessionpath;
sampleClampTrace.sessionName = sessionName;
sampleClampTrace.trialRow = trialRows;
sampleClampTrace.trialNumber = [];
if ~isempty(trialRows)
    sampleClampTrace.trial = trials(trialRows,:);
    sampleClampTrace.trialNumber = trials.TrialNumber(trialRows)';
end

% Deliveries inside the window, rebuilt here so the cell stands on its own
detailEvents = table();
try
    detailEvents = getClampEventTable(behavior,timeNI);
    detailEvents = detailEvents(detailEvents.Time_sec >= sampleStart_sec & ...
                                detailEvents.Time_sec <= sampleEnd_sec,:);
catch ME
    if ~strcmp(ME.identifier,'getClampEventTable:MissingClassification'); rethrow(ME); end
end
sampleClampTrace.events = detailEvents;
detailClampEpochs = getClampEpochs(timeSeries);

% Plot against session time so the interval between deliveries remains visible.
% plotSampleClampTrace marks each delivery with plotEvent (airpuff patches span
% 0.2 s) and colors clamped and unclamped deliveries differently.
[axDA,axClamp,sampleFig] = plotSampleClampTrace(sampleClampTrace, ...
    events=detailEvents,clampEpochs=detailClampEpochs, ...
    trials=trials,trialRows=trialRows,timeNI=timeNI, ...
    xlim=[sampleStart_sec,sampleEnd_sec],airpuffDuration=0.2, ...
    titleText=sprintf('%s | %.3f to %.3f s',sessionName,sampleStart_sec,sampleEnd_sec));
sampleClampTrace.figure = sampleFig;
sampleClampTrace.axes = [axDA,axClamp];
disp('Finished: example DA and clampIndex traces plotted');

%% Find the flattest 30 s windows with equal clamped / unclamped duration (DA only)

% Independent of the cells above: marks when the clamp command is on, scans
% every 30 s window of every photometry channel and keeps the non-overlapping
% windows that spend half their time clamped and whose DA is flattest while
% clamped. Deliveries are drawn but never required, so a window needs no
% particular mix of water / airpuff / tone. Clamp index is not plotted here:
% each tile shows DA alone, one row per photometry channel.

flatWindow_sec = 30;
flatStep_sec = 1;
clampOnThreshold_pct = 1;   % command percent counted as clamping, as in analyzeSessions_clamp
clampGapFill_sec = 0.1;     % fill shorter OFF gaps so one clamp epoch stays one epoch
maxDurationImbalance = 0.1; % |clamped fraction - 0.5| allowed, i.e. 40-60% clamped
minSideSamples = 10;        % DA samples required on each side before scoring a window
nTopWindows = 5;
airpuffDuration_sec = 0.2;

validateattributes(flatWindow_sec,{'numeric'},{'scalar','finite','positive'});
validateattributes(flatStep_sec,{'numeric'},{'scalar','finite','positive'});
validateattributes(clampOnThreshold_pct,{'numeric'},{'scalar','finite','positive'});
validateattributes(clampGapFill_sec,{'numeric'},{'scalar','finite','nonnegative'});
validateattributes(maxDurationImbalance,{'numeric'},{'scalar','finite','nonnegative','<=',0.5});
validateattributes(minSideSamples,{'numeric'},{'scalar','finite','positive','integer'});
validateattributes(nTopWindows,{'numeric'},{'scalar','finite','positive','integer'});

[clampEpochs,clampOn,clampTime_sec] = getClampEpochs(timeSeries, ...
    threshold_pct=clampOnThreshold_pct,gapFill_sec=clampGapFill_sec);
fprintf('Clamp on for %.1f%% of the session, in %d epochs.\n', ...
    100*mean(clampOn),size(clampEpochs,1));

% Deliveries are labels only here; they do not gate the search
flatEvents = getClampEventTable(behavior,timeNI);

% Score every window of every channel: half clamped, and flat while clamped
sessionDuration_sec = timeNI(end)-timeNI(1);
flatWindows = table();
for ch = reshape(photometryIdx,1,[])
    da = double(timeSeries(ch).data(:));
    daFs = double(timeSeries(ch).finalFs);
    daTime_sec = (0:numel(da)-1)'/daFs + double(timeSeries(ch).time_offset);
    % Clamp state on the DA grid, so fractions match the plotted samples
    daClamped = interp1(clampTime_sec,double(clampOn),daTime_sec,'nearest','extrap') > 0.5;
    rows = zeros(0,7);
    for start_sec = 0:flatStep_sec:(sessionDuration_sec-flatWindow_sec)
        end_sec = start_sec+flatWindow_sec;
        first = max(1,ceil((start_sec-daTime_sec(1))*daFs)+1);
        last = min(numel(da),floor((end_sec-daTime_sec(1))*daFs)+1);
        if last-first+1 < 2*minSideSamples; continue; end
        window = da(first:last);
        clampedHere = daClamped(first:last);
        if any(~isfinite(window)); continue; end
        if sum(clampedHere) < minSideSamples || sum(~clampedHere) < minSideSamples; continue; end
        clampFraction = mean(clampedHere);
        if abs(clampFraction-0.5) > maxDurationImbalance; continue; end
        inWindow = flatEvents.Time_sec >= start_sec & flatEvents.Time_sec <= end_sec;
        rows(end+1,:) = [start_sec,end_sec,clampFraction,std(window(clampedHere)), ...
            std(window(~clampedHere)),sum(inWindow & flatEvents.IsClamp), ...
            sum(inWindow & ~flatEvents.IsClamp)]; %#ok<SAGROW>
    end
    if isempty(rows); continue; end
    channelWindows = array2table(rows,'VariableNames',{'Start_sec','End_sec','ClampFraction', ...
        'ClampStd','UnclampStd','NClampDelivery','NUnclampDelivery'});
    channelWindows.ChannelIndex = repmat(ch,height(channelWindows),1);
    flatWindows = [flatWindows;channelWindows]; %#ok<AGROW>
end
if isempty(flatWindows)
    fprintf(['No %g s window is %g-%g%% clamped. Raise maxDurationImbalance, or check that ' ...
             'the clamp epochs are long enough for this window.\n'], ...
        flatWindow_sec,100*(0.5-maxDurationImbalance),100*(0.5+maxDurationImbalance));
    return
end
flatWindows.ChannelName = reshape(string({timeSeries(flatWindows.ChannelIndex).name}),[],1);

% Flattest clamped DA first; neighbouring starts overlap heavily, so keep only
% the best non-overlapping windows that getSampleClampTrace can extract in full.
flatWindows = sortrows(flatWindows,{'ClampStd','Start_sec'},{'ascend','ascend'});
[topWindows,topTraces] = selectTopClampWindows(flatWindows,timeSeries,params,timeNI, ...
    nTop=nTopWindows);
if isempty(topWindows)
    disp('No balanced window has continuous, finite DA.');
    return
end
disp('Flattest windows with equal clamped / unclamped duration:');
disp(topWindows(:,{'ChannelName','Rank','Start_sec','End_sec','ClampFraction','ClampStd', ...
    'UnclampStd','NClampDelivery','NUnclampDelivery'}));

% One row per photometry channel, its top windows across the row
topWindows.Label = string(compose('%.0f%% clamped | std %.3f', ...
    100*topWindows.ClampFraction,topWindows.ClampStd));
flatFig = plotClampWindowGrid(topWindows,topTraces,events=flatEvents, ...
    clampEpochs=clampEpochs,channelIdx=photometryIdx,airpuffDuration=airpuffDuration_sec, ...
    titleText=sprintf('%s | flattest %g s windows, half clamped by duration (DA only)', ...
        sessionName,flatWindow_sec));

% The flattest window becomes the current one, so the detail cell can plot it
[~,bestFlatRow] = min(topWindows.ClampStd);
daIdx = topWindows.ChannelIndex(bestFlatRow);
sampleStart_sec = topWindows.Start_sec(bestFlatRow);
sampleEnd_sec = topWindows.End_sec(bestFlatRow);
referenceIdx = topWindows.ReferenceIdx(bestFlatRow);
timeRange = topWindows.TimeRange(bestFlatRow,:);
trialRows = [];
fprintf(['Current window: %s %g-%g s. Run the detail cell for DA + clamp index, or copy ' ...
         'another window into the opt 1 cell first.\n'], ...
    topWindows.ChannelName(bestFlatRow),sampleStart_sec,sampleEnd_sec);

flatClampWindows.sessionName = sessionName;
flatClampWindows.sessionpath = sessionpath;
flatClampWindows.windows = flatWindows;
flatClampWindows.top = topWindows;
flatClampWindows.traces = topTraces;
flatClampWindows.events = flatEvents;
flatClampWindows.clampOn = clampOn;
flatClampWindows.clampEpochs = clampEpochs;
flatClampWindows.clampTime_sec = clampTime_sec;
flatClampWindows.settings = struct('window_sec',flatWindow_sec,'step_sec',flatStep_sec, ...
    'clampOnThreshold_pct',clampOnThreshold_pct,'clampGapFill_sec',clampGapFill_sec, ...
    'maxDurationImbalance',maxDurationImbalance,'minSideSamples',minSideSamples, ...
    'nTopWindows',nTopWindows);
flatClampWindows.figure = flatFig;
disp('Finished: flattest clamped DA windows plotted');
