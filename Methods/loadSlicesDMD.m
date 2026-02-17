function loadSlicesDMD(epochs,options)

% Create spots.mat for each epoch (i.e. each search)

arguments
    epochs table %full path of the session or epochs.mat

    options.filterSignal logical = true
    options.filterSweeps logical = true % If true, do not analyze sweeps with different Vhold (included == false)

    options.save logical = true
    options.reload logical = false
    options.reloadCells logical = false
    options.reloadCellAnalysis logical = false
    options.calculateQC logical = true
    options.reconstructDMD logical = true % reconstruct DMD if fullSearchTable is not found

    options.outputFs double = 10000
    options.timeRange double = [-20,100] % in ms
    options.analysisWindowLength double = 30 % in ms after stim onset
    options.controlWindowLength double = 30 % in ms before stim onset
    options.eventSample % in sample
    options.nArtifactSamples double = 0 % in sample
    options.rcCheckRecoveryWindow double = 100 % in ms
    options.peakWindow double = 1 % in ms around the peak to average

    options.vholdChannel string = 'AD1'

    options.feature = 'auc'
    options.thresholdFactor double = 3 % 3*std

    options.rawDataPath string
    options.saveDataPath string = 'default'
end

%% General setup

today = char(datetime('today','Format','yyyyMMdd'));
sessionPath = osPathSwitch(epochs{1,"Session"});
dirsplit = split(sessionPath,filesep); expName = dirsplit{end};

% Decide reload if session has already been loaded
if options.reload
    options.reloadCells = true;
    options.reloadCellAnalysis = true;
    resultsFolderName = ['Results-',today];
else
    resultsList = sortrows(struct2cell(dir(fullfile(sessionPath,'Results-*')))',3);
    if isempty(resultsList)
        options.reload = true;
        options.reloadCells = true;
        options.reloadCellAnalysis = true;
        resultsFolderName = ['Results-',today];
    else
        resultsPath = dir(fullfile(resultsList{end,2},resultsList{end,1},'cell*','spots_*.mat'));
        if ~isempty(resultsPath)
            disp('Loading stop: spots file found.');
            resultsFolderName = resultsList{end,1};
            if ~options.reloadCells && ~options.reloadCellAnalysis
                return; 
            elseif options.reloadCells
                options.reloadCellAnalysis = true;
            elseif ~options.reloadCells && options.reloadCellAnalysis
                try load(fullfile(resultsList{end,2},resultsList{end,1},['cells_DMD_',expName,'.mat']));
                catch
                    error('Error: did not find cells_DMD.mat!');
                end
                disp(['Loaded: cells.mat in ',resultsList{end,1}]);
            end
        end
    end
end

% Determine data path
% 1. If saveDataPath == 'default', save files in saveDataPath = rawDataPath
% 2. If saveDataPath ~= 'default', saveDataPath should be specified by the user
% 3. "Sessions" in epochs and cells table should be the same as rawDataPath
if strcmp(options.saveDataPath,'default')
    options.rawDataPath = sessionPath;
    options.saveDataPath = strcat(options.rawDataPath,filesep,resultsFolderName);
    if ~exist(options.saveDataPath, 'dir')
        mkdir(options.saveDataPath);
    end
else
    options.rawDataPath = sessionPath;
    if ~exist(options.saveDataPath, 'dir')
        mkdir(options.saveDataPath);
    end
end
options.saveDataPath = string(options.saveDataPath);

% Turn off some warnings
warning('off','MATLAB:unknownObjectNowStruct');
warning('off','MATLAB:table:RowsAddedExistingVars');

%% Select randomSearch epochs
% randomSearchIdx = cellfun(@(x) contains(x.cycle, 'randomSearch'), epochs.("Protocol"));
% randomSearchIdx = cellfun(@(x) iscell(x), epochs.("Protocol"));
randomSearchIdx = cellfun(@(x) iscell(x) && contains(x{1}.cycle, 'randomSearch', IgnoreCase=true), epochs.("Protocol"));
exp = sortrows(epochs(randomSearchIdx,:),'Cell'); % randomSearchEpochs

%% Build spots.mat and cells.mat
if options.reloadCells
    %% Initialize cells_DMD.mat
    % structure for responseMap
    % For each row, response map is a cell with #searches element
    % within each element, theres a 3-dim matrix (xRange,yRange,depth)
    
    varTypes = {'string','string','string','double','cell','cell','cell',...
                'cell','cell','cell','cell'};
    varNames = {'Session','Animal','Task','Cell','Epochs','Vhold','Protocol',...
                'Response map','Difference map','Stats','Options'};
    cells = table('Size',[max(exp{:,'Cell'}),length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    
    % Initialize cell-level params
    cellResponseMap = {}; isResponseMap_cell = {}; cellHotspotMap = {};
    cellCurrentMap = {}; cellBaselineMap = {};
    cellVhold = []; cellEpochs = {}; searchTotalSpots = {}; 
    cellDepthList = {}; cellProtocols = {}; 
    nullSpotData = []; preStimData = []; baselineData = [];
    cellSpotResponse = {}; cellSpotSequence = {};
    cellMaxResponse = {}; cellMinResponse = {};
    cellMaxTime = {};cellMinTime = {};
    cellAUC = {}; cellEIindex = {};
    cellBaselineAUC = {}; cellBaselineSTD = {};
    cellSpotLocation = {};

    %% Iterate & analyze individual epoch
    for row = 1:size(exp,1)
        clearvars AD*
        cellResultsPath = string(strcat(options.saveDataPath,filesep,['cell',num2str(exp{row,'Cell'})]));
    
        if options.reload
            %% Load basic info
        
            % Initialize fullSearchTable pointer
            % Since I'm iterating sweep in order, the spots are correspond to the
            % order recorded in fullSearchTable
            curFSRow = 0;   % pointer into fullSearchTable rows
            curSpot  = 0;   % number of stored rows in "spots" (main search only)
        
            % Load fullSearchTable
            epoch = exp{row,'Epoch'};
            sweepAcq = exp{row,'Sweep names'}{1};
            epochPath = osPathSwitch(exp{row,'Options'}{1}.rawDataPath);
            
            % Intialize potentia reconstruction params
            reconstructedSearch = false;
            
            try 
                load(fullfile(epochPath,strcat('Epoch',num2str(epoch),'_fullSearchTable.mat')),'fullSearchTable');
            catch
                fullSearchTable = table();
                warning(strcat("Epoch ",num2str(epoch)," (cell ",num2str(exp{row,'Cell'}),") is missing fullSearchTable. ", ...
                    "Attempting reconstruction from *_responseMap.fig..."));
            
                if options.reconstructDMD
                    [reconOK, reconNoise] = reconstructFromResponseMap(epochPath, cellResultsPath, exp(row,:), epoch, options);
                    if reconOK
                        reconstructedSearch = true;
                        % Append to noise buffers so noise_cellX.mat still builds
                        nullSpotData = [nullSpotData, reconNoise.nullSpotData];
                        preStimData  = [preStimData,  reconNoise.preStimData];
                        baselineData = [baselineData, reconNoise.baselineData];
                    else
                        warning(strcat("Reconstruction failed for Epoch ",num2str(epoch),". Skipping this epoch."));
                        continue
                    end
                else
                    warning(strcat("Reconstruction disabled. Skipping Epoch ",num2str(epoch),"."));
                    continue
                end
            end
        
            % Read patching info csv file
            csvOpts = detectImportOptions(fullfile(epochPath,'InfoPatching.xlsx'));
            csvOpts.SelectedVariableNames = 1:9;
            csvOpts.VariableNamesRange = 25;
            patchInfo = readtable(fullfile(epochPath,'InfoPatching.xlsx'),csvOpts);
            patchInfo = rmmissing(patchInfo,DataVariables="acq_");
            patchInfo = rmmissing(patchInfo,DataVariables="epoch");
            
            if iscell(patchInfo.acq_);     patchInfo.acq_     = str2double(patchInfo.acq_); end
            if iscell(patchInfo.epoch);    patchInfo.epoch    = str2double(patchInfo.epoch); end
            if iscell(patchInfo.cyclePos); patchInfo.cyclePos = str2double(patchInfo.cyclePos); end
            if iscell(patchInfo.holding);  patchInfo.holding  = str2double(patchInfo.holding); end

        
            %% Initialize spots table
            
            varTypes = {'string','string','string','double','double','double',...
                        'double','double','string','cell','cell',...
                        'cell','cell',...
                        'cell',...
                        'cell'};
            varNames = {'Session','Animal','Task','Epoch','Cell','Vhold',...
                        'Depth','Repetition','Sweep','Location','Protocol',...
                        'Response','Stats',...
                        'QC',...
                        'Options'};
            spots = table('Size',[height(fullSearchTable),length(varNames)],...
                'VariableTypes',varTypes,'VariableNames',varNames);

            hotOpts = struct( ...
                        'AllowCrossDepth', true, ...
                        'PreferDepthMatch', true, ...
                        'PreferDepthKind', "any", ...      % or "finalDepth"
                        'PulseTolRel', 0.01, ...
                        'PulseTolAbs', NaN, ...
                        'Debug', false);

            % Some sweeps at the end of an epoch (e.g., hotspot validation sweeps)
            % are not represented in Epoch*_fullSearchTable.mat. 
            % We collect them separately using
            % Epoch*_finalHotspots_Depth*.mat or Epoch*_maxSearchHotspots_Depth*.mat
            % and later save them as separate "searches" so analyzeDMDSearch can plot them.
            hotspotSpots = table('Size',[0,length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
            [sweepHotspots, sweepDepths] = localLoadHotspots(epochPath, epoch);
        
            %% Load individual sweeps
            for k = 1:length(sweepAcq)
                % Load sweep traces (.data)
                disp(['Loading ',sweepAcq{k},'.mat for epoch ',num2str(epoch),' of cell ',num2str(exp{row,"Cell"})]);
                try load(fullfile(epochPath,strcat(sweepAcq{k},'.mat'))); 
                catch
                    warning(strcat("Sweep ",num2str(sweepAcq{k}), " not saved, skipping this sweep!"));
                    continue
                end
                
                % Extract experiment protocol from header string
                headerString = eval([sweepAcq{k},'.UserData.headerString']);
                protocol = getCellProtocol(headerString,...
                                           outputFs=options.outputFs,...
                                           rcCheckRecoveryWindow=options.rcCheckRecoveryWindow);
                if k == 1; prevDepth = protocol.depth; end
                
                % ---- Decide which spot geometry to use for this sweep ----
                [sweepSpots, protocol, nPulsesThisSweep, isHotspotSweep, curFSRow] = ...
                            selectSweepSpots(headerString, protocol, fullSearchTable, curFSRow, reconstructedSearch, ...
                                             sweepHotspots, sweepDepths, hotOpts, sweepAcq{k});
                if isempty(sweepSpots); continue; end

                % Extract raw trace
                raw_trace = eval([sweepAcq{k},'.data']);

                % Safety: make stimOnset consistent with what we analyze
                protocol.stimOnset = protocol.stimOnset(1:nPulsesThisSweep);

                % Define baseline window: have two windows, one before pulse and one after pulse
                rcCheckOnset = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='delay') * (options.outputFs/1000);
                rcCheckPulseWidth = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='pulseWidth') * (options.outputFs/1000);
                lastStim = protocol.stimOnset(nPulsesThisSweep);
                postStimStart = lastStim + round(0.200 * options.outputFs); % 200 ms after last pulse

                if rcCheckOnset < protocol.stimOnset(1)
                    rcCheckEnd = rcCheckOnset + rcCheckPulseWidth + (options.rcCheckRecoveryWindow*(options.outputFs/1000));
                    preStimWindow = rcCheckEnd : (protocol.stimOnset(1)-1);
                    postStimWindow = postStimStart:length(raw_trace);
                    baselineWindow = [preStimWindow,postStimWindow];
                else
                    preStimWindow = 1:(protocol.stimOnset(1)-1);
                    postStimWindow = postStimStart:rcCheckOnset;
                    baselineWindow = [preStimWindow,postStimWindow];
                end
        
                % Define time window for plotting
                options.eventSample = protocol.stimOnset;
                timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
                timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
                plotWindowLength = timeRangeEndSample(1) - timeRangeStartSample(1) + 1;
                plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
        
                % Define control window: 50ms before each spot stim onset
                controlWindowSamples = options.controlWindowLength * options.outputFs/1000;
    
                % Define analysis window: 50ms after each spot stim onset
                analysisWindowSamples = options.analysisWindowLength * options.outputFs/1000;
                eventSample = abs(options.outputFs*options.timeRange(1)/1000);
                analysisWindow = eventSample:eventSample+analysisWindowSamples;
        
                % Define time window for peak window analysis
                peakWindowWidth = (options.peakWindow/(2*1000)) * options.outputFs;
    
                % Save time windows to options
                options.baselineWindow = baselineWindow;
                options.baselineWindow_preStim = preStimWindow;
                options.baselineWindow_postStim = postStimWindow;
                options.analysisWindow = analysisWindow;
                options.plotWindowLength = plotWindowLength;
                options.peakWindowWidth = peakWindowWidth;
                options.plotWindowTime = plotWindowTime;
                options.analysisWindowSamples = analysisWindowSamples;
                options.controlWindowSamples = controlWindowSamples;
    
        
                % Extract quality metrics from header string
                qc = getCellQC(headerString,calculate=options.calculateQC,data=raw_trace);
                % Rs = qc.Rs; Rm = qc.Rm; Cm = qc.Cm;


                % Process trace: mean-subtracted, optional LP
                % Mean subtraction
                meanRawBaseline = mean(raw_trace(baselineWindow));
                mean_subtracted = raw_trace - meanRawBaseline;
                if options.filterSignal
                    Fs = options.outputFs; % Sampling frequency  
                    LP = lowpass(mean_subtracted',2000,Fs);
                    % Notch filter
                    d = designfilt('bandstopiir','FilterOrder',2, ...
                                   'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
                                   'DesignMethod','butter','SampleRate',Fs);
                    Notch = filtfilt(d,LP);   
                    % Smooth data using sgolay filter
                    yT = sgolayfilt(LP,5,27); % polynomial order of 5 and framelength of 27
                    y = yT';
                    % Median filter using 0.5ms window
                    y = movmedian(y,6,2);
                    % Subtract mean again (in SeulAh's code)
                    base2 = mean(y(:,baselineWindow),2); baseM2 = repmat(base2,1,size(y,2));
                    processed_trace = y - baseM2;
                else
                    processed_trace = mean_subtracted;
                end
        
                % Find baseline stats
                stats.baseline.avg = mean(processed_trace(baselineWindow));
                stats.baseline.std = std(processed_trace(baselineWindow));
        
                % Determine Vhold
                [vhold,~] = getSweepVhold(epochPath, sweepAcq{k}, meanRawBaseline, options.vholdChannel, infoTable=patchInfo);
        
                % Initialize matrix for noise analysis later
                baselineData = [baselineData, processed_trace(baselineWindow)];
        
                %% Loop through all spots
                for s = 1:nPulsesThisSweep
                    %% Define timeWindow and extract response
                    % Define control window
                    plotWindow = timeRangeStartSample(s):timeRangeEndSample(s);
                    controlWindow = protocol.stimOnset(s)-controlWindowSamples-1 : protocol.stimOnset(s)-1;
                    options.plotWindow = plotWindow;
                    options.controlWindow = controlWindow;
    
                    % Save spot coordinates
                    location = [sweepSpots{s,"xStart"}+1, sweepSpots{s,"xStart"}+sweepSpots{s,"xWidth"},...
                                sweepSpots{s,"yStart"}+1, sweepSpots{s,"yStart"}+sweepSpots{s,"yHeight"}];
        
                    % Save spot responses
                    responses.raw = raw_trace(plotWindow);
                    responses.processed = processed_trace(plotWindow);
                    responses.control = processed_trace(controlWindow);
                    responses.rawTrace = raw_trace;
                    responses.processedTrace = processed_trace;
        
                    %% Identify putative hotspots
    
                    % Seulah's method: 
                    % spots that exceed thresholdFactor scaled MAD from the median for at least 5ms
                    % trace_mad = mad(responses.processed,1);
                    % isResponse = find(sum(isoutlier(trace_mad,median=1,ThresholdFactor=options.thresholdFactor),2)>50); 
    
                    % Alternative method: more than thresholdFactor more than
                    % baseline std for at least 5ms
                    response_threshold =  options.thresholdFactor * stats.baseline.std;
                    responses.isResponse = sweepSpots{s,'response'};
                    if find(sum(abs(responses.processed(analysisWindow)) >= response_threshold)>=50)
                        responses.hotspot = true;
                    else; responses.hotspot = false; 
                    end
                    
                    %% Determine analysis time window per hotspot based on changepoint analysis
                    % if responses.hotspot
                    %     % Changepoint analysis
                    %     cpaWindow= 0; % option to widen the analysis window (unit in 0.1 ms), default is 0
                    %     changeIdx = findchangepts(responses.processed,Statistic='rms',MaxNumChanges=2); % rms works better than std or mean
                    %     % Final analysis time window
                    %     changeIdx(1) = changeIdx(1)-cpaWindow;   
                    %     changeIdx(2) = changeIdx(2)+cpaWindow;
                    % end
    
                    %% Extract baseline and prestim trace for noise analysis
                    % Concat baseline period + prestim period of all spots +
                    % response period for null spots
    
                    preStimData = [preStimData, processed_trace(controlWindow)];
                    if ~responses.hotspot; nullSpotData = [nullSpotData, responses.processed]; end
    
                    %% Calculate statistics
    
                    % Find auc
                    stats.response.auc = sum(responses.processed(analysisWindow)) / options.outputFs;
                    stats.baseline.auc = sum(processed_trace(controlWindow)) / options.outputFs;
        
                    % Find min and max value for stim response
                    trace = responses.processed(analysisWindow);
                    [~,maxIdx] = max(trace); [~,minIdx] = min(trace);
                    % Average around max/min idx to get final value
                    maxWindowStart = max(1,maxIdx-peakWindowWidth);
                    maxWindowEnd = min(maxIdx+peakWindowWidth,length(trace));
                    minWindowStart = max(1,minIdx-peakWindowWidth);
                    minWindowEnd = min(minIdx+peakWindowWidth,length(trace));
                    stats.response.max = mean(trace(maxWindowStart:maxWindowEnd));
                    stats.response.min = mean(trace(minWindowStart:minWindowEnd));
                    stats.response.maxTime = maxIdx * 1000/options.outputFs;
                    stats.response.minTime = minIdx * 1000/options.outputFs;
                    if vhold < -50; stats.response.peak = stats.response.min;
                    elseif vhold > 10; stats.response.peak = stats.response.max; 
                    end
        
                    % Find min and max for control baseline
                    trace = processed_trace(controlWindow);
                    [~,maxIdx] = max(trace); [~,minIdx] = min(trace);
                    % Average around max/min idx to get final value
                    maxWindowStart = max(1,maxIdx-peakWindowWidth);
                    maxWindowEnd = min(maxIdx+peakWindowWidth,length(trace));
                    minWindowStart = max(1,minIdx-peakWindowWidth);
                    minWindowEnd = min(minIdx+peakWindowWidth,length(trace));
                    stats.baseline.max = mean(trace(maxWindowStart:maxWindowEnd));
                    stats.baseline.min = mean(trace(minWindowStart:minWindowEnd));
                    stats.baseline.maxTime = maxIdx * 1000/options.outputFs;
                    stats.baseline.minTime = minIdx * 1000/options.outputFs;
                    if vhold < -50; stats.baseline.peak = stats.baseline.min;
                    elseif vhold > 10; stats.baseline.peak = stats.baseline.max; 
                    end
    
                    % Find E/I index
                    stats.response.EIindex = (abs(stats.response.max)-abs(stats.response.min)) / (abs(stats.response.max)+abs(stats.response.min));
                    stats.baseline.EIindex = (abs(stats.baseline.max)-abs(stats.baseline.min)) / (abs(stats.baseline.max)+abs(stats.baseline.min));
                    %% Store data for the current spot
                    if isHotspotSweep
                        hotIdx = height(hotspotSpots) + 1;
                        hotspotSpots{hotIdx,'Session'} = cellResultsPath;
                        hotspotSpots{hotIdx,'Animal'} = exp{row,'Animal'};
                        hotspotSpots{hotIdx,'Task'} = exp{row,'Task'};
                        hotspotSpots{hotIdx,'Epoch'} = exp{row,'Epoch'};
                        hotspotSpots{hotIdx,'Cell'} = exp{row,'Cell'};
                        hotspotSpots{hotIdx,'Vhold'} = vhold;
                        hotspotSpots{hotIdx,'Depth'} = protocol.depth;
                        hotspotSpots{hotIdx,'Repetition'} = protocol.repetition;
                        hotspotSpots{hotIdx,'Sweep'} = string(sweepAcq{k});
                        hotspotSpots{hotIdx,'Location'} = num2cell(location,[1 2]);
                        hotspotSpots{hotIdx,'Protocol'} = {protocol};
                        hotspotSpots{hotIdx,'Response'} = {responses};
                        hotspotSpots{hotIdx,'Stats'} = {stats};
                        hotspotSpots{hotIdx,'QC'} = {qc};
                        hotspotSpots{hotIdx,'Options'} = {options};
                    else
                        curSpot = curSpot + 1;
                        spots{curSpot,'Session'} = cellResultsPath;
                        spots{curSpot,'Animal'} = exp{row,'Animal'};
                        spots{curSpot,'Task'} = exp{row,'Task'};
                        spots{curSpot,'Epoch'} = exp{row,'Epoch'};
                        spots{curSpot,'Cell'} = exp{row,'Cell'};
                        spots{curSpot,'Vhold'} = vhold;
                        spots{curSpot,'Depth'} = protocol.depth;
                        spots{curSpot,'Repetition'} = protocol.repetition;
                        spots{curSpot,'Sweep'} = string(sweepAcq{k});
                        spots{curSpot,'Location'} = num2cell(location,[1 2]);
                        spots{curSpot,'Protocol'} = {protocol};
                        spots{curSpot,'Response'} = {responses};
                        spots{curSpot,'Stats'} = {stats};
                        spots{curSpot,'QC'} = {qc};
                        spots{curSpot,'Options'} = {options};
                    end
                end
    
                % Save spots.mat for a specific depth
                % If not, spots are too big and matlab can't load later
                searchRowsCount = sum(strcmp(fullSearchTable.sweepStage, 'search'));
                if ~isHotspotSweep && (k == length(sweepAcq) || prevDepth ~= protocol.depth || curFSRow == searchRowsCount)
                    if options.save
                        % Generate spotsAtDepth
                        if k == length(sweepAcq)
                            spotsAtDepth = spots(spots.Depth == protocol.depth,:);
                        else
                            spotsAtDepth = spots(spots.Depth == prevDepth,:);
                            spots(spots.Depth == prevDepth,:) = [];
                            spots = rmmissing(spots,DataVariables='Session');
                        end
                        spotsAtDepth = rmmissing(spotsAtDepth,DataVariables='Session');
    
                        % Save spotsAtDepth
                        filename = strcat('spots_cell',num2str(exp{row,'Cell'}),...
                                          '_epoch',num2str(exp{row,'Epoch'}),...
                                          '_depth',num2str(prevDepth));
                        mkdir(cellResultsPath);
                        save(fullfile(cellResultsPath,filename),'spotsAtDepth','-v7.3');
                        disp(strcat("New spots.mat created & saved: ",filename));
                        
                        prevDepth = protocol.depth; clearvars AD*
                    end
                end
            end

            % after finishing all sweeps for this epoch:
            if options.save && ~isempty(hotspotSpots)
                saveSweepHotspots(hotspotSpots, cellResultsPath, exp{row,'Cell'}, exp{row,'Epoch'});
            end

        end
        
        %% Loop through all depth and add search results to cells.mat

        disp('Ongoing: adding current search to cells_DMD.mat'); clearvars AD*

        % Get a list of spots.mat files of a given cell+epoch.
        % We may have multiple "search keys" for the same epoch (e.g., hotspot sweeps saved as
        % spots_cellX_epochY_hotspot_m70_depth*.mat). Each key is treated as a separate search
        % entry in cells_DMD so analyzeDMDSearch can plot them.
        basePrefix = strcat('spots_cell',num2str(exp{row,'Cell'}),'_epoch',num2str(exp{row,'Epoch'}));
        spotsAll = dir(fullfile(cellResultsPath, strcat(basePrefix,'*_depth*.mat')));

        % Skip if no spots files are found for this epoch
        if isempty(spotsAll)
            filename = basePrefix;
            disp(['Skipped: search epoch ', filename,' not found, skipped instead']);
            if row == size(exp,1)
                responseMap.responseMap = cellResponseMap';
                responseMap.isResponseMap = isResponseMap_cell';
                responseMap.hotspotMap = cellHotspotMap';
                responseMap.currentMap = cellCurrentMap';
                responseMap.baselineMap = cellBaselineMap';
                responseMap.depths = cellDepthList';
                responseMap.spotSequence = cellSpotSequence';
                responseMap.hotspot = cellSpotResponse';
                responseMap.spotLocation = cellSpotLocation';

                cellStats.max = cellMaxResponse';
                cellStats.min = cellMinResponse';
                cellStats.maxTime = cellMaxTime';
                cellStats.minTime = cellMinTime';
                cellStats.auc = cellAUC';
                cellStats.EIindex = cellEIindex';
                cellStats.baseline.auc = cellBaselineAUC';
                cellStats.baseline.std = cellBaselineSTD';

                options.searchTotalSpots = searchTotalSpots;
                options.cellLocation = [spotsAtDepth{1,'Protocol'}{1}.cellX, spotsAtDepth{1,'Protocol'}{1}.cellY];
                options.spotOptions = spotsAtDepth{1,'Options'}{1};

                cells{curCell,'Session'} = options.saveDataPath;
                cells{curCell,'Animal'} = spotsAtDepth{1,'Animal'};
                cells{curCell,'Task'} = spotsAtDepth{1,'Task'};
                cells{curCell,'Cell'} = curCell;
                cells{curCell,'Epochs'} = {cellEpochs'};
                cells{curCell,'Vhold'} = {cellVhold};
                cells{curCell,'Protocol'} = {cellProtocols'};
                cells{curCell,'Response map'} = {responseMap};
                cells{curCell,'Stats'} = {cellStats};
                cells{curCell,'Options'} = {options};

                % Save noise data for this cell
                if options.reload
                    filename = strcat('noise_cell',num2str(curCell));
                    save(fullfile(sessionPath,filename),'nullSpotData','preStimData','baselineData','-v7.3');
                    disp(strcat("Saved noise data for cell: ",num2str(curCell)));
                end

                % Update prevCell and reset cellResponseMap
                cellVhold = []; cellEpochs = {}; searchTotalSpots = {};
                cellDepthList = {}; cellSpotSequence = {}; cellProtocols = {};
                cellResponseMap = {}; isResponseMap_cell = {}; cellHotspotMap = {};
                cellCurrentMap = {}; cellBaselineMap = {};
                nullSpotData = []; preStimData = []; baselineData = [];
                cellSpotResponse = {};
                cellMaxResponse = {}; cellMinResponse = {};
                cellMaxTime = {};cellMinTime = {};
                cellAUC = {}; cellEIindex = {};
                cellBaselineAUC = {}; cellBaselineSTD = {};
                cellSpotLocation = {};
            end
            continue;
        end

        % Derive unique search keys from filenames (strip trailing _depthN.mat)
        names = {spotsAll.name};
        keys = cellfun(@(n) regexprep(n,'_depth\d+\.mat$',''), names, 'UniformOutput', false);
        searchKeys = unique(keys,'stable');

        % Ensure the basePrefix (main search) is processed first if present
        if any(strcmp(searchKeys, basePrefix))
            searchKeys = searchKeys(:);
            searchKeys = [{basePrefix}; searchKeys(~strcmp(searchKeys, basePrefix))];
        end

        % Process each searchKey as a separate "search" entry in cells_DMD
        for kSearch = 1:numel(searchKeys)
            filename = searchKeys{kSearch};
            spotsList = spotsAll(strcmp(keys, filename));

            disp(['Ongoing: adding search epoch ', filename,' to cells_DMD.mat']);
            totalSpots = 0; % total spot for this search

            % Initialize params
            nDepth = length(spotsList);
            searchResponseMap = zeros(684,608,nDepth);
            isResponseMap_search = zeros(684,608,nDepth);
            searchHotspotMap = zeros(684,608,nDepth);
            searchCurrentMap = cell(nDepth,1);
            searchBaselineMap = cell(nDepth,1);
            searchDepthList = zeros(1,nDepth);
            searchProtocols = cell(nDepth,1);
            searchSpotSequence = cell(nDepth,1);
            searchSpotLocation = cell(nDepth,1);

            searchSpotResponse = cell(nDepth,1);
            searchMaxResponse = cell(nDepth,1);
            searchMinResponse = cell(nDepth,1);
            searchMaxTime = cell(nDepth,1);
            searchMinTime = cell(nDepth,1);
            searchAUC = cell(nDepth,1);
            searchEIindex = cell(nDepth,1);
            searchBaselineAUC = cell(nDepth,1);
            searchBaselineSTD = cell(nDepth,1);
    
            % Loop through depth
            for depthIdx = 1:nDepth
                % Initialization
                depthFilePath = fullfile(spotsList(depthIdx).folder,spotsList(depthIdx).name);
                load(depthFilePath,'spotsAtDepth');
                % Get search depth
                d = spotsAtDepth{1,'Depth'};
                % Get search protocol
                searchProtocols{depthIdx} = spotsAtDepth{1,'Protocol'}{1};
                % Get search total spots
                totalSpots = totalSpots + size(spotsAtDepth,1);

                % Get search grid
                locCell = spotsAtDepth.Location;  % cell array, each is [cStart cEnd rStart rEnd]
                colStarts = sort(unique(cellfun(@(L) L(1), locCell)));  % location(1) is col-start (1-indexed)
                rowStarts = sort(unique(cellfun(@(L) L(3), locCell)));  % location(3) is row-start (1-indexed)
                nCol = numel(colStarts);
                nRow = numel(rowStarts);
                nSpotGrid = nCol * nRow;
            
                % Initialize map
                depthResponseMap = zeros(684,608);
                isResponseMap_depth = zeros(684,608);
                depthHotspotMap = zeros(684,608);
                depthCurrentMap  = cell(nSpotGrid,1);
                depthBaselineMap = cell(nSpotGrid,1);
                depthSpotSequence = zeros(size(spotsAtDepth,1),1);
                depthSpotLocation = nan(nSpotGrid,4);
        
                % Build response map
                for s = 1:size(spotsAtDepth,1)
                    % Build depthCurrentMap
                    % depthCurrentMap = cell{nSpotAtDepth,1}
                    % depthCurrentMap{spot1} = zeros(nRep,nSamples)
    
                    % Get spot index
                    location = spotsAtDepth{s,'Location'}{1};
                    % Clip + integerize ranges to valid indices (handles fractional/out-of-bounds locations cleanly)
                    [nX, nY] = size(depthResponseMap);  % xRange indexes rows, yRange indexes columns
                    y1 = max(1, min(nY, floor(min(location(1),location(2)))));
                    y2 = max(1, min(nY,  ceil(max(location(1),location(2)))));
                    x1 = max(1, min(nX, floor(min(location(3),location(4)))));
                    x2 = max(1, min(nX,  ceil(max(location(3),location(4)))));
                    yRange = y1:y2;
                    xRange = x1:x2;

                    % Robust spot index based on matching starts (works with uneven box sizes)
                    colIdx = find(colStarts == location(1), 1);
                    rowIdx = find(rowStarts == location(3), 1);
                
                    % Fallback (should rarely trigger, but avoids hard crash if something is off)
                    if isempty(colIdx)
                        [~, colIdx] = min(abs(colStarts - location(1)));
                    end
                    if isempty(rowIdx)
                        [~, rowIdx] = min(abs(rowStarts - location(3)));
                    end
                
                    spotIdx = (colIdx - 1) + (rowIdx - 1) * nCol + 1;
                    depthSpotLocation(spotIdx,:) = location;

                    % Get response trace
                    trace = spotsAtDepth{s,'Response'}{1}.processed;
                    baseline = spotsAtDepth{s,'Response'}{1}.control;
    
                    % Add to current map
                    depthCurrentMap{spotIdx} = [depthCurrentMap{spotIdx}; trace];
                    depthBaselineMap{spotIdx} = [depthBaselineMap{spotIdx}; baseline];
                    depthSpotSequence(s) = spotIdx;

                    % Get response value
                    originalValue = mean(depthResponseMap(xRange,yRange),'all','omitnan');
                    newValue = spotsAtDepth{s,'Stats'}{1}.response.(options.feature);
                    % Add to depthResponseMap
                    if originalValue==0; depthResponseMap(xRange,yRange) = newValue;
                    else; depthResponseMap(xRange,yRange) = mean([originalValue,newValue]);
                    end
    
                    % Get isResponse value (robust to NaNs)
                    tmp = isResponseMap_depth(xRange,yRange);
                    tmp = tmp(~isnan(tmp));
                    if isempty(tmp); originalValue_isResponse = 0;
                    else; originalValue_isResponse = mode(tmp(:));
                    end
                    newValue_isResponse = spotsAtDepth{s,'Response'}{1}.isResponse;
                    if isempty(newValue_isResponse) || (isnumeric(newValue_isResponse) && isnan(newValue_isResponse))
                        newValue_isResponse = 0;
                    end
                    origBool = (~isnan(originalValue_isResponse)) && (originalValue_isResponse ~= 0);
                    newBool  = (~isnan(newValue_isResponse))      && (newValue_isResponse ~= 0);
                    isResponseMap_depth(xRange,yRange) = origBool | newBool;


                    % Get hotspot value
                    tmp = depthHotspotMap(xRange,yRange);
                    tmp = tmp(~isnan(tmp));
                    if isempty(tmp); originalValue_hotspot = 0;
                    else; originalValue_hotspot = mode(tmp(:));
                    end
                    newValue_hotspot = spotsAtDepth{s,'Response'}{1}.hotspot;
                    if isempty(newValue_hotspot) || (isnumeric(newValue_hotspot) && isnan(newValue_hotspot))
                        newValue_hotspot = 0;
                    end
                    origBool = (~isnan(originalValue_hotspot)) && (originalValue_hotspot ~= 0);
                    newBool  = (~isnan(newValue_hotspot))      && (newValue_hotspot ~= 0);
                    depthHotspotMap(xRange,yRange) = origBool | newBool;
                end        

                % Extract statistics
                hotspotIdx = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Response'}{1}.hotspot,1:height(spotsAtDepth), UniformOutput=false)');
                depthSpotResponse = accumarray(depthSpotSequence(:), hotspotIdx(:), [], @(x) {x});
                maxResponse = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.max,1:height(spotsAtDepth), UniformOutput=false)');
                minResponse = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.min,1:height(spotsAtDepth), UniformOutput=false)');
                depthMaxResponse = accumarray(depthSpotSequence(:), maxResponse(:), [], @(x) {x});
                depthMinResponse = accumarray(depthSpotSequence(:), minResponse(:), [], @(x) {x});
                maxTime = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.maxTime,1:height(spotsAtDepth), UniformOutput=false)');
                minTime = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.minTime,1:height(spotsAtDepth), UniformOutput=false)');
                depthMaxTime = accumarray(depthSpotSequence(:), maxTime(:), [], @(x) {x});
                depthMinTime = accumarray(depthSpotSequence(:), minTime(:), [], @(x) {x});
                auc = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.auc,1:height(spotsAtDepth), UniformOutput=false)');
                EIindex = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.EIindex,1:height(spotsAtDepth), UniformOutput=false)');
                depthAUC = accumarray(depthSpotSequence(:), auc(:), [], @(x) {x});
                depthEIindex = accumarray(depthSpotSequence(:), EIindex(:), [], @(x) {x});
                baselineAUC = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.baseline.auc,1:height(spotsAtDepth), UniformOutput=false)');
                baselineSTD = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.baseline.std,1:height(spotsAtDepth), UniformOutput=false)');
                depthBaselineAUC = accumarray(depthSpotSequence(:), baselineAUC(:), [], @(x) {x});
                depthBaselineSTD = accumarray(depthSpotSequence(:), baselineSTD(:), [], @(x) {x}); 

                % Save statistics
                searchResponseMap(:,:,depthIdx) = depthResponseMap;
                isResponseMap_search(:,:,depthIdx) = isResponseMap_depth;
                searchHotspotMap(:,:,depthIdx) = depthHotspotMap;
                searchCurrentMap{depthIdx} = depthCurrentMap;
                searchBaselineMap{depthIdx} = depthBaselineMap;
                searchDepthList(depthIdx) = d;
                searchSpotSequence{depthIdx} = depthSpotSequence;
                searchSpotLocation{depthIdx} = depthSpotLocation;

                searchSpotResponse{depthIdx} = depthSpotResponse;
                searchMaxResponse{depthIdx} = depthMaxResponse;
                searchMinResponse{depthIdx} = depthMinResponse;
                searchMaxTime{depthIdx} = depthMaxTime;
                searchMinTime{depthIdx} = depthMinTime;
                searchAUC{depthIdx} = depthAUC;
                searchEIindex{depthIdx} = depthEIindex;
                searchBaselineAUC{depthIdx} = depthBaselineAUC;
                searchBaselineSTD{depthIdx} = depthBaselineSTD;
            end
    
            % Save response map
            cellResponseMap{end+1} = searchResponseMap;
            isResponseMap_cell{end+1} = isResponseMap_search;
            cellHotspotMap{end+1} = searchHotspotMap;
            cellCurrentMap{end+1} = searchCurrentMap;
            cellBaselineMap{end+1} = searchBaselineMap;
            cellDepthList{end+1} = searchDepthList;
            cellSpotSequence{end+1} = searchSpotSequence;
            cellSpotLocation{end+1} = searchSpotLocation;

            cellSpotResponse{end+1} = searchSpotResponse;
            cellMaxResponse{end+1} = searchMaxResponse;
            cellMinResponse{end+1} = searchMinResponse;
            cellMaxTime{end+1} = searchMaxTime;
            cellMinTime{end+1} = searchMinTime;
            cellAUC{end+1} = searchAUC;
            cellEIindex{end+1} = searchEIindex;
            cellBaselineAUC{end+1} = searchBaselineAUC;
            cellBaselineSTD{end+1} = searchBaselineSTD;
    
            % Save to cell
            curCell = spotsAtDepth{1,'Cell'};
            cellVhold = [cellVhold; spotsAtDepth{1,'Vhold'}];
            cellEpochs{end+1} = filename;
            searchTotalSpots{end+1} = totalSpots;
            cellProtocols{end+1} = searchProtocols;
        end

        if row == size(exp,1) || curCell ~= exp{row+1,'Cell'}
            responseMap.responseMap = cellResponseMap';
            responseMap.isResponseMap = isResponseMap_cell';
            responseMap.hotspotMap = cellHotspotMap';
            responseMap.currentMap = cellCurrentMap';
            responseMap.baselineMap = cellBaselineMap';
            responseMap.depths = cellDepthList';
            responseMap.spotSequence = cellSpotSequence';
            responseMap.hotspot = cellSpotResponse';
            responseMap.spotLocation = cellSpotLocation';

            cellStats.max = cellMaxResponse';
            cellStats.min = cellMinResponse';
            cellStats.maxTime = cellMaxTime';
            cellStats.minTime = cellMinTime';
            cellStats.auc = cellAUC';
            cellStats.EIindex = cellEIindex';
            cellStats.baseline.auc = cellBaselineAUC';
            cellStats.baseline.std = cellBaselineSTD';

            options.searchTotalSpots = searchTotalSpots;
            options.cellLocation = [spotsAtDepth{1,'Protocol'}{1}.cellX, spotsAtDepth{1,'Protocol'}{1}.cellY];
            options.spotOptions = spotsAtDepth{1,'Options'}{1};
    
            cells{curCell,'Session'} = options.saveDataPath;
            cells{curCell,'Animal'} = spotsAtDepth{1,'Animal'};
            cells{curCell,'Task'} = spotsAtDepth{1,'Task'};
            cells{curCell,'Cell'} = curCell;
            cells{curCell,'Epochs'} = {cellEpochs'};
            cells{curCell,'Vhold'} = {cellVhold};
            cells{curCell,'Protocol'} = {cellProtocols'};
            cells{curCell,'Response map'} = {responseMap};
            cells{curCell,'Stats'} = {cellStats};
            cells{curCell,'Options'} = {options};
    
            % Save noise data for this cell
            if options.reload
                filename = strcat('noise_cell',num2str(curCell));
                save(fullfile(sessionPath,filename),'nullSpotData','preStimData','baselineData','-v7.3');
                disp(strcat("Saved noise data for cell: ",num2str(curCell)));
            end
    
            % Update prevCell and reset cellResponseMap
            cellVhold = []; cellEpochs = {}; searchTotalSpots = {}; 
            cellDepthList = {}; cellSpotSequence = {}; cellProtocols = {};
            cellResponseMap = {}; isResponseMap_cell = {}; cellHotspotMap = {};
            cellCurrentMap = {}; cellBaselineMap = {};
            nullSpotData = []; preStimData = []; baselineData = [];
            cellSpotResponse = {}; 
            cellMaxResponse = {}; cellMinResponse = {};
            cellMaxTime = {};cellMinTime = {};
            cellAUC = {}; cellEIindex = {};
            cellBaselineAUC = {}; cellBaselineSTD = {};
            cellSpotLocation = {};
        end
    end
end

%% Cell specific analysis

if options.reloadCellAnalysis
    % The main goal here is to find common spots and calculate there difference
    % in AUC or other feature
    cells = rmmissing(cells,DataVariables='Session');
    cellList = unique(cells.Cell);
    
    for cellIdx = 1:size(cells,1)
        c = cellList(cellIdx);
        cellData = cells(cells.Cell == c,:);
        searchPerCell = length(cellData.Vhold{1});

        % Load cellStats struct
        cellStats = cells{cellIdx,"Stats"}{1};

        % Initialize new variables
        diff_rmap_cell = {}; common_isResponse_cell = {}; 
        diff_pairs = {}; diff_vholds = []; diff_depthList = {};
    
        %% Build a noise model for a given cell
        % pd1 is a gaussian model fit to null spots
        % pd2 is a gaussian model fit to -30ms to -5ms period prior to
        % photostimulation onset
        % pd3 is a gaussian model fit to pooled datapoints of null spots,
        % prestim period, and baseline period.
    
        disp(['Ongoing: building noise model for cell ',num2str(num2str(c))]);
    
        % Load noise data
        cellResultsPath = string(strcat(options.saveDataPath,filesep,['cell',num2str(c)]));
        filename = strcat('noise_cell',num2str(c),'.mat');
        load(fullfile(sessionPath,filename));
    
        % Initialize noise summary plot
        initializeFig(0.75,1); tiledlayout(1,3);

        % Remove outliers
        outlier1 = isoutlier(nullSpotData);
        outlier2 = isoutlier(preStimData);
        outlier3 = isoutlier(baselineData);

        % Handle empty nullSpotData
        if isempty(nullSpotData)
            warning('Cell %d has NO null spots. Using preStimData as proxy for visualization.', c);
            nullSpotData = preStimData; % Use pre-stim data for the specific "Null Spot" plot
            titleSuffix = ' (Proxy: Pre-stim)';
        else
            titleSuffix = '';
        end
        
        % pd1
        nexttile;
        cleanNull = nullSpotData(~outlier1);
        if numel(cleanNull) > 1
            noise_nullSpot = fitdist(cleanNull(:), 'Normal'); % Use (:) to enforce column vector
            histogram(nullSpotData, 'Normalization', 'pdf'); hold on
            if ~isempty(cleanNull)
                histogram(cleanNull, 'Normalization', 'pdf'); hold on
                xrange = [min(nullSpotData), max(nullSpotData)];
                if diff(xrange)==0; xrange=[xrange(1)-1, xrange(2)+1]; end
                plot(linspace(xrange(1),xrange(2)), pdf(noise_nullSpot, linspace(xrange(1),xrange(2))), 'LineWidth', 2);
            end
        else
            % Fallback distribution if absolutely no data exists (rare)
            noise_nullSpot = makedist('Normal', 'mu', 0, 'sigma', 1);
        end
        title(['Null spot response', titleSuffix]);
        
        % pd2
        nexttile;
        cleanPre = preStimData(~outlier2);
        if numel(cleanPre) > 1
            noise_preStim = fitdist(cleanPre(:), 'Normal'); % Use (:) to enforce column vector
            histogram(preStimData, 'Normalization', 'pdf'); hold on
            if ~isempty(cleanPre)
                histogram(cleanPre, 'Normalization', 'pdf'); hold on
                xrange = [min(preStimData), max(preStimData)];
                if diff(xrange)==0; xrange=[xrange(1)-1, xrange(2)+1]; end
                plot(linspace(xrange(1),xrange(2)), pdf(noise_preStim, linspace(xrange(1),xrange(2))), 'LineWidth', 2);
            end
        else
            % Fallback distribution if absolutely no data exists (rare)
            noise_preStim = makedist('Normal', 'mu', 0, 'sigma', 1);
        end
        title(['Pre-stim response', titleSuffix]);


        % pd3
        nexttile;
        allNullData = [nullSpotData(~outlier1),preStimData(~outlier2),baselineData(~outlier3)];
        if numel(allNullData) > 1
            noise_all = fitdist(allNullData(:), 'Normal'); % Use (:) to enforce column vector
            histogram(allNullData, 'Normalization', 'pdf'); hold on
            xrange = [min(allNullData), max(allNullData)];
            if diff(xrange)==0; xrange=[xrange(1)-1, xrange(2)+1]; end
            plot(linspace(xrange(1),xrange(2)), pdf(noise_all, linspace(xrange(1),xrange(2))), 'LineWidth', 2);
        else
            % Fallback distribution if absolutely no data exists (rare)
            noise_all = makedist('Normal', 'mu', 0, 'sigma', 1);
        end
        title(['Baseline + Pre-stim + Null spot response', titleSuffix]);
    
        % Set threshold based on the noise model (3*standard devation of the
        % symmetric noise)
        cellStats.Ethres = -options.thresholdFactor * noise_all.sigma;
        cellStats.Ithres = options.thresholdFactor * noise_all.sigma;
        cellStats.noiseSD = noise_all.sigma;
        
        % Save noise model
        if options.save
            saveFigures(gcf,['cell',num2str(c),'_noise'],cellResultsPath);
            save(fullfile(sessionPath,filename),'noise_nullSpot','noise_preStim','noise_all','allNullData','-append');
            disp(strcat("Saved noise model for cell: ",num2str(c))); close all;
        end
        
        %% Analyze all pairs of searches
        for search1 = 1:searchPerCell
    
            % Unfinished: build Vhold average map
    
            % Loop through all search pairs
            for search2 = search1+1:searchPerCell
                % Get response map
                search1_rmap = cellData.("Response map"){1}.responseMap{search1};
                search2_rmap = cellData.("Response map"){1}.responseMap{search2};
                search1_isResponse = cellData.("Response map"){1}.isResponseMap{search1};
                search2_isResponse = cellData.("Response map"){1}.isResponseMap{search2};
                search1_depthList = cellData.("Response map"){1}.depths{search1};
                search2_depthList = cellData.("Response map"){1}.depths{search2};
    
                % Check whether Vhold are the same
                diffVhold = cellData.Vhold{1}(search1) ~= cellData.Vhold{1}(search2);
    
                % Check all common depth
                [commonDepthList, search1_commonIdx] = intersect(search1_depthList,search2_depthList);
                [~, search2_commonIdx] = intersect(search2_depthList,search1_depthList);
                search1_rmap = search1_rmap(:,:,search1_commonIdx);
                search2_rmap = search2_rmap(:,:,search2_commonIdx);
                search1_isResponse = search1_isResponse(:,:,search1_commonIdx);
                search2_isResponse = search2_isResponse(:,:,search2_commonIdx);
    
                % Initialize difference map
                diff_rmap = zeros(size(search1_rmap));
                common_isResponse = zeros(size(search1_rmap));
    
                % Calculate difference map for each depth
                for d = 1:length(commonDepthList)
                    if diffVhold
                        % Add two response map together to calculate difference
                        diff_rmap(:,:,d) = search1_rmap(:,:,d) + search2_rmap(:,:,d);
                    else
                        % Subtract two response map to calculate differnece
                        diff_rmap(:,:,d) = search1_rmap(:,:,d) - search2_rmap(:,:,d);
                    end
                    common_isResponse(:,:,d) = search1_isResponse(:,:,d) & search2_isResponse(:,:,d);
                end
    
                % Save difference map
                diff_vholds = [diff_vholds; diffVhold];
                diff_pairs{end+1} = [search1, search2];
                diff_rmap_cell{end+1} = diff_rmap;
                common_isResponse_cell{end+1} = common_isResponse;
                diff_depthList{end+1} = commonDepthList;
            end 
        end
    
        % Save in searchDiff structure
        searchDiff.response = diff_rmap_cell';
        searchDiff.commonSpots = common_isResponse_cell';
        searchDiff.pair = diff_pairs';
        searchDiff.diffVhold = diff_vholds;
        searchDiff.commonDepths = diff_depthList';

        % Save in cells.mat
        cells{cellIdx,'Stats'} = {cellStats};
        cells{cellIdx,'Difference map'} = {searchDiff};
    end
end

%% Save cells.mat

if options.save
    save(strcat(options.saveDataPath,filesep,'cells_DMD_',expName),'cells','-v7.3');
    disp(strcat("Saved cells_DMD.mat: ",expName));
end
close all

end








%% ==================== Helper functions ====================

function [hotspotsByDepth, depths] = localLoadHotspots(epochPath, epochNumber)
    % Load hotspot geometry tables into a map keyed by depth.
    % Supports BOTH:
    %   - Epoch*_finalHotspots_Depth*.mat
    %   - Epoch*_maxSearchHotspots_Depth*.mat
    
    hotspotsByDepth = containers.Map('KeyType','double','ValueType','any');
    depths = [];
    
    % ---- 1) maxDepthHotspots (hotspotSearch stage-1 geometry) ----
    pat1 = fullfile(epochPath, sprintf('Epoch%d_maxSearchHotspots_Depth*.mat', epochNumber));
    files1 = dir(pat1);
    for i = 1:numel(files1)
        fname = files1(i).name;
        tok = regexp(fname, 'Depth(\d+)\.mat$', 'tokens', 'once');
        if isempty(tok); continue; end
        d = str2double(tok{1});
    
        S = load(fullfile(files1(i).folder, fname));
        try
            pack = struct();
            pack.depthKind   = 'maxSearchDepth';
            pack.spotsTable  = S.maxDepthHotspots;
            pack.hotspotMeta = loadHotspotMeta(epochPath, epochNumber, d, 'maxSearchDepth', 'maxDepthHotspots');

            if isKey(hotspotsByDepth, d)
                existing = hotspotsByDepth(d);
                if ~iscell(existing); existing = {existing}; end
                hotspotsByDepth(d) = [existing, {pack}];
            else
                hotspotsByDepth(d) = pack;
            end
            depths(end+1) = d;
            disp(['     [localLoadHotspots]: loaded ', fname]);
        catch
            warning("   maxSearch hotspots loading failed: .mat not found or is empty, skipped");
        end
    end

    % ---- 1) finalHotspots (subdivided geometry) ----
    pat2 = fullfile(epochPath, sprintf('Epoch%d_finalHotspots_Depth*.mat', epochNumber));
    files2 = dir(pat2);
    for i = 1:numel(files2)
        fname = files2(i).name;
        tok = regexp(fname, 'Depth(\d+)\.mat$', 'tokens', 'once');
        if isempty(tok); continue; end
        d = str2double(tok{1});
    
        S = load(fullfile(files2(i).folder, fname));
        try
            pack = struct();
            pack.depthKind   = 'finalDepth';
            pack.spotsTable  = S.finalHotspots;
            pack.hotspotMeta = loadHotspotMeta(epochPath, epochNumber, d, 'finalDepth', 'finalHotspots');

            if isKey(hotspotsByDepth, d)
                existing = hotspotsByDepth(d);
                if ~iscell(existing); existing = {existing}; end
                hotspotsByDepth(d) = [existing, {pack}];
            else
                hotspotsByDepth(d) = pack;
            end
            depths(end+1) = d;
            disp(['     [localLoadHotspots]: loaded ', fname]);
        catch
            warning("   Final hotspots loading failed: .mat not found or is empty, skipped");
        end
    end
    
    depths = unique(depths);
end



function hotspotMeta = loadHotspotMeta(epochPath, epochNum, depthValue, depthKind, stageName)
    % localLoadHotspotMeta
    % Tries to load the hotspotMeta file produced by hotspotSearch_shun.m:
    %   Epoch#_hotspotMeta_maxSearchDepth#.mat  OR  Epoch#_hotspotMeta_finalDepth#.mat
    %
    % Returns:
    %   - the loaded hotspotMeta struct (with hotspotMeta.stage as a struct array), OR
    %   - a compatibility-safe minimal struct if not found / incompatible.

    hotspotMeta = [];

    % Expected filenames from hotspotSearch_shun.m
    % depthKind should be 'maxSearchDepth' or 'finalDepth'
    metaFile = fullfile(epochPath, sprintf('Epoch%d_hotspotMeta_%s.mat', epochNum, depthKind));

    if exist(metaFile, 'file')
        S = load(metaFile);
        if isfield(S, 'hotspotMeta')
            hotspotMeta = S.hotspotMeta;
        else
            % In case variable name differs, attempt best-effort:
            fns = fieldnames(S);
            if ~isempty(fns)
                hotspotMeta = S.(fns{1});
            end
        end
    end

    % ---- Compatibility normalization for downstream analysis ----
    % Downstream code expects:
    %   hotspotMeta.stage is a struct array, and may include stage(ii).acqNums, etc.
    %
    % If loaded file is missing / malformed, create a minimal struct that will NOT break code.
    if isempty(hotspotMeta) || ~isstruct(hotspotMeta)
        hotspotMeta = struct();
    end

    % Some older code patterns mistakenly used hotspotMeta.stage as a string.
    % Normalize so hotspotMeta.stage is always a struct array.
    if ~isfield(hotspotMeta, 'stage') || isempty(hotspotMeta.stage) || ischar(hotspotMeta.stage) || isstring(hotspotMeta.stage)
        % Minimal, analysis-safe stage entry:
        hotspotMeta.stage = struct( ...
            'name',      stageName, ...
            'depthKind', depthKind, ...
            'depth',     depthValue, ...
            'acqNums',   [], ...
            'notes',     '' ...
        );
    else
        % Ensure struct array + ensure required fields exist
        if ~isstruct(hotspotMeta.stage)
            hotspotMeta.stage = struct( ...
                'name',      stageName, ...
                'depthKind', depthKind, ...
                'depth',     depthValue, ...
                'acqNums',   [], ...
                'notes',     '' ...
            );
        else
            for ii = 1:numel(hotspotMeta.stage)
	                % Normalize acquisition-number field naming across versions
	                % (acqNum / acqNumbers / acqNums).
	                if ~isfield(hotspotMeta.stage(ii), 'acqNums')
	                    hotspotMeta.stage(ii).acqNums = [];
	                end
	                if isempty(hotspotMeta.stage(ii).acqNums)
	                    if isfield(hotspotMeta.stage(ii), 'acqNum')
	                        hotspotMeta.stage(ii).acqNums = hotspotMeta.stage(ii).acqNum;
	                    elseif isfield(hotspotMeta.stage(ii), 'acqNumbers')
	                        hotspotMeta.stage(ii).acqNums = hotspotMeta.stage(ii).acqNumbers;
	                    end
	                end
	                % Ensure numeric row vector
	                if ~isempty(hotspotMeta.stage(ii).acqNums)
	                    hotspotMeta.stage(ii).acqNums = double(hotspotMeta.stage(ii).acqNums(:)');
	                end
                if ~isfield(hotspotMeta.stage(ii), 'depth');     hotspotMeta.stage(ii).depth = depthValue; end
                if ~isfield(hotspotMeta.stage(ii), 'name');      hotspotMeta.stage(ii).name = stageName; end
                if ~isfield(hotspotMeta.stage(ii), 'depthKind'); hotspotMeta.stage(ii).depthKind = depthKind; end
                if ~isfield(hotspotMeta.stage(ii), 'notes');     hotspotMeta.stage(ii).notes = ''; end
            end
        end
    end

    % Also store convenience fields (safe even if downstream ignores them)
    hotspotMeta.loadedFrom = metaFile;
end



function [sweepSpots, depthChosen, stageName, depthKind] = ...
    getHotspotSweeps(hotspotsByDepth, depths, desiredDepth, sweepName, expectedNumPulses, depthKindHint)
    % Return a sweepSpots table compatible with fullSearchTable slicing.
    %
    % Handles BOTH hotspot sweep types:
    %   - maxSearchDepth (Epoch*_maxSearchHotspots_Depth*.mat)
    %   - finalDepth     (Epoch*_finalHotspots_Depth*.mat)
    %
    % Selection strategy (robust + backward compatible):
    %   1) If we can parse the acquisition number from sweepName AND any
    %      hotspotMeta file maps that acq number to a stage, prefer that pack.
    %   2) Otherwise, use desiredDepth if present.
    %   3) Otherwise, fall back to the maximum available depth.

    if isempty(depths)
        error('No hotspot geometry files found (Epoch*_finalHotspots_Depth*.mat or Epoch*_maxSearchHotspots_Depth*.mat), but extra sweeps were detected. Cannot map extra sweeps.');
    end

    if nargin < 5; expectedNumPulses = []; end
    if nargin < 6; depthKindHint = ""; end

    % Parse acquisition number (best-effort)
    acqNum = [];
    tok = regexp(sweepName, 'AD0_(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        acqNum = str2double(tok{1});
    end

    % ---- 1) Prefer pack whose hotspotMeta explicitly contains this acqNum ----
    packChosen = [];
    depthChosen = NaN;
    stageName = "hotspot_unknown";
    depthKind = "";

    if ~isempty(acqNum)
        for di = 1:numel(depths)
            d = depths(di);
            packVal = hotspotsByDepth(d);
            if ~iscell(packVal); packVal = {packVal}; end
            for pi = 1:numel(packVal)
                p = packVal{pi};
                if isstruct(p) && isfield(p,'hotspotMeta') && ~isempty(p.hotspotMeta) && isfield(p.hotspotMeta,'stage')
                    st = p.hotspotMeta.stage;
                    for ii = 1:numel(st)
                        if isfield(st(ii),'acqNums') && any(st(ii).acqNums == acqNum)
                            packChosen = p;
                            depthChosen = d;
                            stageName = string(st(ii).name);
                            if isfield(p,'depthKind'); depthKind = string(p.depthKind); end
                            break
                        end
                    end
                end
                if ~isempty(packChosen); break; end
            end
            if ~isempty(packChosen); break; end
        end
    end

    % ---- 2) Fall back to desiredDepth or max(depths) ----
    if isempty(packChosen)
        if ismember(desiredDepth, depths)
            depthChosen = desiredDepth;
        else
            depthChosen = max(depths);
        end

        packVal = hotspotsByDepth(depthChosen);
        if ~iscell(packVal); packVal = {packVal}; end
        
        % Prefer matching depthKindHint (finalDepth vs maxSearchDepth) if provided
        idx = [];
        if depthKindHint ~= ""
            idx = find(cellfun(@(p) isstruct(p) && isfield(p,'depthKind') && string(p.depthKind)==string(depthKindHint), packVal), 1);
        end
        
        % Next prefer matching expectedNumPulses if provided
        if isempty(idx) && ~isempty(expectedNumPulses)
            idx = find(cellfun(@(p) height(getHotspotTable(p)) == expectedNumPulses, packVal), 1);
        end
        
        % Otherwise take the first
        if isempty(idx); idx = 1; end
        packChosen = packVal{idx};
        
        if isstruct(packChosen) && isfield(packChosen,'depthKind')
            depthKind = string(packChosen.depthKind);
        end

        if isstruct(packChosen) && isfield(packChosen,'depthKind')
            depthKind = string(packChosen.depthKind);
        end
    end

    % ---- Extract geometry table from chosen pack (supports old + new formats) ----
    if istable(packChosen)
        hotspotTable = packChosen;
        hotspotMeta = [];
    else
        hotspotMeta = [];
        if isfield(packChosen,'hotspotMeta'); hotspotMeta = packChosen.hotspotMeta; end

        if isfield(packChosen,'spotsTable')
            hotspotTable = packChosen.spotsTable;
        elseif isfield(packChosen,'finalHotspots')
            hotspotTable = packChosen.finalHotspots;
        elseif isfield(packChosen,'maxDepthHotspots')
            hotspotTable = packChosen.maxDepthHotspots;
        elseif isfield(packChosen,'hotspots')
            hotspotTable = packChosen.hotspots;
        else
            error('Unrecognized hotspot pack format; cannot locate geometry table.');
        end
    end

    needVars = {'depth','xStart','yStart','xWidth','yHeight'};
    missing = setdiff(needVars, hotspotTable.Properties.VariableNames);
    if ~isempty(missing)
        error('Hotspot geometry table is missing required variables: %s', strjoin(missing, ', '));
    end

    sweepSpots = hotspotTable(:, needVars);
    sweepSpots.response = false(height(sweepSpots),1);

    % If stageName is still unknown, attempt best-effort from chosen meta
    if stageName == "hotspot_unknown" && ~isempty(acqNum) && ~isempty(hotspotMeta) && isfield(hotspotMeta,'stage')
        for ii = 1:numel(hotspotMeta.stage)
            if isfield(hotspotMeta.stage(ii),'acqNums') && any(hotspotMeta.stage(ii).acqNums == acqNum)
                stageName = string(hotspotMeta.stage(ii).name);
                break
            end
        end
    end
end


function info = detectHotspotSweep(headerString, protocol, hotspotsByDepth, depths, opts)
% Detect and materialize hotspot sweep geometry even when header is missing stage strings.
%
% Returns struct "info" with fields:
%   .isHotspot (logical)
%   .sweepSpots (table)        % geometry table (depth,xStart,yStart,xWidth,yHeight,response)
%   .depthChosen (double)
%   .depthKind (string)        % "finalDepth" / "maxSearchDepth" / "unknown"
%   .stageName (string)        % parsed stage or "hotspot_inferred"
%   .matchMethod (string)      % "header" / "pulseCount" / "none"
%   .why (string)              % debug reason

    arguments
        headerString
        protocol struct
        hotspotsByDepth
        depths double = []
        opts.AllowCrossDepth (1,1) logical = true
        opts.PreferDepthMatch (1,1) logical = true
        opts.PreferDepthKind (1,1) string = "any"  % "finalDepth"|"maxSearchDepth"|"any"
        opts.PulseTolAbs (1,1) double = NaN               % if NaN -> auto
        opts.PulseTolRel (1,1) double = 0.01              % 1% default
        opts.StageRegex (1,1) string = "state\.zDMD\.sweepStage\s*=\s*'([^']+)'"
        opts.RequireGeometryVars (1,:) string = ["depth","xStart","yStart","xWidth","yHeight"]
        opts.Debug (1,1) logical = false
    end

    info = struct( ...
        'isHotspot', false, ...
        'sweepSpots', table(), ...
        'depthChosen', NaN, ...
        'depthKind', "unknown", ...
        'stageName', "", ...
        'matchMethod', "none", ...
        'why', "" );

    if isempty(depths) || isempty(hotspotsByDepth)
        info.why = "No hotspot tables loaded (no Epoch*_maxSearchHotspots_*.mat / finalHotspots_*.mat).";
        return
    end

    % --- Parse stage from header if possible (newer versions)
    stageName = "";
    depthKindHint = "unknown";
    if ~isempty(headerString)
        tok = regexp(string(headerString), opts.StageRegex, "tokens", "once");
        if ~isempty(tok)
            stageName = string(tok{1});
            if contains(stageName,"finalDepth","IgnoreCase",true)
                depthKindHint = "finalDepth"; 
            elseif contains(stageName,"maxDepth","IgnoreCase",true) || contains(stageName,"maxSearch","IgnoreCase",true)
                depthKindHint = "maxSearchDepth";
            else
                disp(['     [detectHotspotSweep] did not find finalDepth or maxSearchDepth in stageName: ', stageName]);
            end
        end
    end

    % --- Get expected pulse count (critical for old headers)
    nP = [];
    if isfield(protocol,'numPulses') && ~isempty(protocol.numPulses) && ~isnan(protocol.numPulses)
        nP = double(protocol.numPulses);
    end

    % --- Candidate packs: collect (depth, depthKind, nBoxes, table)
    cands = collectHotspotCandidates(hotspotsByDepth, depths, opts);
    if isempty(cands)
        info.why = "Hotspot tables exist, but none contain required geometry variables.";
        return
    end

    % --- Filter by depth domain
    d0 = NaN;
    if isfield(protocol,'depth') && ~isempty(protocol.depth) && ~isnan(protocol.depth)
        d0 = double(protocol.depth);
    end
    if opts.PreferDepthMatch && ~isnan(d0)
        % Prefer candidates at same depth; if none, optionally allow cross-depth
        sameDepth = cands([cands.depth] == d0);
        if ~isempty(sameDepth)
            cands = sameDepth;
        elseif ~opts.AllowCrossDepth
            info.why = "No hotspot table at protocol.depth and cross-depth disabled.";
            return
        end
    elseif ~opts.AllowCrossDepth && ~isnan(d0)
        cands = cands([cands.depth] == d0);
        if isempty(cands)
            info.why = "Cross-depth disabled and no candidates at protocol.depth.";
            return
        end
    end

    % --- Filter by depth kind preference
    if opts.PreferDepthKind ~= "any"
        pref = string(opts.PreferDepthKind);
        isPref = strcmpi(string({cands.depthKind}), pref);
        if any(isPref)
            cands = cands(isPref);
        end
    end

    % --- Rank / choose candidate
    chosen = [];
    % 1) If header stage is available, try to honor depthKindHint first
    if stageName ~= ""
        info.stageName = stageName;
        info.matchMethod = "header";

        if depthKindHint ~= "unknown"
            isHint = strcmpi(string({cands.depthKind}), depthKindHint);
            if any(isHint)
                candsH = cands(isHint);
            else
                candsH = cands;
            end
        else
            candsH = cands;
        end

        % If multiple remain, use pulse count if available
        if ~isempty(nP)
            chosen = chooseByPulseCountThenDepth(candsH, nP, d0, opts);
        end

        if isempty(chosen)
            % fallback: choose closest depth (if d0 known), else first
            chosen = chooseByDepthOnly(candsH, d0);
        end
    end

    % 2) If no header stage, use pulse-count matching (old versions)
    if isempty(chosen)
        info.stageName = "hotspot_inferred";
        info.matchMethod = "pulseCount";

        if isempty(nP)
            info.matchMethod = "none";
            info.why = "No sweepStage in header and protocol.numPulses missing -> cannot infer hotspot sweep reliably.";
            return
        end

        chosen = chooseByPulseCountThenDepth(cands, nP, d0, opts);
        if isempty(chosen)
            info.matchMethod = "none";
            info.why = "No hotspot geometry table matches this sweep's numPulses (within tolerance).";
            return
        end
    end

    % --- Materialize output geometry table
    T = chosen.table;
    if ~ismember("response", T.Properties.VariableNames)
        T.response = false(height(T),1);
    end

    info.isHotspot = true;
    info.sweepSpots = T;
    info.depthChosen = chosen.depth;
    info.depthKind = chosen.depthKind;

    if info.matchMethod == "header" && stageName == ""
        info.stageName = "hotspot_header_unknown";
    end

    if opts.Debug
        info.why = sprintf("Chosen: depth=%g kind=%s nBoxes=%d method=%s", ...
            chosen.depth, chosen.depthKind, chosen.nBoxes, info.matchMethod);
    end
end




function cands = collectHotspotCandidates(hotspotsByDepth, depths, opts)
    cands = struct('depth',{},'depthKind',{},'nBoxes',{},'table',{});
    for di = 1:numel(depths)
        d = depths(di);
        pv = hotspotsByDepth(d);
        if ~iscell(pv); pv = {pv}; end
        for pi = 1:numel(pv)
            p = pv{pi};
            T = getHotspotTable(p);
            if isempty(T) || ~istable(T), continue; end
            if ~all(ismember(opts.RequireGeometryVars, string(T.Properties.VariableNames)))
                continue
            end
            cands(end+1).depth = d; %#ok<AGROW>
            cands(end).table = T(:, intersect(T.Properties.VariableNames, cellstr(opts.RequireGeometryVars), 'stable'));
            cands(end).nBoxes = height(T);
            if isstruct(p) && isfield(p,'depthKind')
                cands(end).depthKind = string(p.depthKind);
            else
                cands(end).depthKind = "unknown";
            end
        end
    end
end



function chosen = chooseByPulseCountThenDepth(cands, nP, d0, opts)
    chosen = [];
    tolAbs = opts.PulseTolAbs;
    if isnan(tolAbs)
        tolAbs = max(1, round(opts.PulseTolRel * nP));
    end

    deltas = abs([cands.nBoxes] - nP);
    keep = find(deltas <= tolAbs);
    if isempty(keep), return; end
    c = cands(keep);

    % tie-break: closest depth if d0 known
    if ~isnan(d0)
        [~,ord] = sort(abs([c.depth] - d0), 'ascend');
        c = c(ord);
    end

    % choose first after tie-break
    chosen = c(1);
end



function chosen = chooseByDepthOnly(cands, d0)
    if isempty(cands), chosen = []; return; end
    if ~isnan(d0)
        [~,i] = min(abs([cands.depth] - d0));
        chosen = cands(i);
    else
        chosen = cands(1);
    end
end



function T = getHotspotTable(packChosen)
    T = [];
    try
        if istable(packChosen)
            T = packChosen; return
        end
        if isstruct(packChosen)
            if isfield(packChosen,'spotsTable'),       T = packChosen.spotsTable; return; end
            if isfield(packChosen,'finalHotspots'),    T = packChosen.finalHotspots; return; end
            if isfield(packChosen,'maxDepthHotspots'), T = packChosen.maxDepthHotspots; return; end
        end
    catch
        T = [];
    end
end



function saveSweepHotspots(hotspotSpots, cellResultsPath, cellNum, epochNum)
    % Save hotspot/extra sweeps as separate search keys so analyzeDMDSearch can plot them.
    % Files are saved as:
    %   spots_cellX_epochY_hotspot_<tag>_depthD.mat
    % where tag is based on holding potential (m70 / p10 / other).
    
    if isempty(hotspotSpots); return; end
    if ~exist(cellResultsPath,'dir'); mkdir(cellResultsPath); end
    
    % Tag by holding potential
    % Prefer stage tag from protocol.hotspotStage if present, otherwise fall back to Vhold
    tag = strings(height(hotspotSpots),1);
    kind = strings(height(hotspotSpots),1);

    for i = 1:height(hotspotSpots)
        p = hotspotSpots{i,'Protocol'}{1};
        if isstruct(p)
            % Stage/tag
            if isfield(p,'hotspotStage') && ~isempty(p.hotspotStage)
                st = string(p.hotspotStage);
                if contains(lower(st),'exci') || contains(lower(st),'neg')
                    tag(i) = "exci";
                elseif contains(lower(st),'inhi') || contains(lower(st),'pos') || contains(lower(st),'10')
                    tag(i) = "inhi";
                else
                    tag(i) = "other";
                end

                % Best-effort infer kind from stage name
                if contains(lower(st),'final')
                    kind(i) = "finalDepth";
                elseif contains(lower(st),'max')
                    kind(i) = "maxSearchDepth";
                end
            end

            % Explicit kind from protocol if present (preferred)
            if isfield(p,'hotspotDepthKind') && ~isempty(p.hotspotDepthKind)
                kind(i) = string(p.hotspotDepthKind);
            end
        end
    end
    
    % Fill any still-empty tags using Vhold heuristic
    v = hotspotSpots.Vhold;
    tag(tag=="") = "other";
    tag(tag=="other" & v < -50) = "exci";
    tag(tag=="other" & v > 0)   = "inhi";
    hotspotSpots.Tag = tag;
    hotspotSpots.Kind = kind;

    depths = unique(hotspotSpots.Depth);
    for di = 1:numel(depths)
        d = depths(di);
        subD = hotspotSpots(hotspotSpots.Depth == d, :);
        kinds = unique(subD.Kind);
        for ki = 1:numel(kinds)
            kd = kinds(ki);
            subK = subD(subD.Kind == kd, :);
            tags = unique(subK.Tag);
            for ti = 1:numel(tags)
                t = tags(ti);
                sub = subK(subK.Tag == t, :);

                % Save as spotsAtDepth (expected by downstream loader)
                spotsAtDepth = sub;
                spotsAtDepth = rmmissing(spotsAtDepth,DataVariables='Session');
        
                if kd == "" || kd == "unknown"
                    filename = sprintf('spots_cell%d_epoch%d_hotspot_%s_depth%d', cellNum, epochNum, char(t), d);
                else
                    filename = sprintf('spots_cell%d_epoch%d_hotspot_%s_%s_depth%d', cellNum, epochNum, char(kd), char(t), d);
                end
    
                save(fullfile(cellResultsPath, filename), 'spotsAtDepth', '-v7.3');
                disp(strcat("New spots.mat created & saved (hotspot): ", filename));
            end
        end
    end
end



function [sweepSpots, protocol, nPulsesThisSweep, isHotspotSweep, curFSRow] = ...
    selectSweepSpots(headerString, protocol, fullSearchTable, curFSRow, reconstructedSearch, ...
                     sweepHotspots, sweepDepths, hotOpts, sweepName)

    sweepSpots = table();
    nPulsesThisSweep = 0;
    isHotspotSweep = false;

    if nargin < 9 || isempty(sweepName)
        sweepName = "";
    end

    if ~isfield(protocol,'stimOnset') || isempty(protocol.stimOnset)
        return
    end
    stimCount = numel(protocol.stimOnset);
    nFS = height(fullSearchTable);

    % ---- Case 0: Prefer acqNum matching when available (prevents curFSRow drift) ----
    if ~reconstructedSearch && ...
            ismember('acqNum', fullSearchTable.Properties.VariableNames) && ...
            any(~isnan(fullSearchTable.acqNum))
    
        % Parse acquisition number from sweepName, e.g. "AD0_0123..." -> 123
        thisAcq = nan;
        if ~isempty(sweepName)
            tok = regexp(char(sweepName), 'AD0[_-]?(\d+)', 'tokens', 'once');
            if ~isempty(tok)
                thisAcq = str2double(tok{1});
            end
        end
    
        if ~isnan(thisAcq)
            idx = find(fullSearchTable.acqNum == thisAcq);
            if ~isempty(idx)
    
                sweepSpots = fullSearchTable(idx, :);
    
                % If pulseIndex exists, order spots by pulseIndex (stable sweep ordering)
                if ismember('pulseIndex', sweepSpots.Properties.VariableNames)
                    [~, ord] = sort(sweepSpots.pulseIndex);
                    sweepSpots = sweepSpots(ord, :);
                end
    
                % Respect expected pulses if available; otherwise default to stimCount
                expectedPulses = stimCount;
                if isfield(protocol,'numPulses') && ~isempty(protocol.numPulses) && ~isnan(protocol.numPulses)
                    expectedPulses = double(protocol.numPulses);
                end
    
                nPulsesThisSweep = min([height(sweepSpots), stimCount, expectedPulses]);
                if nPulsesThisSweep <= 0
                    sweepSpots = table();
                    return
                end
                sweepSpots = sweepSpots(1:nPulsesThisSweep, :);
    
                % Detect hotspot stage (same logic as your Case 1)
                if ismember('sweepStage', sweepSpots.Properties.VariableNames)
                    st = unique(string(sweepSpots.sweepStage));
                    st = st(st ~= "" & lower(st) ~= "search");
                    if ~isempty(st)
                        isHotspotSweep = true;
                        protocol.isHotspotSweep = true;
                        protocol.hotspotStage = st(1);
    
                        st1 = lower(st(1));
                        if contains(st1,'final')
                            protocol.hotspotDepthKind = "finalDepth";
                        elseif contains(st1,'max')
                            protocol.hotspotDepthKind = "maxSearchDepth";
                        else
                            protocol.hotspotDepthKind = "";
                        end
                    end
                end
    
                % Trust table depth if consistent
                if ismember('depth', sweepSpots.Properties.VariableNames)
                    dUnique = unique(sweepSpots.depth);
                    dUnique = dUnique(~isnan(dUnique));
                    if numel(dUnique) == 1 && protocol.depth ~= dUnique
                        protocol.depth = dUnique;
                        warning(['Depth mismatch: using fullSearchTable depth = ', num2str(dUnique)]);
                    elseif numel(dUnique) > 1
                        warning('acqNum-matched rows span multiple depths (unexpected).');
                    end
                end
    
                if ~ismember('response', sweepSpots.Properties.VariableNames)
                    sweepSpots.response = false(height(sweepSpots),1);
                end
    
                % IMPORTANT: return early so we do NOT slice by curFSRow
                return
            end
        end
    end


    % ---- Case 1: Use fullSearchTable rows (normal) ----
    if ~reconstructedSearch && curFSRow < nFS

        expectedPulses = stimCount;
        if isfield(protocol,'numPulses') && ~isempty(protocol.numPulses) && ~isnan(protocol.numPulses)
            expectedPulses = double(protocol.numPulses);
        end

        nRemain = nFS - curFSRow;
        nPulsesThisSweep = min([expectedPulses, nRemain, stimCount]);
        if nPulsesThisSweep <= 0
            return
        end

        sweepSpots = fullSearchTable(curFSRow+1 : curFSRow+nPulsesThisSweep, :);
        % ---- Drift check: sliced rows should come from ONE acqNum (or all NaN for old data) ----
        if ismember('acqNum', sweepSpots.Properties.VariableNames)
            a = unique(sweepSpots.acqNum); a = a(~isnan(a));
            if numel(a) > 1
                warning('Likely curFSRow drift: sliced chunk contains multiple acqNum values: %s', mat2str(a'));
            end
        end

        curFSRow = curFSRow + nPulsesThisSweep; % always advance when consuming the table

        % Detect whether this sweep is tagged as a hotspot inside fullSearchTable
        if ismember('sweepStage', sweepSpots.Properties.VariableNames)
            st = unique(string(sweepSpots.sweepStage));
            st = st(st ~= "" & lower(st) ~= "search");
            if ~isempty(st)
                isHotspotSweep = true;
                protocol.isHotspotSweep = true;
                protocol.hotspotStage = st(1);

                st1 = lower(st(1));
                if contains(st1,'final')
                    protocol.hotspotDepthKind = "finalDepth";
                elseif contains(st1,'max')
                    protocol.hotspotDepthKind = "maxSearchDepth";
                else
                    protocol.hotspotDepthKind = "";
                end
            end
        end

        if ismember('depth', sweepSpots.Properties.VariableNames)
            if any(sweepSpots.depth ~= protocol.depth)
                dUnique = unique(sweepSpots.depth);
                if numel(dUnique)==1 && protocol.depth ~= dUnique
                    protocol.depth = dUnique;   % trust geometry table
                    warning(['Selected spots have wrong depth vs protocol.depth. Set to depth as ', num2str(dUnique)]);
                elseif numel(dUnique)>1
                    warning("fullSearchTable slice spans multiple depths -> likely curFSRow drift");
                end
            end
        end

        if ~ismember('response', sweepSpots.Properties.VariableNames)
            sweepSpots.response = false(height(sweepSpots),1);
        end

        return
    end

    % ---- Case 2: reconstructedSearch OR extra sweeps beyond fullSearchTable ----
    hotInfo = detectHotspotSweep(headerString, protocol, sweepHotspots, sweepDepths, ...
                                    PreferDepthMatch = hotOpts.PreferDepthMatch, ...
                                    PreferDepthKind  = hotOpts.PreferDepthKind, ...
                                    Debug            = hotOpts.Debug);

    if ~hotInfo.isHotspot
        % In reconstructedSearch mode: skip non-hotspot sweeps (search already reconstructed from fig)
        % In normal mode: unmappable extra sweep -> skip
        return
    end

    isHotspotSweep = true;
    sweepSpots = hotInfo.sweepSpots;

    protocol.depth = hotInfo.depthChosen;
    protocol.isHotspotSweep = true;
    protocol.hotspotStage = hotInfo.stageName;
    protocol.hotspotDepthKind = hotInfo.depthKind;

    nPulsesThisSweep = min(height(sweepSpots), stimCount);
    sweepSpots = sweepSpots(1:nPulsesThisSweep, :);

    if ~ismember('response', sweepSpots.Properties.VariableNames)
        sweepSpots.response = false(height(sweepSpots),1);
    end
end




%% ======================= Reconstruction functions =======================

function [reconOK, reconNoise] = reconstructFromResponseMap(epochPath, cellResultsPath, expRow, epoch, options)
    % Reconstruct spotsAtDepth depth files from EpochX_DepthY_responseMap.fig
    % so downstream analyzeDMDSearch can run when fullSearchTable is missing.
    
    reconOK = false;
    reconNoise = struct('nullSpotData',[],'preStimData',[],'baselineData',[]);
    
    if ~exist(cellResultsPath,'dir'); mkdir(cellResultsPath); end
    
    figList = dir(fullfile(epochPath, sprintf('Epoch%d_Depth*_responseMap.fig', epoch)));
    if isempty(figList)
        warning(['No responseMap.fig found for Epoch ', num2str(epoch)]);
        return
    end
    
    % Template protocol (try from first sweep)
    protocolTemplate = struct('depth',nan,'pulseWidth',5,'cellX',nan,'cellY',nan,'repetition',1);
    try
        sweepAcq = expRow{1,'Sweep names'}{1};
        if ~isempty(sweepAcq)
            load(fullfile(epochPath, [sweepAcq{1},'.mat']));
            headerString = eval([sweepAcq{1},'.UserData.headerString']);
            p = getCellProtocol(headerString, outputFs=options.outputFs, rcCheckRecoveryWindow=options.rcCheckRecoveryWindow);
            if ~isfield(p,'pulseWidth'); p.pulseWidth = 5; end
            if ~isfield(p,'cellX'); p.cellX = nan; end
            if ~isfield(p,'cellY'); p.cellY = nan; end
            protocolTemplate = p;
        end
    catch
    end
    
    Fs = options.outputFs;
    analysisLen = round(options.analysisWindowLength*Fs/1000)+1;
    ctrlLen     = round(options.controlWindowLength*Fs/1000)+1;
    minHotspotSamples = round(5*Fs/1000); % 5ms
    
    cellNum = expRow{1,'Cell'};
    basePrefix = sprintf('spots_cell%d_epoch%d', cellNum, epoch);
    prefix = basePrefix;
    
    for i = 1:numel(figList)
        figPath = fullfile(figList(i).folder, figList(i).name);
        depthTok = regexp(figList(i).name, 'Depth(\d+)', 'tokens', 'once');
        if isempty(depthTok); continue; end
        depth = str2double(depthTok{1});
        disp(['Reconstructing: ', prefix, ', depth ', num2str(depth)]);

        % 1. Load InfoPatching if needed for Excel lookup (done once per epoch ideally, but okay here)
        infoPatching = [];
        if strcmpi(options.vholdChannel, 'excel')
            try
                csvOpts = detectImportOptions(fullfile(epochPath,'InfoPatching.xlsx'));
                csvOpts.SelectedVariableNames = 1:9; csvOpts.VariableNamesRange = 25;
                infoPatching = readtable(fullfile(epochPath,'InfoPatching.xlsx'), csvOpts);
                infoPatching = rmmissing(infoPatching, DataVariables="acq_");
                if iscell(infoPatching.acq_); infoPatching.acq_ = str2double(infoPatching.acq_); end
                if iscell(infoPatching.holding); infoPatching.holding = str2double(infoPatching.holding); end
            catch
                warning('Could not load InfoPatching.xlsx for reconstruction fallback.');
            end
        end
    
        % 2. Find the sweep file name that corresponds to this depth
        targetSweepName = '';
        try
            sweepNames = expRow{1,'Sweep names'}{1};
            for k = 1:length(sweepNames)
                swName = sweepNames{k};
                swFile = fullfile(epochPath, [swName, '.mat']);
                if ~exist(swFile, 'file')
                    disp(['     [reconstructFromResponseMap] did not find ', swName,' , skipped']);
                    continue; 
                end
                
                % Quick load header to check depth
                S_sw = load(swFile);
                if isfield(S_sw, swName); rec = S_sw.(swName); else; f=fieldnames(S_sw); rec=S_sw.(f{1}); end
                
                % Check protocol depth
                p = getCellProtocol(rec.UserData.headerString, ...
                    outputFs=options.outputFs, rcCheckRecoveryWindow=options.rcCheckRecoveryWindow);
                
                if p.depth == depth
                    targetSweepName = swName;
                    break; 
                end
            end
        catch
        end
    
        nSide = 2^depth;
        nTotal = nSide^2;

        figData = extractFromResponseMap(figPath, nSide);
        if isempty(figData.spotTraces); continue; end
    
        nSpots = numel(figData.spotTraces);
        xEdges0 = round(linspace(0,608,nSide+1));
        yEdges0 = round(linspace(0,684,nSide+1));
    
        plotWindowTime = figData.plotWindowTime;
        [~,eventSample] = min(abs(plotWindowTime));
        analysisWindow = eventSample:min(eventSample+analysisLen-1,numel(plotWindowTime));
    
        depthResponseMap = figData.depthResponseMap;
        if isempty(depthResponseMap); depthResponseMap=zeros(684,608); end
    
        varTypes = {'string','string','string','double','double','double',...
                    'double','double','string','cell','cell','cell','cell','cell','cell'};
        varNames = {'Session','Animal','Task','Epoch','Cell','Vhold',...
                    'Depth','Repetition','Sweep','Location','Protocol',...
                    'Response','Stats','QC','Options'};
        spotsAtDepth = table('Size',[0,length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
    
        for spotIdx = 1:nTotal
            traces = figData.spotTraces{spotIdx};
            if isempty(traces); continue; end
    
            col = mod(spotIdx-1,nSide)+1;
            row = floor((spotIdx-1)/nSide)+1;
            xStart0=xEdges0(col); xWidth=xEdges0(col+1)-xEdges0(col);
            yStart0=yEdges0(row); yHeight=yEdges0(row+1)-yEdges0(row);
            location = [xStart0+1, xStart0+xWidth, yStart0+1, yStart0+yHeight];
    
            yRange = max(1,min(608,location(1):location(2)));
            xRange = max(1,min(684,location(3):location(4)));
            isResponse_img = any(depthResponseMap(xRange,yRange)>0,'all');
    
            for rep = 1:size(traces,1)
                processed = traces(rep,:);
                control   = localMakeSyntheticControl(processed, plotWindowTime, ctrlLen);
    
                [effectiveVhold,~] = getSweepVhold(epochPath, targetSweepName, nan, options.vholdChannel, infoTable=infoPatching);

                baseStd = std(control,0, 'omitnan');
                thr = options.thresholdFactor*baseStd;
                isHotspot = sum(abs(processed(analysisWindow))>=thr)>=minHotspotSamples;
    
                responses = struct('raw',processed,'processed',processed,'control',control,...
                                   'rawTrace',processed,'processedTrace',processed,...
                                   'isResponse',logical(isResponse_img),'hotspot',logical(isHotspot));
    
                stats = localComputeStatsFromProcessed(processed, control, analysisWindow, effectiveVhold, options);
    
                idx = height(spotsAtDepth)+1;
                spotsAtDepth{idx,'Session'}     = string(cellResultsPath);
                spotsAtDepth{idx,'Animal'}      = string(expRow{1,'Animal'});
                spotsAtDepth{idx,'Task'}        = string(expRow{1,'Task'});
                spotsAtDepth{idx,'Epoch'}       = epoch;
                spotsAtDepth{idx,'Cell'}        = cellNum;
                spotsAtDepth{idx,'Vhold'}       = effectiveVhold;
                spotsAtDepth{idx,'Depth'}       = depth;
                spotsAtDepth{idx,'Repetition'}  = rep;
                spotsAtDepth{idx,'Sweep'}       = string(figList(i).name);
                spotsAtDepth{idx,'Location'}    = num2cell(location,[1 2]);
    
                p=protocolTemplate; p.depth=depth; p.repetition=rep; p.pulseWidth=5;
                spotsAtDepth{idx,'Protocol'} = {p};
    
                spotsAtDepth{idx,'Response'}    = {responses};
                spotsAtDepth{idx,'Stats'}       = {stats};
                spotsAtDepth{idx,'QC'}          = {struct('included',true,'reconstructed',true)};
                spotsAtDepth{idx,'Options'}     = {options};
    
                reconNoise.preStimData  = [reconNoise.preStimData, control];
                reconNoise.baselineData = [reconNoise.baselineData, control];
                if ~responses.hotspot
                    reconNoise.nullSpotData = [reconNoise.nullSpotData, processed];
                end
            end
        end
    
        spotsAtDepth = rmmissing(spotsAtDepth,DataVariables='Session');
        if height(spotsAtDepth)==0; continue; end
    
        outName = sprintf('%s_depth%d', prefix, depth);
        save(fullfile(cellResultsPath,outName),'spotsAtDepth','-v7.3');
        disp(['Reconstructed & saved: ', outName]);
        reconOK = true;
    end
end


function figData = extractFromResponseMap(figPath, nSideTotal)
    nTotal = nSideTotal^2;

    figData = struct( ...
        'spotTraces',    {cell(nTotal,1)}, ...
        'plotWindowTime',[], ...
        'depthResponseMap',[] , ...
        'spotRow',[], ...
        'spotCol',[] );

    hFig = openfig(figPath,'invisible');
    cleanupObj = onCleanup(@() close(hFig));

    % --- hotspot response map image (684x608)
    imgs = findall(hFig,'Type','image');
    depthMap = [];
    for i=1:numel(imgs)
        C = imgs(i).CData;
        if isnumeric(C) && ismatrix(C) && all(size(C)==[684 608])
            depthMap = C; break;
        end
    end
    figData.depthResponseMap = depthMap;

    % --- get all axes, exclude any that contains an image (big map panel)
    axAll = findall(hFig,'Type','axes');
    hasImg = false(size(axAll));
    for a = 1:numel(axAll)
        hasImg(a) = ~isempty(findall(axAll(a),'Type','image'));
    end
    axCand = axAll(~hasImg);
    if isempty(axCand); return; end

    % --- choose the tiled spot axes by size (area mode)
    posAll = vertcat(axCand.Position);         % [x y w h] normalized
    area   = posAll(:,3).*posAll(:,4);
    areaQ  = round(area*1e4)/1e4;              % quantize for stable mode
    aMode  = mode(areaQ);

    tol = 0.35;                                % generous layout tolerance
    isSpot = abs(area - aMode) <= tol*aMode;

    axSpot = axCand(isSpot);
    if isempty(axSpot); return; end

    pos = vertcat(axSpot.Position);
    xc = pos(:,1) + pos(:,3)/2;
    yc = pos(:,2) + pos(:,4)/2;

    % --- map axes to grid lanes
    xRef = linspace(min(xc), max(xc), nSideTotal);
    yRef = linspace(max(yc), min(yc), nSideTotal);   % top -> bottom

    [~, col] = min(abs(xc - xRef), [], 2);           % nAxes x 1
    [~, row] = min(abs(yc - yRef), [], 2);           % nAxes x 1

    figData.spotRow = row;
    figData.spotCol = col;

    % --- get plotWindowTime from any axis that actually has a line
    for k = 1:numel(axSpot)
        ln0 = findall(axSpot(k),'Type','line');
        ln0 = ln0(arrayfun(@(h) isnumeric(h.XData) && numel(h.XData)>20, ln0));
        if ~isempty(ln0)
            figData.plotWindowTime = ln0(1).XData(:)';
            break;
        end
    end

    % --- extract traces, place into full nTotal indexing
    for k = 1:numel(axSpot)
        spotIdx = (row(k)-1)*nSideTotal + col(k);
        if spotIdx < 1 || spotIdx > nTotal
            continue;
        end

        ln = findall(axSpot(k),'Type','line');
        ln = ln(arrayfun(@(h) isnumeric(h.YData) && numel(h.YData)>20, ln));
        if isempty(ln)
            % empty spot -> leave figData.spotTraces{spotIdx} = []
            continue;
        end

        % remove thick summary line (same idea as your current code)
        lw = arrayfun(@(h) h.LineWidth, ln);
        [~,imax] = max(lw);
        ln(imax) = [];

        Y = arrayfun(@(h) h.YData(:)', ln, 'UniformOutput', false);
        L = cellfun(@numel, Y);
        if isempty(L); continue; end
        Lmode = mode(L);
        Y = Y(L==Lmode);
        if isempty(Y); continue; end

        figData.spotTraces{spotIdx} = cell2mat(Y(:));   % reps x time
    end
end


function control = localMakeSyntheticControl(processed, plotWindowTime, ctrlLen)
    pre = processed(plotWindowTime<0);
    if isempty(pre); pre = processed(1:min(20,numel(processed))); end
    pre = pre(:)';
    repN = ceil(ctrlLen/numel(pre));
    control = repmat(pre,1,repN);
    control = control(1:ctrlLen);
    end
    
    function stats = localComputeStatsFromProcessed(processed, control, analysisWindow, vhold, options)
    Fs = options.outputFs;
    peakWindowWidth = round(options.peakWindow*Fs/1000);
    
    stats.baseline.avg = mean(control,'omitnan');
    stats.baseline.std = std(control,0,'omitnan');
    
    trace = processed(analysisWindow);
    stats.response.auc = sum(trace,'omitnan')/Fs;
    stats.baseline.auc = sum(control,'omitnan')/Fs;
    
    [~,maxIdx]=max(trace); [~,minIdx]=min(trace);
    mw=max(1,maxIdx-peakWindowWidth):min(numel(trace),maxIdx+peakWindowWidth);
    nw=max(1,minIdx-peakWindowWidth):min(numel(trace),minIdx+peakWindowWidth);
    stats.response.max = mean(trace(mw),'omitnan');
    stats.response.min = mean(trace(nw),'omitnan');
    stats.response.maxTime = maxIdx*1000/Fs;
    stats.response.minTime = minIdx*1000/Fs;
    
    den = abs(stats.response.max)+abs(stats.response.min);
    if den==0; stats.response.EIindex=0;
    else; stats.response.EIindex=(abs(stats.response.max)-abs(stats.response.min))/den;
    end
    
    if ~isfield(stats.response, options.feature)
        stats.response.(options.feature) = stats.response.auc;
    end
end
