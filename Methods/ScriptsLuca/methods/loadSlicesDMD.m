function loadSlicesDMD(epochs,options)

% Create spots.mat for each epoch (i.e. each search)

arguments
    epochs table %full path of the session or epochs.mat

    options.filterSignal logical = true
    options.filterSweeps logical = true % If true, do not analyze sweeps with different Vhold (included == false)

    options.save logical = true
    options.reload logical = false
    options.reloadCells logical = false
    options.calculateQC logical = false

    options.outputFs double = 10000
    options.timeRange double = [-20,100] % in ms
    options.analysisWindowLength double = 50 % in ms after stim onset
    options.controlWindowLength double = 50 % in ms before stim onset
    options.eventSample % in sample
    options.nArtifactSamples double = 0 % in sample
    options.rcCheckRecoveryWindow double = 100 % in ms
    options.peakWindow double = 1 % in ms around the peak to average

    options.feature = 'auc'
    options.thresholdFactor double = 3 % 3*std

    options.rawDataPath string
    options.saveDataPath string = 'default'
end

%% General setup

today = char(datetime('today','Format','yyyyMMdd'));
dirsplit = split(epochs{1,"Session"},filesep); expName = dirsplit{end};

% Decide reload if session has already been loaded
if options.reload
    options.reloadCells = true;
    resultsFolderName = ['Results_',today];
else
    resultsList = sortrows(struct2cell(dir(fullfile(epochs{1,"Session"},'Results_*')))',3);
    resultsPath = fullfile(resultsList{end,2},resultsList{end,1},'cell1','spots_*.mat');
    if ~isempty(resultsPath)
        disp('Loading stop: spots file found.');
        resultsFolderName = resultsList{end,1};
        if ~options.reloadCells; return; end
    end
end

% Determine data path
% 1. If saveDataPath == 'default', save files in saveDataPath = rawDataPath
% 2. If saveDataPath ~= 'default', saveDataPath should be specified by the user
% 3. "Sessions" in epochs and cells table should be the same as rawDataPath
if strcmp(options.saveDataPath,'default')
    options.rawDataPath = epochs{1,"Session"};
    options.saveDataPath = strcat(options.rawDataPath,filesep,resultsFolderName);
else
    options.rawDataPath = epochs{1,"Session"};
end

% Turn off some warnings
warning('off','MATLAB:unknownObjectNowStruct');
warning('off','MATLAB:table:RowsAddedExistingVars');

%% Select randomSearch epochs
randomSearchIdx = cellfun(@(x) contains(x.cycle, 'randomSearch'), epochs.("Protocol"));
exp = epochs(randomSearchIdx,:); % randomSearchEpochs

%% Initialize cells_DMD.mat

% structure for responseMap
% For each row, response map is a cell with #searches element
% within each element, theres a 3-dim matrix (xRange,yRange,depth)

varTypes = {'string','string','string','double','cell','cell',...
            'cell','cell','cell'};
varNames = {'Session','Animal','Task','Cell','Epochs','Vhold',...
            'Response map','Difference map','Options'};
cells = table('Size',[max(exp{:,'Cell'}),length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

% Initialize cell-level params
cellResponseMap = {}; isResponseMap_cell = {}; cellCurrentMap = {};
cellVhold = []; cellEpochs = {}; searchTotalSpots = {};
nullSpotData = []; preStimData = []; baselineData = [];

%% Iterate & analyze individual epoch

for row = 1:size(exp,1)
    clearvars AD*
    cellResultsPath = strcat(options.saveDataPath,filesep,['cell',num2str(exp{row,'Cell'})]);

    if options.reload
        %% Load basic info
    
        % Initialize fullSearchTable pointer
        % Since I'm iterating sweep in order, the spots are correspond to the
        % order recorded in fullSearchTable
        curSpot = 0; prevDepth = 1; totalSpots = 0;
    
        % Load fullSearchTable
        epoch = exp{row,'Epoch'};
        sweepAcq = exp{row,'Sweep names'}{1};
        epochPath = strcat(exp{:,"Session"}{row},filesep,'cell',num2str(exp{row,"Cell"}));
        try load(fullfile(epochPath,strcat('Epoch',num2str(epoch),'_fullSearchTable.mat')),'fullSearchTable')
        catch
            warning(strcat("fullSearchtable for Epoch ",num2str(epoch), " not saved, skipping this epoch!"));
            continue
        end
    
        % Read patching info csv file
        csvOpts = detectImportOptions(fullfile(epochPath,'InfoPatching.xlsx'));
        csvOpts.SelectedVariableNames = 1:9;
        csvOpts.VariableNamesRange = 25;
        info = readtable(fullfile(epochPath,'InfoPatching.xlsx'),csvOpts);
        info = rmmissing(info,DataVariables="acq_");
    
        %% Initialize spots table
        
        varTypes = {'string','string','string','double','double','double',...
                    'double','double','string','cell','cell',...
                    'cell','cell',...
                    'double','double','double',...
                    'cell'};
        varNames = {'Session','Animal','Task','Epoch','Cell','Vhold',...
                    'Depth','Repetition','Sweep','Location','Protocol',...
                    'Response','Stats',...
                    'Rs','Rm','Cm',...
                    'Options'};
        spots = table('Size',[size(exp,1),length(varNames)],...
            'VariableTypes',varTypes,'VariableNames',varNames);
    
        %% Load individual sweeps
        for k = 1:length(sweepAcq)
            % Load sweep traces (.data)
            disp(['Loading ',sweepAcq{k},'.mat for epoch ',num2str(epoch),' of cell ',num2str(exp{row,"Cell"})]);
            try load(fullfile(epochPath,strcat(sweepAcq{k},'.mat'))); 
            catch
                warning(strcat("Sweep ",num2str(sweepAcq{k}), " not saved, skipping this sweep!"));
                continue
            end
    
            % Extract raw trace
            raw_trace = eval([sweepAcq{k},'.data']);
            headerString = eval([sweepAcq{k},'.UserData.headerString']);
            
            % Extract experiment protocol from header string
            protocol = getCellProtocol(headerString,...
                                       outputFs=options.outputFs,...
                                       rcCheckRecoveryWindow=options.rcCheckRecoveryWindow);
    
            % Extract quality metrics from header string
            qc = getCellQC(raw_trace,calculate=options.calculateQC,headerString=headerString);
            Rs = qc.Rs; Rm = qc.Rm; Cm = qc.Cm;
    
            % Find corresponding spot from fullSearchTable
            sweepSpots = fullSearchTable(curSpot+1 : curSpot+protocol.numPulses,:);
            % curSpot = curSpot+protocol.numPulses;
            % Check whether selected rows are of the right depth
            if any(sweepSpots.("depth") ~= protocol.depth)
                error('Error: selected spot are of the wrong depth based on headerString!');
            end
    
    
            % Define baseline window: have two windows, one before pulse and one after pulse
            rcCheckOnset = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='delay') * (options.outputFs/1000);
            rcCheckPulseWidth = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='pulseWidth') * (options.outputFs/1000);
            stimDuration = (protocol.numPulses * protocol.isi) * options.outputFs/1000;
            if rcCheckOnset < protocol.stimOnset(1)
                rcCheckEnd = rcCheckOnset + rcCheckPulseWidth + (options.rcCheckRecoveryWindow*(options.outputFs/1000));
                preStimWindow = rcCheckEnd : (protocol.stimOnset(1)-1);
                postStimWindow = (protocol.stimOnset(1) + stimDuration):length(raw_trace);
                baselineWindow = [preStimWindow,postStimWindow];
            else
                preStimWindow = 1:(protocol.stimOnset(1)-1);
                postStimWindow = (protocol.stimOnset(1) + stimDuration):rcCheckOnset;
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

    
            % Process trace: mean-subtracted, optional LP
            % Mean subtraction
            mean_subtracted = raw_trace - mean(raw_trace(baselineWindow));
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
            % Find in .xlsx file first, if can't find, use the heuristic that
            % cells with negative leak current Vhold=-70, otherwise Vhold = 10;
            acqsplit = split(sweepAcq{k},'_'); acqNum = str2double(acqsplit{end});
            vhold = info{info.acq_ == acqNum,'holding'};
            if isempty(vhold)
                if stats.baseline.avg < 0; vhold = -70;
                else; vhold = 10; end
            end
    
            % Initialize matrix for noise analysis later
            baselineData = [baselineData, processed_trace(baselineWindow)];
    
            %% Loop through all spots
            for s = 1:protocol.numPulses
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
                if find(sum(abs(responses.processed(analysisWindow)) >= response_threshold)>50)
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
                stats.response.EIindex = abs(stats.response.max)-abs(stats.response.min) / abs(stats.response.max)+abs(stats.response.min);
                stats.baseline.EIindex = abs(stats.baseline.max)-abs(stats.baseline.min) / abs(stats.baseline.max)+abs(stats.baseline.min);
    
                %% Store data for the current spot
                curSpot = curSpot + 1;
                totalSpots = totalSpots + 1;
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
                spots{curSpot,'Rs'} = Rs;
                spots{curSpot,'Rm'} = Rm;
                spots{curSpot,'Cm'} = Cm;
                spots{curSpot,'Options'} = {options};
            end

            % Save spots.mat for a specific depth
            % If not, spots are too big and matlab can't load later
            if k == length(sweepAcq) || prevDepth ~= protocol.depth
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
    end
    
    %% Loop through all depth and add search results to cells.mat

    disp('Ongoing: adding current search to cells_DMD.mat'); clearvars AD*

    % Get a list of spots.mat file of a given cell and epoch
    filename = strcat('spots_cell',num2str(exp{row,'Cell'}),'_epoch',num2str(exp{row,'Epoch'}));
    depthfilename = strcat(filename,'_depth*.mat');
    spotsList = dir(fullfile(cellResultsPath,depthfilename));
    disp(['Ongoing: adding search epoch ', filename,' to cells_DMD.mat']);

    % Initialize params
    nDepth = length(spotsList);
    searchResponseMap = zeros(608,684,nDepth);
    isResponseMap_search = zeros(608,684,nDepth);
    searchCurrentMap = cell(nDepth,1);

    % Loop through depth
    for depthIdx = 1:nDepth
        % Initialization
        depthFilePath = fullfile(spotsList(depthIdx).folder,spotsList(depthIdx).name);
        load(depthFilePath,'spotsAtDepth');
        % Get search depth
        d = spotsAtDepth{1,'Depth'};
        
        depthResponseMap = zeros(608,684);
        isResponseMap_depth = zeros(608,684);
        depthCurrentMap = zeros(4^d,size(spotsAtDepth{1,'Response'}{1}.processed,2));
    
        % Build response map
        for s = 1:size(spotsAtDepth,1)
            % Get response value
            location = spotsAtDepth{s,'Location'}{1};
            xRange = location(1):location(2);
            yRange = location(3):location(4);
            originalValue = mean(depthResponseMap(xRange,yRange),'all');
            newValue = spotsAtDepth{s,'Stats'}{1}.response.(options.feature);
    
            % Add to depthResponseMap
            if originalValue==0; depthResponseMap(xRange,yRange) = newValue;
            else; depthResponseMap(xRange,yRange) = mean([originalValue,newValue]);
            end

            % Get isResponse value
            originalValue_isResponse = mode(isResponseMap_depth(xRange,yRange),'all');
            newValue_isResponse = spotsAtDepth{s,'Response'}{1}.isResponse;
    
            % Add to isResponseMap
            if originalValue_isResponse==0; isResponseMap_depth(xRange,yRange) = newValue_isResponse;
            else; isResponseMap_depth(xRange,yRange) = originalValue_isResponse || newValue_isResponse;
            end

            % Get response trace
            square_width = 608/2^d; square_height = 684/2^d;
            x_index = floor(location(1) / square_width);
            y_index = floor(location(3) / square_height);
            spotIdx = y_index + x_index * 2^d + 1;
            trace = spotsAtDepth{s,'Response'}{1}.processed;

            % Add to current map
            originalTrace = depthCurrentMap(spotIdx,:);
            if mean(originalTrace)==0; depthCurrentMap(spotIdx,:) = trace;
            else; depthCurrentMap(spotIdx,:) = (originalTrace + trace)/2;
            end
        end        

        % Save depthResponseMap
        searchResponseMap(:,:,d) = depthResponseMap;
        isResponseMap_search(:,:,d) = isResponseMap_depth;
        searchCurrentMap{d} = depthCurrentMap;
    end

    % Save response map
    cellResponseMap{end+1} = searchResponseMap;
    isResponseMap_cell{end+1} = isResponseMap_search;
    cellCurrentMap{end+1} = searchCurrentMap;

    % Save to cell
    curCell = spotsAtDepth{1,'Cell'};
    cellVhold = [cellVhold; spotsAtDepth{1,'Vhold'}];
    cellEpochs{end+1} = filename;
    searchTotalSpots{end+1} = totalSpots;

    if row == size(exp,1) || curCell ~= exp{row+1,'Cell'}
        responseMap.responseMap = cellResponseMap';
        responseMap.isResponseMap = isResponseMap_cell';
        responseMap.currentMap = cellCurrentMap';
        options.searchTotalSpots = searchTotalSpots;

        cells{curCell,'Session'} = options.saveDataPath;
        cells{curCell,'Animal'} = spotsAtDepth{1,'Animal'};
        cells{curCell,'Task'} = spotsAtDepth{1,'Task'};
        cells{curCell,'Cell'} = curCell;
        cells{curCell,'Epochs'} = {cellEpochs'};
        cells{curCell,'Vhold'} = {cellVhold};
        cells{curCell,'Response map'} = {responseMap};
        cells{curCell,'Options'} = {options};

        % Save noise data for this cell
        filename = strcat('noise_cell',num2str(curCell));
        save(fullfile(cellResultsPath,filename),'nullSpotData','preStimData','baselineData','-v7.3');
        disp(strcat("Saved noise data for cell: ",num2str(curCell)));

        % Update prevCell and reset cellResponseMap
        cellVhold = []; cellEpochs = {}; searchTotalSpots = {};
        cellResponseMap = {}; isResponseMap_cell = {}; cellCurrentMap = {};
        nullSpotData = []; preStimData = []; baselineData = [];
    end
end

%% Cell specific analysis
    
% The main goal here is to find common spots and calculate there difference
% in AUC or other feature

for c = 1:size(cells,1)
    cellData = cells(cells.Cell == c,:);
    searchPerCell = length(cellData.Vhold{1});
    diff_rmap_cell = {}; common_isResponse_cell = {}; 
    diff_pairs = {}; diff_vholds = [];

    %% Build a noise model for a given cell
    % pd1 is a gaussian model fit to null spots
    % pd2 is a gaussian model fit to -30ms to -5ms period prior to
    % photostimulation onset
    % pd3 is a gaussian model fit to pooled datapoints of null spots,
    % prestim period, and baseline period.

    disp(['Ongoing: building noise model for cell ',num2str(num2str(c))]);

    % Load noise data
    cellResultsPath = strcat(options.saveDataPath,filesep,['cell',num2str(c)]);
    filename = strcat('noise_cell',num2str(c),'.mat');
    load(fullfile(cellResultsPath,filename));

    % Initialize noise summary plot
    initializeFig(0.75,1); tiledlayout(1,3);
    
    % pd1
    nexttile; 
    outlier1 = isoutlier(nullSpotData);
    noise_nullSpot = fitdist(nullSpotData(~outlier1)','Normal');
    histogram(nullSpotData,'Normalization','pdf'); hold on
    histogram(nullSpotData(~outlier1),'Normalization','pdf');  hold on
    xrange = [min(nullSpotData),max(nullSpotData)];
    plot(linspace(xrange(1),xrange(2)),pdf(noise_nullSpot,linspace(xrange(1),xrange(2))),'LineWidth',2);
    title('Null spot response');
    
    % pd2
    nexttile; 
    outlier2 = isoutlier(preStimData);
    noise_preStim = fitdist(preStimData(~outlier2)','Normal');
    xrange=[min(preStimData),max(preStimData)];
    histogram(preStimData,'Normalization','pdf'); hold on
    histogram(preStimData(~outlier2),'Normalization','pdf'); hold on
    plot(linspace(xrange(1),xrange(2)),pdf(noise_preStim,linspace(xrange(1),xrange(2))),'LineWidth',2);
    title('Pre-stim period');
    
    % pd3
    nexttile;
    outlier3 = isoutlier(baselineData);
    allNullData = [nullSpotData(~outlier1),preStimData(~outlier2),baselineData(~outlier3)];
    noise_all = fitdist(allNullData','Normal');
    xrange=[min(allNullData),max(allNullData)];
    histogram(allNullData,'Normalization','pdf'); hold on
    plot(linspace(xrange(1),xrange(2)),pdf(noise_all,linspace(xrange(1),xrange(2))),'LineWidth',2);
    title('Baseline + Pre-stim + Null spot response');

    % Set threshold based on the noise model (2*standard devation of the
    % symmetric noise)
    options.Ethres = -2 * noise_all.sigma;
    options.Ithres = 2 * noise_all.sigma;
    
    % Save noise model
    if options.save
        saveFigures(gcf,['cell',num2str(c),'_noise'],cellResultsPath);
        save(fullfile(cellResultsPath,filename),'noise_nullSpot','noise_preStim','noise_all','allNullData','-append');
        disp(strcat("Saved noise model for cell: ",num2str(c)));
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

            % Check max common depth
            maxCommonDepth = min([size(search1_rmap,3),size(search2_rmap,3)]);
            search1_rmap = search1_rmap(:,:,1:maxCommonDepth);
            search2_rmap = search2_rmap(:,:,1:maxCommonDepth);
            search1_isResponse = search1_isResponse(:,:,1:maxCommonDepth);
            search2_isResponse = search2_isResponse(:,:,1:maxCommonDepth);

            % Check whether Vhold are the same
            diffVhold = cellData.Vhold{1}(search1) ~= cellData.Vhold{1}(search2);

            % Initialize difference map
            diff_rmap = zeros(size(search1_rmap));
            common_isResponse = zeros(size(search1_rmap));


            % Calculate difference map for each depth
            for d = 1:maxCommonDepth
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
            diff_rmap_cell{end+1} = diff_rmap;
            common_isResponse_cell{end+1} = common_isResponse;
            diff_pairs{end+1} = [search1, search2];
            diff_vholds = [diff_vholds; diffVhold];
        end 
    end

    % Save in cells.mat
    diff.response = diff_rmap_cell';
    diff.commonSpots = common_isResponse_cell';
    diff.pair = diff_pairs';
    diff.diffVhold = diff_vholds;
    cells{c,'Difference map'} = {diff};
    cells{c,'Options'} = {options};
end

%% Save cells.mat

if options.save
    save(strcat(options.saveDataPath,filesep,'cells_DMD_',expName),'cells','-v7.3');
    disp(strcat("Saved cells_DMD.mat: ",expName));
end
close all

end