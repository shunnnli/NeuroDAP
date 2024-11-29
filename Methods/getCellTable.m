function allCells = getCellTable(epochs,options)

% Combine cell-specific statistics from epochs

arguments
    epochs table
    options.save logical = true
    options.saveDataPath

    options.outputFs double = 10000
    options.timeRange double = [-20,100] % in ms
    options.analysisWindowLength double = 50 % in ms after the stim onset
    options.controlWindowLength double = 50 % in ms before the stim onset
    options.eventSample % in sample
    options.nArtifactSamples double = 0 % in sample
    options.peakWindow double = 2; % in ms around the peak to average

    options.calculate logical = false % if true, re-calculate statistics
end

%% cells table

% Each row is a cell
% Vhold: every epochs of this cell: nEpochs x 1
% Included: all the included sweeps within each epochs: cell, nEpochs x nSweeps
% Sweep names: matlab names for the sweep: cell, nEpochs x nSweeps
% Raw/processed sweeps: cell, each row contains raw sweeps for this epoch
% Peaks/AUCs: nEpochs x nSweeps
% Rs/Rm/Cm: nEpochs x nSweeps

% If epochs contains only one animal, create cells for that animal
% If epochs contains multiple animals:
    % 1. Create cells for each animal, save them into there corresponding
    % folder
    % 2. Concatenate into allCells and return

%% Create cells table by animal

animalList = unique(epochs{:,"Animal"});
allCells = [];

for i = 1:length(animalList)
    % Find epochs for an animal
    animalEpochs = epochs(epochs.Animal == animalList(i),:);
    cellList = unique(animalEpochs{:,"Cell"});
    dirsplit = split(animalEpochs{1,"Session"},filesep); expName = dirsplit{end};

    disp(strcat("Creating cells table for: ",expName));

    varTypes = {'string','string','string','double','cell',...
                'cell','cell','cell','cell',...
                'cell','cell',...
                'cell',...
                'cell','cell'};
    varNames = {'Session','Animal','Task','Cell','Vhold',...
                'Included','Sweep names','Raw sweeps','Processed sweeps',...
                'Protocol','Stats',...
                'QC',...
                'VholdInfo','Options'};
    cells = table('Size',[length(cellList),length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);

    for c = 1:size(cellList,1)
        cellEpochs = animalEpochs(animalEpochs.Cell == cellList(c),:);
        cellProtocols = cellfun(@(v) v.cycle, cellEpochs.Protocol, UniformOutput=false);
        wholeFieldIdx = ~cellfun(@(x) any(contains(x,'randomSearch')),cellProtocols);
        cellEpochs = cellEpochs(wholeFieldIdx,:);
        options.rawDataPath = cellEpochs{1,'Options'}{1}.rawDataPath;

        % Get rows
        EPSCrows = cellEpochs{:,'Vhold'} < -50;
        IPSCrows = cellEpochs{:,'Vhold'} > -10;

        % Get included
        if sum(EPSCrows)>0; EPSCincluded = cellEpochs{EPSCrows,'Included'}{1};
        else; EPSCincluded = []; end
        if sum(IPSCrows)>0; IPSCincluded = cellEpochs{IPSCrows,'Included'}{1};
        else; IPSCincluded = []; end
        if sum(EPSCincluded) == 0; EPSCincluded = ones(length(EPSCincluded),1); end
        if sum(IPSCincluded) == 0; IPSCincluded = ones(length(IPSCincluded),1); end
        EPSCincluded = logical(EPSCincluded);
        IPSCincluded = logical(IPSCincluded);

        % Initialize structure to save cell data
        responses.raw = {}; responses.processed = {}; responses.isResponse = {};
        if options.calculate
            stats.response.auc = {}; stats.response.peak = {}; stats.response.peakTime = {};
            stats.baseline.auc = {}; stats.baseline.peak = {}; stats.baseline.peakTime = {};
            stats.baseline.avg = {}; stats.baseline.std = {};
        else
            epochStats = cellEpochs{:,'Stats'};
            stats = mergeStructs(epochStats,combine=false);
        end

        % Loop over each epoch of this cell
        for epoch = 1:size(cellEpochs,1)
            % Extract epoch info
            protocol = cellEpochs{epoch,'Protocol'}{1};
            raw_trace = cellEpochs{epoch,'Raw sweeps'}{1};
            processed_trace = cellEpochs{epoch,'Processed sweeps'}{1};
            vhold = cellEpochs{epoch,'Vhold'};

            protocol.stimOnset = protocol.stimOnset(1);
            protocol.numPulses = protocol.numPulses(1);

            % Define time window for plotting
            options.eventSample = protocol.stimOnset;
            timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
            timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
            plotWindowLength = timeRangeEndSample(1) - timeRangeStartSample(1) + 1;
            baselineWindow = cellEpochs{epoch,'Options'}{1}.baselineWindow;

            % Initialize response matrix
            responses_raw = zeros(protocol.numPulses, plotWindowLength);
            responses_processed = zeros(protocol.numPulses, plotWindowLength);
            responses_isResponse = zeros(protocol.numPulses,1);

            if options.calculate
                % Define control window: 50ms before each spot stim onset
                controlWindowSamples = options.controlWindowLength * options.outputFs/1000;
    
                % Define analysis window: 50ms after each spot stim onset
                analysisWindowSamples = options.analysisWindowLength * options.outputFs/1000;
                eventSample = abs(options.outputFs*options.timeRange(1)/1000);
                analysisWindow = eventSample:eventSample+analysisWindowSamples;
        
                % Define time window for peak window analysis
                peakWindowWidth = (options.peakWindow/(2*1000)) * options.outputFs;
    
                % Store in options
                options.controlWindowSamples = controlWindowSamples;
                options.analysisWindowSamples = analysisWindowSamples;
                options.baselineWindow = baselineWindow;
                options.analysisWindow = analysisWindow;
                options.peakWindowWidth = peakWindowWidth;
    
                % Initialize stats matrix
                stats_baseline_peak = zeros(protocol.numPulses,1);
                stats_baseline_auc = zeros(protocol.numPulses,1);
                stats_baseline_peakTime = zeros(protocol.numPulses,1);
                stats_response_peak = zeros(protocol.numPulses,1);
                stats_response_auc = zeros(protocol.numPulses,1);
                stats_response_peakTime = zeros(protocol.numPulses,1);
                
                % Find baseline stats
                stats_baseline_avg = mean(processed_trace(baselineWindow));
                stats_baseline_std = std(processed_trace(baselineWindow));
            end

            for s = 1:protocol.numPulses
                % Save spot responses
                responses_raw(s,:) = raw_trace(timeRangeStartSample(s):timeRangeEndSample(s));
                responses_processed(s,:) = processed_trace(timeRangeStartSample(s):timeRangeEndSample(s));

                % Save epoch responses
                responses.raw{end+1} = responses_raw; 
                responses.processed{end+1} = responses_processed; 
                responses.isResponse{end+1} = responses_isResponse; 

                if options.calculate
                    % Define control window
                    controlWindow = protocol.stimOnset(s)-controlWindowSamples-1 : protocol.stimOnset(s)-1;
            
                    % Find auc
                    stats_response_auc(s) = sum(responses_processed(analysisWindow)) / options.outputFs;
                    stats_baseline_auc(s) = sum(processed_trace(controlWindow)) / options.outputFs;
        
                    % Find peak for stim response
                    trace = responses_processed(analysisWindow);
                    if vhold < -30; [~,peakIdx] = max(-trace);
                    else; [~,peakIdx] = max(trace); end
                    peakWindowStart = max(1,peakIdx-peakWindowWidth);
                    peakWindowEnd = min(peakIdx+peakWindowWidth,length(trace));
                    if vhold < -30; stats_response_peak(s) = -mean(trace(peakWindowStart:peakWindowEnd));
                    else; stats_response_peak(s) = mean(trace(peakWindowStart:peakWindowEnd)); end
                    stats_response_peakTime(s) = peakIdx * 1000/options.outputFs;
        
                    % Find peak for control baseline
                    trace = processed_trace(controlWindow);
                    if vhold < -30; [~,peakIdx] = max(-trace);
                    else; [~,peakIdx] = max(trace); end
                    peakWindowStart = max(1,peakIdx-peakWindowWidth);
                    peakWindowEnd = min(peakIdx+peakWindowWidth,length(trace));
                    if vhold < -30; stats_baseline_peak(s) = -mean(trace(peakWindowStart:peakWindowEnd));
                    else; stats_baseline_peak(s) = mean(trace(peakWindowStart:peakWindowEnd)); end
                    stats_baseline_peakTime(s) = peakIdx * 1000/options.outputFs;
        
                    % Determine whether there is a response across threshold
                    response_threshold =  5 * stats_baseline_std;
                    if abs(stats_response_peak) >= response_threshold; responses_isResponse = true;
                    else; responses_isResponse = false; end

                    % Save stats
                    stats.response.auc{end+1} = stats_response_auc;
                    stats.response.peak{end+1} = stats_response_peak;
                    stats.response.peakTime{end+1} = stats_response_peakTime;
                    stats.baseline.auc{end+1} = stats_baseline_auc;
                    stats.baseline.peak{end+1} = stats_baseline_peak;
                    stats.baseline.peakTime{end+1} = stats_baseline_peakTime;
                    stats.baseline.avg{end+1} = stats_baseline_avg;
                    stats.baseline.std{end+1} = stats_baseline_std;
                end
            end
        end

        % Get response statistics
        if options.calculate
            % Not finished: did not consider included
            stats.summary.EPSC.peakAvg = mean(cellfun(@mean, stats.response.peak(EPSCrows)),'all');
            stats.summary.IPSC.peakAvg = mean(cellfun(@mean, stats.response.peak(IPSCrows)),'all');
            stats.summary.EPSC.aucAvg = mean(cellfun(@mean, stats.response.auc(EPSCrows)),'all');
            stats.summary.IPSC.aucAvg = mean(cellfun(@mean, stats.response.auc(IPSCrows)),'all');
        else      
            if sum(EPSCrows) > 1
                warning('More than one epoch of EPSC included. Concatenate all of them.');
            elseif sum(IPSCrows) > 1
                warning('More than one epoch of IPSC included. Concatenate all of them.');
            end
            responseStats = cell2mat(stats.response.min(EPSCrows));
            stats.summary.EPSC.peakAvg = mean(responseStats(EPSCincluded),'all');
            responseStats = cell2mat(stats.response.max(IPSCrows));
            stats.summary.IPSC.peakAvg = mean(responseStats(IPSCincluded),'all');
            responseStats = cell2mat(stats.response.auc(EPSCrows));
            stats.summary.EPSC.aucAvg = mean(responseStats(EPSCincluded),'all');
            responseStats = cell2mat(stats.response.auc(IPSCrows));
            stats.summary.IPSC.aucAvg = mean(responseStats(IPSCincluded),'all');
        end

        cells{c,'Session'} = cellEpochs{1,'Session'};
        cells{c,'Animal'} = cellEpochs{1,'Animal'};
        cells{c,'Task'} = cellEpochs{1,'Task'};
        cells{c,'Cell'} = cellList(c);
        cells{c,'Vhold'} = {cellEpochs{:,'Vhold'}};
        cells{c,'Included'} = {cellEpochs.('Included')};
        cells{c,'Sweep names'} = {cellEpochs.('Sweep names')};
        cells{c,'Raw sweeps'} = {cellEpochs.('Raw sweeps')};
        cells{c,'Processed sweeps'} = {cellEpochs.('Processed sweeps')};
        cells{c,'Protocol'} = {cellEpochs.('Protocol')};
        % cells{c,'Response'} = {responses};
        cells{c,'Stats'} = {stats};
        cells{c,'QC'} = {cellEpochs.('QC')};
        cells{c,'VholdInfo'} = {cellEpochs.('VholdInfo')};
        cells{c,'Options'} = {cellEpochs.('Options')};
    end

    if options.save
        today = char(datetime('today','Format','yyyyMMdd')); 
        save(strcat(animalEpochs{1,'Session'},filesep,'cells_',today),'cells');
    end
    allCells = [allCells; cells];
end

disp("Finished: created cells table");

end
