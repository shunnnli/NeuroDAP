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

    % Initialize cells table
    varTypes = {'string','string','string','double','cell',...
                'cell','cell','cell','cell',...
                'cell','cell','cell',...
                'cell','cell','cell',...
                'cell','cell'};
    varNames = {'Session','Animal','Task','Cell','Vhold',...
                'Included','Sweep names','Raw sweeps','Processed sweeps',...
                'Protocol','Response','Stats',...
                'Rs','Rm','Cm',...
                'VholdInfo','Options'};
    cells = table('Size',[length(cellList),length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);

    for c = 1:size(cellList,1)
        cellEpochs = animalEpochs(animalEpochs.Cell == cellList(c),:);
        wholeFieldIdx = ~cellfun(@(x) contains(x.cycle,'randomSearch'),cellEpochs.Protocol);
        cellEpochs = cellEpochs(wholeFieldIdx,:);
        options.rawDataPath = cellEpochs{1,'Options'}{1}.rawDataPath;

        % Initialize structure to save cell data
        respnoses.raw = {}; respnoses.processed = {}; respnoses.isResponse = {};
        stats.response.auc = {}; stats.response.peak = {}; stats.response.peakTime = {};
        stats.baseline.auc = {}; stats.baseline.peak = {}; stats.baseline.peakTime = {};
        stats.baseline.avg = {}; stats.baseline.std = {};

        for e = 1:size(cellEpochs,1)
            % Extract epoch info
            protocol = cellEpochs{e,'Protocol'}{1};
            raw_trace = cellEpochs{e,'Raw sweeps'}{1};
            processed_trace = cellEpochs{e,'Processed sweeps'}{1};
            vhold = cellEpochs{e,'Vhold'};
            
            % Define time window for plotting
            options.eventSample = protocol.stimOnset;
            timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
            timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
            plotWindowLength = timeRangeEndSample(1) - timeRangeStartSample(1) + 1;
            baselineWindow = cellEpochs{e,'Options'}{1}.baselineWindow;

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

            % Initialize response matrix
            responses_raw = zeros(protocol.numPulses, plotWindowLength);
            responses_processed = zeros(protocol.numPulses, plotWindowLength);
            responses_isResponse = zeros(protocol.numPulses,1);
            stats_baseline_peak = zeros(protocol.numPulses,1);
            stats_baseline_auc = zeros(protocol.numPulses,1);
            stats_baseline_peakTime = zeros(protocol.numPulses,1);
            stats_response_peak = zeros(protocol.numPulses,1);
            stats_response_auc = zeros(protocol.numPulses,1);
            stats_response_peakTime = zeros(protocol.numPulses,1);
            
            % Find baseline stats
            stats_baseline_avg = mean(processed_trace(baselineWindow));
            stats_baseline_std = std(processed_trace(baselineWindow));

            for s = 1:protocol.numPulses
                % Save spot responses
                responses_raw(s,:) = raw_trace(timeRangeStartSample(s):timeRangeEndSample(s));
                responses_processed(s,:) = processed_trace(timeRangeStartSample(s):timeRangeEndSample(s));
    
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
            end

            % Save epoch responses and stats to cellResponses
            respnoses.raw{end+1} = responses_raw; 
            respnoses.processed{end+1} = responses_processed; 
            respnoses.isResponse{end+1} = responses_isResponse; 
            stats.response.auc{end+1} = stats_response_auc;
            stats.response.peak{end+1} = stats_response_peak;
            stats.response.peakTime{end+1} = stats_response_peakTime;
            stats.baseline.auc{end+1} = stats_baseline_auc;
            stats.baseline.peak{end+1} = stats_baseline_peak;
            stats.baseline.peakTime{end+1} = stats_baseline_peakTime;
            stats.baseline.avg{end+1} = stats_baseline_avg;
            stats.baseline.std{end+1} = stats_baseline_std;
        end

        % Get response statistics
        EPSCrows = cellEpochs{:,'Vhold'} < -30;
        IPSCrows = cellEpochs{:,'Vhold'} >= -30;
        stats.EPSC.peakAvg = mean(cellfun(@mean, stats.response.peak(EPSCrows)),'all');
        stats.IPSC.peakAvg = mean(cellfun(@mean, stats.response.peak(IPSCrows)),'all');
        stats.EPSC.aucAvg = mean(cellfun(@mean, stats.response.auc(EPSCrows)),'all');
        stats.IPSC.aucAvg = mean(cellfun(@mean, stats.response.auc(IPSCrows)),'all');

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
        cells{c,'Response'} = {respnoses};
        cells{c,'Stats'} = {stats};
        cells{c,'Rs'} = {cellEpochs.('Rs')};
        cells{c,'Rm'} = {cellEpochs.('Rm')};
        cells{c,'Cm'} = {cellEpochs.('Cm')};
        cells{c,'VholdInfo'} = {cellEpochs.('VholdInfo')};
        cells{c,'Options'} = {options};
    end

    if options.save
        save(strcat(animalEpochs{1,'Session'},filesep,'cells_',expName),'cells');
        disp(strcat("Created cells.mat: cell_",expName));
    end
    allCells = [allCells; cells];
end

end
