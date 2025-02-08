function varargout = loadSlices(exp,options)


arguments
    exp  %full path of the session or epochs.mat

    options.filterSignal logical = false
    options.filterSweeps logical = true % If true, do not analyze sweeps with different Vhold (included == false)
    
    options.animal string
    options.task string = 'random'

    options.outputFs double = 10000
    options.timeRange double = [-20,100] % in ms
    options.eventSample % in sample
    options.nArtifactSamples double = 0 % in sample
    options.analysisWindowLength double = 30 % in ms after stim onset
    options.controlWindowLength double = 30 % in ms
    options.rcCheckRecoveryWindow double = 100 % in ms
    options.peakWindow double = 2 % in ms around the peak to average

    options.defaultStimOnset double

    options.rawDataPath string
    options.saveDataPath string = 'default'

    % edit: properties in epochs to recalculate
    % (currently mostly for included)
    options.edit string
    options.reload logical = false % if true, reload epochs.mat AND cells.mat
    options.reloadCell logical = false % if true, only reload cells.mat
    options.getCellTable logical = true % false for DMD setup
    options.calculateQC logical = true
    options.save logical = true

    options.QCThreshold struct
    options.plot logical = true % plot epoch analysis summary
end

%% General setup

today = char(datetime('today','Format','yyyyMMdd')); 

% If animal is not provided, create a place holder
if ~isfield(options,'animal'); options.animal = "NA"; end

% Determine processing type
% If exp is epochs.mat, then only analyze rows/neurons in epochs.mat 
% (reprocess signal for postQC data)
if ischar(exp); createNew = true;
elseif istable(exp); createNew = false; end

% Determine data path
% 1. If saveDataPath == 'default', save files in saveDataPath = rawDataPath
% 2. If saveDataPath ~= 'default', saveDataPath should be specified by the user
% 3. "Sessions" in epochs and cells table should be the same as rawDataPath
if strcmp(options.saveDataPath,'default')
    if createNew; options.rawDataPath = exp;
    else; options.rawDataPath = exp{1,"Session"}; end
    options.saveDataPath = strcat(options.rawDataPath,filesep,'Epochs-',today);
    mkdir(options.saveDataPath);
    options.sessionPath = options.rawDataPath;
else
    if createNew; options.rawDataPath = exp;
    else; options.rawDataPath = exp{1,"Session"}; end
    options.sessionPath = options.rawDataPath;
end

% Load existing table if reload is false
dirsplit = split(options.sessionPath,filesep); expName = dirsplit{end};
epochsFiles = sortrows(struct2cell(dir(fullfile(options.sessionPath,"epochs_*.mat")))',[1 3]);
cellsFiles = sortrows(struct2cell(dir(fullfile(options.sessionPath,"cells_*.mat")))',[1 3]);
if ~options.reload && (~isempty(epochsFiles) || istable(exp))
    disp(['Loading stop: epochs file found for ',expName]);
    load(strcat(options.sessionPath,filesep,epochsFiles{end,1}));
    if ~exist('epochs','var')
        if istable(exp); epochs = exp;
        elseif exist('epochs_old','var') && istable(epochs_old); epochs = epochs_old; 
        end
    end
    if options.reloadCell
        varargout{1} = epochs;
        options.getCellTable = true;
        options.reload = false;
    else
        varargout{1} = epochs;
        isCellsTable = true;
        if isempty(cellsFiles)
            warning('cells.mat not found. Not loading it!');
            isCellsTable = false;
        else
            load(strcat(exp,filesep,cellsFiles{end,1}));          
        end
        if isCellsTable; varargout{2} = cells; 
        else; varargout{2} = []; end
        return
    end   
else
    options.reload = true;
end

% Turn off some warnings
warning('off','MATLAB:unknownObjectNowStruct');
warning('off','MATLAB:table:RowsAddedExistingVars');

% QC setup
if ~isfield(options,'QCThreshold')
    options.QCThreshold.include = {};
    options.QCThreshold.Rs = 30;
    options.QCThreshold.Verror = 10;
    options.QCThreshold.Ibaseline = -300;
    options.QCThreshold.Ibaseline_std = 20;
else
    if isfield(options.QCThreshold,'include') && isstring(options.QCThreshold.include)
        options.QCThreshold.include = {options.QCThreshold.include};
    end
    if ~isfield(options.QCThreshold,'Rs')
        options.QCThreshold.Rs = 25;
    end
    if ~isfield(options.QCThreshold,'Verror')
        options.QCThreshold.Verror = 10;
    end
    if ~isfield(options.QCThreshold,'Ibaseline')
        options.QCThreshold.Ibaseline = -300;
    end
    if ~isfield(options.QCThreshold,'Ibaseline_std')
        options.QCThreshold.Ibaseline_std = 20;
    end
end

%% Determine recording rig (ie file structure)

items = dir(fullfile(options.sessionPath,'cell*'));
cellFolders = items([items.isdir]);
% Remove '.' and '..' from the list
cellFolders = cellFolders(~ismember({cellFolders.name}, {'.', '..'}));

if isempty(cellFolders); rig = 'Wengang';
else; rig = 'Paolo'; end

options.rig = rig;

%% Generate epochs.mat

if options.reload
    %% List all epochs
    if strcmp(rig,'Wengang')
        if createNew
            epochList = sortrows(struct2cell(dir(fullfile(exp,['AD0_e*','p1avg.mat'])))',3);
            vholdList = sortrows(struct2cell(dir(fullfile(exp,['AD2_e*','p1avg.mat'])))',3);
        else
            epochList = {}; vholdList = {};
            for i = 1:size(exp,1)
                epochPath = strcat(filesep,'AD0_e',num2str(exp{i,"Epoch"}),'p1avg.mat');
                epochList = [epochList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, epochPath))))'];
                vholdPath = strcat(filesep,'AD2_e',num2str(exp{i,"Epoch"}),'p1avg.mat');
                vholdList = [vholdList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, vholdPath))))'];
            end
        end
    
        % other common params that we know will not change
        cellid = 0; % intialize cellid
    
    else
        if createNew
            epochList = {}; vholdList = {};
            for c = 1:length(cellFolders)
                if contains(cellFolders(c).name,'.mat'); continue; end
                cellPath = strcat(cellFolders(c).folder,filesep,cellFolders(c).name);
                cellEpochList = sortrows(struct2cell(dir(fullfile(cellPath,['AD0_e*','p*avg.mat'])))',3);
                cellEpochList(:,end+1) = num2cell(sscanf(cellFolders(c).name,'cell%d'),[1 2]); % store cell number as the last column
                epochList = [epochList; cellEpochList];
                cellVholdList = sortrows(struct2cell(dir(fullfile(cellPath,['AD2_e*','p*avg.mat'])))',3);
                cellVholdList(:,end+1) = num2cell(sscanf(cellFolders(c).name,'cell%d'),[1 2]); % store cell number as the last column
                vholdList = [vholdList; cellVholdList];
            end
        else
            epochList = {}; vholdList = {};
            for i = 1:size(exp,1)
                epochPath = strcat(filesep,'cell',num2str(exp{i,"Cell"}),filesep,'AD0_e',num2str(exp{i,"Epoch"}),'p*avg.mat');
                epochList = [epochList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, epochPath))))'];
                vholdPath = strcat(filesep,'cell',num2str(exp{i,"Cell"}),filesep,'AD2_e',num2str(exp{i,"Epoch"}),'p*avg.mat');
                vholdList = [vholdList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, vholdPath))))'];
            end
        end
    end
    
    %% Initialize epochs table
    
    % Notes:
    % 1. Vhold is calculated based on the heuristic that if the leak current is
    % negative, the Vhold = -70. Otherwise Vhold = 10
    % 2. Vhold epoch/sweep mean/trace are extracted from AD2
    
    varTypes = {'string','string','string','double','double','double',...
                'cell','cell','cell','cell',...
                'cell',...
                'cell',...
                'cell',...
                'cell','cell'};
    varNames = {'Session','Animal','Task','Epoch','Cell','Vhold',...
                'Included','Sweep names','Raw sweeps','Processed sweeps',...
                'Protocol',...
                'Stats',...
                'QC',...
                'VholdInfo','Options'};
    epochs = table('Size',[length(epochList),length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    
    % Check whether there's AD2 to record Vhold
    withVhold = ~isempty(dir(fullfile(epochList{1,2},"AD2*.mat")));
    withVholdAvg = length(epochList)==length(vholdList);
    vholdInfo.vholdEpochMean = nan;
    vholdInfo.vholdSweepsMean = nan;
    vholdInfo.vholdEpochTrace = nan;
    vholdInfo.vholdSweepsTrace = nan;
    
    %% Iterate & analyze individual epoch
    
    for row = 1:size(epochList,1)
        clearvars AD*
    
        % Load epoch file to find individual sweep.mat
        load(fullfile(epochList{row,2},epochList{row,1}));
        namesplit = strsplit(epochList{row,1},{'e','p','avg.mat'}); 
        epoch = str2double(namesplit{2});
        options.rawDataPath = epochList{row,2};
        sweepAcq = eval(['AD0_e',num2str(epoch),'p',namesplit{3},'avg.UserData.Components']);
    
        % Import processed numbers if post QC
        if ~createNew
            options.animal = exp{row,"Animal"};
            options.task = exp{row,"Task"};
            included = exp{row,"Included"}{1};
            cellid = exp{row,"Cell"};
            sweepAcq = exp{row,"Sweep names"}{1};
            vhold = exp{row,"Vhold"};
            vholdInfo = exp{row,"VholdInfo"}{1};
        end
        
        % Initialize some temporary matrix
        sweeps = zeros(length(sweepAcq), size(eval(['AD0_e',num2str(epoch),'p',namesplit{3},'avg.data']),2));
        nAnalyzedSweeps = length(sweepAcq);
        processed = zeros(size(sweeps));
        QCs = cell(length(sweepAcq),1);
        cycles = cell(length(sweepAcq),1); % For detecting whether a sweep uses different cycle within an epoch
        protocols = cell(length(sweepAcq),1);
        statistics = cell(length(sweepAcq),1);

    
        % Load Vhold for epoch avg file (AD2)
        if withVholdAvg && withVhold
            load(fullfile(vholdList{row,2},vholdList{row,1}));
            namesplit = strsplit(vholdList{row,1},{'e','p'}); 
            if epoch ~= str2double(namesplit{2})
                error('Epoch number does not match between AD0 and AD2!!');
                vholdAcq = sweepAcq;
                vholdSweeps = zeros(size(sweeps));
                withVholdAvg = false;
            else
                vholdAcq = eval(['AD2_e',num2str(epoch),'p1avg.UserData.Components']);
                vholdEpoch = eval(['AD2_e',num2str(epoch),'p1avg.data']);
                vholdSweeps = zeros(length(vholdAcq),length(vholdEpoch));
            end
        elseif withVhold && ~withVholdAvg
            vholdAcq = sweepAcq;
            vholdSweeps = zeros(size(sweeps));
        end

        % Load csv file for vhold
        if ~strcmp(rig,'Wengang')
            csvOpts = detectImportOptions(fullfile(epochList{row,2},'InfoPatching.xlsx'));
            csvOpts.SelectedVariableNames = 1:11;%1:9
            csvOpts.VariableNamesRange = 25;
            info = readtable(fullfile(epochList{row,2},'InfoPatching.xlsx'),csvOpts);
            info = rmmissing(info,DataVariables="acq_");
            info = rmmissing(info,DataVariables="epoch");

            if iscell(info.acq_); info.acq_ = str2double(info.acq_); end
            if iscell(info.epoch); info.epoch = str2double(info.epoch); end
            if iscell(info.cyclePos); info.epoch = str2double(info.cyclePos); end
            if iscell(info.holding); info.acq_ = str2double(info.holding); end
        end
    
        %% Load individual sweeps
        for k = 1:length(sweepAcq)
            % Load sweep traces (.data)
            disp(['Loading ',sweepAcq{k},'.mat for epoch ',num2str(epoch)]);
            try load(fullfile(epochList{row,2},strcat(sweepAcq{k},'.mat'))); 
            catch
                warning(strcat("Sweep ",num2str(sweepAcq{k}), " not saved, skipping this sweep!"));
                nAnalyzedSweeps = nAnalyzedSweeps - 1;
                continue
            end
    

            % Extract raw trace
            raw_trace = eval([sweepAcq{k},'.data']);
            headerString = eval([sweepAcq{k},'.UserData.headerString']);
            

            % Extract experiment protocol from header string
            protocol = getCellProtocol(headerString,...
                                       outputFs=options.outputFs,...
                                       rcCheckRecoveryWindow=options.rcCheckRecoveryWindow,...
                                       rig=rig);
            protocols{k} = protocol;
            cycles{k} = protocol.cycle;


            % Skip analysis for some sweeps
            % 1. For whole field cycles, warn user if total length of raw_trace
            % is differnt from avg trace
            % 2. For DMD random search cycles, skip since random searches means
            % some sweeps by default will have different length
            if contains(protocol.cycle,'randomSearch')
                disp(['     Sweep cycle is ',protocol.cycle,', skip epoch-level anlaysis below.']);
                nAnalyzedSweeps = nAnalyzedSweeps - 1;
                continue
            end
            if length(raw_trace) ~= size(sweeps,2)
                warning("Sweep duration is different from epoch avg duration!!");
                nAnalyzedSweeps = nAnalyzedSweeps - 1;
                continue
            end
            sweeps(k,:) = raw_trace;


            % For some recording of wengang, there's no stim, default to be
            % happen at 1sec after the sweep
            if isempty(protocol.stimOnset)
                if isfield(options,'defaultStimOnset')
                    protocol.stimOnset = options.defaultStimOnset;
                    protocol.numPulses = 1;
                    warning('Cannot find stimOnset, set to default!');
                else
                    error('Need to provide stimOnset time if stimOnset is not provided!');
                end
            end
    
    
            % Define time window for baseline & analysis
            % Calculate baseline window: have two windows, one before pulse and one after pulse
            rcCheckOnset = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='delay') * (options.outputFs/1000);
            rcCheckPulseWidth = getHeaderValue(headerString,'state.phys.internal.pulseString_RCCheck',pulseVar='pulseWidth') * (options.outputFs/1000);
            rcCheckEnd = rcCheckOnset + rcCheckPulseWidth + (options.rcCheckRecoveryWindow*(options.outputFs/1000));
            stimDuration = ((protocol.numPulses * protocol.isi)+200) * options.outputFs/1000; % 200ms recovery window after last pulse
            if rcCheckOnset < protocol.stimOnset(1)
                preStimWindow = rcCheckEnd : (protocol.stimOnset(1)-1);
                postStimWindow = (protocol.stimOnset(1) + stimDuration):length(raw_trace);
                baselineWindow = [preStimWindow,postStimWindow];
            else
                preStimWindow = 1:(protocol.stimOnset(1)-1);
                postStimWindow = (protocol.stimOnset(1) + stimDuration):rcCheckOnset;
                baselineWindow = [preStimWindow,postStimWindow];
            end

            % Define time window for plotting
            options.eventSample = protocol.stimOnset(1);
            % timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
            % timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
            % plotWindowLength = timeRangeEndSample(1) - timeRangeStartSample(1) + 1;
            % plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
    
            % Define control window: 50ms before each spot stim onset
            controlWindowSamples = options.controlWindowLength * options.outputFs/1000;
            controlWindow = options.eventSample-controlWindowSamples-1 : options.eventSample-1;
            
            % Define analysis window: 50ms after each spot stim onset
            analysisWindowSamples = options.analysisWindowLength * options.outputFs/1000;
            analysisWindow = options.eventSample:options.eventSample+analysisWindowSamples;
    
            % Define time window for peak window analysis (1ms around peak)
            peakWindowWidth = (options.peakWindow/(2*1000)) * options.outputFs;

            % Save time windows to options and protocols
            options.baselineWindow = baselineWindow;
            options.baselineWindow_preStim = preStimWindow;
            options.baselineWindow_postStim = postStimWindow;
            options.analysisWindow = analysisWindow;
            % options.plotWindowLength = plotWindowLength;
            % options.plotWindowTime = plotWindowTime;
            options.peakWindowWidth = peakWindowWidth;
            options.analysisWindowSamples = analysisWindowSamples;
            options.controlWindowSamples = controlWindowSamples;
            options.controlWindow = controlWindow;
            options.stimDuration = stimDuration;

            % Calcultate QC for all sweeps
            qc = getCellQC(headerString,calculate=options.calculateQC,...
                            data=raw_trace,baselineWindow=baselineWindow,...
                            plot=false,rig=rig);
            QCs{k} = qc;
    

            % Process trace: mean-subtracted, optional LP
            % Mean subtraction
            baselineAvg = qc.Ibaseline;
            mean_subtracted = raw_trace - baselineAvg;
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
                processed(k,:) = processed_trace;
            else
                processed_trace = mean_subtracted;
                processed(k,:) = processed_trace;
            end
    
            
            % Find peak and area during analysis window (for processed)
            if ~contains(protocol.cycle,'randomSearch')
                % Calculate statistics
                % Find auc
                stats.response.auc = sum(processed_trace(analysisWindow)) / options.outputFs;
                stats.baseline.auc = sum(processed_trace(controlWindow)) / options.outputFs;
    
                % Find min and max value for stim response
                trace = processed_trace(analysisWindow);
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

                % Find E/I index
                stats.response.EIindex = abs(stats.response.max)-abs(stats.response.min) / abs(stats.response.max)+abs(stats.response.min);
                stats.baseline.EIindex = abs(stats.baseline.max)-abs(stats.baseline.min) / abs(stats.baseline.max)+abs(stats.baseline.min);

                statistics{k} = stats;
            end


            % Load Vhold traces if needed
            if withVhold
                try load(fullfile(vholdList{row,2},strcat(vholdAcq{k},'.mat'))); 
                catch
                    warning(strcat("Vhold for sweep ",num2str(vholdAcq{k}), " not saved, skipping vhold of this sweep!"));
                    continue
                end
                vhold_trace = eval([vholdAcq{k},'.data']);
                if length(vhold_trace) ~= size(vholdSweeps,2)
                    warning("Sweep duration is different from epoch avg duration!!");
                    continue;
                end
                vholdSweeps(k,:) = vhold_trace;
            end
        end

        %% Merge protocol/qc/stats for each sweep into one for each epoch
        if ~contains(protocol.cycle,{'randomSearch','plasticity'})
            protocols = mergeStructs(protocols);
            statistics = mergeStructs(statistics);
            QC = mergeStructs(QCs);

            %% Calculate peripeak data
    
            % Find min and max value for stim response
            trace = processed(:,analysisWindow);
            avg_trace = mean(trace,1);
            [~,maxIdx] = max(avg_trace); [~,minIdx] = min(avg_trace);
            % Average around max/min idx to get final value
            maxWindowStart = max(1,maxIdx-peakWindowWidth);
            maxWindowEnd = min(maxIdx+peakWindowWidth,length(avg_trace));
            minWindowStart = max(1,minIdx-peakWindowWidth);
            minWindowEnd = min(minIdx+peakWindowWidth,length(avg_trace));
            statistics.response.periMax = mean(trace(:,maxWindowStart:maxWindowEnd),2);
            statistics.response.periMin = mean(trace(:,minWindowStart:minWindowEnd),2);
            statistics.response.periMaxTime = maxIdx * 1000/options.outputFs;
            statistics.response.periMinTime = minIdx * 1000/options.outputFs;
    
            % Find min and max for control baseline
            trace = processed(:,controlWindow);
            avg_trace = mean(trace,1);
            [~,maxIdx] = max(avg_trace); [~,minIdx] = min(avg_trace);
            % Average around max/min idx to get final value
            maxWindowStart = max(1,maxIdx-peakWindowWidth);
            maxWindowEnd = min(maxIdx+peakWindowWidth,length(avg_trace));
            minWindowStart = max(1,minIdx-peakWindowWidth);
            minWindowEnd = min(minIdx+peakWindowWidth,length(avg_trace));
            statistics.baseline.periMax = mean(trace(:,maxWindowStart:maxWindowEnd),2);
            statistics.baseline.periMin = mean(trace(:,minWindowStart:minWindowEnd),2);
            statistics.baseline.periMaxTime = maxIdx * 1000/options.outputFs;
            statistics.baseline.periMinTime = minIdx * 1000/options.outputFs;
        
        
            %% Determine cellID and mean epoch vhold for that cell
            if createNew
                if strcmp(rig,'Wengang')
                    if withVhold
                        if withVholdAvg
                            vholdEpochMean = mean(vholdEpoch(:,baselineWindow),"all");
                        else
                            vholdEpochMean = mean(vholdSweeps(:,baselineWindow),"all"); 
                            vholdEpoch = mean(vholdSweeps(:,baselineWindow),1);
                        end
                        if vholdEpochMean < -50; cellid = cellid + 1; end
                        vholdSweepsMean = mean(vholdSweeps(:,1:20000),2);
                        vholdInfo.vholdSweepsMean = vholdSweepsMean;
                        vholdInfo.vholdEpochTrace = vholdEpoch;
                        vholdInfo.vholdSweepsTrace = vholdSweeps;
                    else
                        cellid = row;
                        vholdEpochMean = 100;
                    end
                    vhold = vholdEpochMean;
                    vholdInfo.vholdEpochMean = vholdEpochMean;
                else
                    cellid = epochList{row,end};
                    
                    % Determine Vhold
                    % Find in .xlsx file first, if can't find, use the heuristic that
                    % cells with negative leak current Vhold=-70, otherwise Vhold = 10;
                    acqsplit = split(sweepAcq{k},'_'); acqNum = str2double(acqsplit{end});
                    vhold = info{info.acq_ == acqNum,'holding'};
                    if isempty(vhold)
                        if baselineAvg < 0; vhold = -70;
                        else; vhold = 10; end
                    end
                end
            end
    
            %% Remove empty/erraneous sweeps
            if createNew || isfield(options,'include')
                included = ones(nAnalyzedSweeps,1);
    
                % Remove sweeps with different Vhold
                if options.filterSweeps && withVhold
                    % rounded_epoch_vhold = roundToTarget(vholdEpochMean,[-70,0,8]);
                    % rounded_sweeps_vhold = roundToTarget(vholdSweepsMean,[-70,0,8]);
                    % included = (rounded_sweeps_vhold == rounded_epoch_vhold);
                    % diffVhold_included = ~isoutlier(vholdInfo.vholdSweepsMean,'mean');
                    % included = all([included,diffVhold_included],2);
                end 
        
                % Remove sweeps with different cycles
                cycles = cycles(~cellfun(@isempty,cycles));
                if ~isempty(cycles)
                    cyclesCount = tabulate(cycles);
                    [~, mostCommonIdx] = max([cyclesCount{:, 2}]);
                    if ~all(strcmp(cycles, cycles{mostCommonIdx}))
                        cycles_included = strcmp(cycles, cycles{mostCommonIdx});
                        % protocol = protocols{mostCommonIdx};
                        included = all([included,cycles_included],2);
                    end
                end
    
                % Remove sweeps based on QC
                if ~isempty(fieldnames(QC))
                    % Remove sweeps from QC based on baseline current
                    if all(isnan(QC.Ibaseline)) || all(isnan(QC.Ibaseline_std))
                        QCIncluded = zeros(length(included),1);
                        QC.included = QCIncluded;
                    else
                        Ibaseline_avg_filter = QC.Ibaseline >= options.QCThreshold.Ibaseline;
                        Ibaseline_std_filter = QC.Ibaseline_std <= options.QCThreshold.Ibaseline_std;
                        QCIncluded = all([Ibaseline_avg_filter,Ibaseline_std_filter],2);
                        QC.included = QCIncluded;
                    end
                    included = all([included,QCIncluded],2);
    
                    % Remove epochs with high Rs or high Verror
                    Rs_final = mean(QC.Rs(included==1));
                    Verror_final = mean(abs(QC.Verror(included==1)));
                    if Rs_final > options.QCThreshold.Rs || Verror_final > options.QCThreshold.Verror
                        included = zeros(length(sweepAcq),1);
                    end
    
                    % Remove sweeps based on other provided criteria
                    if ~isempty(options.QCThreshold.include)
                        for criterion = options.QCThreshold.include
                            if length(eval(criterion{1})) > 1
                                newIncluded = eval(criterion{1});
                            else
                                newIncluded = ones(length(sweepAcq),1) * eval(criterion{1});
                            end
                            included = all([included,newIncluded],2);
                        end
                        if length(included) ~= length(sweepAcq)
                            warning('included have different size than #sweeps! Reset included to true for all!');
                            included = ones(length(sweepAcq),1);
                        end
                    end
                end 
            end
        end
    
        %% Store everything in epochs
        epochs{row,'Session'} = string(options.sessionPath);
        epochs{row,'Animal'} = options.animal;
        epochs{row,'Task'} = options.task;
        epochs{row,'Epoch'} = epoch;
        epochs{row,'Cell'} = cellid;
        epochs{row,'Vhold'} = vhold;
        epochs{row,'Included'} = num2cell(included,[1 2]);
        epochs{row,'Sweep names'} = num2cell(sweepAcq,[1 2]);
        epochs{row,'Raw sweeps'} = num2cell(sweeps,[1 2]);
        epochs{row,'Processed sweeps'} = num2cell(processed,[1 2]);
        epochs{row,'Protocol'} = {protocols};
        epochs{row,'Stats'} = {statistics};
        epochs{row,'QC'} = {QC};
        epochs{row,'VholdInfo'} = {vholdInfo};
        epochs{row,'Options'} = {options};

        %% Plot epoch analysis summary
        if options.plot && ~isempty(fieldnames(QC)) && nAnalyzedSweeps > 0
            close all;
            plotEpochSummary(epochs,row,save=true,saveDataPath=options.saveDataPath);
        end
    end
    
    % If vhold never have -70, plot by each epoch
    if isscalar(unique(epochs{:,'Cell'}))
        epochs{:,'Cell'} = (1:size(epochs,1))';
    end
    
    %% Save epochs.mat
    if options.save
        save(strcat(options.saveDataPath,filesep,'epochs_',today),'epochs','-v7.3');
        save(strcat(options.rawDataPath,filesep,'epochs_',today),'epochs','-v7.3');
        disp(strcat("New epochs.mat created & saved: ",expName));
    end
end

%% Create cells.mat

if options.getCellTable
    cells = getCellTable(epochs,save=options.save,...
                         timeRange=options.timeRange,...
                         outputFs=options.outputFs,...
                         controlWindowLength=options.controlWindowLength,...
                         nArtifactSamples=options.nArtifactSamples,...
                         peakWindow=options.peakWindow);
end

%% Define output
varargout{1} = epochs;
if options.getCellTable; varargout{2} = cells; end

end