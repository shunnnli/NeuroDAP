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
    options.vholdChannel = 'AD2'
    options.sortSweepsByAcq logical = true

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
    if ~isfield(options.QCThreshold,'include')
        options.QCThreshold.include = {};
    elseif isstring(options.QCThreshold.include)
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
            vholdList = sortrows(struct2cell(dir(fullfile(exp,[options.vholdChannel,'_e*','p1avg.mat'])))',3);
        else
            epochList = {}; vholdList = {};
            for i = 1:size(exp,1)
                epochPath = strcat(filesep,'AD0_e',num2str(exp{i,"Epoch"}),'p1avg.mat');
                epochList = [epochList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, epochPath))))'];
                vholdPath = strcat(filesep,'options.vholdChannel_e',num2str(exp{i,"Epoch"}),'p1avg.mat');
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
                cellVholdList = sortrows(struct2cell(dir(fullfile(cellPath,[options.vholdChannel,'_e*','p*avg.mat'])))',3);
                cellVholdList(:,end+1) = num2cell(sscanf(cellFolders(c).name,'cell%d'),[1 2]); % store cell number as the last column
                vholdList = [vholdList; cellVholdList];
            end
        else
            epochList = {}; vholdList = {};
            for i = 1:size(exp,1)
                epochPath = strcat(filesep,'cell',num2str(exp{i,"Cell"}),filesep,'AD0_e',num2str(exp{i,"Epoch"}),'p*avg.mat');
                epochList = [epochList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, epochPath))))'];
                vholdPath = strcat(filesep,'cell',num2str(exp{i,"Cell"}),filesep,options.vholdChannel,'_e',num2str(exp{i,"Epoch"}),'p*avg.mat');
                vholdList = [vholdList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, vholdPath))))'];
            end
        end
    end
    
    %% Group multiple p*avg files into one epoch entry
    n = size(epochList,1);
    epochNum = nan(n,1);
    pulseNum = nan(n,1);
    cellNum  = zeros(n,1);
    
    for i = 1:n
        fn = epochList{i,1};                 % e.g. AD0_e140p19avg.mat
        ns = strsplit(fn, {'e','p','avg.mat'});
        epochNum(i) = str2double(ns{2});
        pulseNum(i) = str2double(ns{3});
        if ~strcmp(rig,'Wengang')
            cellNum(i) = epochList{i,end};   % Paolo: cell number stored in last col
        end
    end
    
    key = strcat(string(epochList(:,2)), "|c", string(cellNum), "|e", string(epochNum));
    [ukey,~,g] = unique(key,'stable');
    nGroups = numel(ukey);
    
    grpFolder   = cell(nGroups,1);
    grpEpoch    = nan(nGroups,1);
    grpCell     = nan(nGroups,1);
    grpAvgFiles = cell(nGroups,1);   % cell array of filenames
    grpPulses   = cell(nGroups,1);   % numeric array per group (same order as files)
    
    for gi = 1:nGroups
        idx = find(g==gi);
        grpFolder{gi} = epochList{idx(1),2};
        grpEpoch(gi)  = epochNum(idx(1));
        grpCell(gi)   = cellNum(idx(1));
    
        % sort p-files by PulseNum so things are deterministic
        [~,ord] = sort(pulseNum(idx));
        grpAvgFiles{gi} = epochList(idx(ord),1);
        grpPulses{gi}   = pulseNum(idx(ord));
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
    epochs = table('Size',[nGroups,length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    
    % Check whether vhold source exists (AD channel files or Excel)
    if strcmpi(string(options.vholdChannel),'excel')
        withVhold = exist(fullfile(grpFolder{1},'InfoPatching.xlsx'),'file') == 2;
    else
        withVhold = ~isempty(dir(fullfile(grpFolder{1}, options.vholdChannel + "_*.mat")));
        if ~withVhold
            % fallback if your files don't use the underscore convention
            withVhold = ~isempty(dir(fullfile(grpFolder{1}, options.vholdChannel + "*.mat")));
        end
    end
    % Wengang recordings do not have cell folders, so retain an explicit
    % session-level counter rather than implicitly reusing loop variables.
    wengangCellId = 0;
    
    %% Iterate & analyze individual epoch
    
    for row = 1:nGroups
        clearvars AD*
    
        % --- NEW: load all avg files for this epoch and combine acquisitions ---
        epoch = grpEpoch(row);
        options.rawDataPath = grpFolder{row};
        avgFilesThisEpoch = grpAvgFiles{row};
        pulseNumsThisEpoch = grpPulses{row};
        
        sweepAcq = {};
        pulseNumByAcq = [];     % same length as sweepAcq
        avgLen = [];
        for f = 1:numel(avgFilesThisEpoch)
            fn = avgFilesThisEpoch{f};
            p  = pulseNumsThisEpoch(f);
            load(fullfile(options.rawDataPath, fn));
            comps = eval(sprintf('AD0_e%dp%davg.UserData.Components', epoch, p));
            sweepAcq = [sweepAcq; comps(:)];
            pulseNumByAcq = [pulseNumByAcq; repmat(p, numel(comps), 1)];
            thisLen = size(eval(sprintf('AD0_e%dp%davg.data', epoch, p)), 2);
            if isempty(avgLen)
                avgLen = thisLen;
            elseif thisLen ~= avgLen
                warning('Epoch %d has p%davg with different avg length (%d vs %d). Those sweeps may be skipped later.', ...
                        epoch, p, thisLen, avgLen);
                % (you can decide to handle this more strictly if you want)
            end
        end
    
        % Import processed numbers if post QC
        if ~createNew
            options.animal = exp{row,"Animal"};
            options.task = exp{row,"Task"};
            included = exp{row,"Included"}{1};
            cellid = exp{row,"Cell"};
            sweepAcq = exp{row,"Sweep names"}{1};
            vhold = exp{row,"Vhold"};
            vholdInfo = exp{row,"VholdInfo"}{1};
            prot = exp{row,"Protocol"}{1};
            if isstruct(prot) && isfield(prot,'PulseNumByAcq')
                pulseNumByAcq = prot.PulseNumByAcq;
            end
        else
            included = true(numel(sweepAcq),1);
            vhold = NaN;
            if strcmp(rig,'Paolo')
                cellid = grpCell(row);
            else
                % A DMD search can precede the first voltage-clamp epoch.
                % Use cell 1 until the holding-voltage heuristic detects a
                % later cell transition.
                cellid = max(wengangCellId,1);
            end
        end

        
        % Initialize some temporary matrix
        % sweeps = zeros(length(sweepAcq), size(eval(['AD0_e',num2str(epoch),'p',namesplit{3},'avg.data']),2));
        sweeps = zeros(length(sweepAcq), avgLen);
        nAnalyzedSweeps = length(sweepAcq);
        processed = zeros(size(sweeps));
        QCs = cell(length(sweepAcq),1);
        cycles = cell(length(sweepAcq),1); % For detecting whether a sweep uses different cycle within an epoch
        protocols = cell(length(sweepAcq),1);
        statistics = cell(length(sweepAcq),1);
        loadedSweep = false(length(sweepAcq),1);
        analyzedSweep = false(length(sweepAcq),1);
        QC = struct();
        

        % --- Vhold bookkeeping (per-sweep) via getSweepVhold ---
        % Keep vholdInfo fields the same as before:
        %   vholdInfo.vholdEpochMean   : most common vhold across sweeps in this epoch
        %   vholdInfo.vholdSweepsMean  : per-sweep scalar vhold values
        %   vholdInfo.vholdEpochTrace  : mean vhold trace (if available) or constant vholdEpochMean
        %   vholdInfo.vholdSweepsTrace : per-sweep vhold traces (if available; NaNs otherwise)
        vholdSweeps = nan(length(sweepAcq), avgLen);
        vholdSweepsMean = nan(length(sweepAcq), 1);
        vholdInfo.vholdEpochMean = NaN;
        vholdInfo.vholdSweepsMean = vholdSweepsMean;
        vholdInfo.vholdEpochTrace = nan(1,avgLen);
        vholdInfo.vholdSweepsTrace = vholdSweeps;

        % Load InfoPatching.xlsx once (if present) and pass into getSweepVhold
        infoTable = [];
        xlsxPath = fullfile(grpFolder{row}, 'InfoPatching.xlsx');
        if exist(xlsxPath,'file')
            infoTable = readInfoPatchingTable(xlsxPath, 11);
            infoTable = rmmissing(infoTable, DataVariables="acq_");
            infoTable = rmmissing(infoTable, DataVariables="epoch");

            if iscell(infoTable.acq_); infoTable.acq_ = str2double(infoTable.acq_); end
            if iscell(infoTable.epoch); infoTable.epoch = str2double(infoTable.epoch); end
            if iscell(infoTable.cyclePos); infoTable.cyclePos = str2double(infoTable.cyclePos); end
            if iscell(infoTable.holding); infoTable.holding = str2double(infoTable.holding); end
        end


        % --- OPTIONAL: sort sweep names by acquisition number (AD0_###) ---
        if options.sortSweepsByAcq && ~isempty(sweepAcq)
            % Extract numeric acquisition number from strings like "AD0_342"
            acqNums = nan(numel(sweepAcq),1);
            for i = 1:numel(sweepAcq)
                s = string(sweepAcq{i});
                tok = regexp(s, '_(\d+)$', 'tokens', 'once');
                if ~isempty(tok)
                    acqNums(i) = str2double(tok{1});
                end
            end
            % Sort, putting any non-parsable names at the end (NaNs last)
            [~, ord] = sort(acqNums, 'ascend', 'MissingPlacement','last');
            sweepAcq = sweepAcq(ord);
            % Keep any per-sweep vectors aligned with sweepAcq
            if exist('pulseNumByAcq','var') && numel(pulseNumByAcq) == numel(ord)
                pulseNumByAcq = pulseNumByAcq(ord);
            end
            if exist('included','var') && numel(included) == numel(ord)
                included = included(ord);
            end
            if exist('vholdInfo','var') && numel(vholdInfo) == numel(ord)
                vholdInfo = vholdInfo(ord);
            end
        end
    
        %% Load individual sweeps
        for k = 1:length(sweepAcq)
            % Load sweep traces (.data)
            disp(['Loading ',sweepAcq{k},'.mat for epoch ',num2str(epoch)]);
            try load(fullfile(grpFolder{row},strcat(sweepAcq{k},'.mat'))); 
            catch
                warning(strcat("Sweep ",num2str(sweepAcq{k}), " not saved, skipping this sweep!"));
                nAnalyzedSweeps = nAnalyzedSweeps - 1;
                continue
            end
            loadedSweep(k) = true;
    

            % Extract raw trace
            raw_trace = eval([sweepAcq{k},'.data']);
            headerString = eval([sweepAcq{k},'.UserData.headerString']);
            

            % Extract experiment protocol from header string
            protocol = getCellProtocol(headerString,...
                                       outputFs=options.outputFs,...
                                       rcCheckRecoveryWindow=options.rcCheckRecoveryWindow,...
                                       rig=rig);
            protocol.PulseNum = pulseNumByAcq(k);  % the p# that this acq came from
            protocols{k} = protocol;
            cycles{k} = protocol.cycle;


            % Skip analysis for some sweeps
            % 1. For whole field cycles, warn user if total length of raw_trace
            % is differnt from avg trace
            % 2. For DMD random search cycles, skip since random searches means
            % some sweeps by default will have different length
            if isRandomSearchCycle(protocol.cycle)
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


            % --- Vhold (per-sweep) via getSweepVhold ---
            try
                [vh, vhMeta] = getSweepVhold(grpFolder{row}, sweepAcq{k}, baselineAvg, options.vholdChannel, infoTable=infoTable);
                vholdSweepsMean(k) = vh;
                if isfield(vhMeta,'trace') && ~isempty(vhMeta.trace)
                    tr = vhMeta.trace(:)';
                    if numel(tr) == avgLen
                        vholdSweeps(k,:) = tr;
                    else
                        warning('Vhold trace length mismatch for %s (got %d, expected %d).', sweepAcq{k}, numel(tr), avgLen);
                    end
                end
            catch ME
                warning('getSweepVhold failed for %s: %s', sweepAcq{k}, ME.message);
            end
    
            
            % Find peak and area during analysis window (for processed)
            if ~isRandomSearchCycle(protocol.cycle)
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
                analyzedSweep(k) = true;
            end
        end

        %% Classify every sweep before deciding how to store the epoch
        cycleNames = getSweepCycleNames(protocols,numel(sweepAcq));
        validProtocol = loadedSweep & strlength(cycleNames) > 0;
        isRandomSweep = validProtocol & isRandomSearchCycle(cycleNames);
        isPlasticitySweep = validProtocol & contains(cycleNames,"plasticity",IgnoreCase=true);
        hasRandomSearch = any(isRandomSweep);
        hasPlasticity = any(isPlasticitySweep);
        hasFullField = any(validProtocol & ~isRandomSweep & ~isPlasticitySweep);

        % Calculate epoch Vhold and cell identity independently of protocol
        % representation. DMD-only epochs legitimately have NaN here because
        % loadSlicesDMD calculates Vhold from each raw sweep.
        if createNew
            vh = vholdSweepsMean(~isnan(vholdSweepsMean));
            if isempty(vh)
                vholdEpochMean = NaN;
                vholdEpochTrace = nan(1,avgLen);
            else
                vholdEpochMean = mode(round(vh));
                if any(~isnan(vholdSweeps(:)))
                    vholdEpochTrace = mean(vholdSweeps,1,'omitnan');
                else
                    vholdEpochTrace = vholdEpochMean * ones(1,avgLen);
                end
            end

            if strcmp(rig,'Wengang')
                if ~isnan(vholdEpochMean) && vholdEpochMean < -50
                    wengangCellId = wengangCellId + 1;
                end
                cellid = max(wengangCellId,1);
            else
                cellid = grpCell(row);
            end
            vhold = vholdEpochMean;
        else
            vholdEpochMean = vhold;
            if any(~isnan(vholdSweeps(:)))
                vholdEpochTrace = mean(vholdSweeps,1,'omitnan');
            else
                vholdEpochTrace = vhold * ones(1,avgLen);
            end
        end

        vholdInfo.vholdEpochMean = vholdEpochMean;
        vholdInfo.vholdSweepsMean = vholdSweepsMean;
        vholdInfo.vholdEpochTrace = vholdEpochTrace;
        vholdInfo.vholdSweepsTrace = vholdSweeps;

        if hasRandomSearch && hasFullField
            warning(['Cell %g, epoch %g contains both full-field and random-search sweeps. ', ...
                     'Sweep-level protocols are preserved; loadSlicesDMD will analyze only the searches.'], ...
                    cellid,epoch);
        end

        %% Merge only homogeneous full-field epochs
        % Search and plasticity protocols remain sweep-level cell arrays.
        % This also preserves the individual cycles in a mixed acquisition
        % group so loadSlicesDMD can select only confirmed search sweeps.
        if any(analyzedSweep) && ~hasRandomSearch && ~hasPlasticity
            analysisSweep = find(analyzedSweep);
            protocols = mergeStructs(protocols(analysisSweep));
            protocols.PulseNumByAcq = pulseNumByAcq(analysisSweep);
            protocols.PulseNumsInEpoch = unique(pulseNumByAcq(analysisSweep))';
            protocols.AcqPulseMap = table(string(sweepAcq(analysisSweep)), ...
                                          pulseNumByAcq(analysisSweep), ...
                                          'VariableNames',{'Acq','PulseNum'});
            statistics = mergeStructs(statistics(analysisSweep));
            QC = mergeStructs(QCs(analysisSweep));

            %% Calculate peripeak data from successfully analyzed sweeps
            trace = processed(analysisSweep,analysisWindow);
            avg_trace = mean(trace,1);
            [~,maxIdx] = max(avg_trace); [~,minIdx] = min(avg_trace);
            maxWindowStart = max(1,maxIdx-peakWindowWidth);
            maxWindowEnd = min(maxIdx+peakWindowWidth,length(avg_trace));
            minWindowStart = max(1,minIdx-peakWindowWidth);
            minWindowEnd = min(minIdx+peakWindowWidth,length(avg_trace));
            statistics.response.periMax = mean(trace(:,maxWindowStart:maxWindowEnd),2);
            statistics.response.periMin = mean(trace(:,minWindowStart:minWindowEnd),2);
            statistics.response.periMaxTime = maxIdx * 1000/options.outputFs;
            statistics.response.periMinTime = minIdx * 1000/options.outputFs;

            trace = processed(analysisSweep,controlWindow);
            avg_trace = mean(trace,1);
            [~,maxIdx] = max(avg_trace); [~,minIdx] = min(avg_trace);
            maxWindowStart = max(1,maxIdx-peakWindowWidth);
            maxWindowEnd = min(maxIdx+peakWindowWidth,length(avg_trace));
            minWindowStart = max(1,minIdx-peakWindowWidth);
            minWindowEnd = min(minIdx+peakWindowWidth,length(avg_trace));
            statistics.baseline.periMax = mean(trace(:,maxWindowStart:maxWindowEnd),2);
            statistics.baseline.periMin = mean(trace(:,minWindowStart:minWindowEnd),2);
            statistics.baseline.periMaxTime = maxIdx * 1000/options.outputFs;
            statistics.baseline.periMinTime = minIdx * 1000/options.outputFs;

            %% Apply inclusion criteria while preserving sweep alignment
            localIncluded = true(numel(analysisSweep),1);
            analysisCycles = cycleNames(analysisSweep);
            uniqueCycles = unique(analysisCycles,'stable');
            cycleCounts = zeros(numel(uniqueCycles),1);
            for cycleIdx = 1:numel(uniqueCycles)
                cycleCounts(cycleIdx) = sum(analysisCycles == uniqueCycles(cycleIdx));
            end
            [~,mostCommonIdx] = max(cycleCounts);
            localIncluded = localIncluded & analysisCycles == uniqueCycles(mostCommonIdx);

            if ~isempty(fieldnames(QC))
                if all(isnan(QC.Ibaseline)) || all(isnan(QC.Ibaseline_std))
                    QCIncluded = false(numel(analysisSweep),1);
                else
                    QCIncluded = QC.Ibaseline >= options.QCThreshold.Ibaseline & ...
                                 QC.Ibaseline_std <= options.QCThreshold.Ibaseline_std;
                end
                QCIncluded = QCIncluded(:);
                QC.included = QCIncluded;
                localIncluded = localIncluded & QCIncluded;

                Rs_final = mean(QC.Rs(localIncluded),'omitnan');
                Verror_final = mean(abs(QC.Verror(localIncluded)),'omitnan');
                if Rs_final > options.QCThreshold.Rs || Verror_final > options.QCThreshold.Verror
                    localIncluded(:) = false;
                end
            end

            for criterion = options.QCThreshold.include
                criterionIncluded = eval(criterion{1});
                if isscalar(criterionIncluded)
                    criterionIncluded = repmat(logical(criterionIncluded),numel(analysisSweep),1);
                elseif numel(criterionIncluded) == numel(sweepAcq)
                    criterionIncluded = logical(criterionIncluded(analysisSweep));
                elseif numel(criterionIncluded) == numel(analysisSweep)
                    criterionIncluded = logical(criterionIncluded(:));
                else
                    warning('Ignoring QC criterion with %d values; expected 1, %d, or %d.', ...
                            numel(criterionIncluded),numel(analysisSweep),numel(sweepAcq));
                    continue
                end
                localIncluded = localIncluded & criterionIncluded;
            end

            included = false(numel(sweepAcq),1);
            included(analysisSweep) = localIncluded;
        else
            % Search/mixed rows are not subjected to epoch-level full-field
            % QC. Keep successful acquisitions aligned; loadSlicesDMD will
            % apply its own QC after selecting search sweeps.
            included = loadedSweep;
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
