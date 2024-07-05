function varargout = loadSlices(exp,options)


arguments
    exp  %full path of the session or epochs.mat

    options.filterSignal logical = false
    options.filterSweeps logical = true % If true, do not analyze sweeps with different Vhold (included == false)

    options.save logical = true
    options.reload logical = false
    options.getCellTable logical = true % false for DMD setup
    options.calculateQC logical = false
    
    options.animal string
    options.task string = 'random'

    options.outputFs double = 10000
    options.timeRange double = [-20,100] % in ms
    options.eventSample % in sample
    options.nArtifactSamples double = 0 % in sample
    options.rcCheckRecoveryWindow double = 100 % in ms
    options.peakWindow double = 2; % in ms around the peak to average

    options.rawDataPath string
    options.saveDataPath string = 'default'
end

%% General setup

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
    options.saveDataPath = options.rawDataPath;
else
    if createNew; options.rawDataPath = exp;
    else; options.rawDataPath = exp{1,"Session"}; end
end

if createNew
    % Decide reload if session has already been loaded
    dirsplit = split(exp,filesep); expName = dirsplit{end};
    if ~isempty(dir(fullfile(exp,"epochs_*.mat")))
        if ~options.reload
            disp('Loading stop: epochs file found.');
            load(strcat(exp,filesep,'epochs_',expName,'.mat'));
            load(strcat(exp,filesep,'cells_',expName,'.mat'));
            if ~exist('epochs','var')
                if istable(exp); epochs = exp;
                elseif exist('epochs_old','var') && istable(epochs_old); epochs = epochs_old; end
            end
            varargout{1} = epochs;
            varargout{2} = cells;
            return
        end
    else
        options.reload = true;
    end
else
    dirsplit = split(exp{1,"Session"},filesep); expName = dirsplit{end};
end

%% Determine recording rig (ie file structure)

cellFolders = dir(fullfile(exp,'cell*'));

if isempty(cellFolders); rig = 'Wengang';
else; rig = 'Paolo'; end

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
    % timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000; % for plotWindow
    timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
    analysisWindow = (options.eventSample+options.nArtifactSamples) : timeRangeEndSample;
    peakWindowWidth = (options.peakWindow/(2*1000)) * options.outputFs;

else
    if createNew
        epochList = {}; vholdList = {};
        for c = 1:length(cellFolders)
            cellPath = strcat(cellFolders(c).folder,filesep,cellFolders(c).name);
            cellEpochList = sortrows(struct2cell(dir(fullfile(cellPath,['AD0_e*','p*avg.mat'])))',3);
            cellEpochList(:,end+1) = num2cell(str2double(cellFolders(c).name(5)),[1 2]); % store cell number as the last column
            epochList = [epochList; cellEpochList];
            cellVholdList = sortrows(struct2cell(dir(fullfile(cellPath,['AD2_e*','p*avg.mat'])))',3);
            cellVholdList(:,end+1) = num2cell(str2double(cellFolders(c).name(5)),[1 2]); % store cell number as the last column
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
            'cell','cell',...
            'cell','cell','cell',...
            'cell','cell'};
varNames = {'Session','Animal','Task','Epoch','Cell','Vhold',...
            'Included','Sweep names','Raw sweeps','Processed sweeps',...
            'Protocol',...
            'Peaks','AUCs',...
            'Rs','Rm','Cm',...
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
    sweepAcq = eval(['AD0_e',num2str(epoch),'p',namesplit{3},'avg.UserData.Components']);

    % Import processed numbers if post QC
    if ~createNew
        included = exp{row,"Included"}{1};
        cellid = exp{row,"Cell"};
        sweepAcq = exp{row,"Sweep names"}{1};
        Vhold = exp{row,"Vhold"};
        protocol = exp{row,"Protocol"}{1};
    end
    
    % Initialize some temporary matrix
    sweeps = zeros(length(sweepAcq), size(eval(['AD0_e',num2str(epoch),'p',namesplit{3},'avg.data']),2));
    processed = zeros(size(sweeps));
    vholds = zeros(length(sweepAcq),1);
    AUCs = zeros(length(sweepAcq),1);
    peaks = zeros(length(sweepAcq),1);
    Rss = zeros(length(sweepAcq),1);
    Rms = zeros(length(sweepAcq),1);
    Cms = zeros(length(sweepAcq),1);

    % For detecting whether a sweep uses different cycle within an epoch
    cycles = cell(length(sweepAcq),1);
    protocols = cell(length(sweepAcq),1);

    % Load Vhold for epoch avg file (AD2)
    if withVholdAvg && withVhold
        load(fullfile(vholdList{row,2},vholdList{row,1}));
        namesplit = strsplit(vholdList{row,1},{'e','p'}); 
        if epoch ~= str2double(namesplit{2})
            error('Epoch number does not match between AD0 and AD2!!');
        end
        vholdAcq = eval(['AD2_e',num2str(epoch),'p1avg.UserData.Components']);
        vholdEpoch = eval(['AD2_e',num2str(epoch),'p1avg.data']);
        vholdSweeps = zeros(length(vholdAcq),length(vholdEpoch));
    elseif withVhold && ~withVholdAvg
        vholdAcq = sweepAcq;
        vholdSweeps = zeros(size(sweeps));
    end

    %% Load individual sweeps
    for k = 1:length(sweepAcq)
        % Load sweep traces (.data)
        disp(['Loading ',sweepAcq{k},'.mat for epoch ',num2str(epoch)]);
        try load(fullfile(epochList{row,2},strcat(sweepAcq{k},'.mat'))); 
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
        protocols{k} = {protocol};
        cycles{k} = protocol.cycle;
        vholds(k) = mean(raw_trace(baselineWindow));

        % Skip analysis for some sweeps
        % 1. For whole field cycles, warn user if total length of raw_trace
        % is differnt from avg trace
        % 2. For DMD random search cycles, skip since random searches means
        % some sweeps by default will have different length
        if contains(protocol.cycle,'randomSearch')
            disp(['     Sweep cycle is ',protocol.cycle,', skip epoch-level anlaysis below.']);
            continue
        end
        if length(raw_trace) ~= size(sweeps,2)
            warning("Sweep duration is different from epoch avg duration!!");
            continue
        end
        sweeps(k,:) = raw_trace;

        % Extract quality metrics from header string
        qc = getCellQC(raw_trace,calculate=options.calculateQC,headerString=headerString);
        Rss(k) = qc.Rs; Rms(k) = qc.Rm; Cms(k) = qc.Cm;


        % Define time window for baseline & analysis
        % Calculate baseline window: have two windows, one before pulse and one after pulse
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
        % Calculate analysis window
        if ~strcmp(rig,'Wengang')
            options.eventSample = protocol.stimOnset;
            timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
            analysisWindow = (options.eventSample+options.nArtifactSamples) : timeRangeEndSample;
            peakWindowWidth = (options.peakWindow/(2*1000)) * options.outputFs;
        end

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
            processed(k,:) = processed_trace;
        else
            processed_trace = mean_subtracted;
            processed(k,:) = processed_trace;
        end


        % Load Vhold traces if needed
        if withVhold
            try load(fullfile(vholdList{row,2},strcat(vholdAcq{k},'.mat'))); 
            catch
                warning(strcat("Vhold for sweep ",num2str(vholdAcq{k}), " not saved, skipping this sweep!"));
                continue
            end
            vhold_trace = eval([vholdAcq{k},'.data']);
            if length(vhold_trace) ~= size(vholdSweeps,2)
                warning("Sweep duration is different from epoch avg duration!!");
                continue;
            end
            vholdSweeps(k,:) = vhold_trace;
        end


        % Find peak and area during analysis window (for processed)
        if ~contains(protocol.cycle,'randomSearch')
            AUCs(k) = sum(processed_trace(analysisWindow)) / options.outputFs;
            [~,peakIdx] = max(abs(processed_trace(analysisWindow)));
            peakWindowStart = max(1,peakIdx-peakWindowWidth) + analysisWindow(1);
            peakWindowEnd = min(analysisWindow(end),peakIdx+peakWindowWidth) + analysisWindow(1);
            peaks(k) = mean(processed_trace(peakWindowStart:peakWindowEnd));
        end
    end


    % Determine cellID and mean epoch vhold for that cell
    if strcmp(rig,'Wengang')
        if withVhold
            if withVholdAvg
                vholdEpochMean = mean(vholdEpoch(:,baselineWindow),"all");
            else
                vholdEpochMean = mean(vholdSweeps(:,baselineWindow),"all"); 
                vholdEpoch = mean(vholdSweeps(:,baselineWindow),1);
            end
    
            if vholdEpochMean < -50 && createNew; cellid = cellid + 1; end
    
            vholdSweepsMean = mean(vholdSweeps(:,1:20000),2);
        else
            cellid = row;
            vholdEpochMean = 100;
        end
    else
        cellid = epochList{row,end};
    end

    if mode(vholds) < 0; vhold = -70;
    else; vhold = 10; end


    % Remove empty/erraneous sweeps
    if createNew
        included = ones(length(sweepAcq),1);

        % Remove sweeps with different Vhold
        if options.filterSweeps && withVhold
            % rounded_epoch_vhold = roundToTarget(vholdEpochMean,[-70,0,8]);
            % rounded_sweeps_vhold = roundToTarget(vholdSweepsMean,[-70,0,8]);
            % included = (rounded_sweeps_vhold == rounded_epoch_vhold);
            included = ~isoutlier(vholdSweepsMean,'mean');
        end 

        % Remove sweeps with different cycles
        cycles = cycles(~cellfun(@isempty,cycles));
        cyclesCount = tabulate(cycles);
        [~, mostCommonIdx] = max([cyclesCount{:, 2}]);
        if ~all(strcmp(cycles, cycles{mostCommonIdx}))
            included = strcmp(cycles, cycles{mostCommonIdx});
            protocol = protocols{mostCommonIdx};
        end
    end

    %% Store everything in epochs
    epochs{row,'Session'} = string(options.rawDataPath);
    epochs{row,'Animal'} = options.animal;
    epochs{row,'Task'} = options.task;
    epochs{row,'Epoch'} = epoch; %string(epochList{row,1});
    epochs{row,'Cell'} = cellid;
    epochs{row,'Vhold'} = vhold;
    epochs{row,'Included'} = num2cell(included,[1 2]);
    epochs{row,'Sweep names'} = num2cell(sweepAcq,[1 2]);
    epochs{row,'Raw sweeps'} = num2cell(sweeps,[1 2]);
    epochs{row,'Processed sweeps'} = num2cell(processed,[1 2]);
    epochs{row,'Protocol'} = {protocol};
    epochs{row,'Peaks'} = num2cell(peaks,[1 2]);
    epochs{row,'AUCs'} = num2cell(AUCs,[1 2]);
    epochs{row,'Rs'} = num2cell(Rss,[1 2]);
    epochs{row,'Rm'} = num2cell(Rms,[1 2]);
    epochs{row,'Cm'} = num2cell(Cms,[1 2]);

    if withVhold
        vholdInfo.vholdEpochMean = vholdEpochMean;
        vholdInfo.vholdSweepsMean = vholdSweepsMean;
        vholdInfo.vholdEpochTrace = vholdEpoch;
        vholdInfo.vholdSweepsTrace = vholdSweeps;
    end
    epochs{row,'VholdInfo'} = {vholdInfo};
    epochs{row,'Options'} = {options};
end

% epochs = rmmissing(epochs);

% If vhold never have -70, plot by each epoch
if isscalar(unique(epochs{:,'Cell'}))
    epochs{:,'Cell'} = (1:size(epochs,1))';
end

%% Save epochs.mat
if options.save
    if createNew
        save(strcat(options.saveDataPath,filesep,'epochs_',expName),'epochs','-v7.3');
        disp(strcat("New epochs.mat created & saved: ",expName));
    else
        % Save old epochs
        today = char(datetime('today','Format','yyyyMMdd')); epochs_old = exp;
        save(strcat(options.saveDataPath,filesep,'epochs_',today),'epochs_old','-v7.3');
        save(strcat(options.saveDataPath,filesep,'epochs_',expName),'epochs','-append');
        disp(strcat("Existing epochs.mat reprocessed: ",expName));
    end
end

%% Create cells.mat

if options.getCellTable
    cells = getCellTable(epochs,save=options.save);
end

%% Define output
varargout{1} = epochs;
if options.getCellTable; varargout{2} = cells; end

end