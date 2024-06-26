function varargout = loadSlices(exp,options)


arguments
    exp  %full path of the session or epochs.mat

    options.filterSignal logical = false
    options.filterSweeps logical = true % If true, do not analyze sweeps with different Vhold (included == false)

    options.save logical = true
    options.reload logical = false
    options.calculateQC logical = false
    
    options.animal string
    options.task string = 'random'

    options.outputFs double = 10000
    options.timeRange double = [-20,100] % in ms
    options.eventSample double = 10000 % in sample
    options.nArtifactSamples double = 0 % in sample
    options.baselineWindow double = 1:10000;
    options.peakWindow double = 2; % in ms

    options.rawDataPath string = 'default'
    options.saveDataPath string
end

%% General setup

timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
analysisWindow = (options.eventSample+options.nArtifactSamples) : timeRangeEndSample;
peakWindowWidth = (options.peakWindow/(2*1000)) * options.outputFs;

% Determine processing type
% If exp is epochs.mat, then only analyze rows/neurons in epochs.mat 
% (reprocess signal for postQC data)
if ischar(exp); createNew = true;
elseif istable(exp); createNew = false; end

if strcmp(options.rawDataPath,'default')
    if createNew; options.rawDataPath = exp;
    else; options.rawDataPath = exp{1,"Session"}; end
    options.saveDataPath = options.rawDataPath;
else
    if createNew; options.saveDataPath = exp;
    else; options.saveDataPath = exp{1,"Session"}; end
end

if createNew
    % Decide reload if session has already been loaded
    dirsplit = split(exp,filesep); expName = dirsplit{end};
    if ~isempty(dir(fullfile(exp,"epochs_*.mat")))
        if ~options.reload
            disp('Loading stop: epochs file found.');
            load(strcat(exp,filesep,'epochs_',expName,'.mat'));
            if ~exist('epochs','var')
                if istable(exp); epochs = exp;
                elseif exist('epochs_old','var') && istable(epochs_old); epochs = epochs_old; end
            end
            return
        end
    else
        options.reload = true;
    end
else
    dirsplit = split(exp{1,"Session"},filesep); expName = dirsplit{end};
end

%% Grab corresponding recordings through epoch average files

if createNew
    epochList = sortrows(struct2cell(dir(fullfile(exp,['AD0_e*','p1avg.mat'])))',3);
    vholdList = sortrows(struct2cell(dir(fullfile(options.rawDataPath,['AD2_e*','p1avg.mat'])))',3);
else
    epochList = {};
    for i = 1:size(exp,1)
        epochPath = strcat(filesep,'AD0_e',num2str(exp{i,"Epoch"}),'p1avg.mat');
        epochList = [epochList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, epochPath))))'];
    end
    vholdList = {};
    for i = 1:size(exp,1)
        vholdPath = strcat(filesep,'AD2_e',num2str(exp{i,"Epoch"}),'p1avg.mat');
        vholdList = [vholdList; struct2cell(dir(fullfile(strcat(exp{:,"Session"}{i}, vholdPath))))'];
    end
end

% Initialize epochs.mat
varTypes = {'string','string','string','double','double',...
            'cell','cell','cell','cell',...
            'cell','cell',...
            'cell','cell','cell',...
            'double','cell','cell','cell'};
varNames = {'Session','Animal','Task','Epoch','Cell',...
            'Included','Sweep names','Raw sweeps','Processed sweeps',...
            'Peaks','AUCs',...
            'Rs','Rm','Cm',...
            'Vhold epoch mean','Vhold sweep mean','Vhold epoch trace','Vhold sweep trace'};
epochs = table('Size',[length(epochList),length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);


% Check whether there's AD2 to record Vhold
withVhold = ~isempty(dir(fullfile(options.rawDataPath,"AD2*.mat")));
withVholdAvg = length(epochList)==length(vholdList);

%% Load individual epoch
cellid = 0;
for row = 1:size(epochList,1)

    % Load epoch file to find individual sweep.mat
    load(fullfile(options.rawDataPath,epochList{row,1}));
    namesplit = strsplit(epochList{row,1},{'e','p1avg'}); 
    epoch = str2double(namesplit{2});
    sweepAcq = eval(['AD0_e',num2str(epoch),'p1avg.UserData.Components']);

    % Import processed numbers if post QC
    if ~createNew
        included = exp{row,"Included"}{1};
        cellid = exp{row,"Cell"};
        sweepAcq = exp{row,"Sweep names"}{1};
    end
    
    % Initialize some temporary matrix
    sweeps = zeros(length(sweepAcq), size(eval(['AD0_e',num2str(epoch),'p1avg.data']),2));
    processed = zeros(size(sweeps));
    AUCs = zeros(length(sweepAcq),1);
    peaks = zeros(length(sweepAcq),1);
    Rss = zeros(length(sweepAcq),1);
    Rms = zeros(length(sweepAcq),1);
    Cms = zeros(length(sweepAcq),1);

    % Load Vhold for epoch avg file
    if withVholdAvg && withVhold
        load(fullfile(options.rawDataPath,vholdList{row,1}));
        namesplit = strsplit(epochList{row,1},{'e','p1avg'}); 
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
        try load(fullfile(options.rawDataPath,strcat(sweepAcq{k},'.mat'))); 
        catch
            warning(strcat("Sweep ",num2str(sweepAcq{k}), " not saved, skipping this sweep!"));
            continue
        end
        raw_trace = eval([sweepAcq{k},'.data']);
        
        if length(raw_trace) ~= size(sweeps,2)
            warning("Sweep duration is different from epoch avg duration!!");
            continue;
        end
        sweeps(k,:) = raw_trace;


        % Extract quality metrics from header string
        qc = getCellQC(raw_trace,calculate=options.calculateQC,...
                    headerString=eval([sweepAcq{k},'.UserData.headerString']));
        Rss(k) = qc.Rs;
        Rms(k) = qc.Rm;
        Cms(k) = qc.Cm;
        

        % Process trace: mean-subtracted, optional LP
        % Mean subtraction
        mean_subtracted = raw_trace - mean(raw_trace(options.baselineWindow));
        if options.filterSignal
            Fs = 10000; % Sampling frequency  
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
            base2 = mean(y(:,options.baselineWindow),2); baseM2 = repmat(base2,1,size(y,2));
            processed_trace = y - baseM2;
            processed(k,:) = processed_trace;
        else
            processed_trace = mean_subtracted;
            processed(k,:) = processed_trace;
        end


        % Load Vhold traces if needed
        if withVhold
            try load(fullfile(options.rawDataPath,strcat(vholdAcq{k},'.mat'))); 
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
        AUCs(k) = sum(processed_trace(analysisWindow)) / options.outputFs;
        [~,peakIdx] = max(abs(processed_trace(analysisWindow)));
        peakWindowStart = max(1,peakIdx-peakWindowWidth) + analysisWindow(1);
        peakWindowEnd = min(analysisWindow(end),peakIdx+peakWindowWidth) + analysisWindow(1);
        peaks(k) = mean(processed_trace(peakWindowStart:peakWindowEnd));
    end


    % Determine cellID and mean epoch vhold for that cell
    if withVhold
        if withVholdAvg
            vholdEpochMean = mean(vholdEpoch(:,options.baselineWindow),"all");
        else
            vholdEpochMean = mean(vholdSweeps(:,options.baselineWindow),"all"); 
            vholdEpoch = mean(vholdSweeps(:,options.baselineWindow),1);
        end

        if vholdEpochMean < -50 && createNew; cellid = cellid + 1; end

        vholdSweepsMean = mean(vholdSweeps(:,1:20000),2);
    end


    % Remove empty/erraneous sweeps (optional)
    if createNew
        if options.filterSweeps && withVhold
            % rounded_epoch_vhold = roundToTarget(vholdEpochMean,[-70,0,8]);
            % rounded_sweeps_vhold = roundToTarget(vholdSweepsMean,[-70,0,8]);
            % included = (rounded_sweeps_vhold == rounded_epoch_vhold);
            included = ~isoutlier(vholdSweepsMean,'mean');
        else
            included = ones(length(sweepAcq),1);
        end 
    end

    %% Store everything in epochs
    epochs{row,'Session'} = options.rawDataPath;
    epochs{row,'Animal'} = options.animal;
    epochs{row,'Task'} = options.task;
    epochs{row,'Epoch'} = epoch; %string(epochList{row,1});
    epochs{row,'Included'} = num2cell(included,[1 2]);
    epochs{row,'Sweep names'} = num2cell(sweepAcq,[1 2]);
    epochs{row,'Raw sweeps'} = num2cell(sweeps,[1 2]);
    epochs{row,'Processed sweeps'} = num2cell(processed,[1 2]);
    epochs{row,'Peaks'} = num2cell(peaks,[1 2]);
    epochs{row,'AUCs'} = num2cell(AUCs,[1 2]);
    epochs{row,'Rs'} = num2cell(Rss,[1 2]);
    epochs{row,'Rm'} = num2cell(Rms,[1 2]);
    epochs{row,'Cm'} = num2cell(Cms,[1 2]);

    if withVhold
        epochs{row,'Cell'} = cellid;
        epochs{row,'Vhold epoch mean'} = vholdEpochMean;
        epochs{row,'Vhold sweep mean'} = num2cell(vholdSweepsMean,[1 2]);
        epochs{row,'Vhold epoch trace'} = num2cell(vholdEpoch,[1 2]);
        epochs{row,'Vhold sweep trace'} = num2cell(vholdSweeps,[1 2]);
    else
        if createNew; epochs{row,'Cell'} = epoch;
        else; epochs{row,'Cell'} = cellid; end
        epochs{row,'Vhold epoch mean'} = 100;
        epochs{row,'Vhold sweep mean'} = num2cell(100,[1 2]);
        epochs{row,'Vhold epoch trace'} = num2cell(100,[1 2]);
        epochs{row,'Vhold sweep trace'} = num2cell(100,[1 2]);
    end
end

epochs = rmmissing(epochs);

% If vhold never have -70, plot by each epoch
if isscalar(unique(epochs{:,'Cell'}))
    epochs{:,3} = (1:size(epochs,1))';
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

cells = getCellTable(epochs,save=options.save);

%% Define output

varargout(1) = epochs;
varargout(2) = cells;

end