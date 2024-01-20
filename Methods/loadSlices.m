function epochs = loadSlices(expPath,options)


arguments
    expPath string %full path of the session

    options.filterSignal logical = false;
    options.filterSweeps logical = true; % remove sweeps that have different vhold

    options.save logical = true;
    options.reload logical = false;

    options.outputFs double = 10000
    options.timeRange double = [-20,50] % in ms
    options.eventSample double = 10000 % in sample
    options.nArtifactSamples double = 0 % in sample
    options.baselineWindow double = 1:10000;
    options.peakWindow double = 2; % in ms

    options.rawDataPath string
end

%% Setup

if ~isfield(options,'rawDataPath')
    options.rawDataPath = expPath;
end

% Decide reload if session has already been loaded
dirsplit = split(expPath,filesep); expName = dirsplit{end};
if ~isempty(dir(fullfile(expPath,"epochs_*.mat")))
    if ~options.reload
        disp('Loading stop: epochs file found.');
        load(strcat(expPath,filesep,'epochs_',expName,'.mat'));
        return; 
    end
else
    options.reload = true;
end

timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
analysisWindow = (options.eventSample+options.nArtifactSamples) : timeRangeEndSample;
peakWindowWidth = (options.peakWindow/(2*1000)) * options.outputFs;

%% Grab corresponding recordings through epoch average files
epochList = sortrows(struct2cell(dir(fullfile(expPath,['AD0_e*','p1avg.mat'])))',3);
vholdList = sortrows(struct2cell(dir(fullfile(options.rawDataPath,['AD2_e*','p1avg.mat'])))',3);

% Check whether there's AD2 to record Vhold
withVhold = ~isempty(dir(fullfile(options.rawDataPath,"AD2*.mat")));
withVholdAvg = length(epochList)==length(vholdList);

if withVhold
    varTypes = {'string','string','double','cell','cell','cell','cell',...
                'double','cell','cell','cell',...
                'cell','cell',...
                'cell','cell','cell'};
    varNames = {'Session','Epoch','Cell','Included','Sweep names','Raw sweeps','Processed sweeps',...
        'Vhold epoch mean','Vhold sweep mean',...
        'Vhold epoch trace','Vhold sweep trace',...
        'Peaks','AUCs',...
        'Rin','Rs','Cm'};
    epochs = table('Size',[length(epochList),length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
else
    varTypes = {'string','string','double','cell','cell','cell','cell',...
                'double',...
                'cell','cell',...
                'cell','cell','cell'};
    varNames = {'Session','Epoch','Cell','Included','Sweep names','Raw sweeps','Processed sweeps',...
                'Vhold epoch mean',...
                'Peaks','AUCs',...
                'Rin','Rs','Cm'};
    epochs = table('Size',[length(epochList),length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
end

%% Load individual epoch
cellid = 0;
for row = 1:size(epochList,1)

    % Load epoch file to find individual sweep.mat
    load(fullfile(options.rawDataPath,epochList{row,1}));
    namesplit = strsplit(epochList{row,1},{'e','p1avg'}); 
    epoch = str2double(namesplit{2});
    sweepAcq = eval(['AD0_e',num2str(epoch),'p1avg.UserData.Components']);
    
    % Initialize some temporary matrix
    sweeps = zeros(length(sweepAcq), size(eval(['AD0_e',num2str(epoch),'p1avg.data']),2));
    processed = zeros(size(sweeps));
    AUCs = zeros(length(sweepAcq),1);
    peaks = zeros(length(sweepAcq),1);
    Rins = zeros(length(sweepAcq),1);
    Rss = zeros(length(sweepAcq),1);
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


        % Calculate quality metrics
        % Rin
        Rins(k) = getRin(raw_trace,headerString=eval([sweepAcq{k},'.UserData.headerString']));
        % Rs
        % Rss(k) = getRs();
        % Cm
        % Cms(k) = getCm();

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
        AUCs(k) = sum(processed_trace(analysisWindow));
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

        if vholdEpochMean < -50; cellid = cellid + 1; end

        vholdSweepsMean = mean(vholdSweeps(:,1:20000),2);
    end


    % Remove empty/erraneous sweeps (optional)
    if options.filterSweeps && withVhold
        rounded_epoch_vhold = roundToTarget(vholdEpochMean,[-70,0,8]);
        rounded_sweeps_vhold = roundToTarget(vholdSweepsMean,[-70,0,8]);
        included = (rounded_sweeps_vhold == rounded_epoch_vhold);
    else
        included = ones(length(sweepAcq),1);
    end

    %% Store everything in epochs
    epochs{epoch,'Session'} = string(options.rawDataPath);
    epochs{epoch,'Epoch'} = string(epochList{row,1});
    epochs{epoch,'Included'} = num2cell(included,[1 2]);
    epochs{epoch,'Sweep names'} = num2cell(sweepAcq,[1 2]);
    epochs{epoch,'Raw sweeps'} = num2cell(sweeps,[1 2]);
    epochs{epoch,'Processed sweeps'} = num2cell(processed,[1 2]);
    epochs{epoch,'Peaks'} = num2cell(peaks,[1 2]);
    epochs{epoch,'AUCs'} = num2cell(AUCs,[1 2]);
    epochs{epoch,'Rin'} = num2cell(Rins,[1 2]);
    epochs{epoch,'Rs'} = num2cell(Rss,[1 2]);
    epochs{epoch,'Cm'} = num2cell(Cms,[1 2]);

    if withVhold
        epochs{epoch,'Cell'} = cellid;
        epochs{epoch,'Vhold epoch mean'} = vholdEpochMean;
        epochs{epoch,'Vhold sweep mean'} = num2cell(vholdSweepsMean,[1 2]);
        epochs{epoch,'Vhold epoch trace'} = num2cell(vholdEpoch,[1 2]);
        epochs{epoch,'Vhold sweep trace'} = num2cell(vholdSweeps,[1 2]);
    else
        epochs{epoch,'Cell'} = epoch;
        epochs{epoch,'Vhold epoch mean'} = 100;
    end
end

epochs = rmmissing(epochs);

% If vhold never have -70, plot by each epoch
if length(unique(epochs{:,'Cell'})) == 1
    epochs{:,3} = (1:size(epochs,1))';
end

if options.save
    save(strcat(expPath,filesep,'epochs_',expName),'epochs','-v7.3');
    disp(strcat("Saved: ",expName));
end

end