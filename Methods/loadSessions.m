function [] = loadSessions(sessionpath,options)

arguments
   sessionpath string

   % Record params
   options.NISetup string = 'Shun'
   options.withPhotometryNI logical = false
   options.photometryNI_mod logical = false
   options.photometryNI_modFreq double = 0
   
   % Labjack concat options
   options.labjackSetup string = 'Shun'
   options.followOriginal logical = true

   % Neuropixel recording options
   options.noSyncPulse logical = false
   
   % Reload options (unfinished)
   options.reloadAll logical = false    % If false, skip the whole function if sync_.mat is found
   options.reloadNI logical = false     % If false, skip loading NI and use previous loaded ones
   options.reloadLJ logical = false     % If false, skip loading LJ and use previous loaded ones
   options.reloadCam logical = false    % If false, skip loading Cam and use previous loaded ones
   options.reloadImec logical = false   % If false, skip loading Imec and use previous loaded ones

   % Stimulation digital line params
   options.invertStim logical = true % Invert signal from arduino for stimulation
   options.getConsecutive logical = true
   
   % Photometry related params
   options.recordLJ double = [1,1,0]
   options.nSampPerDemodBin double = 1 % for labjack demod (originally specturalWindowNew)
   options.rollingWindowTime double = 60 % in seconds
   options.LPFreq double = 0; % low pass freq (0: no LP)
   options.downsampleFs double = 50 % downsample NI photometry
   options.dsMethod string = 'resample' % check downsampleSignal.m for explaination
   options.movingAvergeFilter logical = false % moving average filter after downsample
   options.movingAverageWindowSize double = 2 % moving average window size after downsample
   options.removeTwoEnds logical = false % see demodulateSignal.m
   options.plotPhotometry logical = true % Plot photometry signals or not
   options.modFreq double

   % Sync related params
   options.syncPulseWindow double = 5000 % #of sync pulse to xcorr

   % Output file name
   options.outputName string % use the same name as the session folder if emtpy
   options.outputPath string % should be a folder, outputs to the session folder if empty
   
end

%% Notes
% 2023/09/05
%   1. previous oscillation in LJ is due to LP filter, disabled it by default
%   2. add default downsample for LJ if without modulation (default to 50Hz)
%   3. Can partially analyzed by selecting which system to reload

% 2023/10/13
%   1. package code for concat LJ into a function (concatLabjack)

%% Load file

% Select session
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; 
if ispc; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
elseif isunix; session.projectPath = strcat('/',fullfile(dirsplit{2:end-1}));
end
dirsplit = strsplit(sessionName,{'-','_'});
date = dirsplit{1}; animal = dirsplit{2}; sessionTask = dirsplit{3};
clear dirsplit

disp(strcat('**********',sessionName,'**********'));
session.name = sessionName; session.path = sessionpath;
session.date = date;        session.animal = animal;
session.sessionTask = sessionTask;

withNI = ~isempty(dir(fullfile(session.path,'*.nidq.bin')));
withRecording = ~isempty(dir(fullfile(session.path,'*_imec0')));
withCamera = ~isempty(dir(fullfile(session.path,'times_cam1*.csv')));
withPhotometry = isfolder(fullfile(session.path,'Photometry'));
session.nSystems = sum([withRecording,withCamera,withPhotometry,withNI]);

% Decide reload if session has already been loaded and synced
if ~isempty(dir(fullfile(session.path,"sync_*.mat"))) && ~any([options.reloadAll,options.reloadNI,options.reloadLJ,options.reloadCam,options.reloadImec])
    disp('Loading stop: sync file found.'); 
    return; 
else
    % Name all output files
    if ~isfield(options,'outputName'); options.outputName = sessionName; end
    if ~isfield(options,'outputPath'); options.outputPath = sessionpath; end
    syncOutputName = strcat(options.outputPath,filesep,'sync_',options.outputName);
    timeseriesOutputName = strcat(options.outputPath,filesep,'timeseries_',options.outputName);
    dataOutputName = strcat(options.outputPath,filesep,'data_',options.outputName);
    behaviorOutputName = strcat(options.outputPath,filesep,'behavior_',options.outputName);

    % Create folder if neccessary
    if ~exist(options.outputPath,'dir'); mkdir(options.outputPath); end

    % Initialize all .mat files
    options.reloadAll = true;
    save(behaviorOutputName,'sessionName','session','-v7.3');
    save(syncOutputName,'sessionName','session','-v7.3');
    if withPhotometry; save(timeseriesOutputName,'sessionName','session','-v7.3'); end
    save(dataOutputName,'sessionName','session','-v7.3');
end

% Check whether there's multiple recordings in a session
if withCamera
    camerapath = dir(fullfile(session.path,'times_cam1*.csv'));
    if size(camerapath,1) > 1
        warning('More than 1 camera recording found. Skipped!');
        withCamera = 0;
    end
end

% Check event systems
% Either or both Labjack or NIDAQ has to be present
if ~withNI && ~withPhotometry
    error('No NI and LJ recordings found!');
elseif withNI
    session.baselineSystem = 'NI';
elseif withPhotometry && ~withNI
    session.baselineSystem = 'LJ';
    options.withPhotometryNI = false;
end

%% Set up params

if withNI
    %% Load session paths
    session.pathNidq = strcat(session.path,filesep);
    session.nidqBin = strcat(sessionName,'_t0.nidq.bin');
    nidq.meta = ReadMeta(session.nidqBin, session.pathNidq);
    nidq.Fs = str2double(nidq.meta.niSampRate);
    if withPhotometry; session.pathPhotometry = strcat(session.path,filesep,'Photometry',filesep); end 
    
    % NI settings
    params.sync.behaviorFs = nidq.Fs;
    %{ Set time blocks for NIDAQ data
    % This is needed cause data is too big to load as a whole
    % blockTime = 60; %s
    % nBlocks = 20; % total 1200 sec
    % ap.nSampPerBlock = floor(blockTime * ap.Fs);
    % ap.totalSampIncluded = nBlocks * ap.nSampPerBlock;
    % % lfp.nSampPerBlock = floor(blockTime * lfp.Fs);
    % % lfp.totalSampIncluded = nBlocks * lfp.nSampPerBlock;
    % nidq.nSampPerBlock = floor(blockTime * nidq.Fs);
    % nidq.totalSampIncluded = nBlocks * nidq.nSampPerBlock; 
    %}
    
    % Use the whole session
    nidqTotalSecs = str2double(nidq.meta.fileSizeBytes) / nidq.Fs / str2double(nidq.meta.nSavedChans) / 2;
    nidq.totalSampIncluded = floor(nidqTotalSecs * nidq.Fs);
    blockTime = nidqTotalSecs;
    nBlocks = 1;
    session.blockTime = blockTime;
    session.nBlocks = nBlocks;
    session.totalBehaviorSamp = nidq.totalSampIncluded;

    %% Read NIDAQ data
    
    nidqData = ReadBin(0, nidq.totalSampIncluded, nidq.meta, session.nidqBin, session.pathNidq);
    
    % Read Analog Channels in NIDAQ data
    % For an analog channel: gain correct saved channel ch (1-based for MATLAB)
    a_ch_start = str2double(nidq.meta.niXAChans1(1));
    a_ch_end = str2double(nidq.meta.niXAChans1(3));
    a_ch_list = a_ch_start:a_ch_end;
    
    analogNI = GainCorrectNI(nidqData, a_ch_list+1, nidq.meta); % remove +1 if analyzing session before 0512
    
    % for i = 1:numel(a_ch_list)
    %     figure(10+i);
    %     plot(analogNI(a_ch_list(i)+1,:)); % remove +1 if analyzing session before 0512
    %     title(['NIDAQ analog channel ',num2str(a_ch_list(i))]);
    % end
    % autoArrangeFigures();

    % for i = 1:size(analogNI,1)
    %     figure;
    %     plot(analogNI(i,:));
    %     title(i);
    % end
    % autoArrangeFigures();
    
    % Read Digital Channels in NIDAQ data
    d_ch_start = str2double(nidq.meta.niXDChans1(1));
    d_ch_end = str2double(nidq.meta.niXDChans1(3));
    d_ch_list = d_ch_start:d_ch_end;
    
    % For a digital channel: read this digital word dw in the saved file
    % (1-based). For imec data there is never more than one saved digital word.
    dw = 1;
    digitalNI = ExtractDigital(nidqData, nidq.meta, dw, d_ch_list);
    % syncNI = digitalNI(2,:); % Save sync channel from NIDAQ data separately
    
    close all;

    %% Read behavior/photometry (NI) data
    
    if options.reloadAll || options.reloadNI
        channels = assignNIDAQChannel(analogNI, digitalNI,user=options.NISetup,...
            invertStim=options.invertStim,getConsecutive=options.getConsecutive);
        
        % Save assigned channels separately for later analysis
        if strcmpi(options.NISetup,'Shun')
            leftLick = channels{1}; rightLick = channels{2};
            blink = channels{3}; gyro = channels{4}; photometry_raw = channels{5};

            airpuff = channels{6}; allTones = channels{13};
            leftSolenoid = channels{7}; rightSolenoid = channels{8};
            leftTone = channels{9}; rightTone = channels{10};
            redLaser = channels{11}; blueLaser = channels{12};
            syncNI = channels{14};
        elseif strcmpi(options.NISetup,'Kevin')
            syncNI = channels{1}; leftLick = channels{2};
            rightLick_analog = channels{3}; leftNose = channels{4};

            rightLED = channels{5}; rightSolenoid = channels{6};
            rightLick = channels{7}; rightNose = channels{8};
            centerLED = channels{9}; centerNose = channels{10};
            leftLED = channels{11}; leftSolenoid = channels{12};
        end
        
        
        %% Save behavioral data
        if strcmpi(options.NISetup,'Shun')
            if options.withPhotometryNI
                save(dataOutputName,...
                    'airpuff','leftLick','rightLick','leftTone','rightTone',....
                    'leftSolenoid','rightSolenoid','allTones','blink','gyro',...
                    'photometry_raw','blueLaser','redLaser','-append');
            else
                save(dataOutputName,...
                    'airpuff','leftLick','rightLick','leftTone','rightTone',....
                    'leftSolenoid','rightSolenoid','allTones','blink','gyro',...
                    'blueLaser','redLaser','-append');
            end
        elseif strcmpi(options.NISetup,'Kevin')
            save(dataOutputName,...
                    'leftLick','rightLick','leftNose','rightNose',....
                    'leftSolenoid','rightSolenoid','leftLED','rightLED','centerLED',...
                    'centerNose','rightLick_analog','-append');
        end
        disp("Finished: NIDAQ data saved in data_.mat");
    end
end

%% Read ephys (imec) data
if withRecording && (options.reloadAll || options.reloadImec)
    %% Load imec data
    if ~isempty(dir(fullfile(session.path,'catgt_*')))
        % For old catgt recordings
        session.pathImec = strcat(session.path,filesep,'catgt_', sessionName,filesep);
        session.apBin = strcat(sessionName,'_tcat.imec0.ap.bin');
    else
        % For new spikeinterface recordings
        session.pathImec = strcat(session.path,filesep,sessionName,'_imec0');
        session.apBin = strcat(sessionName,'_t0.imec0.ap.bin');
    end
    
    ap.meta = ReadMeta(session.apBin, session.pathImec);
    params.sync.apFs = str2double(ap.meta.imSampRate);

    apTotalSecs = str2double(ap.meta.fileSizeBytes) / params.sync.apFs / str2double(ap.meta.nSavedChans) / 2;
    ap.totalSampIncluded = floor(apTotalSecs * params.sync.apFs);
    session.totalBehaviorSamp = nidq.totalSampIncluded;
    session.totalImecSamp = ap.totalSampIncluded;

    %% Load LFP data
    session.pathLFP = strcat(session.path,filesep,sessionName,'_imec0',filesep);
    session.lfpBin = strcat(sessionName,'_t0.imec0.lf.bin');
    lfp.meta = ReadMeta(session.lfpBin,session.pathLFP);
    params.sync.lfpFs = str2double(lfp.meta.imSampRate);
    params.ephys.lfpFs = params.sync.lfpFs;

    lfpTotalSecs = str2double(lfp.meta.fileSizeBytes) / params.sync.lfpFs / str2double(lfp.meta.nSavedChans) / 2;
    lfp.totalSampIncluded = floor(lfpTotalSecs * params.sync.lfpFs);
    session.totalLFPSamp = lfp.totalSampIncluded;

    %% Read in sync channel
    % Normally takes ~10 seconds for 60s data; 230s for 1200s data
    if ~options.noSyncPulse
        tic
        disp('Ongoing: reading sync channel of Imec...');
        syncImec = ReadBinByCh(0, ap.totalSampIncluded, ap.meta, session.apBin, session.pathImec, 385);
        disp('Ongoing: reading sync channel of LFP...');
        syncLFP = ReadBinByCh(0, lfp.totalSampIncluded, lfp.meta, session.lfpBin, session.pathLFP, 385);
        disp('time for reading sync channel from Imec&LFP data');
        toc
    end
    
    %% Save ephys data
    save(dataOutputName,'ap','lfp','-append');
    disp('Finished: Neuropixel data saved in data_.mat');
end

%% Read photometry data
if (withPhotometry || options.withPhotometryNI) && (options.reloadAll || options.reloadLJ)
    %% Concatenate raw.mat files

    if withPhotometry
        disp("Ongoing: concatenate and save raw photometry data");
        if isempty(dir(fullfile(session.path,filesep,'data_labjack.mat')))
            if strcmpi(options.labjackSetup,"Shun")
                concatLabjack(session.path,save=true,plot=false,record=options.recordLJ,followOriginal=options.followOriginal,...
                              outputPath=options.outputPath);
            elseif strcmpi(options.labjackSetup,"Shijia")
                concatLabjack_shijia(session.path,save=true,plot=false,record=options.recordLJ);
            end
        end
        load(strcat(session.path,filesep,'data_labjack.mat'));
        
        % Reload concatLabjack if labjack.record does not agree
        updateLabjack = false;
        if ~isfield(labjack,'record')
            disp('     Did not find labjack.record, use options.reloadLJ instead');
            labjack.record = options.recordLJ; 
            labjack.nSignals = sum(labjack.record);
            % Remove non-recorded channels
            labjack.raw(find(~labjack.record),:) = [];
            labjack.modulation(find(~labjack.record),:) = [];
            labjack.name(find(~labjack.record)) = [];
            updateLabjack = true;
        end
        if sum(labjack.record == options.recordLJ) ~= length(labjack.record)
            disp(['labjack.record: ',labjack.record]);
            disp(['options.recordLJ: ',options.recordLJ]);
            warning(['labjack.record does not agree with recordLJ, use labjack.record = ',num2str(labjack.record)]); 
        end
        disp('Finished: concatenate and saved raw photometry data in data_labjack.mat');
        disp(labjack);

        % Update modFreq to user input if neccessary
        if ~isfield(options,'modFreq'); options.modFreq = labjack.modFreq; end
    
        % Store relevant info
        params.photometry.params = labjack;
        params.sync.labjackFs = labjack.samplerate;
        totalDuration_LJ = length(sync_labjack) / params.sync.labjackFs;
        params.photometry.totalDuration = totalDuration_LJ;
        params.sync.photometryFs = [];

        % Remove PMT if neccessary
        if ~options.withPhotometryNI
            PMT_idx = find(cellfun(@(c)strcmp(c,"PMT"),labjack.name));
            if labjack.record(PMT_idx)
                labjack.nSignals = labjack.nSignals-1; 
                labjack.raw(3,:) = [];
                labjack.modulation(3,:) = [];
                updateLabjack = true;
            end
        else
            params.sync.ni_photometryFs = [];
        end

        % Read behavior from Labjack if neccessary
        if isfield(labjack,'recordBehavior') && labjack.recordBehavior
            structvars(labjack.behavior);
            save(dataOutputName,...
                'airpuff','leftLick','rightLick','leftTone','rightTone',....
                'leftSolenoid','rightSolenoid','allTones','blink','gyro',...
                'blueLaser','redLaser','-append');
        end

        % Save updated labjack
        if updateLabjack; save(strcat(session.path,filesep,'data_labjack'),'labjack','-append'); end
        save(dataOutputName,'labjack','-append');
    end

    %% Loop through all signals
    if withPhotometry
        for i = 1:size(labjack.raw,1)
            
            if labjack.mod(i)
                processed = demodulateSignal(labjack.raw(i,:),...
                                targetFs=options.downsampleFs,...
                                originalFs=params.sync.labjackFs,...
                                modFreq=options.modFreq(i),...
                                removeTwoEnds=options.removeTwoEnds,...
                                rollingWindowTime=options.rollingWindowTime);
    
                % Check final Fs
                finalFs = length(processed.demodData) / totalDuration_LJ;
                if options.downsampleFs ~= finalFs
                    disp(['finalFs: ',num2str(finalFs)]);
                    disp(['targetFs: ',num2str(options.downsampleFs)]);
                    warning([labjack.name{i},': Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!']);
                end
                params.sync.photometryFs = [finalFs, params.sync.photometryFs];
    
                % Store process information
                timeSeries(i).name = labjack.name{i};
                timeSeries(i).data = processed.demodData;
                timeSeries(i).finalFs = finalFs;
                timeSeries(i).system = 'LJ';
                timeSeries(i).time_offset = NaN;
                timeSeries(i).demux = true;
                timeSeries(i).demux_freq = processed.options.modFreq;
                timeSeries(i).detrend = true;
                timeSeries(i).detrend_type = 'rolling-z';
                timeSeries(i).detrend_window = processed.options.rollingWindowTime;
                timeSeries(i).options = processed.options;
        
            else
                disp("Ongoing: Process unmodulated photometry data");
                % Downsample
                processed = downsampleSignal(labjack.raw(i,:),...
                            targetFs=options.downsampleFs,...
                            originalFs=params.sync.labjackFs,...
                            rollingZ=true,rollingWindowTime=options.rollingWindowTime,...
                            dsMethod='resample');
                
                % Check final Fs
                finalFs = length(processed.dsData) / totalDuration_LJ;
                if options.downsampleFs ~= finalFs
                    disp(['finalFs: ',num2str(finalFs)]);
                    disp(['targetFs: ',num2str(options.downsampleFs)]);
                    warning([labjack.name{i},': Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!']);
                end
                params.sync.photometryFs = [finalFs, params.sync.photometryFs];
                
                % Store process information
                timeSeries(i).name = labjack.name{i};
                timeSeries(i).data = processed.dsData;
                timeSeries(i).finalFs = finalFs;
                timeSeries(i).system = 'LJ';
                timeSeries(i).time_offset = NaN;
                timeSeries(i).demux = true;
                timeSeries(i).demux_freq = NaN;
                timeSeries(i).detrend = true;
                timeSeries(i).detrend_type = 'rolling-z';
                timeSeries(i).detrend_window = processed.options.rollingWindowTime;
                timeSeries(i).options = processed.options;
                disp("Finished: no modulation, skip demodulation and downsampled");
            
                % Low pass filter data (skip by default)
                if options.LPFreq > 0
                    params.photometry.nSampPerDemodBin = 41;
                    params.photometry.filtCut = options.LPFreq;
                    lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',params.photometry.filtCut,...
                        'PassbandRipple',0.01, 'Samplerate',params.sync.labjackFs/params.photometry.nSampPerDemodBin);
                    % Lowpass filter
                    dsGreenLP = filtfilt(lpFilt,double(dsGreen));
                    rollingGreenLP = rollingZ(dsGreenLP,options.rollingWindowTime);
            
                    save(dataOutputName,'dsGreenLP','rollingGreenLP','-append');
                end
            end
        end
    end

    % Process NIDAQ photometry data
    if options.withPhotometryNI == true
        params.sync.ni_photometryFs = [];
        totalDuration_NI = length(photometry_raw) / params.sync.behaviorFs;
        disp("Ongoing: Process NI photometry data");

        % Find row to store
        if ~withPhotometry; i = 1;
        else; i = size(timeSeries,2) + 1; end

        % Process
        if options.photometryNI_mod
            processed = demodulateSignal(photometry_raw,...
                            targetFs=options.downsampleFs,...
                            originalFs=nidq.Fs,...
                            modFreq=options.photometryNI_modFreq,...
                            removeTwoEnds=options.removeTwoEnds,...
                            rollingWindowTime=options.rollingWindowTime);

            % Check final Fs
            finalFs = length(processed.demodData) / totalDuration_NI;
            if options.downsampleFs ~= finalFs
                disp(['finalFs: ',num2str(finalFs)]);
                disp(['targetFs: ',num2str(options.downsampleFs)]);
                warning('PMT: Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!');
            end
            params.sync.ni_photometryFs = [finalFs, params.sync.ni_photometryFs];

            % Store params
            timeSeries(i).name = 'PMT';
            timeSeries(i).data = processed.demodData;
            timeSeries(i).finalFs = finalFs;
            timeSeries(i).system = 'NI';
            timeSeries(i).time_offset = 0;
            timeSeries(i).demux = true;
            timeSeries(i).demux_freq = processed.options.modFreq;
            timeSeries(i).detrend = true;
            timeSeries(i).detrend_type = 'rolling-z';
            timeSeries(i).detrend_window = processed.options.rollingWindowTime;
            timeSeries(i).options = processed.options;
            disp('Finished: NI photometry demodulation');

        else
            % Downsample
            processed = downsampleSignal(photometry_raw,...
                                targetFs=options.downsampleFs,originalFs=nidq.Fs,...
                                rollingWindowTime=options.rollingWindowTime,...
                                dsMethod=options.dsMethod);

            % Check final Fs
            finalFs = length(processed.dsData) / totalDuration_NI;
            if options.downsampleFs ~= finalFs
                disp(['finalFs: ',num2str(finalFs)]);
                disp(['targetFs: ',num2str(options.downsampleFs)]);
                warning("PMT: Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!");
            end
            params.sync.ni_photometryFs = [finalFs, params.sync.ni_photometryFs];
            
            % Store params
            timeSeries(i).name = 'PMT';
            timeSeries(i).data = processed.dsData;
            timeSeries(i).finalFs = finalFs;
            timeSeries(i).system = 'NI';
            timeSeries(i).time_offset = 0;
            timeSeries(i).demux = false;
            timeSeries(i).demux_freq = NaN;
            timeSeries(i).detrend = true;
            timeSeries(i).detrend_type = 'rolling-z';
            timeSeries(i).detrend_window = processed.options.rollingWindowTime;
            timeSeries(i).options = processed.options;
            disp('Finished: NI photometry downsampling');
        end
    end

    %% Plot signals and save channels

    if options.plotPhotometry
        initializeFig(1,1); tiledlayout(size(timeSeries,2),2);

        if options.withPhotometryNI && withPhotometry
            totalDurationMin = min([totalDuration_NI,totalDuration_LJ]);
        elseif withPhotometry; totalDurationMin = totalDuration_LJ; 
        else; totalDurationMin = totalDuration_NI; 
        end

        if totalDurationMin < 30; plotDuration = totalDurationMin; 
        else; plotDuration = 30; end % in sec

        for i = 1:size(timeSeries,2)
            if strcmp(timeSeries(i).system,'NI')
                rawTrace = photometry_raw; 
                rawFs = params.sync.behaviorFs;
                totalDuration = totalDuration_NI;
            else
                rawTrace = labjack.raw(i,:);
                rawFs = params.sync.labjackFs;
                totalDuration = totalDuration_LJ;
            end
            meanTrace = movmean(rawTrace,60*rawFs);

            % Plot 30s sample
            nexttile;
            yyaxis left % Plot raw green
            plotSamples = rawFs * plotDuration;
            plot(linspace(1,plotDuration,plotSamples),rawTrace(1:plotSamples),Color=[0.3010 0.7450 0.9330]); hold on
            yyaxis right % Plot processed green
            plotSamples = floor(timeSeries(i).finalFs * plotDuration);
            plot(linspace(1,plotDuration,plotSamples),timeSeries(i).data(1:plotSamples),...
                LineWidth=2); hold on
            legend({'Raw','Processed'}); xlabel('Time (s)');
            title(timeSeries(i).name);
            % Plot all session
            nexttile;
            yyaxis left % Plot raw all session
            % plotSamples = rawFs * totalDuration;
            plot(linspace(1,totalDuration,length(rawTrace)),rawTrace,Color=[0.3010 0.7450 0.9330]); hold on
            plot(linspace(1,totalDuration,length(meanTrace)),meanTrace,Color='b',LineWidth=2); hold on
            yyaxis right % Plot processed all session
            % plotSamples = floor(timeSeries(i).finalFs * totalDuration);
            plot(linspace(1,totalDuration,length(timeSeries(i).data)),timeSeries(i).data,Color='#D95319',LineWidth=2); hold on
            legend({'Raw','Mean','Processed'}); xlabel('Time (s)');
            title(timeSeries(i).name);
        end

        % Save channels
        saveas(gcf,strcat(options.outputPath,filesep,'Summary_photometry_processing.png')); 
        % saveas(gcf,strcat(options.outputPath,filesep,'Summary_photometry_processing.fig'));
    end
    disp("Finished: processed and saved photometry data in data_.mat");
end

%% Read camera data
if withCamera && (options.reloadAll || options.reloadCam)
    camerapath = dir(fullfile(session.path,'times_cam1*.csv'));
    csvpath = [camerapath.folder,filesep,camerapath.name];

    % Read camera table
    opts = detectImportOptions(csvpath);

    % Determine whether there's eye tracking built in
    if length(opts.VariableNames) == 8
        session.withEyeTracking = true;
        opts.VariableNames = {'Sync','FrameCounter','Time','Date',...
                              'PupilArea','EyeArea','Blink','EyePixelIntensity'};
        opts.VariableTypes{end} = 'logical';
        opts.VariableTypes{4} = 'datetime';
        opts.VariableOptions(1,4).DatetimeFormat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSxxxxx';
        opts.VariableOptions(1,4).InputFormat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSxxxxx';
        opts.VariableOptions(1,4).TimeZone = 'America/New_York';
    
        camera = readtable(csvpath,opts);
        camera{:,end+1} = zeros(size(camera,1),1);
        camera.Properties.VariableNames{end} = 'SkippedFrame';
        disp("Finished: read camera csv file");
    
        % Calculate frame rate 
        camFs = 1/mean(seconds(diff(camera{:,'Date'})));
        params.sync.camFs = camFs;
    
        % Check skipped frames
        skipped_frame = find(diff(camera{:,2}) > 1);
        if ~isempty(skipped_frame) && length(skipped_frame) <= 1000
            disp(['Ongoing: filling ', num2str(length(skipped_frame)), ' skipped frames']);
        
            % Fill in skip frames with the previous frame
            fill_in = [];
            for i = 1:length(skipped_frame)
                prev_frame = camera{skipped_frame(i),2};
                post_frame = camera{skipped_frame(i)+1,2};
                frame_diff = post_frame - prev_frame;
        
                fill = repelem(camera(skipped_frame(i),:),frame_diff-1,1);
                fill{:,'SkippedFrame'} = 1;
                fill{:,'FrameCounter'} = linspace(prev_frame+1,post_frame-1,frame_diff-1)';
                fill_in = [fill_in;fill];
            end
        
            % Re-check again
            camera = sortrows([camera;fill_in],2);
            skipped_frame = find(diff(camera{:,2}) > 1);
            if ~isempty(skipped_frame)
                disp(['Found ', num2str(length(skipped_frame)), ' skipped frames after filling!']);
            end
        end
    
        if isempty(skipped_frame)
            % Process eye, pupil area data
            pupilArea = camera{:,5}; eyeArea = camera{:,6}; 
            eyePixelIntensity = -camera{:,8}; % convert sign to make blinking going down

            traces = {pupilArea; eyeArea; eyePixelIntensity};
            trace_names = {'pupilArea';'eyeArea';'eyePixelIntensity'};
            for i = 1:size(traces,1)
                processed_cam = downsampleSignal(traces{i},...
                                targetFs=options.downsampleFs,...
                                originalFs=params.sync.camFs,...
                                rollingZ=true,rollingWindowTime=options.rollingWindowTime,...
                                dsMethod='cloestInteger');

                % Save to timeseries struct
                row = size(timeSeries,2) + 1;
                timeSeries(row).name = trace_names{i};
                timeSeries(row).data = processed_cam.dsData;
                timeSeries(row).finalFs = options.downsampleFs;
                timeSeries(row).system = 'Cam';
                timeSeries(row).time_offset = 0;
                timeSeries(row).demux = false;
                timeSeries(row).demux_freq = NaN;
                timeSeries(row).detrend = true;
                timeSeries(row).detrend_type = 'rolling-z';
                timeSeries(row).detrend_window = processed_cam.options.rollingWindowTime;
                timeSeries(row).options = processed_cam.options;
            end
            disp("Finished: detrend eye, pupil area");

            % Save camera table
            save(dataOutputName,'camera',...
                'pupilArea','eyeArea','-append');
            disp("Finished: saved camera csv file in data_.mat");
        else
            disp(['Warning: Too much skipped frames (', num2str(length(skipped_frame)), ') were found, skipped camera anlysis!']);
            withCamera = 0;
        end


    elseif length(opts.VariableNames) == 7
        session.withEyeTracking = true;
        opts.VariableNames = {'Sync','FrameCounter','Time','Date',...
                              'PupilArea','EyeArea','Blink'};
        opts.VariableTypes{end} = 'logical';
        opts.VariableTypes{4} = 'datetime';
        opts.VariableOptions(1,4).DatetimeFormat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSxxxxx';
        opts.VariableOptions(1,4).InputFormat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSxxxxx';
        opts.VariableOptions(1,4).TimeZone = 'America/New_York';
    
        camera = readtable(csvpath,opts);
        camera{:,end+1} = zeros(size(camera,1),1);
        camera.Properties.VariableNames{end} = 'SkippedFrame';
        disp("Finished: read camera csv file");
    
        % Calculate frame rate 
        camFs = 1/mean(seconds(diff(camera{:,'Date'})));
        params.sync.camFs = camFs;
    
        % Check skipped frames
        skipped_frame = find(diff(camera{:,2}) > 1);
        if ~isempty(skipped_frame) && length(skipped_frame) <= 1000
            disp(['Ongoing: filling ', num2str(length(skipped_frame)), ' skipped frames']);
        
            % Fill in skip frames with the previous frame
            fill_in = [];
            for i = 1:length(skipped_frame)
                prev_frame = camera{skipped_frame(i),2};
                post_frame = camera{skipped_frame(i)+1,2};
                frame_diff = post_frame - prev_frame;
        
                fill = repelem(camera(skipped_frame(i),:),frame_diff-1,1);
                fill{:,'SkippedFrame'} = 1;
                fill{:,'FrameCounter'} = linspace(prev_frame+1,post_frame-1,frame_diff-1)';
                fill_in = [fill_in;fill];
            end
        
            % Re-check again
            camera = sortrows([camera;fill_in],2);
            skipped_frame = find(diff(camera{:,2}) > 1);
            if ~isempty(skipped_frame)
                disp(['Found ', num2str(length(skipped_frame)), ' skipped frames after filling!']);
            end
        end
    
        if isempty(skipped_frame)
            % Process eye, pupil area data
            pupilArea = camera{:,5}; eyeArea = camera{:,6};

            traces = {pupilArea; eyeArea};
            trace_names = {'pupilArea';'eyeArea'};
            for i = 1:size(traces,1)
                processed_cam = downsampleSignal(traces{i},...
                                targetFs=options.downsampleFs,...
                                originalFs=params.sync.camFs,...
                                rollingZ=true,rollingWindowTime=options.rollingWindowTime,...
                                dsMethod='cloestInteger');

                % Save to timeseries struct
                row = size(timeSeries,2) + 1;
                timeSeries(row).name = trace_names{i};
                timeSeries(row).data = processed_cam.dsData;
                timeSeries(row).finalFs = options.downsampleFs;
                timeSeries(row).system = 'Cam';
                timeSeries(row).time_offset = 0;
                timeSeries(row).demux = false;
                timeSeries(row).demux_freq = NaN;
                timeSeries(row).detrend = true;
                timeSeries(row).detrend_type = 'rolling-z';
                timeSeries(row).detrend_window = processed_cam.options.rollingWindowTime;
                timeSeries(row).options = processed_cam.options;
            end
            disp("Finished: detrend eye, pupil area");

            % Save camera table
            save(dataOutputName,'camera',...
                'pupilArea','eyeArea','-append');
            disp("Finished: saved camera csv file in data_.mat");
        else
            disp(['Warning: Too much skipped frames (', num2str(length(skipped_frame)), ') were found, skipped camera anlysis!']);
            withCamera = 0;
        end

    elseif length(opts.VariableNames) == 4
        session.withEyeTracking = false;
        opts.VariableNames = {'Sync','FrameCounter','Time','Date'};
        opts.VariableTypes{4} = 'datetime';
        opts.VariableOptions(1,4).DatetimeFormat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS';
        opts.VariableOptions(1,4).InputFormat = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSS';
        opts.VariableOptions(1,4).TimeZone = 'America/New_York';
    
        camera = readtable(csvpath,opts);
        camera{:,end+1} = zeros(size(camera,1),1);
        camera.Properties.VariableNames{end} = 'SkippedFrame';
        disp("Finished: read camera csv file");
    
        % Calculate frame rate 
        camFs = 1/mean(seconds(diff(camera{:,'Date'})));
        params.sync.camFs = camFs;
    
        % Check skipped frames
        skipped_frame = find(diff(camera{:,2}) > 1);
        if ~isempty(skipped_frame) && length(skipped_frame) <= 1000
            disp(['Ongoing: filling ', num2str(length(skipped_frame)), ' skipped frames']);
        
            % Fill in skip frames with the previous frame
            fill_in = [];
            for i = 1:length(skipped_frame)
                prev_frame = camera{skipped_frame(i),2};
                post_frame = camera{skipped_frame(i)+1,2};
                frame_diff = post_frame - prev_frame;
        
                fill = repelem(camera(skipped_frame(i),:),frame_diff-1,1);
                fill{:,'SkippedFrame'} = 1;
                fill{:,'FrameCounter'} = linspace(prev_frame+1,post_frame-1,frame_diff-1)';
                fill_in = [fill_in;fill];
            end
        
            % Re-check again
            camera = sortrows([camera;fill_in],2);
            skipped_frame = find(diff(camera{:,2}) > 1);
            if ~isempty(skipped_frame)
                disp(['Found ', num2str(length(skipped_frame)), ' skipped frames after filling!']);
            end
        end
    end
end

%% (Ver5) Sync data using xcorr

% Convert sync signals to digital signal (boolean);
if withNI
    syncNI = (syncNI==1);
    % Extract location of rising edge
    temp = [false,diff(syncNI)]; % diff(syncImec) records time of rising edge (1) and falling edge (0)
    syncNI_diff = (temp==1);
end

if withRecording && ~options.noSyncPulse
    syncImec = (syncImec > (max(syncImec)/2));
    temp = [false,diff(syncImec)];
    syncImec_diff = (temp==1);

    syncLFP = (syncLFP > (max(syncLFP)/2));
    temp = [false,diff(syncLFP)];
    syncLFP_diff = (temp==1);

    if ~any(find(syncImec_diff, 1) & find(syncLFP_diff,1), "all") 
        options.noSyncPulse = true;
    end
end

if withCamera
    syncCam = (camera{:,1} > max(camera{:,1}/2))';
    temp = [false,diff(syncCam)];
    syncCam_diff = (temp==1);
end

if withPhotometry
    sync_labjack = (sync_labjack > (max(sync_labjack)/2));
    temp = [false,diff(sync_labjack)];
    syncLJ_diff = (temp==1);
end

% Extract start of sync pulse
% 1. Loop over first n pulse
syncPulseWindow = options.syncPulseWindow;
ni_idx = []; lj_idx = []; cam_idx = []; imec_idx = []; lfp_idx = [];
if withNI; ni_idx = find(syncNI_diff>0,syncPulseWindow); end
if withPhotometry; lj_idx = find(syncLJ_diff>0,syncPulseWindow); end
if withCamera; cam_idx = find(syncCam_diff>0,syncPulseWindow); end
if withRecording && ~options.noSyncPulse
    imec_idx = find(syncImec_diff>0,syncPulseWindow); 
    lfp_idx = find(syncLFP_diff>0,syncPulseWindow); 
end

minSyncPulse = min(nonzeros([length(ni_idx), length(lj_idx), length(cam_idx), length(imec_idx), length(lfp_idx)]));
if syncPulseWindow > minSyncPulse
    syncPulseWindow = minSyncPulse;
    ni_idx = []; lj_idx = []; cam_idx = []; imec_idx = []; lfp_idx = [];
    if withNI; ni_idx = find(syncNI_diff>0,syncPulseWindow); end
    if withPhotometry; lj_idx = find(syncLJ_diff>0,syncPulseWindow); end
    if withCamera; cam_idx = find(syncCam_diff>0,syncPulseWindow); end
    if withRecording && ~options.noSyncPulse
        imec_idx = find(syncImec_diff>0,syncPulseWindow); 
        lfp_idx = find(syncLFP_diff>0,syncPulseWindow); 
    end
end

% 2. Calculate ISI
if withNI; ISI_ni = (ni_idx(2:end) - ni_idx(1:end-1)) / nidq.Fs; end
if withCamera; ISI_cam = (cam_idx(2:end) - cam_idx(1:end-1)) / params.sync.camFs; end
if withPhotometry; ISI_lj = (lj_idx(2:end) - lj_idx(1:end-1)) / params.sync.labjackFs; end
if withRecording && ~options.noSyncPulse
    ISI_imec = (imec_idx(2:end) - imec_idx(1:end-1)) / params.sync.apFs; 
    ISI_lfp = (lfp_idx(2:end) - lfp_idx(1:end-1)) / params.sync.lfpFs; 
end

% 3. Cross correlation
if strcmpi(session.baselineSystem,'NI')
    if withCamera; [r_cam, l_cam] = xcorr(zscore(ISI_ni), zscore(ISI_cam),'normalized'); end
    if withPhotometry; [r_lj, l_lj] = xcorr(zscore(ISI_ni), zscore(ISI_lj),'normalized'); end
    if withRecording && ~options.noSyncPulse
        [r_imec, l_imec] = xcorr(zscore(ISI_ni), zscore(ISI_imec),'normalized'); 
        [r_lfp,l_lfp] = xcorr(zscore(ISI_ni), zscore(ISI_lfp),'normalized');
    end
elseif strcmpi(session.baselineSystem,'LJ')
    if withCamera; [r_cam, l_cam] = xcorr(zscore(ISI_lj), zscore(ISI_cam),'normalized'); end
end

% 4. Find the max of xcorr
initializeFig(0.67,0.67); tiledlayout('flow');
if strcmpi(session.baselineSystem,'NI')
    Cam_inBaseline = 0; LJ_inBaseline = 0; Imec_inBaseline = 0; LFP_inBaseline = 0;
    Cam_inBaselineIdx = 0; LJ_inBaselineIdx = 0; Imec_inBaselineIdx = 0; LFP_inBaselineIdx = 0;

    if withCamera
        [~,maxIdx] = max(r_cam); lags_cam = l_cam(maxIdx); 
        if lags_cam >= 0; firstPulse_inCam = cam_idx(1); Cam_inBaseline = ni_idx(1+lags_cam);
        else; firstPulse_inCam = cam_idx(1-lags_cam); Cam_inBaseline = ni_idx(1); end
        Cam_inBaselineIdx = max([1+lags_cam,1]); firstPulse_inCamIdx = max([1-lags_cam,1]);
        nexttile; stem(l_cam,r_cam); title('NI vs Cam xcor');
    end
    if withPhotometry
        [~,maxIdx] = max(r_lj); lags_lj = l_lj(maxIdx); 
        if lags_lj >= 0; firstPulse_inLJ = lj_idx(1); LJ_inBaseline = ni_idx(1+lags_lj);
        else; firstPulse_inLJ = lj_idx(1-lags_lj); LJ_inBaseline = ni_idx(1); end
        LJ_inBaselineIdx = max([1+lags_lj,1]); firstPulse_inLJIdx = max([1-lags_lj,1]);
        nexttile; stem(l_lj,r_lj); title('NI vs Labjack xcor');
    end
    if withRecording
        if ~options.noSyncPulse
            [~,maxIdx] = max(r_imec); lags_imec = l_imec(maxIdx); 
            [~,maxIdx] = max(r_lfp); lags_lfp = l_lfp(maxIdx);
            if lags_imec >= 0
                firstPulse_inImec = imec_idx(1); Imec_inBaseline = ni_idx(1+lags_imec);
                firstPulse_inLFP = lfp_idx(1); LFP_inBaseline = ni_idx(1+lags_lfp);
            else
                firstPulse_inImec = imec_idx(1-lags_imec); Imec_inBaseline = ni_idx(1);
                firstPulse_inLFP = lfp_idx(1-lags_lfp); LFP_inBaseline = ni_idx(1);
            end
            Imec_inBaselineIdx = max([1+lags_imec,1]); firstPulse_inImecIdx = max([1-lags_imec,1]);
            LFP_inBaselineIdx = max([1+lags_lfp,1]); firstPulse_inLFPIdx = max([1-lags_lfp,1]);
            nexttile; stem(l_imec,r_imec); title('NI vs Imec xcor');
            nexttile; stem(l_lfp,r_lfp); title('NI vs LFP xcor');
        end
    end
elseif strcmpi(session.baselineSystem,'LJ')
    Cam_inBaseline = 0; LJ_inBaseline = 0; Imec_inBaseline = 0; LFP_inBaseline = 0;
    Cam_inBaselineIdx = 0; LJ_inBaselineIdx = 0; Imec_inBaselineIdx = 0; LFP_inBaselineIdx = 0;

    if withCamera
        [~,maxIdx] = max(r_cam); lags_cam = l_cam(maxIdx); 
        if lags_cam >= 0; firstPulse_inCam = cam_idx(1); Cam_inBaseline = lj_idx(1+lags_cam);
        else; firstPulse_inCam = cam_idx(1-lags_cam); Cam_inBaseline = lj_idx(1); end
        Cam_inBaselineIdx = max([1+lags_cam,1]); firstPulse_inCamIdx = max([1-lags_cam,1]);
        nexttile; stem(l_cam,r_cam); title('LJ vs Cam xcor');
    end
end

% 5. Find the first common sync pulse of all system
if ~(withRecording || withPhotometry || withCamera) % Only nidq
    firstPulse_inNI = find(syncNI_diff>0,1);
elseif ~(withRecording || withNI || withCamera) % Only Labjack
    firstPulse_inLJ = find(syncLJ_diff>0,1);
    if isempty(firstPulse_inLJ); firstPulse_inLJ = 1; end
else
    firstPulse_inNI = find(syncNI_diff>0,1);
    firstSamp = [firstPulse_inNI,Imec_inBaseline,Cam_inBaseline,LJ_inBaseline,LFP_inBaseline];
    firstIdx = [firstPulse_inNI,Imec_inBaselineIdx,Cam_inBaselineIdx,LJ_inBaselineIdx,LFP_inBaselineIdx];
    if strcmpi(session.baselineSystem,'NI')
        [firstPulse_inNI,system] = max(firstSamp);
        if system ~= 1 % if some system comes online after NI
            if withCamera
                diffIdx = firstIdx(system) - Cam_inBaselineIdx; 
                firstPulse_inCam = cam_idx(firstPulse_inCamIdx + diffIdx);
            end
            if withPhotometry
                diffIdx = firstIdx(system) - LJ_inBaselineIdx; 
                firstPulse_inLJ = lj_idx(firstPulse_inLJIdx + diffIdx);
            end
            if withRecording && ~options.noSyncPulse
                diffIdx_imec = firstIdx(system) - Imec_inBaselineIdx; 
                firstPulse_inImec = imec_idx(firstPulse_inImecIdx + diffIdx_imec);
                diffIdx_lfp = firstIdx(system) - LFP_inBaselineIdx; 
                firstPulse_inLFP = lfp_idx(firstPulse_inLFPIdx + diffIdx_lfp);
            end
        end
    elseif strcmpi(session.baselineSystem,'LJ')
        [firstPulse_inLJ,system] = max(firstSamp);
        if system ~= 1 % if some system comes online after baseline
            if withCamera
                diffIdx = firstIdx(system) - Cam_inBaselineIdx; 
                firstPulse_inCam = cam_idx(firstPulse_inCamIdx + diffIdx);
            end
        end
    end
end

% 5.1. Save all sync related data
if withCamera; params.sync.commonStartCamera = firstPulse_inCam; end
if withPhotometry; params.sync.commonStartPhotometry = firstPulse_inLJ; end
if withRecording
    if options.noSyncPulse
        firstPulse_inNI = 0;
        firstPulse_inImec = firstPulse_inNI * (params.sync.apFs/params.sync.behaviorFs); 
        firstPulse_inLFP = firstPulse_inNI * (params.sync.lfpFs/params.sync.behaviorFs);
    end
    params.sync.commonStartAP = firstPulse_inImec; 
    params.sync.commonStartLFP = firstPulse_inLFP;
end
if withNI; params.sync.commonStartNI = firstPulse_inNI; end

disp("Finished: find common start using xcorr");

%% (Ver5) Assign common time stamp
% 6.1. Initialize time stamp array
if withNI; timeNI = zeros(1,length(syncNI)); end
if withCamera; timeCamera = zeros(1,length(syncCam)); end
if withPhotometry; timePhotometry = zeros(1,length(sync_labjack)); end
if withRecording
    if options.noSyncPulse
        timeImec = zeros(1,session.totalImecSamp); 
        timeLFP = zeros(1,session.totalLFPSamp); 
    else
        timeImec = zeros(1,length(syncImec)); 
        timeLFP = zeros(1,length(syncLFP)); 
    end
end

% 6.2. Filling in time
% Fill in the Imec and NI time
if withRecording
    if options.noSyncPulse
        for i=1:session.totalImecSamp
            timeImec(i) = (i-firstPulse_inImec)/params.sync.apFs; % Time of each timebin of Imec in sec
        end
        for i=1:session.totalLFPSamp
            timeLFP(i) = (i-firstPulse_inLFP)/params.sync.lfpFs; % Time of each timebin of Imec in sec
        end
        t0 = 0;
        disp("Finished: Imec&LFP timestamp assigned");
    else
        for i=1:length(syncImec)
            timeImec(i) = (i-firstPulse_inImec)/params.sync.apFs; % Time of each timebin of Imec in sec
        end
        idx_Imec = find(syncImec_diff(firstPulse_inImec:end)==1)+firstPulse_inImec-1;
        t0 = timeImec(firstPulse_inImec);
        disp("Finished: Imec/NI timestamp assigned");
    end
elseif withNI
    for i=1:length(syncNI)
        timeNI(i) = (i-firstPulse_inNI)/params.sync.behaviorFs; % Time of each timebin of NI in sec
    end
    idx_NI = find(syncNI_diff(firstPulse_inNI:end)==1)+firstPulse_inNI-1;
    t0 = timeNI(firstPulse_inNI);
    disp("Finished: Imec/NI timestamp assigned");
end

% Fill in photometry time
if withPhotometry
    if strcmpi(session.baselineSystem,'LJ')
        for i=1:length(sync_labjack)
            timePhotometry(i) = (i-firstPulse_inLJ)/params.sync.labjackFs; % Time of each timebin of labjack in sec
        end
        idx_LJ = find(syncLJ_diff(firstPulse_inLJ:end)==1)+firstPulse_inLJ-1;
        t0 = timePhotometry(firstPulse_inLJ);
        disp("Finished: Labjack timestamp assigned");
    elseif strcmpi(session.baselineSystem,'NI') 
        if withRecording && ~options.noSyncPulse
            % Index of all sync pulses starting from the first common pulse
            idx_LJ = find(syncLJ_diff(firstPulse_inLJ:end)==1)+firstPulse_inLJ-1;
            [timePhotometry,l_LJImec] = assignTimeStamp(timePhotometry,timeImec,idx_LJ,idx_Imec,params.sync.labjackFs);
            timePhotometry = timePhotometry - t0; % align timePhotometry to timeNI
        else
            % Index of all sync pulses starting from the first common pulse
            idx_LJ = find(syncLJ_diff(firstPulse_inLJ:end)==1)+firstPulse_inLJ-1;
            [timePhotometry,l_LJNI] = assignTimeStamp(timePhotometry,timeNI,idx_LJ,idx_NI,params.sync.labjackFs);
            timePhotometry = timePhotometry - t0; % align timePhotometry to timeNI
        end
    end
    disp("Finished: labjack timestamp assigned");
end

% Fill in camera
if withCamera
    if withRecording && ~options.noSyncPulse
        % Index of all sync pulses starting from the first common pulse
        idx_Cam = find(syncCam_diff(firstPulse_inCam:end)==1)+firstPulse_inCam-1;
        [timeCamera,l_CamImec] = assignTimeStamp(timeCamera,timeImec,idx_Cam,idx_Imec,params.sync.camFs);
        timeCamera = timeCamera - t0; % align timeCamera to timeImec
    else
        if strcmpi(session.baselineSystem,'NI')
            % Index of all sync pulses starting from the first common pulse
            idx_Cam = find(syncCam_diff(firstPulse_inCam:end)==1)+firstPulse_inCam-1;
            [timeCamera,l_CamNI] = assignTimeStamp(timeCamera,timeNI,idx_Cam,idx_NI,params.sync.camFs);
            timeCamera = timeCamera - t0; % align timeCamera to timeNI
        elseif strcmpi(session.baselineSystem,'LJ')
            % Index of all sync pulses starting from the first common pulse
            idx_Cam = find(syncCam_diff(firstPulse_inCam:end)==1)+firstPulse_inCam-1;
            [timeCamera,l_CamLJ] = assignTimeStamp(timeCamera,timePhotometry,idx_Cam,idx_LJ,params.sync.camFs);
            timeCamera = timeCamera - t0; % align timeCamera to timeNI
        end
    end
    disp("Finished: camera timestamp assigned");
end

% Fill in NI if withRecording == true
if withRecording && ~options.noSyncPulse
    % Index of all sync pulses starting from the first common pulse
    % idx_Imec = find(syncImec_diff(firstPulse_inImec:end)==1)+firstPulse_inImec-1;
    % [timeImec,l_ImecNI] = assignTimeStamp(timeImec,timeNI,idx_Imec,idx_NI,params.sync.apFs);
    idx_NI = find(syncNI_diff(firstPulse_inNI:end)==1)+firstPulse_inNI-1;
    [timeNI,l_NIImec] = assignTimeStamp(timeNI,timeImec,idx_NI,idx_Imec,params.sync.behaviorFs);
    timeNI = timeNI - t0; % align timePhotometry to timeNI

    idx_LFP = find(syncLFP_diff(firstPulse_inLFP:end)==1)+firstPulse_inLFP-1;
    [timeLFP,l_LFPImec] = assignTimeStamp(timeLFP,timeImec,idx_LFP,idx_Imec,params.sync.lfpFs);
    timeLFP = timeLFP - t0; % align timePhotometry to timeNI
    disp("Finished: LFP timestamp assigned");
end

% Align to common start pulse
if withRecording && ~options.noSyncPulse; timeImec = timeImec - t0;
elseif strcmpi(session.baselineSystem,'NI'); timeNI = timeNI - t0; 
elseif strcmpi(session.baselineSystem,'LJ'); timePhotometry = timePhotometry-t0;
end

% 7. Plot summary plot
if withRecording && ~options.noSyncPulse
    title('LFP vs Imec');
    nexttile; plot(timeImec); title('Imec time in sec');
    nexttile; plot(timeLFP); title('LFP time in sec');
    nexttile; plot(diff(timeImec)); title('d(timeImec)');
    nexttile; plot(diff(timeLFP)); title('d(timeLFP)');
   
    nexttile; title('NI vs Imec');
    nexttile; plot(timeNI); title('NI time in sec');
    nexttile; plot(diff(timeNI)); title('d(timeNI)');
    nexttile; plot(timeImec(idx_Imec(1:l_NIImec))-timeNI(idx_NI(1:l_NIImec))); title('Imec-NI');
    nexttile; plot(timeLFP(idx_LFP(1:l_LFPImec))-timeLFP(idx_LFP(1:l_LFPImec))); title('LFP-NI');
    disp("Finished: Imec/NI timestamp plotted");
elseif withNI
    nexttile; plot(timeNI); title('NI time in sec');
    nexttile; plot(diff(timeNI)); title('d(timeNI)');
    disp("Finished: Imec/NI timestamp plotted");  
end

if withCamera
    title('Camera vs Imec/NI');
    nexttile; plot(timeCamera); title('Camera time in sec'); 
    nexttile; plot(diff(timeCamera)); title('d(timeCamera)');
    if withRecording && ~options.noSyncPulse
        nexttile; plot(timeImec(idx_Imec(1:l_CamImec))-timeCamera(idx_Cam(1:l_CamImec))); title('Imec-cam');
    elseif strcmpi(session.baselineSystem,'NI')
        nexttile; plot(timeNI(idx_NI(1:l_CamNI))-timeCamera(idx_Cam(1:l_CamNI))); title('NI-cam');
    elseif strcmpi(session.baselineSystem,'LJ')
        nexttile; plot(timePhotometry(idx_LJ(1:l_CamLJ))-timeCamera(idx_Cam(1:l_CamLJ))); title('LJ-cam');
    end
    disp("Finished: camera timestamp plotted");
end

if withPhotometry
    if strcmpi(session.baselineSystem,'NI')
        title('Labjack vs Imec/NI');
        nexttile; plot(timePhotometry); title('Photometry time in sec'); 
        nexttile; plot(diff(timePhotometry)); title('d(timePhotometry)');
        if withRecording && ~options.noSyncPulse
            nexttile; plot(timeImec(idx_Imec(1:l_LJImec))-timePhotometry(idx_LJ(1:l_LJImec))); title('Imec-photometry');
        else
            nexttile; plot(timeNI(idx_NI(1:l_LJNI))-timePhotometry(idx_LJ(1:l_LJNI))); title('NI-photometry');
        end
        disp("Finished: Labjack timestamp plotted");
    elseif strcmpi(session.baselineSystem,'LJ')
        nexttile; plot(timePhotometry); title('LJ time in sec');
        nexttile; plot(diff(timePhotometry)); title('d(timePhotometry)');
        disp("Finished: LJ timestamp plotted");
    end
end
saveas(gcf,strcat(options.outputPath,filesep,'Summary_syncing_timestamp.png'));

%% (Ver5) Save sync and other data

disp("Ongoing: saving common time stamp");

if withNI; params.sync.timeNI = timeNI; end

if withRecording
    params.sync.timeImec = timeImec;
    params.sync.timeLFP = timeLFP;
    params.sync.behaviorOffset.Imec = timeImec(1) - timeNI(1);
    params.sync.behaviorOffset.LFP = timeImec(1) - timeNI(1);
end

if withCamera
    params.sync.timeCamera = timeCamera;
    if strcmpi(session.baselineSystem,'NI')
        params.sync.behaviorOffset.Camera = timeCamera(1) - timeNI(1);
    elseif strcmpi(session.baselineSystem,'LJ')
        params.sync.behaviorOffset.Camera = timeCamera(1) - timePhotometry(1);
    end
    % Add time_offset to timeSeries data
    for i = 1:size(timeSeries,2)
        if strcmp(timeSeries(i).system,'Cam')
            timeSeries(i).time_offset = params.sync.behaviorOffset.Camera;
        end
    end
end

if withPhotometry
    params.sync.timePhotometry = timePhotometry;
    if withNI
        params.sync.behaviorOffset.Photometry = timePhotometry(1) - timeNI(1);
    else
        params.sync.behaviorFs = params.sync.labjackFs;
        params.sync.behaviorOffset.Photometry = 0;
    end
    % Add time_offset to timeSeries data
    for i = 1:size(timeSeries,2)
        if strcmp(timeSeries(i).system,'LJ')
            timeSeries(i).time_offset = params.sync.behaviorOffset.Photometry;
        end
    end
end

% Save recording systems involved in this session
session.withNI = withNI;
session.withRecording = withRecording;
session.withPhotometry = withPhotometry;
session.withCamera = withCamera;
session.withPhotometryNI = options.withPhotometryNI;
params.session = session;

% Create toml data
% toml.subject_id = animal;
% toml.fiber.light_source = "Thorlabs LED";
% toml.fiber.fiber_type = "Doric";
% toml.fiber.fiber_diameter = "200 um";
% toml.processing_parameters.behavior_offset = params.sync.behaviorOffset.Photometry;
% toml.processing_parameters.final_z = "true";
% toml.processing_parameters.z_window = options.rollingWindowTime;
% toml.processing_parameters.bandpass_bandwidth = processed.options.bandwidth;
% toml.processing_parameters.sampling_frequency = labjack.samplerate;
% toml.processing_parameters.downsample_frequency = options.downsampleFs;
% toml.processing_parameters.transform = "spectrogram";
% toml.processing_parameters.no_per_segment = processed.options.spectralWindow;
% toml.processing_parameters.no_overlap = processed.options.spectralWindowOverlap;
% toml.processing_parameters.left.carrier_frequency_g = 
% toml.processing_parameters.left.carrier_frequency_r = 
% toml.processing_parameters.right.carrier_frequency_g = 
% toml.processing_parameters.right.carrier_frequency_r = 
% toml.signal_indices.total_channels = labjack.nSignals;


save(syncOutputName,'params','-append');
if withPhotometry; save(timeseriesOutputName,'timeSeries','-append'); end
disp('Finished: struct params, session saved in sync_.mat');

return

end
