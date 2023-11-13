function [] = loadSessions(sessionpath,options)

arguments
   sessionpath string
   
   % Reload options (unfinished)
   options.reloadAll logical = false       % If false, skip the whole function if sync_.mat is found
   options.reloadNI logical = false     % If false, skip loading NI and use previous loaded ones
   options.reloadLJ logical = false     % If false, skip loading LJ and use previous loaded ones
   options.reloadCam logical = false    % If false, skip loading Cam and use previous loaded ones
   options.reloadImec logical = false   % If false, skip loading Imec and use previous loaded ones

   % Stimulation digital line params
   options.invertStim logical = true % Invert signal from arduino for stimulation
   
   % Photometry related params
   options.recordLJ double = [1,1,0]
   options.nSampPerDemodBin double = 1 % for labjack demod (originally specturalWindowNew)
   options.rollingWindowTime double = 60 % in seconds
   options.LPFreq double = 0; % low pass freq (0: no LP)
   options.downsampleFs double = 50 % downsample NI photometry
   options.dsMethod string = 'resample' % check downsamplePhotometry.m for explaination
   options.movingAvergeFilter logical = false % moving average filter after downsample
   options.movingAverageWindowSize double = 2 % moving average window size after downsample
   options.removeTwoEnds logical = false % see demodulatePhotometry.m
   options.plotPhotometry logical = true % Plot photometry signals or not

   options.withPhotometryNI logical = false
   options.photometryNI_mod logical = false
   options.photometryNI_modFreq double = 0

   % Sync related params
   options.syncPulseWindow double = 5000 % #of sync pulse to xcorr
   
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
clear dirsplit

disp(strcat('**********',sessionName,'**********'));
session.name = sessionName; session.path = sessionpath;

% Decide reload if session has already been loaded and synced
if ~isempty(dir(fullfile(session.path,"sync_*.mat")))
    if ~any([options.reloadAll,options.reloadNI,options.reloadLJ,options.reloadCam,options.reloadImec])
        disp('Loading stop: sync file found.'); 
        return; 
    end
else
    options.reloadAll = true;
    % Initialize all .mat files
    save(strcat(session.path,filesep,'sync_',sessionName),'sessionName','session','-v7.3');
    save(strcat(session.path,filesep,'timeseries_',sessionName),'sessionName','session','-v7.3');
    save(strcat(session.path,filesep,'data_',sessionName),'sessionName','session','-v7.3');
    save(strcat(session.path,filesep,'behavior_',sessionName),'sessionName','session','-v7.3');
end

withRecording = ~isempty(dir(fullfile(session.path,'catgt_*\*imec*.ap.bin')));
withCamera = ~isempty(dir(fullfile(session.path,'times_cam1*.csv')));
withPhotometry = isfolder(fullfile(session.path,'Photometry'));
session.nSystems = sum([withRecording,withCamera,withPhotometry,1]);

%% Set up params

% Load session paths
session.pathNidq = strcat(session.path,filesep);
session.nidqBin = strcat(sessionName,'_t0.nidq.bin');
%session.nidqBin = dir(fullfile(session.path,'*.nidq.bin'));
nidq.meta = ReadMeta(session.nidqBin, session.pathNidq);
nidq.Fs = str2double(nidq.meta.niSampRate);
if withPhotometry; session.pathPhotometry = strcat(session.path,filesep,'Photometry',filesep); end % dir(fullfile(session.path,'\Photometry\')); end % strcat(session.path,'\Photometry','\')

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

% Read Digital Channels in NIDAQ data
d_ch_start = str2double(nidq.meta.niXDChans1(1));
d_ch_end = str2double(nidq.meta.niXDChans1(3));
d_ch_list = d_ch_start:d_ch_end;

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;
digitalNI = ExtractDigital(nidqData, nidq.meta, dw, d_ch_list);
syncNI = digitalNI(2,:); % Save sync channel from NIDAQ data separately

close all;

%% Read behavior/photometry (NI) data

if options.reloadAll || options.reloadNI

    % Convert other analog signals to digital 
    % Rising edge: temp == 1
    % Falling edge: temp == -1
    leftLick = (analogNI(1,:)> (max(analogNI(1,:))/2));
    temp = [false, diff(leftLick)];
    leftLick = (temp==-1);
    
    rightLick = (analogNI(2,:)> (max(analogNI(2,:))/2));
    temp = [false, diff(rightLick)];
    rightLick = (temp==-1);
    
    blink = (analogNI(3,:)> (max(analogNI(3,:))/2));
    temp = [false, diff(blink)];
    blinkRaw = (temp==1);
    blink = find(blinkRaw);
    
    % temperature = analogNI(3,:);
    gyro = single(analogNI(4:6,:));
    % Photometry PMT
    photometry_raw = analogNI(8,:);
    
    % Digital channels
    temp = false(size(digitalNI));
    for i=1:size(digitalNI,1)
        temp2 = [false, diff(digitalNI(i,:))];
        temp(i,:) = (temp2==1);
    end
    airpuffRaw = digitalNI(1,:);
    leftSolenoidRaw = digitalNI(3,:);
    rightSolenoidRaw = digitalNI(4,:);
    leftToneRaw = digitalNI(5,:);
    rightToneRaw = digitalNI(6,:);
    if options.invertStim
        redLaserRaw = ~digitalNI(7,:); % or ENL
        blueLaserRaw = ~digitalNI(8,:);
    else
        redLaserRaw = digitalNI(7,:); % or ENL
        blueLaserRaw = digitalNI(8,:);
    end
    
    % Calculate length of each pulse
    airpuff = temp(1,:) .* getConsecutive(airpuffRaw)./nidq.Fs;
    leftSolenoid = temp(3,:) .* getConsecutive(leftSolenoidRaw)./nidq.Fs;
    rightSolenoid = temp(4,:) .* getConsecutive(rightSolenoidRaw)./nidq.Fs;
    leftTone = temp(5,:) .* getConsecutive(leftToneRaw)./nidq.Fs;
    rightTone = temp(6,:) .* getConsecutive(rightToneRaw)./nidq.Fs;
    redLaser = temp(7,:) .* getConsecutive(redLaserRaw)./nidq.Fs;
    blueLaser = temp(8,:) .* getConsecutive(blueLaserRaw)./nidq.Fs;
    
    
    % Set up allTrials: not including omissionTrials
    % 1 is left trial, 2 is right trial
    allTones = zeros(size(leftTone));
    allTones(leftTone ~= 0) = 1;
    allTones(rightTone ~= 0) = 2;
    
    %% Save behavioral data
    if options.withPhotometryNI
        save(strcat(session.path,filesep,'data_',sessionName),...
            'airpuff','leftLick','rightLick','leftTone','rightTone',....
            'leftSolenoid','rightSolenoid','allTones','blink','gyro',...
            'photometry_raw','blueLaser','redLaser','-append');
    else
        save(strcat(session.path,filesep,'data_',sessionName),...
            'airpuff','leftLick','rightLick','leftTone','rightTone',....
            'leftSolenoid','rightSolenoid','allTones','blink','gyro',...
            'blueLaser','redLaser','-append');
    end
    disp("Finished: NIDAQ data saved in data_.mat");
end

%% Read ephys (imec) data
if withRecording && (options.reloadAll || options.reloadImec)
    %% Load imec data
    session.pathImec = strcat(session.path,filesep,'catgt_', sessionName,filesep);
    session.apBin = strcat(sessionName,'_tcat.imec0.ap.bin');
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
    
    %% Load kilosort data
%     [clusterLabel,spike_times,spike_clusters,cluster_info] = readNPYData(session.pathImec);
% 
%     % Spike specific data
%     ap.Fs = params.sync.apFs; params.ephys.apFs = params.sync.apFs;
%     ap.goodClusters = clusterLabel.cluster_id(find(clusterLabel.KSLabel=='good')); % Good units
%     % Remove unit with cluster_id == 0
%     ap.goodClusters(ap.goodClusters == 0) = [];
%     ap.nGoodClusters = length(ap.goodClusters);
%     ap.clusterToGoodClusterIndex = zeros(max(ap.goodClusters), 1);  % Index column is cluster_id, first column is goodClusters id -> converst cluster_id to goodClusters id
%     for counter=1:length(ap.goodClusters)
%         ap.clusterToGoodClusterIndex(ap.goodClusters(counter)) = counter;
%     end
% 
%     % Create separate array for cluster_id of good cluster spikes
%     goodClusterSpikeIndices = find(ismember(uint64(spike_clusters), uint64(ap.goodClusters)));
%     % Spike times in samples
%     ap.goodSpikeTimes = spike_times(goodClusterSpikeIndices);
%     % Spike cluster_id in the order of spike occurance
%     ap.goodSpikeClusters = spike_clusters(goodClusterSpikeIndices);
%     ap.cluster_info = cluster_info; 
%     clear goodClusterSpikeIndices cluster_info
% 
%     % Generate spike sparse matrix (downsample to 500Hz)
%     spikeIdx = uint64(ap.clusterToGoodClusterIndex(ap.goodSpikeClusters)); % number of good spikes
%     params.ephys.spikeIdx = spikeIdx;
%     params.ephys.ap = ap;
%     params.ephys.finalFs = 500;
%     spikeDownSample = params.sync.apFs/params.ephys.finalFs; % recording is downsampled by 600x to get a final sample rate of 50 Hz
%     nDownSamples = floor(ap.totalSampIncluded/spikeDownSample);
%     
%     spikeTime_downsampled = floor(ap.goodSpikeTimes/spikeDownSample)+1;
%     % Remove spikes outside of totalSampIncluded
%     outsideSpikes = find(spikeTime_downsampled > nDownSamples);
%     if ~isempty(outsideSpikes); spikeIdx(outsideSpikes) = []; spikeTime_downsampled(outsideSpikes) = []; end
%     nSpikes = length(spikeIdx); val = ones(nSpikes, 1);
%     
%     spikes = sparse(spikeIdx, spikeTime_downsampled, val, ap.nGoodClusters, nDownSamples, nSpikes);
%     params.ephys.nSamples = nDownSamples;
%     params.ephys.nAPSamplePerBin = spikeDownSample;
%     params.ephys.nSpikes = nSpikes;
    

    %% Read in sync channel
    % Normally takes ~10 seconds for 60s data; 230s for 1200s data
    tic
    disp('Ongoing: reading sync channel of Imec...');
    syncImec = ReadBinByCh(0, ap.totalSampIncluded, ap.meta, session.apBin, session.pathImec, 385);
    syncLFP = ReadBinByCh(0, lfp.totalSampIncluded, lfp.meta, session.lfpBin, session.pathLFP, 385);
    disp('time for reading sync channel from Imec data');
    toc
    
    %% Save ephys data
    save(strcat(session.path,filesep,'data_',sessionName),'ap','lfp','-append');
    disp('Finished: Neuropixel data saved in data_.mat');
end

%% Read photometry data
if (withPhotometry || options.withPhotometryNI) && (options.reloadAll || options.reloadLJ)
    %% Concatenate raw.mat files

    if withPhotometry
        disp("Ongoing: concatenate and save raw photometry data");
        if isempty(dir(fullfile(session.path,filesep,'data_labjack.mat')))
            concatLabjack(session.path,save=true,plot=false,record=options.recordLJ);
        end
        load(strcat(session.path,filesep,'data_labjack.mat'));
        
        % Reload concatLabjack if labjack.record does not agree
        if ~isfield(labjack,'record')
            labjack.record = options.recordLJ; 
            labjack.nSignals = sum(labjack.record);
            % Remove non-recorded channels
            labjack.raw(find(~labjack.record),:) = [];
            labjack.modulation(find(~labjack.record),:) = [];
            labjack.name(find(~labjack.record)) = [];
        end
        if sum(labjack.record == options.recordLJ) ~= 3
            disp(['labjack.record: ',labjack.record]);
            disp(['options.recordLJ: ',options.recordLJ]);
            warning("labjack.record does not agree with recordLJ, reload using recordLJ"); 
            concatLabjack(session.path,save=true,plot=false,record=options.recordLJ);
        end
        disp('Finished: concatenate and saved raw photometry data in data_labjack.mat');
        disp(labjack);
    
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
            end
        else
            params.sync.ni_photometryFs = [];
        end

        % Save updated labjack
        save(strcat(sessionpath,filesep,'data_labjack'),'labjack','-append');
    end

    %% Loop through all signals
    i = 1;
    if withPhotometry
        for i = 1:labjack.nSignals
            
            if labjack.mod(i)
                demodulated = demodulatePhotometry(labjack.raw(i,:),...
                                targetFs=options.downsampleFs,...
                                originalFs=params.sync.labjackFs,...
                                modFreq=labjack.modFreq(i),...
                                removeTwoEnds=options.removeTwoEnds,...
                                rollingWindowTime=options.rollingWindowTime);
    
                % Check final Fs
                finalFs = length(demodulated.demodData) / totalDuration_LJ;
                if options.downsampleFs ~= finalFs
                    disp(['finalFs: ',num2str(finalFs)]);
                    disp(['targetFs: ',num2str(options.downsampleFs)]);
                    warning([labjack.name{i},': Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!']);
                end
                params.sync.photometryFs = [finalFs, params.sync.photometryFs];
    
                % Store process information
                timeSeries(i).name = labjack.name{i};
                timeSeries(i).data = demodulated.demodData;
                timeSeries(i).finalFs = finalFs;
                timeSeries(i).system = 'LJ';
                timeSeries(i).time_offset = NaN;
                timeSeries(i).demux = true;
                timeSeries(i).demux_freq = demodulated.options.modFreq;
                timeSeries(i).detrend = true;
                timeSeries(i).detrend_type = 'rolling-z';
                timeSeries(i).detrend_window = demodulated.options.rollingWindowTime;
                timeSeries(i).options = demodulated.options;
        
            else
                disp("Ongoing: Process unmodulated photometry data");
                % Downsample
                downsampled = downsamplePhotometry(labjack.raw(i,:),...
                            targetFs=options.downsampleFs,...
                            originalFs=params.sync.labjackFs,...
                            rollingZ=true,rollingWindowTime=options.rollingWindowTime,...
                            dsMethod='resample');
                
                % Check final Fs
                finalFs = length(downsampled.dsData) / totalDuration_LJ;
                if options.downsampleFs ~= finalFs
                    disp(['finalFs: ',num2str(finalFs)]);
                    disp(['targetFs: ',num2str(options.downsampleFs)]);
                    warning([labjack.name{i},': Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!']);
                end
                params.sync.photometryFs = [finalFs, params.sync.photometryFs];
                
                % Store process information
                timeSeries(i).name = labjack.name{i};
                timeSeries(i).data = downsampled.dsData;
                timeSeries(i).finalFs = finalFs;
                timeSeries(i).system = 'LJ';
                timeSeries(i).time_offset = NaN;
                timeSeries(i).demux = true;
                timeSeries(i).demux_freq = NaN;
                timeSeries(i).detrend = true;
                timeSeries(i).detrend_type = 'rolling-z';
                timeSeries(i).detrend_window = downsampled.options.rollingWindowTime;
                timeSeries(i).options = downsampled.options;
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
            
                    save(strcat(session.path,filesep,'data_',sessionName),'dsGreenLP','rollingGreenLP','-append');
                end
            end
        end
    end

    % Process NIDAQ photometry data
    if options.withPhotometryNI == true
        totalDuration_NI = length(photometry_raw) / params.sync.behaviorFs;
        disp("Ongoing: Process NI photometry data");
        % Process
        if options.photometryNI_mod
            demodulate_NI = demodulatePhotometry(photometry_raw,...
                            targetFs=options.downsampleFs,...
                            originalFs=nidq.Fs,...
                            modFreq=options.photometryNI_modFreq,...
                            removeTwoEnds=options.removeTwoEnds,...
                            rollingWindowTime=options.rollingWindowTime);

            % Check final Fs
            finalFs = length(demodulate_NI.demodData) / totalDuration_NI;
            if options.downsampleFs ~= finalFs
                disp(['finalFs: ',num2str(finalFs)]);
                disp(['targetFs: ',num2str(options.downsampleFs)]);
                warning('PMT: Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!');
            end
            params.sync.ni_photometryFs = [finalFs, params.sync.ni_photometryFs];

            % Store params
            timeSeries(i).name = 'PMT';
            timeSeries(i).data = demodulate_NI.demodData;
            timeSeries(i).finalFs = finalFs;
            timeSeries(i).system = 'NI';
            timeSeries(i).time_offset = 0;
            timeSeries(i).demux = true;
            timeSeries(i).demux_freq = demodulate_NI.options.modFreq;
            timeSeries(i).detrend = true;
            timeSeries(i).detrend_type = 'rolling-z';
            timeSeries(i).detrend_window = demodulate_NI.options.rollingWindowTime;
            timeSeries(i).options = demodulate_NI.options;
            disp('Finished: NI photometry demodulation');

        else
            % Downsample
            downsampled_NI = downsamplePhotometry(photometry_raw,...
                                targetFs=options.downsampleFs,originalFs=nidq.Fs,...
                                rollingWindowTime=options.rollingWindowTime,...
                                dsMethod=options.dsMethod);

            % Check final Fs
            finalFs = length(downsampled_NI.dsData) / totalDuration_NI;
            if options.downsampleFs ~= finalFs
                disp(['finalFs: ',num2str(finalFs)]);
                disp(['targetFs: ',num2str(options.downsampleFs)]);
                warning("PMT: Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!");
            end
            params.sync.ni_photometryFs = [finalFs, params.sync.ni_photometryFs];
            
            % Store params
            timeSeries(i).name = 'PMT';
            timeSeries(i).data = downsampled_NI.dsData;
            timeSeries(i).finalFs = finalFs;
            timeSeries(i).system = 'NI';
            timeSeries(i).time_offset = 0;
            timeSeries(i).demux = false;
            timeSeries(i).demux_freq = NaN;
            timeSeries(i).detrend = true;
            timeSeries(i).detrend_type = 'rolling-z';
            timeSeries(i).detrend_window = downsampled_NI.options.rollingWindowTime;
            timeSeries(i).options = downsampled_NI.options;
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
            plotSamples = rawFs * totalDuration;
            plot(linspace(1,totalDuration,plotSamples),rawTrace,Color=[0.3010 0.7450 0.9330]); hold on
            plot(linspace(1,totalDuration,length(meanTrace)),meanTrace,Color='b',LineWidth=2); hold on
            yyaxis right % Plot processed all session
            plotSamples = floor(timeSeries(i).finalFs * totalDuration);
            plot(linspace(1,totalDuration,plotSamples),timeSeries(i).data,Color='#D95319',LineWidth=2); hold on
            legend({'Raw','Mean','Processed'}); xlabel('Time (s)');
            title(timeSeries(i).name);
        end

        % Save channels
        saveas(gcf,strcat(session.path,filesep,'Summary_photometry_processing.png')); 
        saveas(gcf,strcat(session.path,filesep,'Summary_photometry_processing.fig'));
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
                processed = downsamplePhotometry(traces{i},...
                        targetFs=options.downsampleFs,...
                        originalFs=params.sync.camFs,...
                        rollingZ=true,rollingWindowTime=options.rollingWindowTime,...
                        dsMethod='cloestInteger');

                % Save to timeseries struct
                row = size(timeSeries,2) + 1;
                timeSeries(row).name = trace_names{i};
                timeSeries(row).data = processed.dsData;
                timeSeries(row).finalFs = options.downsampleFs;
                timeSeries(row).system = 'Cam';
                timeSeries(row).time_offset = 0;
                timeSeries(row).demux = false;
                timeSeries(row).demux_freq = NaN;
                timeSeries(row).detrend = true;
                timeSeries(row).detrend_type = 'rolling-z';
                timeSeries(row).detrend_window = processed.options.rollingWindowTime;
                timeSeries(row).options = processed.options;
            end
            disp("Finished: detrend eye, pupil area");

            % Save camera table
            save(strcat(session.path,filesep,'data_',sessionName),'camera',...
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
                processed = downsamplePhotometry(traces{i},...
                        targetFs=options.downsampleFs,...
                        originalFs=params.sync.camFs,...
                        rollingZ=true,rollingWindowTime=options.rollingWindowTime,...
                        dsMethod='cloestInteger');

                % Save to timeseries struct
                row = size(timeSeries,2) + 1;
                timeSeries(row).name = trace_names{i};
                timeSeries(row).data = processed.dsData;
                timeSeries(row).finalFs = options.downsampleFs;
                timeSeries(row).system = 'Cam';
                timeSeries(row).time_offset = 0;
                timeSeries(row).demux = false;
                timeSeries(row).demux_freq = NaN;
                timeSeries(row).detrend = true;
                timeSeries(row).detrend_type = 'rolling-z';
                timeSeries(row).detrend_window = processed.options.rollingWindowTime;
                timeSeries(row).options = processed.options;
            end
            disp("Finished: detrend eye, pupil area");

            % Save camera table
            save(strcat(session.path,filesep,'data_',sessionName),'camera',...
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
syncNI = (syncNI==1);
% Extract location of rising edge
temp = [false,diff(syncNI)]; % diff(syncImec) records time of rising edge (1) and falling edge (0)
syncNI_diff = (temp==1);

if withRecording
    syncImec = (syncImec > (max(syncImec)/2));
    temp = [false,diff(syncImec)];
    syncImec_diff = (temp==1);
    % syncImec_diff = syncImec_downsampled_risingEdge;

    syncLFP = (syncLFP > (max(syncLFP)/2));
    temp = [false,diff(syncLFP)];
    syncLFP_diff = (temp==1);
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
ni_idx = find(syncNI_diff>0,syncPulseWindow);
lj_idx = []; cam_idx = []; imec_idx = []; lfp_idx = [];
if withPhotometry; lj_idx = find(syncLJ_diff>0,syncPulseWindow); end
if withCamera; cam_idx = find(syncCam_diff>0,syncPulseWindow); end
if withRecording
    imec_idx = find(syncImec_diff>0,syncPulseWindow); 
    lfp_idx = find(syncLFP_diff>0,syncPulseWindow); 
end

minSyncPulse = min(nonzeros([length(ni_idx), length(lj_idx), length(cam_idx), length(imec_idx), length(lfp_idx)]));
if syncPulseWindow > min(nonzeros([length(ni_idx), length(lj_idx), length(cam_idx), length(imec_idx), length(lfp_idx)]))
    syncPulseWindow = minSyncPulse;
    ni_idx = find(syncNI_diff>0,syncPulseWindow);
    lj_idx = []; cam_idx = []; imec_idx = []; lfp_idx = [];
    if withPhotometry; lj_idx = find(syncLJ_diff>0,syncPulseWindow); end
    if withCamera; cam_idx = find(syncCam_diff>0,syncPulseWindow); end
    if withRecording
        imec_idx = find(syncImec_diff>0,syncPulseWindow); 
        lfp_idx = find(syncLFP_diff>0,syncPulseWindow); 
    end
end

% 2. Calculate ISI
ISI_ni = (ni_idx(2:end) - ni_idx(1:end-1)) / nidq.Fs;
if withCamera; ISI_cam = (cam_idx(2:end) - cam_idx(1:end-1)) / params.sync.camFs; end
if withPhotometry; ISI_lj = (lj_idx(2:end) - lj_idx(1:end-1)) / params.sync.labjackFs; end
if withRecording
    ISI_imec = (imec_idx(2:end) - imec_idx(1:end-1)) / params.sync.apFs; 
    ISI_lfp = (lfp_idx(2:end) - lfp_idx(1:end-1)) / params.sync.lfpFs; 
end

% 3. Cross correlation
if withCamera; [r_cam, l_cam] = xcorr(zscore(ISI_ni), zscore(ISI_cam),'normalized'); end
if withPhotometry; [r_lj, l_lj] = xcorr(zscore(ISI_ni), zscore(ISI_lj),'normalized'); end
if withRecording
    [r_imec, l_imec] = xcorr(zscore(ISI_ni), zscore(ISI_imec),'normalized'); 
    [r_lfp,l_lfp] = xcorr(zscore(ISI_ni), zscore(ISI_lfp),'normalized');
end

% 4. Find the max of xcorr
initializeFig(0.67,0.67); tiledlayout('flow');
CamNI_ni = 0; LJNI_ni = 0; ImecNI_ni = 0; LFPNI_ni = 0;
CamNI_niIdx = 0; LJNI_niIdx = 0; ImecNI_niIdx = 0; LFPNI_niIdx = 0;
if withCamera
    [~,maxIdx] = max(r_cam); lags_cam = l_cam(maxIdx); 
    if lags_cam >= 0; CamNI_cam = cam_idx(1); CamNI_ni = ni_idx(1+lags_cam);
    else; CamNI_cam = cam_idx(1-lags_cam); CamNI_ni = ni_idx(1); end
    CamNI_niIdx = max([1+lags_cam,1]); CamNI_camIdx = max([1-lags_cam,1]);
    nexttile; stem(l_cam,r_cam); title('NI vs Cam xcor');
end
if withPhotometry
    [~,maxIdx] = max(r_lj); lags_lj = l_lj(maxIdx); 
    if lags_lj >= 0; LJNI_lj = lj_idx(1); LJNI_ni = ni_idx(1+lags_lj);
    else; LJNI_lj = lj_idx(1-lags_lj); LJNI_ni = ni_idx(1); end
    LJNI_niIdx = max([1+lags_lj,1]); LJNI_ljIdx = max([1-lags_lj,1]);
    nexttile; stem(l_lj,r_lj); title('NI vs Labjack xcor');
end
if withRecording
    [~,maxIdx] = max(r_imec); lags_imec = l_imec(maxIdx); 
    [~,maxIdx] = max(r_lfp); lags_lfp = l_lfp(maxIdx);
    if lags_imec >= 0
        ImecNI_imec = imec_idx(1); ImecNI_ni = ni_idx(1+lags_imec);
        LFPNI_lfp = lfp_idx(1); LFPNI_ni = ni_idx(1+lags_lfp);
    else
        ImecNI_imec = imec_idx(1-lags_imec); ImecNI_ni = ni_idx(1);
        LFPNI_lfp = lfp_idx(1-lags_lfp); LFPNI_ni = ni_idx(1);
    end
    ImecNI_niIdx = max([1+lags_imec,1]); ImecNI_imecIdx = max([1-lags_imec,1]);
    LFPNI_niIdx = max([1+lags_lfp,1]); LFPNI_lfpIdx = max([1-lags_lfp,1]);
    nexttile; stem(l_imec,r_imec); title('NI vs Imec xcor');
    nexttile; stem(l_lfp,r_lfp); title('NI vs LFP xcor');
end

% 5. Find the first common sync pulse of all system
if ~(withRecording || withPhotometry || withCamera) % Only nidq
    syncNI_first = find(syncNI_diff>0,1);
else
    firstSamp = [ImecNI_ni,CamNI_ni,LJNI_ni,LFPNI_ni];
    firstIdx = [ImecNI_niIdx,CamNI_niIdx,LJNI_niIdx,LFPNI_niIdx];
    [syncNI_first,system] = max(firstSamp);
    if system ~= 1 % if some system comes online after NI
        if withCamera
            diffIdx = firstIdx(system) - CamNI_niIdx; 
            CamNI_cam = cam_idx(CamNI_camIdx + diffIdx);
        end
        if withPhotometry
            diffIdx = firstIdx(system) - LJNI_niIdx; 
            LJNI_lj = lj_idx(LJNI_ljIdx + diffIdx);
        end
        if withRecording
            diffIdx_imec = firstIdx(system) - ImecNI_niIdx; 
            ImecNI_imec = imec_idx(ImecNI_imecIdx + diffIdx_imec);
            diffIdx_lfp = firstIdx(system) - LFPNI_niIdx; 
            LFPNI_lfp = lfp_idx(LFPNI_lfpIdx + diffIdx_lfp);
        end
    end
end

% 5.1. Save all sync related data
params.sync.commonStartNI = syncNI_first;
if withCamera; params.sync.commonStartCamera = CamNI_cam; end
if withPhotometry; params.sync.commonStartPhotometry = LJNI_lj; end
if withRecording
    params.sync.commonStartAP = ImecNI_imec; 
    params.sync.commonStartLFP = LFPNI_lfp;
end

disp("Finished: find common start using xcorr");

%% (Ver5) Assign common time stamp
% 6.1. Initialize time stamp array
timeNI = zeros(1,length(syncNI));
if withCamera; timeCamera = zeros(1,length(syncCam)); end
if withPhotometry; timePhotometry = zeros(1,length(sync_labjack)); end
if withRecording
    timeImec = zeros(1,length(syncImec)); 
    timeLFP = zeros(1,length(syncLFP)); 
end

% 6.2. Filling in time
% 6.2.1. Only nidq
if ~(withRecording || withPhotometry || withCamera)
    syncNI_first = find(syncNI_diff>0,1);
    disp("Finished: NI timestamp assigned");
end

% 6.2.2. Fill in the baseline time
if withRecording
    for i=1:length(syncImec)
        timeImec(i) = (i-ImecNI_imec)/params.sync.apFs; % Time of each timebin of Imec in sec
    end
    idx_Imec = find(syncImec_diff(ImecNI_imec:end)==1)+ImecNI_imec-1;
    t0 = timeImec(ImecNI_imec);
    disp("Finished: Imec/NI timestamp assigned");
else
    for i=1:length(syncNI)
        timeNI(i) = (i-syncNI_first)/params.sync.behaviorFs; % Time of each timebin of NI in sec
    end
    idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;
    t0 = timeNI(syncNI_first);
    disp("Finished: Imec/NI timestamp assigned");
end


% 6.2.3. Fill in camera
if withCamera
    if withRecording
        % Index of all sync pulses starting from the first common pulse
        idx_Cam = find(syncCam_diff(CamNI_cam:end)==1)+CamNI_cam-1;
        [timeCamera,l_CamImec] = assignTimeStamp(timeCamera,timeImec,idx_Cam,idx_Imec,params.sync.camFs);
        timeCamera = timeCamera - t0; % align timeCamera to timeImec
    else
        % Index of all sync pulses starting from the first common pulse
        idx_Cam = find(syncCam_diff(CamNI_cam:end)==1)+CamNI_cam-1;
        [timeCamera,l_CamNI] = assignTimeStamp(timeCamera,timeNI,idx_Cam,idx_NI,params.sync.camFs);
        timeCamera = timeCamera - t0; % align timeCamera to timeNI
    end
    disp("Finished: camera timestamp assigned");
end

% 6.2.4. Fill in photometry
if withPhotometry
    if withRecording
        % Index of all sync pulses starting from the first common pulse
        idx_LJ = find(syncLJ_diff(LJNI_lj:end)==1)+LJNI_lj-1;
        [timePhotometry,l_LJImec] = assignTimeStamp(timePhotometry,timeImec,idx_LJ,idx_Imec,params.sync.labjackFs);
        timePhotometry = timePhotometry - t0; % align timePhotometry to timeNI
    else
        % Index of all sync pulses starting from the first common pulse
        idx_LJ = find(syncLJ_diff(LJNI_lj:end)==1)+LJNI_lj-1;
        [timePhotometry,l_LJNI] = assignTimeStamp(timePhotometry,timeNI,idx_LJ,idx_NI,params.sync.labjackFs);
        timePhotometry = timePhotometry - t0; % align timePhotometry to timeNI
    end
    disp("Finished: labjack timestamp assigned");
end

% 6.2.5. Fill in NI if withRecording == true
if withRecording
    % Index of all sync pulses starting from the first common pulse
    % idx_Imec = find(syncImec_diff(ImecNI_imec:end)==1)+ImecNI_imec-1;
    % [timeImec,l_ImecNI] = assignTimeStamp(timeImec,timeNI,idx_Imec,idx_NI,params.sync.apFs);
    idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;
    [timeNI,l_NIImec] = assignTimeStamp(timeNI,timeImec,idx_NI,idx_Imec,params.sync.behaviorFs);
    timeNI = timeNI - t0; % align timePhotometry to timeNI

    idx_LFP = find(syncLFP_diff(LFPNI_lfp:end)==1)+LFPNI_lfp-1;
    [timeLFP,l_LFPImec] = assignTimeStamp(timeLFP,timeImec,idx_LFP,idx_Imec,params.sync.lfpFs);
    timeLFP = timeLFP - t0; % align timePhotometry to timeNI
    disp("Finished: LFP timestamp assigned");
end

% Align to common start pulse
if withRecording; timeImec = timeImec - t0;
else; timeNI = timeNI - t0; end

% 7. Plot summary plot
if withRecording
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
else
    nexttile; plot(timeNI); title('NI time in sec');
    nexttile; plot(diff(timeNI)); title('d(timeNI)');
    disp("Finished: Imec/NI timestamp plotted");
end

if withCamera
    title('Camera vs Imec/NI');
    nexttile; plot(timeCamera); title('Camera time in sec'); 
    nexttile; plot(diff(timeCamera)); title('d(timeCamera)');
    if withRecording
        nexttile; plot(timeImec(idx_Imec(1:l_CamImec))-timeCamera(idx_Cam(1:l_CamImec))); title('Imec-cam');
    else
        nexttile; plot(timeNI(idx_NI(1:l_CamNI))-timeCamera(idx_Cam(1:l_CamNI))); title('NI-cam');
    end
    disp("Finished: camera timestamp plotted");
end
if withPhotometry
    title('Labjack vs Imec/NI');
    nexttile; plot(timePhotometry); title('Photometry time in sec'); 
    nexttile; plot(diff(timePhotometry)); title('d(timePhotometry)');
    if withRecording
        nexttile; plot(timeImec(idx_Imec(1:l_LJImec))-timePhotometry(idx_LJ(1:l_LJImec))); title('Imec-photometry');
    else
        nexttile; plot(timeNI(idx_NI(1:l_LJNI))-timePhotometry(idx_LJ(1:l_LJNI))); title('NI-photometry');
    end
    disp("Finished: Labjack timestamp plotted");
end
saveas(gcf,strcat(session.path,filesep,'Summary_syncing_timestamp.png'));

%% (Ver5) Save sync and other data

disp("Ongoing: saving common time stamp");
if withRecording
    params.sync.timeImec = timeImec;
    params.sync.timeLFP = timeLFP;
    params.sync.behaviorOffset.Imec = timeImec(1) - timeNI(1);
    params.sync.behaviorOffset.LFP = timeImec(1) - timeNI(1);
end

if withCamera
    params.sync.timeCamera = timeCamera;
    params.sync.behaviorOffset.Camera = timeCamera(1) - timeNI(1);
end

if withPhotometry
    params.sync.timePhotometry = timePhotometry;
    params.sync.behaviorOffset.Photometry = timePhotometry(1) - timeNI(1);
    % Add time_offset to timeSeries data
    for i = 1:size(timeSeries,2)
        if strcmp(timeSeries(i).system,'LJ')
            timeSeries(i).time_offset = params.sync.behaviorOffset.Photometry;
        end
    end
end

% Save recording systems involved in this session
session.withRecording = withRecording;
session.withPhotometry = withPhotometry;
session.withCamera = withCamera;
session.withPhotometryNI = options.withPhotometryNI;

params.session = session;
params.sync.timeNI = timeNI;
save(strcat(session.path,filesep,'sync_',sessionName),'params','timeNI','-append');
save(strcat(session.path,filesep,'sync_',sessionName),'params','session','-append');
save(strcat(session.path,filesep,'timeseries_',sessionName),'timeSeries','-append');
disp('Finished: struct params, session saved in sync_.mat');

return

end
