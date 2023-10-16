function [] = loadSessions(sessionpath,options)

arguments
   sessionpath string
   options.ni_photometry logical = true
   
   % Reload options (unfinished)
   options.reloadAll logical = false       % If false, skip the whole function if sync_.mat is found
   options.reloadNI logical = false     % If false, skip loading NI and use previous loaded ones
   options.reloadLJ logical = true     % If false, skip loading LJ and use previous loaded ones
   options.reloadCam logical = false    % If false, skip loading Cam and use previous loaded ones
   options.reloadImec logical = false   % If false, skip loading Imec and use previous loaded ones
   
   % Photometry related params
   options.nSampPerDemodBin double = 1 % for labjack demod (originally specturalWindowNew)
   options.rollingWindowTime double = 60 % in seconds
   options.LPFreq double = 0; % low pass freq (0: no LP)
   options.downsampleFs double = 50 % downsample NI photometry
   options.dsMethod string = 'resample' % check downsamplePhotometry.m for explaination
   options.movingAvergeFilter logical = false % moving average filter after downsample
   options.movingAverageWindowSize double = 2 % moving average window size after downsample
   options.plotLJPhotometry logical = true; % Plot photometry signals or not

   % Sync related params
   options.syncPulseWindow double = 200 % #of sync pulse to xcorr
   
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

% Save options params
ni_photometryON = options.ni_photometry;

% Decide reload if session has already been loaded and synced
if ~isempty(dir(fullfile(session.path,"sync_*.mat")))
    if ~any([options.reloadAll,options.reloadNI,options.reloadLJ,options.reloadCam,options.reloadImec])
        disp('Loading stop: sync file found.'); 
        return; 
    end
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


%% Save retrieved data

if options.reloadAll || isempty(dir(fullfile(session.path,"sync_*.mat")))
    save(strcat(session.path,filesep,'sync_',sessionName),...
        'sessionName','nidq','session','syncNI','-v7.3');
    disp("Finished: Session data saved");
end

%% Read behavior/photometry (NI) data
close all;

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
    blink = (temp==1);
    
    % temperature = analogNI(3,:);
    gyro = single(analogNI(4:6,:));
    
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
    redLaserRaw = digitalNI(7,:); % or ENL
    blueLaserRaw = digitalNI(8,:);
    
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
    
    % Process NIDAQ photometry data
    if ni_photometryON == true
        % Downsample
        ds_photometry = downsamplePhotometry(analogNI(8,:),options.downsampleFs,nidq.Fs,...
                                     movingAvergeFilter=options.movingAvergeFilter,...
                                     movingAverageWindowSize=options.movingAverageWindowSize,...
                                     dsMethod=options.dsMethod);
        % Rolling z score
        rollingSize = options.rollingWindowTime; % rolling window in sec
        photometryNI = rollingZ(ds_photometry,rollingSize);
        
        % Store params
        params.sync.ni_photometryFs = options.downsampleFs;
        totalDuration_NI = length(allTones) / params.sync.behaviorFs;
        disp('Finished: NI photometry processing');

        % Plot processed vs raw (30sec)
        if options.plotLJPhotometry
            initializeFig(0.5,0.5); t = tiledlayout('flow');
            if totalDuration_NI < 30; plotDuration = totalDuration_NI; 
            else; plotDuration = 30; end % in sec

            % Plot raw green
            yyaxis left
            plotSamples = params.sync.behaviorFs * plotDuration;
            plot(linspace(1,plotDuration,plotSamples),photometry_raw(1:plotSamples),'-k'); hold on
            % Plot processed green
            yyaxis right
            ax2 = axes(t); ax2.Layout.Tile = 1;
            plotSamples = params.sync.ni_photometryFs * plotDuration;
            plot(linspace(1,plotDuration,plotSamples),photometryNI(1:plotSamples),'-r'); hold on
            legend({'Raw green','Processed green'});
            saveas(gcf,strcat(session.path,'\Photometry_NI_processed',sessionName,'.png'))
        end 
    end
    
    %% (Unused) Integrate to Bernardo's code
    %{
        Ch1: allTrials      Ch2: leftTone;      Ch3: rightTone
        Ch4: leftLick       Ch5: rightLick
        Ch6: leftSolenoid   Ch7: rightSolenoid  Ch8: airpuff
        Ch9: blueLaser      Ch10: redLaser
        Ch11: gyroX         Ch12: gyroY         Ch13: gyroZ
        Ch14: photometry    Ch15: blink
    %}
    % raw = double([allTones; leftTone; rightTone;...
    %     leftLick; rightLick; leftSolenoid; rightSolenoid; airpuff;...
    %     blueLaser; redLaser; gyro(1,:); gyro(2,:); gyro(3,:); ...
    %     photometry_raw; blink]);
    
    % Plotting stuffs
    % figure(1);
    % plot(syncNI_diff); title('sync pulse from Imec/NI');
    % 
    % figure(2);
    % plot(leftLick); title('left lick'); hold on
    % plot(leftSolenoid,'LineWidth',2); hold on
    % plot(leftTone,'LineWidth',2);
    %
    % figure(3);
    % plot(rightLick); title('right lick'); hold on
    % plot(rightSolenoid,'LineWidth',2); title('right reward'); hold on
    % plot(rightTone,'LineWidth',2); title('right tone');
    %
    % figure; timewindow = 1:400*nidq.Fs;
    % plot(airpuff(timewindow)); hold on
    % plot(rightSolenoid(timewindow)); hold on
    
    %% Save behavioral data
    if ni_photometryON
        save(strcat(session.path,filesep,'sync_',sessionName),...
            'airpuff','leftLick','rightLick','leftTone','rightTone',....
            'leftSolenoid','rightSolenoid','allTones','blink','gyro',...
            'photometry_raw','blueLaser','redLaser',...
            'photometryNI','-append');
    else
        save(strcat(session.path,filesep,'sync_',sessionName),...
            'airpuff','leftLick','rightLick','leftTone','rightTone',....
            'leftSolenoid','rightSolenoid','allTones','blink','gyro',...
            'blueLaser','redLaser','-append');
    end
    
    % save(strcat('sync_',sessionName),'omissionTrials','-append');
    disp("Finished: NIDAQ data saved");
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

%     signalRange = 1:spikeDownSample*nDownSamples;
%     params.ephys.signalRange = signalRange;
%     processed.ephys.signals = spikes;
%     syncImec_downsampled_risingEdge = ...
%         (squeeze(sum(...
%         reshape([diff(syncImec(signalRange), 1, 2) zeros(1, 1)]==1, ...
%         1, spikeDownSample, nDownSamples), ...
%         2))>0);
% 
%     syncLFP_downsampled_risingEdge = ...
%         (squeeze(sum(...
%         reshape([diff(syncLFP(signalRange), 1, 2) zeros(1, 1)]==1, ...
%         1, spikeDownSample, nDownSamples), ...
%         2))>0);
    
    
    %% Save ephys data
    save(strcat(session.path,filesep,'sync_',sessionName),'ap','lfp','-append');
    disp('Finished: Neuropixel data loaded');
end

%% Read photometry data
if withPhotometry && (options.reloadAll || options.reloadLJ)
    % rawGreen: raw green signals
    % modGreen: green modulation signal of blue LED
    % demodGreen: demodulation or downsampled detrendGreen
    % rollingGreen: rollingZ(demodGreen)
    % rollingGreenLP: rollingZ(demodGreenLP)

    %% Concatenate raw.mat files

    concatLabjack(session.path,save=true);
    load(strcat(session.path,filesep,'photometryLJ_raw.mat'));
    disp('Finished: concatenate and saved raw photometry data');

    % Store relevant info
    params.sync.labjackFs = samplerate;
    params.photometry.freqMod = freqMod;
    if params.photometry.freqMod
        params.photometry.modFreqGreen = modFreqGreen;
        params.photometry.modFreqRed = modFreqRed;
    end
    totalDuration_LJ = length(sync_labjack) / params.sync.labjackFs;
    params.photometry.totalDuration = totalDuration_LJ;

    %% Downsample (no freq mod)

    if ~params.photometry.freqMod
        % detrendGreen = rollingZ(rawGreen,options.rollingWindowTime);
        % disp('Finished: detrend raw photometry data');

        % Downsample
        demodGreen = downsamplePhotometry(rawGreen,options.downsampleFs,params.sync.labjackFs,...
                         dsMethod='resample');
        if options.downsampleFs ~= length(demodGreen) / totalDuration_LJ
            warning("Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!");
        end
        params.sync.photometryFs = length(demodGreen) / totalDuration_LJ;
        disp("Finished: no modulation, skip demodulation and downsampled");
    
        % Rolling zscore for demodGreen and demodRed
        params.photometry.signalDetrendWindow = floor(options.rollingWindowTime*params.sync.photometryFs);
        rollingGreen = rollingZ(demodGreen,params.photometry.signalDetrendWindow);
    
        % Low pass filter data (skip by default)
        if options.LPFreq > 0
            params.photometry.nSampPerDemodBin = 41;
            params.photometry.filtCut = options.LPFreq;
            lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',params.photometry.filtCut,...
                'PassbandRipple',0.01, 'Samplerate',params.sync.labjackFs/params.photometry.nSampPerDemodBin);
            % Lowpass filter
            demodGreenLP = filtfilt(lpFilt,double(demodGreen));
            rollingGreenLP = rollingZ(demodGreenLP,params.photometry.signalDetrendWindow);
    
            save(strcat(session.path,filesep,'sync_',sessionName),'demodGreenLP','rollingGreenLP','-append');
        end

        % Plot processed vs raw (30sec)
        if options.plotLJPhotometry
            initializeFig(0.5,0.5); t = tiledlayout('flow');
            if totalDuration_LJ < 30; plotDuration = totalDuration_LJ; 
            else; plotDuration = 30; end % in sec

            % Plot raw green
            yyaxis left
            plotSamples = params.sync.labjackFs * plotDuration;
            plot(linspace(1,plotDuration,plotSamples),rawGreen(1:plotSamples),'-k'); hold on
            % Plot processed green
            yyaxis right
            plotSamples = params.sync.photometryFs * plotDuration;
            plot(linspace(1,plotDuration,plotSamples),rollingGreen(1:plotSamples),'-r'); hold on
            legend({'Raw green','Processed green'});
            saveas(gcf,strcat(session.path,'\Photometry_LJ_processed',sessionName,'.png'))
        end

        % Save green channels
        save(strcat(session.path,filesep,'sync_',sessionName), ...
            'demodGreen','rollingGreen','-append');
    end

    %% Demodulation (with freq mod)

    if params.photometry.freqMod
        demodulate_green = demodulatePhotometry(rawGreen,options.downsampleFs,params,...
                            modFreq=params.photometry.modFreqGreen);

        % Store process information
        rollingGreen = demodulate_green.demodData;
        if options.downsampleFs ~= length(demodulate_green.demodData) / totalDuration_LJ
            warning("Desired downsampleFs is different from calculated Fs! Used calculated Fs instead!");
        end
        params.sync.photometryFs = length(demodulate_green.demodData) / totalDuration_LJ;
        params.photometry.demodParams_LJGreen = demodulate_green.options;

        if options.plotLJPhotometry
            initializeFig(0.5,0.5); t = tiledlayout('flow');
            if totalDuration_LJ < 30; plotDuration = totalDuration_LJ; 
            else; plotDuration = 30; end % in sec

            % Plot green
            nexttile;
            % Raw
            yyaxis left
            plotSamples = params.sync.labjackFs * plotDuration;
            plot(ax1,linspace(1,plotDuration,plotSamples),rawGreen(1:plotSamples),'-k'); hold on
            legend('Raw green');
            % Processed
            yyaxis right
            plotSamples = params.sync.photometryFs * plotDuration;
            plot(ax2,linspace(1,plotDuration,plotSamples),rollingGreen(1:plotSamples),'-r'); hold on
            legend({'Raw green','Processed green'});

            saveas(gcf,strcat(session.path,'\Photometry_LJ_processed',sessionName,'.png'));
        end

        % Save green channels
        save(strcat(session.path,filesep,'sync_',sessionName),'rollingGreen','-append');
    end
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
        
            % Rolling z score (60s window)
            rollingSize = 60; % in sec
            rollingmean = movmean(eyeArea,rollingSize*camFs);
            rollingstd = movstd(eyeArea,rollingSize*camFs);
            eyeArea_detrend = (eyeArea - rollingmean)./rollingstd;
            % Rolling z score (60s window)
            rollingmean = movmean(pupilArea,rollingSize*camFs);
            rollingstd = movstd(pupilArea,rollingSize*camFs);
            pupilArea_detrend = (pupilArea - rollingmean)./rollingstd;
            % Rolling z score (60s window)
            rollingmean = movmean(eyePixelIntensity,rollingSize*camFs);
            rollingstd = movstd(eyePixelIntensity,rollingSize*camFs);
            eyePixelIntensity_detrend = (eyePixelIntensity - rollingmean)./rollingstd;
            disp("Finished: detrend eye, pupil area");

            % Save camera table
            save(strcat(session.path,filesep,'sync_',sessionName),'camera',...
                'pupilArea','eyeArea','eyeArea_detrend','pupilArea_detrend',...
                'eyePixelIntensity_detrend','-append');
            disp("Finished: saved camera csv file");
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
            % Rolling z score (60s window)
            rollingSize = 60; % in sec
            rollingmean = movmean(eyeArea,rollingSize*camFs);
            rollingstd = movstd(eyeArea,rollingSize*camFs);
            eyeArea_detrend = (eyeArea - rollingmean)./rollingstd;
            % Rolling z score (60s window)
            rollingmean = movmean(pupilArea,rollingSize*camFs);
            rollingstd = movstd(pupilArea,rollingSize*camFs);
            pupilArea_detrend = (pupilArea - rollingmean)./rollingstd;
            disp("Finished: detrend eye, pupil area");

            % Save camera table
            save(strcat(session.path,filesep,'sync_',sessionName),'camera',...
                'pupilArea','eyeArea','eyeArea_detrend','pupilArea_detrend','-append');
            disp("Finished: saved camera csv file");
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
CamNI_ni = 0; LJNI_ni = 0; ImecNI_ni = 0; LFPNI_ni = 0;
CamNI_niIdx = 0; LJNI_niIdx = 0; ImecNI_niIdx = 0; LFPNI_niIdx = 0;
if withCamera
    [~,maxIdx] = max(r_cam); lags_cam = l_cam(maxIdx); 
    if lags_cam >= 0; CamNI_cam = cam_idx(1); CamNI_ni = ni_idx(1+lags_cam);
    else; CamNI_cam = cam_idx(1-lags_cam); CamNI_ni = ni_idx(1); end
    CamNI_niIdx = max([1+lags_cam,1]); CamNI_camIdx = max([1-lags_cam,1]);
    figure; stem(l_cam,r_cam)
end
if withPhotometry
    [~,maxIdx] = max(r_lj); lags_lj = l_lj(maxIdx); 
    if lags_lj >= 0; LJNI_lj = lj_idx(1); LJNI_ni = ni_idx(1+lags_lj);
    else; LJNI_lj = lj_idx(1-lags_lj); LJNI_ni = ni_idx(1); end
    LJNI_niIdx = max([1+lags_lj,1]); LJNI_ljIdx = max([1-lags_lj,1]);
    figure; stem(l_lj,r_lj)
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
    initializeFig(0.66,0.5); tiledlayout(1,2);
    nexttile; stem(l_imec,r_imec); title('Imec time xcor');
    nexttile; stem(l_lfp,r_lfp); title('LFP time xcor');
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
end

% 6.2.2. Fill in the baseline time
if withRecording
    for i=1:length(syncImec)
        timeImec(i) = (i-ImecNI_imec)/params.sync.apFs; % Time of each timebin of Imec in sec
    end
    idx_Imec = find(syncImec_diff(ImecNI_imec:end)==1)+ImecNI_imec-1;
    t0 = timeImec(ImecNI_imec);
else
    for i=1:length(syncNI)
        timeNI(i) = (i-syncNI_first)/params.sync.behaviorFs; % Time of each timebin of NI in sec
    end
    idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;
    t0 = timeNI(syncNI_first);
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
end

% Align to common start pulse
if withRecording; timeImec = timeImec - t0;
else; timeNI = timeNI - t0; end

% 7. Plot summary plot
if withRecording
    figure; tiledlayout(2,2); title('LFP vs Imec');
    nexttile; plot(timeImec); title('Imec time in sec');
    nexttile; plot(timeLFP); title('LFP time in sec');
    nexttile; plot(diff(timeImec)); title('d(timeImec)');
    nexttile; plot(diff(timeLFP)); title('d(timeLFP)');
    
    
    figure; title('NI vs Imec');
    subplot(1,4,1); plot(timeNI); title('NI time in sec');
    subplot(1,4,2); plot(diff(timeNI)); title('d(timeNI)');
    subplot(1,4,3); plot(timeImec(idx_Imec(1:l_NIImec))-timeNI(idx_NI(1:l_NIImec))); title('Imec-NI');
    subplot(1,4,4); plot(timeLFP(idx_LFP(1:l_LFPImec))-timeLFP(idx_LFP(1:l_LFPImec))); title('LFP-NI');
else
    figure; 
    subplot(1,2,1); plot(timeNI); title('NI time in sec');
    subplot(1,2,2); plot(diff(timeNI)); title('d(timeNI)');
end

if withCamera
    figure; title('Camera vs Imec/NI');
    subplot(1,3,1); plot(timeCamera); title('Camera time in sec'); 
    subplot(1,3,2); plot(diff(timeCamera)); title('d(timeCamera)');
    if withRecording
        subplot(1,3,3); plot(timeImec(idx_Imec(1:l_CamImec))-timeCamera(idx_Cam(1:l_CamImec))); title('Imec-cam');
    else
        subplot(1,3,3); plot(timeNI(idx_NI(1:l_CamNI))-timeCamera(idx_Cam(1:l_CamNI))); title('NI-cam');
    end
end
if withPhotometry
    figure; title('Labjack vs Imec/NI');
    subplot(1,3,1); plot(timePhotometry); title('Photometry time in sec'); 
    subplot(1,3,2); plot(diff(timePhotometry)); title('d(timePhotometry)');
    if withRecording
        subplot(1,3,3); plot(timeImec(idx_Imec(1:l_LJImec))-timePhotometry(idx_LJ(1:l_LJImec))); title('Imec-photometry');
    else
        subplot(1,3,3); plot(timeNI(idx_NI(1:l_LJNI))-timePhotometry(idx_LJ(1:l_LJNI))); title('NI-photometry');
    end
end
autoArrangeFigures();

%% (Ver5) Save sync and other data

if withRecording
    params.sync.timeImec = timeImec;
    params.sync.timeLFP = timeLFP;
    save(strcat(session.path,filesep,'sync_',sessionName),'timeImec','timeLFP','-append');
end

if withCamera
    params.sync.timeCamera = timeCamera;
    save(strcat(session.path,filesep,'sync_',sessionName),'timeCamera','-append');
end

if withPhotometry
    params.sync.timePhotometry = timePhotometry;
    save(strcat(session.path,filesep,'sync_',sessionName),'timePhotometry','-append');
end

% Save recording systems involved in this session
session.withRecording = withRecording;
session.withPhotometry = withPhotometry;
session.withCamera = withCamera;
session.ni_photometryON = ni_photometryON;

params.session = session;
params.sync.timeNI = timeNI;
save(strcat(session.path,filesep,'sync_',sessionName),'params','timeNI','-append');
save(strcat(session.path,filesep,'sync_',sessionName),'params','session','-append');

disp('Finished: struct params, session saved');
return

end
