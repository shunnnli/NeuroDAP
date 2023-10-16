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
   options.movingAvergeFilter logical = true % moving average filter after downsample
   options.movingAverageWindowSize double = 2 % moving average window size after downsample
   options.plotLJPhotometry logical = false; % Plot photometry signals or not

   % Sync related params
   options.syncPulseWindow double = 200 % #of sync pulse to xcorr
   
end

%% Notes
% 2023/09/05
%   1. previous oscillation in LJ is due to LP filter, disabled it by default
%   2. add default downsample for LJ if without modulation (default to 50Hz)
%   3. Can partially analyzed by selecting which system to reload

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
    if ~options.reloadAll; disp('Loading stop: sync file found.'); return; end
else
    options.reloadAll = true;
end

withRecording = ~isempty(dir(fullfile(session.path,'catgt_*\*imec*.ap.bin')));
withCamera = ~isempty(dir(fullfile(session.path,'times_cam1*.csv')));
withPhotometry = isfolder(fullfile(session.path,'Photometry'));
session.nSystems = sum([withRecording,withCamera,withPhotometry,1]);

% Save recording systems involved in this session
session.withRecording = withRecording;
session.withPhotometry = withPhotometry;
session.withCamera = withCamera;
session.ni_photometryON = ni_photometryON;

% % Set up analysis options
% commonDownSample = false; % if true, downsample all data into common freq (will misalign)

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

if ~options.reloadAll
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
    
    % Preprocess NIDAQ photometry data
    if ni_photometryON == true
        photometry_raw = analogNI(8,:);
        [photometryNI,photometry_detrended] = detrendAndDownsample(photometry_raw,...
                                                behaviorFs=nidq.Fs,targetFs=options.downsampleFs,...
                                                rollingWindowTime=options.rollingWindowTime,...
                                                movingAvergeFilter=options.movingAvergeFilter,...
                                                movingAverageWindowSize=options.movingAverageWindowSize,...
                                                dsMethod=options.dsMethod);
        params.sync.ni_photometryFs = options.downsampleFs;
        disp('Finished: NI photometry preprocessing');
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
            'photometryNI','photometry_detrended','-append');
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
    %% Load photometry params
    params.sync.labjackFs = 2052; 
    params.photometry.inclFreqWin = 4; % Number of frequency bins to average (on either side of peak freq) (not used)
    params.photometry.freqStep = 1; % Step sfize in Hz for freqency calculations
    params.photometry.freqStepWidth = 7;
    params.photometry.rollingWindowTime = options.rollingWindowTime; % in seconds
    params.photometry.dropFirstDetrendWindow = 1; % drop the signals before the end of the first detrend time window
    params.photometry.lowPassCorner = 100;
    params.photometry.ptsKeep_before = 40;
    params.photometry.ptsKeep_after = 60;
    
    % params.photometry.freqRange1 = [163.5:5:178.5]; % frequencies for channel 1 spectrogram (target 171 Hz) 
    % params.photometry.freqRange2 = [220.5:5:235.5]; % frequencies for channel 2 spectrogram (target 228 Hz)
    params.photometry.freqRange1 = [166,171,176]; % frequencies for channel 1 spectrogram (target 171 Hz) 
    params.photometry.freq1 = 171;
    params.photometry.freqRange2 = [223,228,233]; % frequencies for channel 2 spectrogram (target 228 Hz)
    params.photometry.freq2 = 228;
    params.photometry.useFreqRange = 0:5:500; %frequencies for all spectrograms in figure
    
    % Load info.mat
    load(strcat(session.pathPhotometry,'info.mat'));
    
    %% Concatenate raw.mat files
    D = dir(strcat(session.pathPhotometry,'Raw_*.mat')); 
    filename = {D.name}; load(strcat(session.pathPhotometry,filename{1}));
    numChannels = length(temp)/params.sync.labjackFs;
    output = zeros(1,(length(D)*length(temp)));
    for i = 1:length(D)
        load(strcat(session.pathPhotometry,filename{i}));
        output(((i-1)*length(temp)+1):(i*length(temp))) = temp;
    end
    totalLen = length(output);
    % timeInSec = (1/params.sync.labjackFs)*(1:round(totalLen/numChannels)); % time in seconds for x axis of figure

    % In the order of scanning
    % Ch0 = find(mod(1:totalLen,numChannels)==0);   % sync pulse
    % Ch1 = find(mod(1:totalLen,numChannels)==1);   % photodetector #1 -green
    % Ch2 = find(mod(1:totalLen,numChannels)==2);   % photodetector #2 -red
    % Ch3 = find(mod(1:totalLen,numChannels)==3);   % copy of 488 modulation
    % Ch4 = find(mod(1:totalLen,numChannels)==4);   % copy of 560 modulation

    % Make green and red arrays and demodulate photodiode signals  
    sync_labjack = output(mod(1:totalLen,numChannels)==0);  % sync pulse
    rawGreen = output(mod(1:totalLen,numChannels)==1);      % photodetector #1 -green
    rawRed = output(mod(1:totalLen,numChannels)==2);        % photodetector #2 -red
    modGreen = output(mod(1:totalLen,numChannels)==3);      % copy of 488 modulation
    modRed = output(mod(1:totalLen,numChannels)==4);        % copy of 560 modulation

    % Plot photometry summary plot (skipped)
    if options.plotLJPhotometry
        figure; 
        subplot(5,1,1); plot(rawGreen); title('green raw'); box off
        subplot(5,1,2); plot(rawRed); title('red raw'); box off
        subplot(5,1,3); plot(modGreen); title('green modulation'); box off
        subplot(5,1,4); plot(modRed); title('red modulation'); box off
        subplot(5,1,5); plot(sync_labjack); title('sync'); box off
        saveas(gcf,strcat(session.path,filesep,'Photometry_raw_',sessionName,'.fig'));
    end
    disp('Finished: reading photometry data');

    %% Detrend (rolling zscore)
    Carriers = lcm(floor(params.sync.labjackFs/params.photometry.freq1),floor(params.sync.labjackFs/params.photometry.freq2));
    rawDetrendWindow = Carriers * floor(options.rollingWindowTime*params.sync.labjackFs/Carriers); % in seconds
                                
    green_mean = movmean(rawGreen, rawDetrendWindow);
    green_std = movstd(rawGreen, rawDetrendWindow);
    detrendGreen = (rawGreen-green_mean)./green_std;
    g_stdZeros = find(green_std==0);
    if ~isempty(g_stdZeros)
        disp('WARNING: Found zeros in standard deviation. Setting Infs to 0');
        detrendGreen(g_stdZeros)=0;
    end

    red_mean = movmean(rawRed, rawDetrendWindow);
    red_std = movstd(rawRed, rawDetrendWindow);
    detrendRed = (rawRed-red_mean)./red_std;
    r_stdZeros = find(red_std==0);
    if ~isempty(r_stdZeros)
        disp('WARNING: Found zeros in standard deviation. Setting Infs to 0');
        detrendRed(r_stdZeros)=0;
    end
    
    if options.plotLJPhotometry
        figure; 
        subplot(2,1,1); plot(detrendGreen,'g'); title("Green detrended"); box off
        subplot(2,1,2); plot(detrendRed,'r'); title("Red detrended"); box off
    end
    
    % % Save all data
    % save(strcat(session.path,filesep,'Photometry_raw.mat'),...
    %     'rawGreen','rawRed','modGreen','modRed',...
    %     'detrendGreen','detrendRed');
    disp('Finished: detrend raw photometry data');

    %% Demodulation or downsample 
    if (~exist('freqMod','var') || freqMod == true)
        
        params.photometry.spectralWindow = 2*9*12; % window size for frequency calculation in spectrogram (number of samples)
        totalDuration_LJ = length(sync_labjack) / params.sync.labjackFs;
        % params.photometry.spectSample = 0.01; % Step size for spectrogram (sec)
        
        params.photometry.nSampPerDemodBin = options.nSampPerDemodBin; % originally spectrualWindowNew
        params.photometry.spectralWindowOverlap = params.photometry.spectralWindow-params.photometry.nSampPerDemodBin; % overlap between windows in spectrogram
        % params.photometry.spectralWindowOverlap = round((length(sync_labjack)-params.photometry.spectralWindow)/(totalDuration_LJ * params.finalFs));

        params.photometry.filtCut = options.LPFreq; % Cut off frequency of 5Hz for low pass filter of processed data
        params.photometry.dsRate = 0.05; % Time steps for down-sampling (seconds) (average of every 100 samples)

        % SAK added 3.11.22 to try detrend before demod
        green = detrendGreen; red = detrendRed;

        % Convert spectrogram window size and overlap from time to samples
        disp(['Spectrum window ', num2str(params.photometry.spectralWindow./params.sync.labjackFs), ' sec; ',...
            num2str(params.photometry.spectralWindow), ' samples at ', ...
            num2str(params.sync.labjackFs), ' Hz with ',...
            num2str(params.photometry.spectralWindowOverlap),' overlap samples']);
        % Create low pass filter for final data
        lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',params.photometry.filtCut,...
            'PassbandRipple',0.01, 'Samplerate',params.sync.labjackFs/params.photometry.nSampPerDemodBin);


        figure; colormap winter; set(gcf,'color','w');
        % Calculate spectrogram channel 1
        [spectVals1,spectFreqs1,filtTimes] = spectrogram(green,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange1,params.sync.labjackFs);
        demodGreen = mean(abs(spectVals1),1); % Avg power of freqRange1 at each time point
        demodGreenLP = filtfilt(lpFilt,double(demodGreen)); % Low pass filter the signals (5Hz)
        subplot(2,3,1); 
        spectrogram(green,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.labjackFs,'yaxis');
        title('green'); box off
        disp("Finished: green modulated with blue freq");

        % Calculate spectrogram channel 2 modulated by green LED
        [spectVals2,spectFreqs2,~] = spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange2,params.sync.labjackFs);
        demodRed = mean(abs(spectVals2),1);
        demodRedLP = filtfilt(lpFilt,double(demodRed)); % Low pass filter the signals
        subplot(2,3,2); 
        spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.labjackFs,'yaxis');
        title('red'); box off
        disp("Finished: red modulated with green freq");

        % Calculate spectrogram channel 2 modulated by blue LED
        [spectVals3,spectFreqs3,~]=spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange1,params.sync.labjackFs);
        rawSig3 = mean(abs(spectVals3),1);
        filtSig3 = filtfilt(lpFilt,double(rawSig3));% Low pass filter the signals
        subplot(2,3,3); 
        spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.labjackFs,'yaxis');
        title('red'); box off
        disp("Finished: red modulated with blue freq");

        % Calculate spectrogram channel 3
        [spectVals4,spectFreqs4,~] = spectrogram(modGreen,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange1,params.sync.labjackFs);
        rawSig4 = mean(abs(spectVals4),1);
        subplot(2,2,3); 
        spectrogram(modGreen,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.labjackFs,'yaxis');
        title('modGreen'); box off
        disp("Finished: modGreen channel");

        % Calculate spectrogram channel 4
        [spectVals5,spectFreqs5,~] = spectrogram(modRed,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange2,params.sync.labjackFs);
        rawSig5 = mean(abs(spectVals5),1);
        subplot(2,2,4); 
        spectrogram(modRed,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.labjackFs,'yaxis');
        title('modRed'); box off
        disp("Finished: modRed channel");
        %saveas(gcf,strcat(session.path,'\Photometry_spectrogram_',sessionName,'.fig'));

        % Store and plot results
        params.photometry.t0_demod = filtTimes(1);
        
        if options.plotLJPhotometry
            figure; set(gcf,'color','w')
            subplot(3,2,1);
            plot(filtTimes,demodGreen,'g'); title('green demod with blue freq'); box off
            subplot(3,2,2);
            plot(filtTimes,demodRed,'r'); title('red demod with green freq'); box off
            subplot(3,2,3);
            plot(filtTimes,rawSig3,'r'); title('red demod with blue freq'); box off
            subplot(3,2,4);
            plot(filtTimes,demodGreenLP,'g'); title('LP filtered: green demod with blue freq'); box off
            subplot(3,2,5);
            plot(filtTimes,demodRedLP,'r'); title('LP filtered: red demod with green freq'); box off
            subplot(3,2,6);
            plot(filtTimes,filtSig3,'r'); title('LP filtered: red demod with blue freq'); box off
            saveas(gcf,strcat(session.path,'\Photometry_demod-LPfiltered_',sessionName,'.fig'))
            autoArrangeFigures;
        end
    
    % Downsample if no modulation
    else        
        totalDuration_LJ = length(sync_labjack) / params.sync.labjackFs;

        % Downsample (default using closestInteger nSampPerBin)
        demodGreen = downsamplePhotometry(detrendGreen,options.downsampleFs,params.sync.labjackFs,...
                                     movingAvergeFilter=options.movingAvergeFilter,...
                                     movingAverageWindowSize=options.movingAverageWindowSize,...
                                     dsMethod='resample');

        disp("Finished: no modulated LED, skip demodulation and downsampled");
    end

    %% Rolling z score agian and save results
    
    % Rolling zscore for demodGreen and demodRed
    params.sync.photometryFs = length(demodGreen) / totalDuration_LJ;
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

        save(strcat(session.path,filesep,'sync_',sessionName),'rollingGreenLP','-append');
    end
    
    if (~exist('freqMod','var') || freqMod == true)
        rollingRed = rollingZ(demodRed,params.photometry.signalDetrendWindow);
        % rollingRedLP = rollingZ(demodRedLP,params.photometry.signalDetrendWindow);
    end

    % greenForSync = detrendGreen; params.photometry.nSampPerDemodBin = 1;
    % greenForSync = demodGreen; params.photometry.nSampPerDemodBin = 41;
    
    % signalRange = 1:length(greenForSync)*params.photometry.nSampPerDemodBin;
    % params.photometry.nSamples_demod = length(demodGreen);
    % params.photometry.finalTimeStep = 1 / params.sync.photometryFs;

    %% Downsample sync channel (not needed)
    % sync_labjack_downsampled = ...
    %     squeeze(sum(...
    %     reshape(sync_labjack(signalRange), ...
    %     1, params.photometry.nSampPerDemodBin, length(greenForSync)), ...
    %     2));

    % sync_labjack_downsampled_risingEdge = ...
    %     (squeeze(sum(...
    %     reshape([diff(sync_labjack(signalRange), 1, 2) zeros(1, 1)]==1, ...
    %     1, params.photometry.nSampPerDemodBin, length(greenForSync)), ...
    %     2))>0);

    %% Save photometry data to sync_.mat

    % rawGreen: raw green signals
    % modGreen: green modulation signal of blue LED
    % detrendGreen: rollingZ(rawGreen)
    % demodGreen: demodulation or downsampled detrendGreen
    % demodGreenLP: LP(demodGreen)
    % rollingGreen: rollingZ(demodGreen)
    % rollingGreenLP: rollingZ(demodGreenLP)

    % Save green channels
    save(strcat(session.path,filesep,'sync_',sessionName),...
        'rawGreen','modGreen',...
        'detrendGreen','demodGreen','demodGreenLP',...
        'rollingGreen','-append');

    % Save red channel if neccessary
    if (~exist('freqMod','var') || freqMod == true)
        save(strcat(session.path,filesep,'sync_',sessionName),...
            'rawRed','modRed',...
            'detrendRed','demodRed','demodRedLP',...
            'rollingRed','-append');
    end
    
end

%% Read camera data
if withCamera && (options.reloadAll || options.reloadCam)
    camerapath = dir(fullfile(session.path,'times_cam1*.csv'));
    csvpath = [camerapath.folder,filesep,camerapath.name];

    % Read camera table
    opts = detectImportOptions(csvpath);
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
    if ~isempty(skipped_frame)
        disp(['Found ', num2str(length(skipped_frame)), ' skipped frames!']);
    end
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
        disp(['Found ', num2str(length(skipped_frame)), ' skipped frames!']);
    end

    % Process eye, pupil area data
    pupilArea = camera{:,5}; 
    eyeArea = camera{:,6};
    % Detrend and downsample
    % [photometryNI,photometry_detrended] = detrendAndDownsample(photometry_raw,...
    %     behaviorFs=nidq.Fs,targetFs=options.downsampleFs,rollingWindowTime=options.rollingWindowTime);
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
    figure; tiledlayout(2,2);
    nexttile; plot(timeImec); title('Imec time in sec');
    nexttile; plot(timeLFP); title('LFP time in sec');
    
    nexttile; plot(diff(timeImec)); title('d(timeImec)');
    nexttile; plot(diff(timeLFP)); title('d(timeLFP)');
    
    
    figure; 
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
    figure; 
    subplot(1,3,1); plot(timeCamera); title('Camera time in sec'); 
    subplot(1,3,2); plot(diff(timeCamera)); title('d(timeCamera)');
    if withRecording
        subplot(1,3,3); plot(timeImec(idx_Imec(1:l_CamImec))-timeCamera(idx_Cam(1:l_CamImec))); title('Imec-cam');
    else
        subplot(1,3,3); plot(timeNI(idx_NI(1:l_CamNI))-timeCamera(idx_Cam(1:l_CamNI))); title('NI-cam');
    end
end
if withPhotometry
    figure; 
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

params.session = session;
params.sync.timeNI = timeNI;
save(strcat(session.path,filesep,'sync_',sessionName),'params','timeNI','-append');

disp('Finished: Sync data saved');

% processed.params = params; processed.session = session;
save(strcat(session.path,filesep,'sync_',sessionName),'params','session','-append');

disp('Finished: struct params, session saved');
return

end

%% Integrate with Bernardo's code
%{
Data structure overview:
    processed.params = params
    processed.session = session
    processed.photometry:
        processed.photometry.signals: cell, each cell contains time series of one channel (rolling->detrend->rolling)
        processed.photometry.signals_demod: cell (rolling->detrend)
        processed.photometry.signals_aligned: cell (rolling->detrend->rolling->start from common sync)
        processed.photometry.signalMoments: moments of each photometry channel (from rollingGreen)
    processed.behavior:
        processed.behavior.downsampled: downsampled behavior channel (each row is one channel)
        processed.behavior.risingEdge: downsampled rising edge
        processed.behavior.fallingEdge: downsampled falling edge
        processed.behavior.occupancy: (not needed)
    processed.trialTable
    processed.timeNI/Imec/Cam/Photometry
%}

% furtherDS = false; % if true, downsample photometry even more
% 
% % 1. Downsample to target frequency
% if withPhotometry && ~furtherDS
%     params.finalFs = params.sync.photometryFs;
%     params.finalTimeStep = params.photometry.finalTimeStep;
% else
%     params.finalFs = 100; %Hz
%     params.finalTimeStep = 1/params.finalFs;
% end
% 
% % 2. Behavior
% % 2.1. Aligned to common start
% nSamples_alingedStartToEnd = session.totalBehaviorSamp - params.sync.commonStartNI+1;
% alignedDuration_NI = nSamples_alingedStartToEnd / params.sync.behaviorFs;
% params.behavior.nSamples_aligned = floor(alignedDuration_NI * params.finalFs); % calculate maximum number of timepoints allowed
% 
% nRawSampPerBin = floor(params.sync.behaviorFs / params.finalFs); %97.48... (always not integer, therefore the time drifts away)
% signalRange = params.sync.commonStartNI : params.sync.commonStartNI+params.behavior.nSamples_aligned*nRawSampPerBin-1;
% params.behavior.alignedDuration_NI = length(signalRange) / params.sync.behaviorFs;
% % processed.behavior.aligned = raw(:,signalRange);
% processed.behavior.signalRange = signalRange;
% 
% % 2.2. Downsample behavior
% % if commonDownSample
% 
%     % Opt 1: downsample directly using resample
%     % [p,q] = rat(params.sync.photometryFs / params.sync.behaviorFs);
%     % n = 10; beta = 5; % n: length of filter window (default 10); beta: smoothing (default 5)
%     % processed.behavior.downSampled = resample(processed.behavior.aligned,p,q,n,beta,'Dimension',2);
%     % % Peak index
%     % processed.behavior.peak = logical(size(processed.behavior.downSampled));
%     % for i = 1:size(processed.behavior.downSampled,1)
%     %     [~,locs] = findpeaks(processed.behavior.downSampled(i,:),'MinPeakDistance',params.sync.photometryFs);
%     %     processed.behavior.peak(i,locs) = 1;
%     % end
% 
% %     % Opt 2: upsample first to integer mulplicates of params.sync.photometryFs
% %     [p,q] = rat(params.finalFs * 100 / params.sync.behaviorFs);
% %     n = 10; beta = 5; % n: length of filter window (default 10); beta: smoothing (default 5)
% %     upsampled = resample(processed.behavior.aligned,p,q,n,beta,'Dimension',2);
% % 
% %     processed.behavior.downSampled = ...
% %         squeeze(sum(...
% %         reshape(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), ...
% %         size(upsampled,1), 100, []), ...
% %         2));
% % 
% %     processed.behavior.risingEdge = ...
% %         squeeze(sum(...
% %         reshape([diff(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), 1, 2) zeros(size(upsampled,1), 1)]==1, ...
% %         size(upsampled,1), 100, []), ...
% %         2))>0;
% % 
% %     processed.behavior.fallingEdge = ...
% %         squeeze(sum(...
% %         reshape([diff(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), 1, 2) zeros(size(upsampled,1), 1)]==-1, ...
% %         size(upsampled,1), 100, []), ...
% %         2))>0;
% %     
% %     processed.behavior.occupance = ...
% %         squeeze(sum(...
% %         reshape(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), ...
% %         size(upsampled,1), 100, []), ...
% %         2)) > 0;
% % end
% 
% % 3. Photometry
% if withPhotometry
%     % 3.1. Aligned photometry signal to common start
%     signalRange = params.sync.commonStartPhotometry : length(demodGreen);
%     params.photometry.t0_aligned = LJNI_lj; % not sure
%     params.photometry.nSamples_aligned = length(signalRange);
%     params.photometry.signalRange = signalRange;
% 
%     % 3.2. Store aligned data
%     processed.photometry.signals{1} = rollingGreen;
%     processed.photometry.signals{2} = rollingRed;
%     processed.photometry.signals_demod{1} = demodGreen;
%     processed.photometry.signals_demod{2} = demodRed;
%     processed.photometry.signals_aligned{1} = rollingGreen;
%     processed.photometry.signals_aligned{2} = rollingRed;
%     processed.photometry.signals_rollingLP{1} = rollingGreenLP;
%     processed.photometry.signals_rollingLP{2} = rollingRedLP;
% 
%     % 3.3. Calculate photometry moments
%     for ccc=1:length(processed.photometry.signals)
%         if ~isempty(processed.photometry.signals{ccc})
%             processed.photometry.signalMoments(ccc, 1) = mean(processed.photometry.signals{ccc});
%             processed.photometry.signalMoments(ccc, 2) = var(processed.photometry.signals{ccc});
%             processed.photometry.signalMoments(ccc, 3) = skewness(processed.photometry.signals{ccc}, 1);
%             processed.photometry.signalMoments(ccc, 4) = kurtosis(processed.photometry.signals{ccc}, 1);
%             processed.photometry.signalMoments(ccc, 5) = skewness(processed.photometry.signals{ccc}, 0); % correct bias
%             processed.photometry.signalMoments(ccc, 6) = kurtosis(processed.photometry.signals{ccc}, 0); % correct bias
%         end
%     end
% 
%     % 3.4. Downsample further if needed
%     if furtherDS
%         % params.photometry.alignedDuration_LJ = length(signalRange)/params.sync.photometryFs;
%         % further downsample
%     end
% end
% 
% % 4. Camera
% if withCamera
%     processed.camera.eyeIntensityRaw = eye_pixel_raw;
%     processed.camera.eyeIntensityDetrend = eye_pixel_detrend;
% 
%     if commonDownSample
%         % Downsample 
%     end
% end
% 
% % 5. Recording
% if withRecording
%     params.ap = ap;
% end
% 
% % 6. Save params & session into processed as well
% % 6.1. Save timestamp data
% processed.timeNI = timeNI;
% if withPhotometry; processed.timePhotometry = timePhotometry; end
% if withCamera; processed.timeCamera = timeCamera; end
% if withRecording; processed.timeImec = timeImec; end
% 
% % 6.2. Save params and session to processed
% processed.params = params; processed.session = session;
% save(strcat(session.path,filesep,'sync_',sessionName),...
%     'processed','params','session','-append');
% 
% disp('Finished: struct processed, params, session saved');
