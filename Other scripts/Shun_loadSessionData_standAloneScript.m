% Shun_loadSessionData_standAloneScript
% Shun Li, 2023/07/25
% Stand alone script for aligning signal from different systems in my rig

%% Load files
clear; close all;
%addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
% addpath(genpath('/Volumes/neurobio/MICROSCOPE/Shun/Neuropixel analysis/Methods'));
% addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
addpath(genpath(fullfile('\\','research.files.med.harvard.edu','neurobio','MICROSCOPE','Shun','Analysis','Methods')));

% Select session via uigetdir
sessionpath = uigetdir('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1}));
clear dirsplit

disp(strcat('**********',sessionName,'**********'));
session.name = sessionName; session.path = sessionpath;

% Decide reload if session has already been loaded and synced
if ~isempty(dir(fullfile(session.path,"sync_*.mat")))
    answer = questdlg("This session has been synced and loaded, do you want to reload again?",...
        "Sync file found","Yes","No",'No');
    switch answer
        case 'Yes'; reload = true;
        case 'No'; reload = false;
    end
    if ~reload; disp('Loading stop: sync file found.'); return; end
end

% Decide whether NI photometry is recorded
answer = questdlg("Is there photometry recorded in NIDAQ?",...
    "Photometry systems","Yes","No",'No');
switch answer
    case 'Yes'; ni_photometryON = true; 
    case 'No'; ni_photometryON = false;
end

withRecording = ~isempty(dir(fullfile(session.path,'catgt_*\*imec*.ap.bin')));
withCamera = 0;%~isempty(dir(fullfile(session.path,'cam*.avi')));
withPhotometry = isfolder(fullfile(session.path,'Photometry'));
session.nSystems = sum([withRecording,withCamera,withPhotometry,1]);
% Save recording systems involved in this session
session.withRecording = withRecording;
session.withPhotometry = withPhotometry;
session.withCamera = withCamera;
session.ni_photometryON = ni_photometryON;

% Set up analysis options
commonDownSample = true; % if true, downsample all data into common freq (will misalign)

%% Set up params

% Load session paths
session.pathNidq = strcat(session.path, '\');
session.nidqBin = strcat(sessionName,'_t0.nidq.bin');
%session.nidqBin = dir(fullfile(session.path,'*.nidq.bin'));
nidq.meta = ReadMeta(session.nidqBin, session.pathNidq);
nidq.Fs = str2double(nidq.meta.niSampRate);
if withPhotometry; session.pathPhotometry = strcat(session.path,'\Photometry','\'); end % dir(fullfile(session.path,'\Photometry\')); end % strcat(session.path,'\Photometry','\')

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

for i = 1:numel(a_ch_list)
    figure(10+i);
    plot(analogNI(a_ch_list(i)+1,:)); % remove +1 if analyzing session before 0512
    title(['NIDAQ analog channel ',num2str(a_ch_list(i))]);
end

% Read Digital Channels in NIDAQ data
d_ch_start = str2double(nidq.meta.niXDChans1(1));
d_ch_end = str2double(nidq.meta.niXDChans1(3));
d_ch_list = d_ch_start:d_ch_end;

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;
digitalNI = ExtractDigital(nidqData, nidq.meta, dw, d_ch_list);

syncNI = digitalNI(2,:); % Save sync channel from NIDAQ data separately

autoArrangeFigures();

%% Save retrieved data
save(strcat(session.path,'/','sync_',sessionName),...
    'sessionName','nidq','session','syncNI','-v7.3');
disp("Finished: Session data saved");

%% Extract behavior from NIDAQ (temperture)
close all;

% Convert other analog signals to rising edge
leftLick = (analogNI(1,:)> (max(analogNI(1,:))/2));
temp = [false, diff(leftLick)];
leftLick = (temp==-1);

rightLick = (analogNI(2,:)> (max(analogNI(2,:))/2));
temp = [false, diff(rightLick)];
rightLick = (temp==-1);

temperature = analogNI(3,:);
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

% Integrate to Bernardo's code
%{
    Ch1: allTrials      Ch2: leftTone;      Ch3: rightTone
    Ch4: leftLick       Ch5: rightLick
    Ch6: leftSolenoid   Ch7: rightSolenoid  Ch8: airpuff
    Ch9: blueLaser      Ch10: redLaser
    Ch11: gyroX         Ch12: gyroY         Ch13: gyroZ
    Ch14: photometry    Ch15: temperature
%}
% raw = double([allTones; leftTone; rightTone;...
%     leftLick; rightLick; leftSolenoid; rightSolenoid; airpuff;...
%     blueLaser; redLaser; gyro(1,:); gyro(2,:); gyro(3,:); ...
%     photometry_raw; temperature]);

% Preprocess NIDAQ photometry data
if ni_photometryON == true
    photometry_raw = analogNI(8,:);
    [photometryNI,ni_photometryFs,photometry_detrended] = preprocessNIPhotometry(photometry_raw,params);
    params.sync.ni_photometryFs = ni_photometryFs;
end

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
    save(strcat(session.path,'\','sync_',sessionName),...
        'airpuff','leftLick','rightLick','leftTone','rightTone',....
        'leftSolenoid','rightSolenoid','allTones','temperature','gyro',...
        'photometry_raw','blueLaser','redLaser',...
        'photometryNI','photometry_detrended','-append');
else
    save(strcat(session.path,'\','sync_',sessionName),...
        'airpuff','leftLick','rightLick','leftTone','rightTone',....
        'leftSolenoid','rightSolenoid','allTones','temperature','gyro',...
        'blueLaser','redLaser','-append');
end

% save(strcat('sync_',sessionName),'omissionTrials','-append');
disp("Finished: NIDAQ data saved");

%% Read ephys data
if withRecording
    %% Load imec data
    session.pathImec = strcat(session.path, '\catgt_', sessionName,'\');
    session.apBin = strcat(sessionName,'_tcat.imec0.ap.bin');
    ap.meta = ReadMeta(session.apBin, session.pathImec);
    params.sync.apFs = str2double(ap.meta.imSampRate);

    apTotalSecs = str2double(ap.meta.fileSizeBytes) / params.sync.apFs / str2double(ap.meta.nSavedChans) / 2;
    ap.totalSampIncluded = floor(apTotalSecs * params.sync.apFs);
    session.totalBehaviorSamp = nidq.totalSampIncluded;
    session.totalImecSamp = ap.totalSampIncluded;

    %% Load LFP data
    session.pathLFP = strcat(session.path,'\',sessionName,'_imec0\');
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
    save(strcat(session.path,'\','sync_',sessionName),'ap','lfp','-append');
    disp('Finished: Neuropixel data loaded');
end

%% Read photometry data
if withPhotometry
    %% Load photometry params
    params.sync.photometryFs = 2052; 
    params.photometry.inclFreqWin = 4; % Number of frequency bins to average (on either side of peak freq) (not used)
    params.photometry.freqStep = 1; % Step sfize in Hz for freqency calculations
    params.photometry.freqStepWidth = 7;
    params.photometry.detrendWindowTime = 60; % in seconds
    params.photometry.dropFirstDetrendWindow = 1; % drop the signals before the end of the first detrend time window
    params.photometry.lowPassCorner = 100;
    params.photometry.ptsKeep_before = 40;
    params.photometry.ptsKeep_after = 60;
    
    % params.finalFs = 20; %params.rawSampleFreq/(12*9);
    % params.finalTimeStep = 1/params.finalFs;
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
    numChannels = length(temp)/params.sync.photometryFs;
    output = zeros(1,(length(D)*length(temp)));
    for i = 1:length(D)
        load(strcat(session.pathPhotometry,filename{i}));
        output(((i-1)*length(temp)+1):(i*length(temp))) = temp;
    end
    totalLen = length(output);
    timeInSec = (1/params.sync.photometryFs)*(1:round(totalLen/numChannels)); % time in seconds for x axis of figure

    % In the order of scanning
    Ch1 = find(mod(1:totalLen,numChannels)==1);   % photodetector #1 -green
    Ch2 = find(mod(1:totalLen,numChannels)==2);   % photodetector #2 -red
    Ch3 = find(mod(1:totalLen,numChannels)==3);   % copy of 488 modulation
    Ch4 = find(mod(1:totalLen,numChannels)==4);   % copy of 560 modulation
    Ch5 = find(mod(1:totalLen,numChannels)==0);   % sync pulse

    % Make green and red arrays and demodulate photodiode signals  
    rawgreen = output(Ch1);
    rawred = output(Ch2);
    modgreen = output(Ch3);
    modred = output(Ch4);
    sync_labjack = output(Ch5);

    % Remove artifacts (caused by stimulation)

    figure; 
    subplot(5,1,1); plot(rawgreen); title('green raw'); box off
    subplot(5,1,2); plot(rawred); title('red raw'); box off
    subplot(5,1,3); plot(modgreen); title('green modulation'); box off
    subplot(5,1,4); plot(modred); title('red modulation'); box off
    subplot(5,1,5); plot(sync_labjack); title('sync'); box off
    saveas(gcf,strcat(session.path,'\Photometry_raw_',sessionName,'.fig'));
    disp('Finished: reading photometry data');

    %% Detrend (60s rolling zscore)
    Carriers = lcm(floor(params.sync.photometryFs/params.photometry.freq1),floor(params.sync.photometryFs/params.photometry.freq2));
    rawDetrendWindow = Carriers * floor(params.photometry.detrendWindowTime*params.sync.photometryFs/Carriers); % in seconds
                                
    green_mean = movmean(rawgreen, rawDetrendWindow);
    green_std = movstd(rawgreen, rawDetrendWindow);
    detrendGreen = (rawgreen-green_mean)./green_std;
    g_stdZeros = find(green_std==0);
    if ~isempty(g_stdZeros)
        disp('WARNING: Found zeros in standard deviation. Setting Infs to 0');
        detrendGreen(g_stdZeros)=0;
    end

    red_mean = movmean(rawred, rawDetrendWindow);
    red_std = movstd(rawred, rawDetrendWindow);
    detrendRed = (rawred-red_mean)./red_std;
    r_stdZeros = find(red_std==0);
    if ~isempty(r_stdZeros)
        disp('WARNING: Found zeros in standard deviation. Setting Infs to 0');
        detrendRed(r_stdZeros)=0;
    end
    
    figure; 
    subplot(2,1,1); plot(detrendGreen,'g'); title("Green detrended"); box off
    subplot(2,1,2); plot(detrendRed,'r'); title("Red detrended"); box off
    
    % Save all data
    save(strcat(session.path,'\Photometry_raw.mat'),...
        'rawgreen','rawred','modgreen','modred',...
        'detrendGreen','detrendRed');
    disp('Finished: detrend raw photometry data');

    %% Demodulation    
    if (~exist('freqMod','var') || freqMod == true)
        % params.photometry.spectSample = 0.01; % Step size for spectrogram (sec)
        params.photometry.spectralWindow = 2*9*12; % window size for frequency calculation in spectrogram (number of samples)
        totalDuration_LJ = length(sync_labjack) / params.sync.photometryFs;
        % params.photometry.spectralWindowOverlap = round((length(sync_labjack)-params.photometry.spectralWindow)/(totalDuration_LJ * params.finalFs));
        params.photometry.spectralWindowNew = 1;
        params.photometry.spectralWindowOverlap = params.photometry.spectralWindow-params.photometry.spectralWindowNew; % overlap between windows in spectrogram

        params.photometry.filtCut = 100/params.photometry.spectralWindowNew; % Cut off frequency of 5Hz for low pass filter of processed data
        params.photometry.dsRate = 0.05; % Time steps for down-sampling (seconds) (average of every 100 samples)

        % SAK added 3.11.22 to try detrend before demod
        green = detrendGreen; red = detrendRed;

        % Convert spectrogram window size and overlap from time to samples
        disp(['Spectrum window ', num2str(params.photometry.spectralWindow./params.sync.photometryFs), ' sec; ',...
            num2str(params.photometry.spectralWindow), ' samples at ', ...
            num2str(params.sync.photometryFs), ' Hz with ',...
            num2str(params.photometry.spectralWindowOverlap),' overlap samples']);
        % Create low pass filter for final data
        lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',params.photometry.filtCut,...
            'PassbandRipple',0.01, 'Samplerate',params.sync.photometryFs/params.photometry.spectralWindowNew);


        figure; colormap winter; set(gcf,'color','w');
        % Calculate spectrogram channel 1
        [spectVals1,spectFreqs1,filtTimes] = spectrogram(green,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange1,params.sync.photometryFs);
        demodGreen = mean(abs(spectVals1),1); % Avg power of freqRange1 at each time point
        demodGreenLP = filtfilt(lpFilt,double(demodGreen)); % Low pass filter the signals (5Hz)
        subplot(2,3,1); 
        spectrogram(green,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.photometryFs,'yaxis');
        title('green'); box off
        disp("Finished: green modulated with blue freq");

        % Calculate spectrogram channel 2 modulated by green LED
%         [spectVals2,spectFreqs2,~] = spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange2,params.sync.photometryFs);
%         demodRed = mean(abs(spectVals2),1);
%         demodRedLP = filtfilt(lpFilt,double(demodRed)); % Low pass filter the signals
%         subplot(2,3,2); 
%         spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.photometryFs,'yaxis');
%         title('red'); box off
%         disp("Finished: red modulated with green freq");

        % Calculate spectrogram channel 2 modulated by blue LED
%         [spectVals3,spectFreqs3,~]=spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange1,params.sync.photometryFs);
%         rawSig3 = mean(abs(spectVals3),1);
%         filtSig3 = filtfilt(lpFilt,double(rawSig3));% Low pass filter the signals
%         subplot(2,3,3); 
%         spectrogram(red,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.photometryFs,'yaxis');
%         title('red'); box off
%         disp("Finished: red modulated with blue freq");

        % Calculate spectrogram channel 3
        [spectVals4,spectFreqs4,~] = spectrogram(modgreen,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange1,params.sync.photometryFs);
        rawSig4 = mean(abs(spectVals4),1);
        subplot(2,2,3); 
        spectrogram(modgreen,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.photometryFs,'yaxis');
        title('modgreen'); box off
        disp("Finished: modgreen channel");

        % Calculate spectrogram channel 4
%         [spectVals5,spectFreqs5,~] = spectrogram(modred,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.freqRange2,params.sync.photometryFs);
%         rawSig5 = mean(abs(spectVals5),1);
%         subplot(2,2,4); 
%         spectrogram(modred,params.photometry.spectralWindow,params.photometry.spectralWindowOverlap,params.photometry.useFreqRange,params.sync.photometryFs,'yaxis');
%         title('modred'); box off
%         disp("Finished: modred channel");
%         %saveas(gcf,strcat(session.path,'\Photometry_spectrogram_',sessionName,'.fig'));

        % Store and plot results
        params.photometry.t0_demod = filtTimes(1);
        
        figure; set(gcf,'color','w')
        subplot(3,2,1);
        plot(filtTimes,demodGreen,'g'); title('green demod with blue freq'); box off
%         subplot(3,2,2);
%         plot(filtTimes,demodRed,'r'); title('red demod with green freq'); box off
%         subplot(3,2,3);
%         plot(filtTimes,rawSig3,'r'); title('red demod with blue freq'); box off
        subplot(3,2,4);
        plot(filtTimes,demodGreenLP,'g'); title('LP filtered: green demod with blue freq'); box off
%         subplot(3,2,5);
%         plot(filtTimes,demodRedLP,'r'); title('LP filtered: red demod with green freq'); box off
%         subplot(3,2,6);
%         plot(filtTimes,filtSig3,'r'); title('LP filtered: red demod with blue freq'); box off
        % saveas(gcf,strcat(session.path,'\Photometry_demod-LPfiltered_',sessionName,'.fig'))
        autoArrangeFigures;
    else        
        totalDuration_LJ = length(sync_labjack) / params.sync.photometryFs;
        % Create low pass filter for final data
        params.photometry.spectralWindowNew = 1;
        params.photometry.filtCut = 5;
        lpFilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',params.photometry.filtCut,...
            'PassbandRipple',0.01, 'Samplerate',params.sync.photometryFs/params.photometry.spectralWindowNew);
        
        green = detrendGreen; %red = detrendRed;
        demodGreen = detrendGreen; %demodRed = detrendRed;
        demodGreenLP = filtfilt(lpFilt,double(detrendGreen));
        %demodRedLP = filtfilt(lpFilt,double(detrendRed));
        disp("Finished: no modulated LED, skip demodulation");
    end

    %% Rolling z score agian and save results
    params.photometry.nSamples_demod = length(demodGreen);
    params.photometry.finalFs = length(demodGreen) / totalDuration_LJ;
    params.sync.photometryFs = params.photometry.finalFs;
    params.photometry.finalTimeStep = 1 / params.photometry.finalFs;
    signalRange = 1: length(demodGreen)*params.photometry.spectralWindowNew;
    
    % Rolling zscore for demodGreen and demodRed
    params.photometry.signalDetrendWindow = floor(params.photometry.detrendWindowTime*params.photometry.finalFs);
    rollingGreen = rollingZ(demodGreen,params.photometry.signalDetrendWindow);
    %rollingRed = rollingZ(demodRed,params.photometry.signalDetrendWindow);
    rollingGreenLP = rollingZ(demodGreenLP,params.photometry.signalDetrendWindow);
    %rollingRedLP = rollingZ(demodRedLP,params.photometry.signalDetrendWindow);

    %% Downsample sync channel
    sync_labjack_downsampled = ...
        squeeze(sum(...
        reshape(sync_labjack(signalRange), ...
        1, params.photometry.spectralWindowNew, length(demodGreen)), ...
        2));

    sync_labjack_downsampled_risingEdge = ...
        (squeeze(sum(...
        reshape([diff(sync_labjack(signalRange), 1, 2) zeros(1, 1)]==1, ...
        1, params.photometry.spectralWindowNew, length(demodGreen)), ...
        2))>0);

    sync_labjack_downsampled_fallingEdge = ...
        (squeeze(sum(...
        reshape([diff(sync_labjack(signalRange), 1, 2) zeros(1, 1)]==-1, ...
        1, params.photometry.spectralWindowNew, length(demodGreen)), ...
        2))>0);

    %% Save photometry data to sync_.mat
%     if (~exist('freqMod','var') || freqMod == true)
%         save(strcat(session.path,'\','sync_',sessionName),...
%         'green','red','modgreen','modred',...
%         'demodGreen','demodRed','demodGreenLP','demodRedLP',...
%         'rollingGreen','rollingRed','rollingGreenLP','rollingRedLP','-append');
%     else
%         save(strcat(session.path,'\','sync_',sessionName),...
%         'rawgreen','rawred','modgreen','modred','green','red',...
%         'rollingGreen','rollingRed','rollingGreenLP','rollingRedLP','-append');
%     end
    
    if (~exist('freqMod','var') || freqMod == true)
        save(strcat(session.path,'\','sync_',sessionName),...
        'green','modgreen',...
        'demodGreen','demodGreenLP',...
        'rollingGreen','rollingGreenLP','-append');
    else
        save(strcat(session.path,'\','sync_',sessionName),...
        'rawgreen','modgreen','modred','green',...
        'rollingGreen','rollingGreenLP','-append');
    end
    
end

%% Read camera data
if withCamera
    camerapath = dir(strcat(session.path,'\cam_analyzed.csv'));
    camera = load([camerapath.folder,'\',camerapath.name]);
    params.sync.camFs = 143.93;
    params.camera.raw = camera;

    % Check skipped frames
    skipped_frame = find(diff(camera(:,2)) > 1);
    if ~isempty(skipped_frame)
        disp(['Found ', num2str(length(skipped_frame)), ' skipped frames!']);
    end
    % What should I do next???????????
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
    syncCam = (camera(:,1) > max(camera(:,1)/2))';
    temp = [false,diff(syncCam)];
    syncCam_diff = (temp==1);
end

if withPhotometry
    syncLJ_diff = sync_labjack_downsampled_risingEdge;
end

% Extract start of sync pulse
% 1. Loop over first n pulse
syncPulseWindow = 200;
ni_idx = find(syncNI_diff>0,syncPulseWindow);
if withCamera; cam_idx = find(syncCam_diff>0,syncPulseWindow); end
if withPhotometry; lj_idx = find(syncLJ_diff>0,syncPulseWindow); end
if withRecording
    imec_idx = find(syncImec_diff>0,syncPulseWindow); 
    lfp_idx = find(syncLFP_diff>0,syncPulseWindow); 
end

% 2. Calculate ISI
ISI_ni = (ni_idx(2:end) - ni_idx(1:end-1)) / nidq.Fs;
if withCamera; ISI_cam = (cam_idx(2:end) - cam_idx(1:end-1)) / params.sync.camFs; end
if withPhotometry; ISI_lj = (lj_idx(2:end) - lj_idx(1:end-1)) / params.photometry.finalFs; end
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

%% (Ver5) Assign common time stamp
% 6.1. Initialize time stamp array
timeNI = zeros(1,length(syncNI));
if withCamera; timeCamera = zeros(1,length(syncCam)); end
if withPhotometry; timePhotometry = zeros(1,length(sync_labjack_downsampled)); end
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
        [timePhotometry,l_LJImec] = assignTimeStamp(timePhotometry,timeImec,idx_LJ,idx_Imec,params.photometry.finalFs);
        timePhotometry = timePhotometry - t0; % align timePhotometry to timeNI
    else
        % Index of all sync pulses starting from the first common pulse
        idx_LJ = find(syncLJ_diff(LJNI_lj:end)==1)+LJNI_lj-1;
        [timePhotometry,l_LJNI] = assignTimeStamp(timePhotometry,timeNI,idx_LJ,idx_NI,params.photometry.finalFs);
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
    save(strcat(session.path,'\','sync_',sessionName),'timeImec','timeLFP','-append');
end

if withCamera
    params.sync.timeCamera = timeCamera;
    eye_pixel_raw = camera(:,5); eye_pixel_detrend = camera(:,6);
    save(strcat(session.path,'\','sync_',sessionName),...
    'eye_pixel_raw','eye_pixel_detrend','timeCamera','-append');
end

if withPhotometry
    params.sync.timePhotometry = timePhotometry;
    save(strcat(session.path,'\','sync_',sessionName),'timePhotometry','-append');
end

params.session = session;
params.sync.timeNI = timeNI;
save(strcat(session.path,'\','sync_',sessionName),'params','timeNI','-append');

disp('Finished: Sync data saved');

processed.params = params; processed.session = session;
save(strcat(session.path,'\','sync_',sessionName),...
    'processed','params','session','-append');

disp('Finished: struct processed, params, session saved');
return

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

furtherDS = false; % if true, downsample photometry even more

% 1. Downsample to target frequency
if withPhotometry && ~furtherDS
    params.finalFs = params.photometry.finalFs;
    params.finalTimeStep = params.photometry.finalTimeStep;
else
    params.finalFs = 100; %Hz
    params.finalTimeStep = 1/params.finalFs;
end

% 2. Behavior
% 2.1. Aligned to common start
nSamples_alingedStartToEnd = session.totalBehaviorSamp - params.sync.commonStartNI+1;
alignedDuration_NI = nSamples_alingedStartToEnd / params.sync.behaviorFs;
params.behavior.nSamples_aligned = floor(alignedDuration_NI * params.finalFs); % calculate maximum number of timepoints allowed

nRawSampPerBin = floor(params.sync.behaviorFs / params.finalFs); %97.48... (always not integer, therefore the time drifts away)
signalRange = params.sync.commonStartNI : params.sync.commonStartNI+params.behavior.nSamples_aligned*nRawSampPerBin-1;
params.behavior.alignedDuration_NI = length(signalRange) / params.sync.behaviorFs;
% processed.behavior.aligned = raw(:,signalRange);
processed.behavior.signalRange = signalRange;

% 2.2. Downsample behavior
% if commonDownSample
    
    % Opt 1: downsample directly using resample
    % [p,q] = rat(params.photometry.finalFs / params.sync.behaviorFs);
    % n = 10; beta = 5; % n: length of filter window (default 10); beta: smoothing (default 5)
    % processed.behavior.downSampled = resample(processed.behavior.aligned,p,q,n,beta,'Dimension',2);
    % % Peak index
    % processed.behavior.peak = logical(size(processed.behavior.downSampled));
    % for i = 1:size(processed.behavior.downSampled,1)
    %     [~,locs] = findpeaks(processed.behavior.downSampled(i,:),'MinPeakDistance',params.photometry.finalFs);
    %     processed.behavior.peak(i,locs) = 1;
    % end

%     % Opt 2: upsample first to integer mulplicates of params.photometry.finalFs
%     [p,q] = rat(params.finalFs * 100 / params.sync.behaviorFs);
%     n = 10; beta = 5; % n: length of filter window (default 10); beta: smoothing (default 5)
%     upsampled = resample(processed.behavior.aligned,p,q,n,beta,'Dimension',2);
% 
%     processed.behavior.downSampled = ...
%         squeeze(sum(...
%         reshape(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), ...
%         size(upsampled,1), 100, []), ...
%         2));
% 
%     processed.behavior.risingEdge = ...
%         squeeze(sum(...
%         reshape([diff(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), 1, 2) zeros(size(upsampled,1), 1)]==1, ...
%         size(upsampled,1), 100, []), ...
%         2))>0;
% 
%     processed.behavior.fallingEdge = ...
%         squeeze(sum(...
%         reshape([diff(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), 1, 2) zeros(size(upsampled,1), 1)]==-1, ...
%         size(upsampled,1), 100, []), ...
%         2))>0;
%     
%     processed.behavior.occupance = ...
%         squeeze(sum(...
%         reshape(upsampled(:,1:size(upsampled,2)-mod(size(upsampled,2),100)), ...
%         size(upsampled,1), 100, []), ...
%         2)) > 0;
% end

% 3. Photometry
if withPhotometry
    % 3.1. Aligned photometry signal to common start
    signalRange = params.sync.commonStartPhotometry : length(demodGreen);
    params.photometry.t0_aligned = LJNI_lj; % not sure
    params.photometry.nSamples_aligned = length(signalRange);
    params.photometry.signalRange = signalRange;
    
    % 3.2. Store aligned data
    processed.photometry.signals{1} = rollingGreen;
    processed.photometry.signals{2} = rollingRed;
    processed.photometry.signals_demod{1} = demodGreen;
    processed.photometry.signals_demod{2} = demodRed;
    processed.photometry.signals_aligned{1} = rollingGreen;
    processed.photometry.signals_aligned{2} = rollingRed;
    processed.photometry.signals_rollingLP{1} = rollingGreenLP;
    processed.photometry.signals_rollingLP{2} = rollingRedLP;
    
    % 3.3. Calculate photometry moments
    for ccc=1:length(processed.photometry.signals)
        if ~isempty(processed.photometry.signals{ccc})
            processed.photometry.signalMoments(ccc, 1) = mean(processed.photometry.signals{ccc});
            processed.photometry.signalMoments(ccc, 2) = var(processed.photometry.signals{ccc});
            processed.photometry.signalMoments(ccc, 3) = skewness(processed.photometry.signals{ccc}, 1);
            processed.photometry.signalMoments(ccc, 4) = kurtosis(processed.photometry.signals{ccc}, 1);
            processed.photometry.signalMoments(ccc, 5) = skewness(processed.photometry.signals{ccc}, 0); % correct bias
            processed.photometry.signalMoments(ccc, 6) = kurtosis(processed.photometry.signals{ccc}, 0); % correct bias
        end
    end
    
    % 3.4. Downsample further if needed
    if furtherDS
        % params.photometry.alignedDuration_LJ = length(signalRange)/params.photometry.finalFs;
        % further downsample
    end
end

% 4. Camera
if withCamera
    processed.camera.eyeIntensityRaw = eye_pixel_raw;
    processed.camera.eyeIntensityDetrend = eye_pixel_detrend;

    if commonDownSample
        % Downsample 
    end
end

% 5. Recording
if withRecording
    params.ap = ap;
end

% 6. Save params & session into processed as well
% 6.1. Save timestamp data
processed.timeNI = timeNI;
if withPhotometry; processed.timePhotometry = timePhotometry; end
if withCamera; processed.timeCamera = timeCamera; end
if withRecording; processed.timeImec = timeImec; end

% 6.2. Save params and session to processed
processed.params = params; processed.session = session;
save(strcat(session.path,'\','sync_',sessionName),...
    'processed','params','session','-append');

disp('Finished: struct processed, params, session saved');
return

%% (Not neccesary) Read All Channels in IMEC data

temp = ReadBin(0, ap.totalSampIncluded, ap.meta, session.apBin, session.pathImec);
L = size(temp,2);

% common average referencing
temp = temp - repmat(median(temp,1),size(temp,1),1);

for i=1:10 %pick channels to plot
    Y = fft(temp(i,:));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    figure(1);
    plot([1:1:L]*1/Fs,temp(i,:));
    title(['raw AP trace ch',num2str(i)]);
    xlabel('s');
    ylabel('raw voltage');
    
    figure(2);
    loglog(f,P1); hold on;
    title('Single-Sided Amplitude Spectrum of X(t)');
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
    
    % autoArrangeFigures();
end

%% (Ver4) Extract rising edge and match sample freq

% Convert sync signals to digital signal (boolean);
syncNI = (syncNI==1);
% Extract location of rising edge
temp = [false,diff(syncNI)]; % diff(syncImec) records time of rising edge (1) and falling edge (0)
syncNI_diff = (temp==1);

if withRecording
    syncImec = (syncImec > (max(syncImec)/2));
    temp = [false,diff(syncImec)];
    syncImec_diff = (temp==1);
end

if withCamera
    syncCam = (camera(:,1) > max(camera(:,1)/2))';
    temp = [false,diff(syncCam)];
    syncCam_diff = (temp==1);
end

if withPhotometry
    syncLJ = (sync_labjack == 1);
    temp = [false,diff(syncLJ)];
    syncLJ_diff = (temp==1);
end

% Extract start of sync pulse
% 1. Loop over first n pulse
ni_idx = find(syncNI_diff>0,20);
if withCamera; cam_idx = find(syncCam_diff>0,20); end
if withPhotometry; lj_idx = find(syncLJ_diff>0,20); end
if withRecording; imec_idx = find(syncImec_diff>0,20); end

% 2. First 60s of sync channel signal
if 60*nidq.Fs > nidq.totalSampIncluded
    niSample = syncNI;
    if withCamera; camSample = syncCam; end
    if withPhotometry; ljSample = syncLJ; end
    if withRecording; imecSample = syncImec; end
else
    niSample = syncNI(1:floor(nidq.Fs*60));
    if withCamera; camSample = syncCam(1:floor(params.sync.camFs*60)); end
    if withPhotometry; ljSample = syncLJ(1:floor(params.photometry.finalFs*60)); end
    if withRecording; imecSample = syncImec(1:floor(ap.Fs*60)); end
end

% Upsample to match sampling rate
if withRecording
    [scalarNI,ni_upsample] = upsampleToMatch(nidq.Fs,ap.Fs,niSample);
    if withCamera; [scalarCam,cam_upsample] = upsampleToMatch(params.sync.camFs,ap.Fs,camSample); end
    if withPhotometry; [scalarLJ,lj_upsample] = upsampleToMatch(params.sync.photometryFs,ap.Fs,ljSample); end
else
    if withCamera; [scalarCam,cam_upsample] = upsampleToMatch(params.sync.camFs,nidq.Fs,camSample); end
    if withPhotometry; [scalarLJ,lj_upsample] = upsampleToMatch(params.sync.photometryFs,nidq.Fs,ljSample); end
end

figure;
plot(niSample); hold on
plot(lj_upsample);
legend;

%% (Ver4) Find common first pulse

if withRecording
    [NIImec_ni,NIImec_imec,NIImec_niIdx,NIImec_imecIdx] = ...
        findCommonFirstPulse(ni_upsample,imecSample,ni_idx,imec_idx,scalarNI,'NI','Imec');
    
    % Initialize first sync pulse
    CamImec_imec = 0; LJImec_imec = 0; CamImec_imecIdx = 0; LJImec_imecIdx = 0;
    
    % Find first sync pulse
    if withCamera
        [CamImec_cam,CamImec_imec,CamImec_camIdx,CamImec_imecIdx] = ...
        findCommonFirstPulse(cam_upsample,imecSample,cam_idx,imec_idx,scalarCam,'Cam','Imec');
    end
    if withPhotometry
        [LJImec_lj,LJImec_imec,LJImec_ljIdx,LJImec_imecIdx] = ...
        findCommonFirstPulse(lj_upsample,imecSample,lj_idx,imec_idx,scalarLJ,'LJ','Imec');
    end

    % Find the first common sync pulse of all system
    firstSamp = [NIImec_imec,CamImec_imec,LJImec_imec];
    firstIdx = [NIImec_imecIdx,CamImec_imecIdx,LJImec_imecIdx];
    [syncImec_first,system] = max(firstSamp);
    if system ~= 1 % if some system comes online later
        if withCamera
            diff = firstIdx(system) - CamImec_imecIdx; 
            CamImec_cam = cam_idx(CamImec_camIdx + diff);
        end
        if withPhotometry
            diff = firstIdx(system) - LJImec_imecIdx; 
            LJImec_lj = lj_idx(LJImec_ljIdx + diff);
        end
    end

else
    % Initialize first sync pulse
    CamNI_ni = 0; LJNI_ni = 0; CamNI_niIdx = 0; LJNI_niIdx = 0;

    if withCamera
        [CamNI_cam,CamNI_ni,CamNI_camIdx,CamNI_niIdx] = findCommonFirstPulse(cam_upsample,niSample, ...
            cam_idx,ni_idx,scalarCam,'Cam','NI');
    end
    if withPhotometry
        [LJNI_lj,LJNI_ni,LJNI_ljIdx,LJNI_niIdx] = findCommonFirstPulse(lj_upsample,niSample, ...
            lj_idx,ni_idx,scalarLJ,'LJ','NI');
    end
    
    % Find the first common sync pulse of all system
    firstSamp = [CamNI_ni,LJNI_ni];
    firstIdx = [CamNI_niIdx,LJNI_niIdx];
    [syncNI_first,system] = max(firstSamp);
    if system ~= 1 % if some system comes online later
        if withCamera
            diff = firstIdx(system) - CamNI_niIdx; 
            CamNI_cam = cam_idx(CamNI_camIdx + diff);
        end
        if withPhotometry
            diff = firstIdx(system) - LJNI_niIdx; 
            LJNI_lj = lj_idx(LJNI_ljIdx + diff);
        end
    end
end

%% (Ver4) Assign each system with a common time stamp

% Initialize time stamp array
timeNI = zeros(1,length(syncNI));
if withCamera; timeCamera = zeros(1,length(syncCam)); end
if withPhotometry; timePhotometry = zeros(1,length(syncLJ)); end
if withRecording; timeImec = zeros(1,length(syncImec)); end

% Filling in the time
if withRecording
    % Fill in the baseline time
    for i=1:length(syncImec)
        timeImec(i)=(i-syncImec_first)*(1/ap.Fs); % Time of each timebin of IMEC in sec
    end
    idx_Imec = find(syncImec_diff(syncImec_first:end)==1)+syncImec_first-1;
    t0 = timeImec(1);

    % Fill in NI
    idx_NI = find(syncNI_diff(NIImec_ni:end)==1)+NIImec_ni-1;
    [timeNI,l_NIImec] = assignTimeStamp(timeNI,timeImec,idx_NI,idx_Imec,nidq.Fs);
    timeNI = timeNI - t0; % align timeNI to timeImec

    % Fill in camera
    if withCamera
        % Index of all sync pulses starting from the first common pulse
        idx_Cam = find(syncCam_diff(CamImec_cam:end)==1)+CamImec_cam-1;
        [timeCamera,l_CamImec] = assignTimeStamp(timeCamera,timeImec,idx_Cam,idx_Imec,params.sync.camFs);
        timeCamera = timeCamera - t0; % align timeCamera to timeImec
    end

    % Fill in photometry
    if withPhotometry
        % Index of all sync pulses starting from the first common pulse
        idx_LJ = find(syncLJ_diff(LJImec_lj:end)==1)+LJImec_lj-1;
        [timePhotometry,l_LJImec] = assignTimeStamp(timePhotometry,timeImec,idx_LJ,idx_Imec,params.sync.photometryFs);
        timePhotometry = timePhotometry - t0; % align timePhotometry to timeImec
    end
    

else
    % Only nidq
    if ~(withRecording || withPhotometry || withCamera)
        syncNI_first = find(syncNI_diff>0,1);
    end

    % Fill in the baseline time
    for i=1:length(syncNI)
        timeNI(i)=(i-syncNI_first)*(1/nidq.Fs); % Time of each timebin of NI in sec
    end
    idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;
    t0 = timeNI(1);

    % Fill in camera
    if withCamera
        % Index of all sync pulses starting from the first common pulse
        idx_Cam = find(syncCam_diff(CamNI_cam:end)==1)+CamNI_cam-1;
        [timeCamera,l_CamNI] = assignTimeStamp(timeCamera,timeNI,idx_Cam,idx_NI,params.sync.camFs);
        timeCamera = timeCamera - t0; % align timeCamera to timeNI
    end

    % Fill in photometry
    if withPhotometry
        % Index of all sync pulses starting from the first common pulse
        idx_LJ = find(syncLJ_diff(LJNI_lj:end)==1)+LJNI_lj-1;
        [timePhotometry,l_LJNI] = assignTimeStamp(timePhotometry,timeNI,idx_LJ,idx_NI,params.sync.photometryFs);
        timePhotometry = timePhotometry - t0; % align timePhotometry to timeNI
    end
end

timeNI = timeNI - t0; % align timeNI to timeImec

figure; 
subplot(1,2,1); plot(timeNI); title('NI time in sec');
subplot(1,2,2); plot(diff(timeNI)); title('d(timeNI)');
if withRecording
    subplot(1,3,1); plot(timeImec); title('Imec time in sec');
    subplot(1,3,2); plot(diff(timeImec)); title('d(timeImec)');
    subplot(1,3,3); plot(timeImec(idx_Imec(1:l_NIImec))-timeNI(idx_NI(1:l_NIImec))); title('Imec-NI');

    if withCamera
        figure; 
        subplot(1,3,1); plot(timeCamera); title('Camera time in sec'); 
        subplot(1,3,2); plot(diff(timeCamera)); title('d(timeCamera)');
        subplot(1,3,3); plot(timeImec(idx_Imec(1:l_CamImec))-timeCamera(idx_Cam(1:l_CamImec))); title('Imec-cam');
    end
    if withPhotometry
        figure; 
        subplot(1,3,1); plot(timePhotometry); title('Photometry time in sec'); 
        subplot(1,3,2); plot(diff(timePhotometry)); title('d(timePhotometry)');
        subplot(1,3,3); plot(timeImec(idx_Imec(1:l_LJImec))-timePhotometry(idx_LJ(1:l_LJImec))); title('Imec-photometry');
    end
else
    if withCamera
        figure; 
        subplot(1,3,1); plot(timeCamera); title('Camera time in sec'); 
        subplot(1,3,2); plot(diff(timeCamera)); title('d(timeCamera)');
        subplot(1,3,3); plot(timeNI(idx_NI(1:l_CamNI))-timeCamera(idx_Cam(1:l_CamNI))); title('NI-cam');
    end
    if withPhotometry
        figure; 
        subplot(1,3,1); plot(timePhotometry); title('Photometry time in sec'); 
        subplot(1,3,2); plot(diff(timePhotometry)); title('d(timePhotometry)');
        subplot(1,3,3); plot(timeNI(idx_NI(1:l_LJNI))-timePhotometry(idx_LJ(1:l_LJNI))); title('NI-photometry');
    end
end
autoArrangeFigures();

%% (Ver4) Save sync and camera & photometry data

sync.nidqFs = nidq.Fs; 
save(strcat(session.path,'\','sync_',sessionName),'sync','timeNI','-append');

if withRecording
    sync.apFs = ap.Fs;
    save(strcat(session.path,'\','sync_',sessionName),'sync','timeImec','-append');
end

if withCamera
    sync.camFs = params.sync.camFs;
    eye_pixel_raw = camera(:,5); eye_pixel_detrend = camera(:,6);
    save(strcat(session.path,'\','sync_',sessionName),'sync',...
    'eye_pixel_raw','eye_pixel_detrend','timeCamera','-append');
end

if withPhotometry
    sync.photometryFs = params.sync.photometryFs;
    save(strcat(session.path,'\','sync_',sessionName),'sync',...
    'rawSig1','rawSig2','filtSig1','filtSig2',...
    'sync','timePhotometry','-append');
end

disp('Finished: Video, photometry, sync data saved');
return

%% (Ver3: without camera) Save sync time
if ~withCamera
    % Extract location of rising edge
    temp = [false,diff(syncNI)];
    syncNI_diff = (temp==1);
    syncNI_first = find(syncNI_diff>0,1);
    
    % Filling in the time for IMEC timebin using the sampling rate
    timeNI = zeros(1,length(syncNI));
    for i=1:length(syncNI)
        % Time of each timebin of IMEC in sec
        timeNI(i)=(i-syncNI_first)*(1/nidq.Fs);
    end
    
    % Align timeNI and timeImec to start at 0
    t0 = timeNI(1);
    timeNI = timeNI - t0; % align timeNI to timeImec
    
    figure; plot(timeNI); title('NI time in sec');
    
    save(strcat(session.path,'\','sync_',sessionName),'timeNI','-append');
    disp("Time data saved");
    return
end

%% (Ver3: With camera) Load data and extract rising edge

% Convert sync signals to digital signal (boolean);
syncNI = (syncNI==1);
syncCam = (camera(:,1) > max(camera(:,1)/2))';
% plot(camera(:,1));

% Extract location of rising edge
temp = [false,diff(syncNI)];
syncNI_diff = (temp==1);
temp = [false,diff(syncCam)];
syncCam_diff = (temp==1);

% **Extract start and match sample freq**

% Loop over first n pulse
y_idx = find(syncNI_diff>0,20);
z_idx = find(syncCam_diff>0,20);

% First 60s of sync channel signal
if 60*nidq.Fs > nidq.totalSampIncluded
    y = syncNI;
    z = syncCam;
else
    y = syncNI(1:floor(nidq.Fs*60));
    z = syncCam(1:floor(params.sync.camFs*60));
end

% Upsample camera to match sampling rate of nidq
c2 = floor(nidq.Fs/params.sync.camFs);
z_upsample = upsample(z,c2);
% Fill upsampled steps with 1 if during sync pulse
for i = 1:length(z_upsample)-c2
    if z_upsample(i) == 1 && z_upsample(i+c2) == 1
        for j = 1:c2-1
            z_upsample(i+j) = 1;
        end
    end
end

plot(y); hold on
plot(z_upsample);
legend;

%% (Ver3: With camera) Find common first pulse of NI and camera

maxdot = 0;
for i = 1:length(y_idx)
    for j = 1:length(z_idx)
        l = min(length(y)-y_idx(i),length(z_upsample)-c2*z_idx(j));
        dotprod = dot(y(y_idx(i):y_idx(i)+l),z_upsample(c2*z_idx(j):c2*z_idx(j)+l));
        disp(['dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot),...
            ', NI pulse: ',num2str(i), ', Camera pulse: ',num2str(j)]);

        if dotprod > maxdot
            y_first = y_idx(i);
            z_first = z_idx(j);
            maxdot = dotprod;

            figure;
            plot(y(y_idx(i):y_idx(i)+l)); hold on
            plot(2*z_upsample(c2*z_idx(j):c2*z_idx(j)+l));
            title(['Imec pulse=',num2str(i), ', Camera pulse=',num2str(j),...
                ', dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot)]);
        end
    end
end

syncNI_first = y_first;
syncCam_first = z_first;
autoArrangeFigures();

%% (Ver3: With camera) Assign each timestep of NI and camera with common time stamp

% Initialize time stamp array
timeNI = zeros(1,length(syncNI));
if withCamera; timeCamera = zeros(1,length(syncCam)); end
if withPhotometry; timePhotometry = zeros(1,length(syncLJ)); end
if withRecording; timeImec = zeros(1,length(syncImec)); end

% Filling in the time for IMEC timebin using the sampling rate
if withRecording
    % Fill in the baseline time
    for i=1:length(syncImec)
        timeImec(i)=(i-syncImec_first)*(1/ap.Fs); % Time of each timebin of IMEC in sec
    end
    idx_Imec = find(syncImec_diff(syncImec_first:end)==1)+syncImec_first-1;
else
    % Fill in the baseline time
    for i=1:length(syncNI)
        timeNI(i)=(i-syncNI_first)*(1/nidq.Fs); % Time of each timebin of NI in sec
    end
    idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;
end

% Filling in the time for camera using the sampling rate and interpolation
% reference = IMEC
timeCamera=zeros(1,length(syncCam));
% Index of all sync pulses starting from the first pulse
idx_Cam = find(syncCam_diff(syncCam_first:end)==1)+syncCam_first-1;


% Fill in time
idx_previous_cam = 0;
l = min([length(idx_NI),length(idx_Cam)]);
for i=1:l
    % Sync the time using the rising edge of the current sync pulse
    timeCamera(idx_Cam(i)) = timeNI(idx_NI(i));
    % Filled in time in between two pulses
    if idx_previous_cam < idx_Cam(i)
        timeCamera(idx_previous_cam+1 : idx_Cam(i)) = ...
            timeCamera(idx_Cam(i)) + (1/params.sync.camFs)*[idx_previous_cam-idx_Cam(i)+1 :1: 0];
    end
    % Filled in time after the last syncNI pulse
    if i == l
        timeCamera(idx_Cam(i):end) = ...
            timeCamera(idx_Cam(i)) + (1/params.sync.camFs)*[0 :1: length(timeCamera)-idx_Cam(i)];
    end

    idx_previous_cam = idx_Cam(i);
end

% Align timeNI and timeImec to start at 0
t0 = timeNI(1);
timeNI = timeNI - t0; % align timeNI to timeImec
timeCamera = timeCamera - t0; % align timeCamera to timeImec

figure;
subplot(3,2,1); plot(timeNI); title('NI time in sec');
subplot(3,2,2); plot(timeCamera); title('Camera time in sec');
subplot(3,2,3); plot(diff(timeNI)); title('d(timeNI)');
subplot(3,2,4); plot(diff(timeCamera)); title('d(timeCamera)');
subplot(3,2,5); plot(timeNI(idx_NI(1:l))-timeCamera(idx_Cam(1:l))); title('NI-cam');
autoArrangeFigures();

%% (Ver3: With camera) Save sync time and analyzed behavior

sync.nidqFs = nidq.Fs;
sync.camFs = params.sync.camFs;
eye_pixel = camera(:,6);

save(strcat(session.path,'\','sync_',sessionName),...
    'eye_pixel','sync','timeNI','timeCamera','-append');

disp('Video and sync data saved');