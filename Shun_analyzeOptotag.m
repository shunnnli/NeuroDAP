%% Load data
clear; close all;
% addpath(genpath('D:\Shun\Analysis\Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
setenv('NEUROPIXEL_MAP_FILE', which('neuropixPhase3B2_kilosortChanMap.mat'));
% setenv('NEUROPIXEL_DATAROOT', 'D:\Shun\Analysis\Result-Neuropixel utils');
[~,~,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed] = loadColors;

% Select session via uigetdir
sessionpath = uigetdir('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1})); clear dirsplit

disp(strcat('********** ',sessionName,'**********'));
load(strcat(sessionpath,'\','sync_',sessionName,'.mat'));
disp(['Session ',sessionName,' loaded']);

% Load waveforms
gwfparams.fileName = session.apBin; gwfparams.nCh = 385;
% gwfparams.fileName = 'temp_wh.dat'; gwfparams.nCh = 383;
gwfparams.dataType = 'int16';           % Data type of .dat file (this should be BP filtered)
gwfparams.nWf = 385;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes = spike_times;         % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = spike_clusters;   % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

fileName = fullfile(session.pathImec,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % Determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'},'Writable',true);
chMap = readNPY(fullfile(session.pathImec, 'channel_map.npy'))+1; % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);
disp(['Finished: ', gwfparams.fileName,' for session ',sessionName,' loaded']);

% Load kilosort
[clusterLabel,spike_times,spike_clusters,cluster_info] = readNPYData(session.pathImec);
ap.goodClusters = clusterLabel.cluster_id(find(clusterLabel.KSLabel=='good' & clusterLabel.group == 'unsorted')); % Good units
% Remove unit with cluster_id == 0
ap.goodClusters(ap.goodClusters == 0) = [];
ap.nGoodClusters = length(ap.goodClusters);
ap.clusterToGoodClusterIndex = zeros(max(ap.goodClusters), 1);  % Index column is cluster_id, first column is goodClusters id -> converst cluster_id to goodClusters id
for counter=1:length(ap.goodClusters)
    ap.clusterToGoodClusterIndex(ap.goodClusters(counter)) = counter;
end

% Create separate array for cluster_id of good cluster spikes
goodClusterSpikeIndices = find(ismember(uint64(spike_clusters), uint64(ap.goodClusters)));
% Spike times in samples
ap.goodSpikeTimes = spike_times(goodClusterSpikeIndices);
% Spike cluster_id in the order of spike occurance
ap.goodSpikeClusters = spike_clusters(goodClusterSpikeIndices);
ap.cluster_info = cluster_info; 
clear goodClusterSpikeIndices cluster_info

% Generate spike sparse matrix-related
spikeIdx = uint64(ap.clusterToGoodClusterIndex(ap.goodSpikeClusters)); % number of good spikes
params.ephys.spikeIdx = spikeIdx;
params.ephys.ap = ap;

params.ephys.sparseMatrixFs = 500;
spikeDownSample = params.sync.apFs/params.ephys.sparseMatrixFs; % recording is downsampled by 600x to get a final sample rate of 50 Hz
nDownSamples = floor(ap.totalSampIncluded/spikeDownSample);

spikeTime_downsampled = floor(ap.goodSpikeTimes/spikeDownSample)+1;
% Remove spikes outside of totalSampIncluded
outsideSpikes = find(spikeTime_downsampled > nDownSamples);
if ~isempty(outsideSpikes); spikeIdx(outsideSpikes) = []; spikeTime_downsampled(outsideSpikes) = []; end
nSpikes = length(spikeIdx); val = ones(nSpikes, 1);

spikes = sparse(spikeIdx, spikeTime_downsampled, val, ap.nGoodClusters, nDownSamples, nSpikes);
params.ephys.nSamples = nDownSamples;
params.ephys.nAPSamplePerBin = spikeDownSample;
params.ephys.nSpikes = nSpikes;
signalRange = 1:spikeDownSample*nDownSamples;
params.ephys.signalRange = signalRange;
processed.ephys.signals = spikes;

% Neuropixels util
session.fullImecPath = strcat(session.pathImec, session.apBin);
imec = Neuropixel.ImecDataset(session.fullImecPath);
fprintf('Duration of recording %s is %g minutes\n', imec.fileStem, imec.nSamplesAP / imec.fsAP / 60);
ks = Neuropixel.KilosortDataset(session.pathImec);
% ks.load()

disp(['Finished: kilosort data for session ',sessionName,' loaded']);

%% Generate laser-triggered firing rate for each good cluster

% 20220523
% laser_pulse_duration(1:300) = 0.01;
% laser_pulse_duration(301:600) = 0.05;
% laser_pulse_duration(601:900) = 0.02;
% stim_per_pattern = 100; nPatterns = 9; % npulse_per_stim = 1; 

% 20220525
% stim_per_pattern = 200; nPatterns = 5; % npulse_per_stim = 1; 

% 20220908
%stim_per_pattern = 50; nPatterns = 3;

% 20220612
% 25uW - 300uW at 25uW increments, repeat 50x, 20ms pulse @ 1Hz

%% Load peri-stim spike times

% % For ally
% optoStimIdx = find(blueLaser); % Find blue laser pulses
% % Define params
% timeRangeLong = [-1, 5];
% timeRangeShort = [-0.05, 0.2];
% binSizeLong = 0.05; % in sec
% binSizeShort = 0.001;

% For shun
blueLaserON = find(blueLaser); % Find blue laser pulses
optoStimIdx = blueLaserON(1:end);
redLaserON = find(redLaser);
% Define params
timeRangeLong = [-0.5, 1];
timeRangeShort = [-0.05, 0.2];
binSizeLong = 0.01; % in sec
binSizeShort = 0.001;

% Find spikes around each laser pulse
[~,optoSpikeRateLong,opto_long_params] = getSpikes(timeRangeLong,binSizeLong,optoStimIdx,params);
[optoSpikeShort,optoSpikeRateShort,opto_short_params] = getSpikes(timeRangeShort,binSizeShort,optoStimIdx,params);

return

%% Plot individual graph for a single neuron
clusterList = [53]; % ally
initializeFig(0.5,0.5);
plotSpikes(optoSpikeRateLong,clusterList,timeRangeLong,bluePurpleRed,average=false,textOn=true,text_source=ap,unit='ms');
plotEvent('Stim',0.02,[7 136 225]);

%% Plot different timeRange
initializeFig(0.5,0.5);
clusterList = [70]; textOn = true;

% plot raster
subplot(4,1,1);
plotSpikeRaster(optoSpikeRateShort,clusterList,[-0.005,0.01],unit='ms',originalTimeRange=timeRangeShort);
plotEvent('Stim',0.02,[7 136 225],unit='ms');

subplot(4,1,2);
plotSpikeRaster(optoSpikeRateShort,clusterList,timeRangeShort,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');

subplot(4,1,3);
plotSpikes(optoSpikeRateShort,clusterList,timeRangeShort,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');

subplot(4,1,4);
plotSpikes(optoSpikeRateLong,clusterList,timeRangeLong,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');

%% Plot average of all neurons

% mean_optoSpikeRateLong = squeeze(mean(optoSpikeRateLong,2));
% mean_optoSpikeRateShort = squeeze(mean(optoSpikeRateShort,2));

initializeFig(0.5,0.5);

subplot(2,1,1)
plotSpikes(optoSpikeRateShort,1:ap.nGoodClusters,timeRangeShort,bluePurpleRed,average=true,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');
xlabel('Time (s)'); ylabel('Spikes/s'); box off

subplot(2,1,2)
plotSpikes(optoSpikeRateLong,1:ap.nGoodClusters,timeRangeLong,bluePurpleRed,average=true,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');
xlabel('Time (s)'); ylabel('Spikes/s'); box off

% save figure
saveas(gcf,strcat(sessionpath,'\PSTH-Opto-afterPhy\psth_stim_population_average.png'));

%% Find spikes for SALT
% Find baseline segment
baselineIdx = round(rand([1 length(optoStimIdx)]) * length(airpuff(1:optoStimIdx(1))));
[spt_baseline,~,spike_params_baseline] = getSpikes([0,2],0.001,baselineIdx,params);

% Find optostim segment
[spt_test,~,spt_test_params] = getSpikes([0,0.015],binSizeShort,optoStimIdx,params);

%% Run salt for each laser pulse for each cell

salt_p = nan(ap.nGoodClusters,1);
salt_i = nan(ap.nGoodClusters,1);
salt_hist = nan(ap.nGoodClusters,11,301);

for neuron = 1:ap.nGoodClusters
    neuron_spt_baseline = logical(squeeze(spt_baseline(neuron,:,:)));
    neuron_spt_test = logical(squeeze(optoSpikeShort(neuron,:,:)));
    [salt_p(neuron),salt_i(neuron),salt_hist(neuron,:,:)] = salt(neuron_spt_baseline,neuron_spt_test,opto_short_params.binSize);
    disp(['Finished: SALT for neuron ',num2str(neuron)]);

    
end

% initializeFig(0.5,0.5);
% subplot(1,2,1); histogram(salt_p,50); box off
% subplot(1,2,2); histogram(salt_i,50); box off

% Plot latency histogram
neuron =5;
initializeFig(0.5,0.5);
hist_bl = mean(squeeze(salt_hist(neuron,:,1:end-1)),2);
bar(hist_bl'); hold on
bar(squeeze(salt_hist(neuron,:,end))');
title(['Cluster: ',num2str(neuron)]);
autoArrangeFigures();

optotag_candidates = find(salt_p <= 0.01);
box off

%% Get waveform (calculate correlation between waveforms)

timeRange = [-0.001,0.001];  % Time (sec) before and after spiketime to include in waveform

optotag_candidates = 1:ap.nGoodClusters;
[optotag_wvf,control_wvf] = getSpikeWaveforms(mmf.Data.x,optoStimIdx,optotag_candidates,timeRange,params);

%% Plot waveform
t = linspace(timeRange(1),timeRange(2),size(control_wvf,3)) * 1000;

for i = 1:length(optotag_candidates)
    initializeFig(0.5,0.5);
    optotag_mean = mean(squeeze(optotag_wvf(i,:,:)));
    baseline_mean = mean(squeeze(control_wvf(i,:,:)));

%     for j = 1:size(optotag_wvf,2)
%         plot(t,squeeze(optotag_wvf(i,j,:)),"Color",blueWhiteRed(400,:)); hold on
%         plot(t,squeeze(control_wvf(i,j,:)),"Color",blueWhiteRed(100,:)); hold on
%     end
%     plot(t,optotag_mean,"Color",blueWhiteRed(500,:),'LineWidth',1.5); hold on
%     plot(t,baseline_mean,"Color",blueWhiteRed(1,:),'LineWidth',1.5); hold on

    plotSEM(t,squeeze(optotag_wvf(i,:,:)),blueWhiteRed(500,:)); hold on
    plotSEM(t,squeeze(control_wvf(i,:,:)),blueWhiteRed(1,:)); hold on

    xlabel('Time (ms)'); ylabel('Voltage (mV)'); title(['Cluster: ', num2str(optotag_candidates(i))]);
    legendLines = [optotag_mean(1) baseline_mean(1)];
    legend(legendLines,'Optotag waveform','Baseline waveform');
    box off
end
autoArrangeFigures;

%% Light-evoked energy vs waveform correlation (Cohen et al., 2012)

optotag_mean = mean(squeeze(optotag_wvf(1,:,:)));
baseline_mean = mean(squeeze(control_wvf(1,:,:)));

% Calculate waveform correlation
r = xcorr(optotag_mean,baseline_mean,0,'normalized');

% Calculate light-evoked energy (integral of the squared voltage values)
opto_energy = sum(optotag_mean.^2);

%% Plot every single neuron and save in a folder

for idx = 1:ap.nGoodClusters
    i = idx;
    initializeFig(1,1); textOn = true;

    % plot raster
    % subplot(4,1,1);
    subplot(4,3,[1 2]);
    plotSpikeRaster(optoSpikeRateShort,i,[-0.005,0.01],unit='ms',originalTimeRange=timeRangeShort);
    plotEvent('Stim',0.02,[7 136 225],unit='ms');
    
    % subplot(4,1,2);
    subplot(4,3,[4 5]);
    plotSpikeRaster(optoSpikeRateShort,i,timeRangeShort,unit='ms');
    plotEvent('Stim',0.02,[7 136 225],unit='ms');
    
    % subplot(4,1,3);
    subplot(4,3,[7 8]);
    plotSpikes(optoSpikeRateShort,i,timeRangeShort,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
    plotEvent('Stim',0.02,[7 136 225],unit='ms');
    
    subplot(4,1,4);
    subplot(4,3,[10 11]);
    plotSpikes(optoSpikeRateLong,i,timeRangeLong,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
    plotEvent('Stim',0.02,[7 136 225],unit='ms');

    % plot raster

    % Get waveform
    timeRange = [-0.001,0.001];
    [optotag_wvf,control_wvf] = getSpikeWaveforms(mmf.Data.x,optoStimIdx,i,timeRange,params);
    t = linspace(timeRange(1),timeRange(2),size(control_wvf,3)) * 1000;

    subplot(4,3,[3 6])
    plotSEM(t,squeeze(optotag_wvf(1,:,:)),blueWhiteRed(500,:)); hold on
    plotSEM(t,squeeze(control_wvf(1,:,:)),blueWhiteRed(1,:)); hold on
    xlabel('Time (ms)'); ylabel('Voltage (mV)'); title(['Cluster: ', num2str(optotag_candidates(i))]);
    legend('Optotag waveform','Baseline waveform','location','southeast');
    box off
    optotag_mean = mean(squeeze(optotag_wvf(1,:,:)));
    baseline_mean = mean(squeeze(control_wvf(1,:,:)));
    r = xcorr(optotag_mean,baseline_mean,0,'normalized');

    [val_max(1),x_max(1)] = max(abs(optotag_mean));
    [val_max(2),x_max(2)] = max(abs(baseline_mean));
    [~,max_trace] = max([val_max(1),val_max(2)]);
    xtxtpos = t(x_max(max_trace));
    ytxtpos = val_max(max_trace);
    text(xtxtpos,ytxtpos,strcat('r = ',num2str(r))); hold on

    subplot(4,3,[9 12])
    hist_bl = mean(squeeze(salt_hist(i,:,1:end-1)),2);
    bar(hist_bl'); hold on
    bar(squeeze(salt_hist(neuron,:,end))');
    xlabel('Time (10ms)'); ylabel('Normalized occurance'); legend('Baseline','Opto');
    title(['Cluster: ',num2str(neuron)]); box off

    
    % save figure
    saveas(gcf,strcat(sessionpath,'\PSTH-Opto-afterPhy\psth_stim_cluster',num2str(i),'.png'));
    close;
end
return

%% (Old) Plot distribution of first significant bin after stim 

%{
Rationale: only calculating first spike latency biases units with high
firing rates. Since striatal unit is quiet, code below perform following
three steps:

1. Calculate number of spikes for time bins before and after stimulation.
There's a spike_in_bin matrix where row is neuron and column
is spike count for each time bin accumulated across all stim trials.

    - Note: baseline is calculated as the averaged time bin (accumulated 
    across multiple time bin before baseline and calculate the mean)
    - Each time bin creates an approx Poisson? distribution, as basically
    the neuron can be fire or not fire during the time bin, across n stim
    trials

2. After calculating spike count, thereby generating a distribution of
spikes for each time bin, I compare whether two distribution is
statiscially significant or not
    
    - 

3. Plot the histogram of the first time bin where there's statistical
significance, which should be a two modal peak
%}

timeRange = [-0.2 0.2]; binSize = 0.001;
nBaselineBins = abs(timeRange(1)/binSize);
blueLaserON = find(blueLaser);
optoStimIdx = blueLaserON(1:200);
redLaserON = find(redLaser);
optotag_window = [5,25]; % in ms; time of valid optotag spikes

% Intialize first_spike_latency: row is neurons, column is stim num
spikes_in_bin = zeros(ap.nGoodClusters,timeRange(2)/binSize + 1);

% Looping over each stim
for i = 1:length(optoStimIdx)
    spikes = getSpikes(timeRange,binSize,laserOnset(i),ap,nidq,timeImec,timeNI);
    baseline = sum(spikes(:,1:nBaselineBins),2)/nBaselineBins;
    
    % Add to spikes_in_bin
    %disp(spikes_in_bin(388-12,2)); disp(spikes(388-12,5));pause
    spikes_in_bin(:,1) = spikes_in_bin(:,1) + baseline;
    spikes_in_bin(:,2:end) = spikes_in_bin(:,2:end) + spikes(:,nBaselineBins+1:end);
end

% Plot histogram distribution of first_spike_latency
neuron = 46;
%baseline_dist = fitdist(spikes_in_bin(neuron,1),'Binomial');
%x = 0:0.001:5; y = pdf(baseline_dist,x); plot(x,y); hold on
histogram(spikes_in_bin(neuron,1),size(spikes_in_bin,2)); hold on
% histfit(spikes_in_bin(:,1),50,'poisson');

%histogram(spikes_in_bin(:,5),50); hold on
histogram(spikes_in_bin(neuron,10),50); hold on
%histogram(spikes_in_bin(neuron,end),50); hold on

%% (Not used) calculate first spike latency

timeRange = [-0.2 0.2]; binSize = 0.005;
blueLaserON = find(blueLaser);
% laserOnset = find(redLaser);
optoStimIdx = blueLaserON;
optotag_window = [5,25]; % in ms; time of valid optotag spikes

% Intialize first_spike_latency: row is neurons, column is stim num
first_spike_latency = nan(ap.nGoodClusters,length(optoStimIdx));
nBins = (timeRange(2)-timeRange(1)) / (1/binSize);
% first_spike_time = zeros(ap.nGoodClusters,length(optoStimIdx));

% Looping over each stim
for i = 1:length(optoStimIdx)
    spikes = getSpikes(timeRange,binSize,blueLaserON,ap,nidq,timeImec,timeNI);
    for neuron = 1:ap.nGoodClusters
        first_spike_bin = find(spikes(neuron,:),1);
        if isempty(first_spike_bin); continue; end
        first_spike_latency(neuron,i) = binSize * first_spike_bin;
    end
end

% Convert to ms
first_spike_latency = first_spike_latency .* 1000;

% Plot histogram distribution of first_spike_latency
first_spike_latency_mean = mean(first_spike_latency,2,'omitnan');
optotag_candidates = find(first_spike_latency_mean >= optotag_window(1) & first_spike_latency_mean <= optotag_window(2));
figure(2); histogram(mean(first_spike_latency,2,'omitnan'),100,'FaceColor',colors(1));
xlabel('First spike latency (ms)'); ylabel('Number of neurons');
box off

%% (Old) Retrieve spike waveform of optotag unit (for only one stim pattern)

timeRange = [-0.0005,0.001];  % Time (sec) before and after spiketime to include in waveform
timesteps = floor(ap.Fs*(timeRange(2)-timeRange(1)));

optotag_candidates = [113 22 233]; optoStimIdx = blueLaserON;
% optotag_wvf = number of candidates x number of optostim x timesteps
optotag_wvf = zeros(length(optotag_candidates),length(optoStimIdx),timesteps);
control_wvf = zeros(length(optotag_candidates),length(optoStimIdx),timesteps);

% Extract control trials start time
%control_stim = [blueLaserOnIdx(1),blueLaserOnIdx(50)];

for i = 1:length(optotag_candidates)
    % Get best channel of the candidate unit
    cluster_id = ap.goodClusters(optotag_candidates(i));
    best_channel = ap.cluster_info(find(ap.cluster_info(:,1) == cluster_id),6);

    % Get opto-triggered waveform
    for j = 1:length(optoStimIdx)
        % Find Imec stim index
        niStimIdx = optoStimIdx(j) + floor(timeRange(1)*nidq.Fs);
        [~, imecStimIdx] = min(abs(timeImec-timeNI(niStimIdx)));

        % Find first spike
        neuron_spikes_idx = find(ap.goodSpikeClusters == cluster_id);
        neuron_spikes_time = ap.goodSpikeTimes(neuron_spikes_idx);
        first_spikes_time = neuron_spikes_time(find(neuron_spikes_time > imecStimIdx,1));
        % disp(1000*(first_spikes_time - imecStimIdx)/ap.Fs);

        % Find first and last imec index
        imecFirstIdx = first_spikes_time + floor(ap.Fs*timeRange(1));
        imecLastIdx = imecFirstIdx + timesteps - 1;

        % Read waveform within range
        wf = mmf.Data.x(best_channel, imecFirstIdx:imecLastIdx);

        % Spike waveform for each stim
        optotag_wvf(i,j,:) = wf; 

        disp(['Stim ',num2str(j),'/',num2str(length(optoStimIdx)),...
            ' for neuron ',num2str(optotag_candidates(i)),' is analyzed']);
    end

    % Randomly select non optotim spikes from this unit
    random_spikes_time = neuron_spikes_time(randi(length(optoStimIdx),[length(optoStimIdx),1]));
    for s = 1:length(random_spikes_time)
        % Find first and last imec index
        imecFirstIdx = random_spikes_time(s) + floor(ap.Fs*timeRange(1));
        imecLastIdx = imecFirstIdx + timesteps - 1;
        
        % Read waveform within range
        wf = mmf.Data.x(best_channel, imecFirstIdx:imecLastIdx);

        % Spike waveform for each stim
        control_wvf(i,s,:) = wf; 

        disp(['Baseline ',num2str(s),'/',num2str(length(optoStimIdx)),...
            ' for neuron ',num2str(optotag_candidates(i)),' is analyzed']);
    end
end

%% (Old) Plot spike waveform

t = linspace(timeRange(1),timeRange(2),timesteps) * 1000;

for i = 1:length(optotag_candidates)
    initializeFig(0.5,0.5);
    for j = 1:size(optotag_wvf,2)
        plot(t,squeeze(optotag_wvf(i,j,:)),"Color",blueWhiteRed(400,:)); hold on
        plot(t,squeeze(control_wvf(i,j,:)),"Color",blueWhiteRed(100,:)); hold on
    end
    optotag_mean = plot(t,mean(squeeze(optotag_wvf(i,:,:))),"Color",blueWhiteRed(500,:),'LineWidth',1.5); hold on
    baseline_mean = plot(t,mean(squeeze(control_wvf(i,:,:))),"Color",blueWhiteRed(1,:),'LineWidth',1.5); hold on

    %plotSEM(t,squeeze(optotag_wvf(i,:,:)),blueWhiteRed(500,:)); hold on
    %plotSEM(t,squeeze(control_wvf(i,:,:)),blueWhiteRed(1,:)); hold on

    xlabel('Time (ms)'); ylabel('Signal'); title(['Cluster: ', num2str(optotag_candidates(i))]);
    legendLines = [optotag_mean(1) baseline_mean(1)];
    legend(legendLines,'Optotag waveform','Baseline waveform');
    box off
end
autoArrangeFigures;

%% (Old) Plot firing rate difference

% Set up parameters
timeRange = [-1, 50];
binSize = 0.005; % in sec
nStim = nPatterns * stim_per_pattern;
blueLaserOnIdx = blueLaserON(101);
laser_pulse_duration = zeros(1,nStim);
% Generate PSTH
[optoSpikes,optoSpikeRate] = getSpikes(timeRange,binSize,blueLaserOnIdx,...
                                    ap,nidq,timeImec,timeNI);

% Plot PSTH for candidate neuron across trials
clusterList = optotag_candidates; textOn = true;
figure; drawSpikeRate('Stim',timeRange,optoSpikeRate,clusterList,textOn,ap,colors);

clusterList = 123; textOn = true;
figure; drawSpikeRate('Stim',timeRange,optoSpikeRate,clusterList,textOn,ap,colors);
for i = 1:round(timeRange(2))
    xline(i,'-' ,'Stim','Color','r','LineWidth',1.5);
end

% Example artifact unit
% 20220523: artifact: 56, 318
% 20220525: artifact: 741, 64; others: 225, 23 
%clusterList = ap.clusterToGoodClusterIndex(64); textOn = true;
%figure; drawSpikeRate('Stim',timeRange,optoSpikeRate,clusterList,textOn,ap,colors);
