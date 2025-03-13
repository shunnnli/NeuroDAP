% Shun_analyzeOptotag
% Shun Li, 12/16/2022

%% Load data
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

setenv('NEUROPIXEL_MAP_FILE', which('neuropixPhase3B2_kilosortChanMap.mat'));
% setenv('NEUROPIXEL_DATAROOT', 'D:\Shun\Analysis\Result-Neuropixel utils');
[~,~,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed] = loadColors;

% Select session via uigetdir
sessionpath = uigetdir(osPathSwitch('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun'));
dirsplit = strsplit(sessionpath,filesep); 
sessionName = dirsplit{end}; session.projectPath = strcat('\\',fullfile(dirsplit{2:end-1})); clear dirsplit

disp(strcat('********** ',sessionName,'**********'));
load(strcat(sessionpath,filesep,'data_',sessionName,'.mat'));
load(strcat(sessionpath,filesep,'sync_',sessionName,'.mat'));
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

% chMap = readNPY(fullfile(session.pathImec, 'channel_map.npy'))+1; % Order in which data was streamed to disk; must be 1-indexed for Matlab
%For kevin
chMap = readNPY(fullfile(session.sorterOutput, 'channel_map.npy'))+1;

nChInMap = numel(chMap);
disp(['Finished: ', gwfparams.fileName,' for session ',sessionName,' loaded']);

%% Load kilosort
[clusterLabel,spike_times,spike_clusters,cluster_info] = readNPYData(session.sorterOutput); % originally pathImec

% Decide which label is used for good unit
answer = questdlg("Which label criteria will be used for good unit?",...
    "Select label criteria","Kilosort good","Manual good","Manual good & unsorted","Manual good & unsorted");
switch answer
    case 'Kilosort'
        ap.goodClusters = clusterLabel.cluster_id(find(clusterLabel.KSLabel == 'good')); % Good units;
    case 'Manual good'
        ap.goodClusters = clusterLabel.cluster_id(find(clusterLabel.KSLabel == 'good' & clusterLabel.group == 'good')); % Good units;
    case 'Manual good & unsorted'
        ap.goodClusters = clusterLabel.cluster_id(find((clusterLabel.group == 'good' | clusterLabel.group == 'unsorted') & clusterLabel.KSLabel == 'good')); % Good units;
end
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
% session.fullImecPath = strcat(session.pathImec, session.apBin);
% imec = Neuropixel.ImecDataset(session.fullImecPath);
% fprintf('Duration of recording %s is %g minutes\n', imec.fileStem, imec.nSamplesAP / imec.fsAP / 60);
% ks = Neuropixel.KilosortDataset(session.pathImec);
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
% binSizeShort = 0.0001;

% For shun
blueLaserON = find(blueLaser); % Find blue laser pulses
optoStimIdx = blueLaserON(1:end);
redLaserON = find(redLaser);
% Define params
timeRangeLong = [-0.5, 1];
timeRangeShort = [-0.05, 0.2];
binSizeLong = 0.01; % in sec
binSizeShort = 0.0001;

% Find spikes around each laser pulse
[~,optoSpikeRateLong,opto_long_params] = getSpikes(timeRangeLong,binSizeLong,optoStimIdx,params);
[optoSpikeShort,optoSpikeRateShort,opto_short_params] = getSpikes(timeRangeShort,binSizeShort,optoStimIdx,params);

%% Find peri-water spikes
waterON = find(leftSolenoid);
airpuffON = find(airpuff);
leftToneON = find(leftTone);
rightToneON = find(rightTone);
timeRange = [-1.5,2];

[~,waterSpikeRate,water_params] = getSpikes(timeRange,0.01,waterON,params);
[~,airpuffSpikeRate,~] = getSpikes(timeRange,0.01,airpuffON,params);
% [~,leftToneSpikeRate,~] = getSpikes(timeRange,0.01,leftToneON,params);
% [~,rightToneSpikeRate,~] = getSpikes(timeRange,0.01,rightToneON,params);

return
%% Calculate modulation index

timeWindowPre = 100:150; timeWindowPost = 151:200;

water_pre = waterSpikeRate(:, :, timeWindowPre);
water_post = waterSpikeRate(:, :, timeWindowPost);
[significance_water, pvalue_water] = ranksumSignificance(water_pre,water_post,1:ap.nGoodClusters);
water_modIdx = getModulationIdx(water_pre,water_post,1:ap.nGoodClusters);
water_significant_unit = find(abs(water_modIdx) > 0.5);
water_unit_modIdx = water_modIdx(water_significant_unit);

airpuff_pre = airpuffSpikeRate(:, :, timeWindowPre);
airpuff_post = airpuffSpikeRate(:, :, timeWindowPost);
[significance_airpuff, pvalue_airpuff] = ranksumSignificance(airpuff_pre,airpuff_post,1:ap.nGoodClusters);
airpuff_modIdx = getModulationIdx(airpuff_pre,airpuff_post,1:ap.nGoodClusters);
airpuff_significant_unit = find(abs(airpuff_modIdx) > 0.5);
airpuff_unit_modIdx = airpuff_modIdx(airpuff_significant_unit);

% leftTone_pre = leftToneSpikeRate(:, :, timeWindowPre);
% leftTone_post = leftToneSpikeRate(:, :, timeWindowPost);
% [significance_leftTone, pvalue_leftTone] = ranksumSignificance(leftTone_pre,leftTone_post,1:ap.nGoodClusters);
% leftTone_modIdx = getModulationIdx(leftTone_pre,leftTone_post,1:ap.nGoodClusters);
% leftTone_significant_unit = find(abs(leftTone_modIdx) > 0.5);
% leftTone_unit_modIdx = leftTone_modIdx(leftTone_significant_unit);
% 
% rightTone_pre = rightToneSpikeRate(:, :, timeWindowPre);
% rightTone_post = rightToneSpikeRate(:, :, timeWindowPost);
% [significance_rightTone, pvalue_rightTone] = ranksumSignificance(rightTone_pre,rightTone_post,1:ap.nGoodClusters);
% rightTone_modIdx = getModulationIdx(rightTone_pre,rightTone_post,1:ap.nGoodClusters);
% rightTone_significant_unit = find(abs(rightTone_modIdx) > 0.5);
% rightTone_unit_modIdx = rightTone_modIdx(rightTone_significant_unit);
% 
% % all tone
% allToneSpikeRate = cat(2,leftToneSpikeRate,rightToneSpikeRate);
% allTone_pre = allToneSpikeRate(:, :, timeWindowPre);
% allTone_post = allToneSpikeRate(:, :, timeWindowPost);
% [significance_allTone, pvalue_allTone] = ranksumSignificance(allTone_pre,allTone_post,1:ap.nGoodClusters);
% allTone_modIdx = getModulationIdx(allTone_pre,allTone_post,1:ap.nGoodClusters);
% allTone_significant_unit = find(abs(allTone_modIdx) > 0.5);
% allTone_unit_modIdx = allTone_modIdx(allTone_significant_unit);

%% Peri-water spike rates
clusterList = optotagged;

initializeFig(0.5,0.5);
plotSpikes(waterSpikeRate,clusterList,timeRange,bluePurpleRed,average=false,unit='ms',smooth=true);
plotEvent('Water',0,[7 136 225],unit='ms');
xlabel('Time (s)'); ylabel('Spikes/s'); box off

%% Summary plot of all modulated units
clusterList = [55 162];

initializeFig(1,1);
tiledlayout(2,1);
nexttile;
plotSpikes(waterSpikeRate,clusterList,timeRange,bluePurpleRed,average=false,unit='ms',smooth=true);
plotEvent('Water',0,[7 136 225],unit='ms'); plotEvent('Tone',0,[7 136 225],unit='ms',x=-1000); box off

nexttile;
plotSpikes(airpuffSpikeRate,clusterList,timeRange,bluePurpleRed,average=false,unit='ms',smooth=true);
plotEvent('Airpuff',0,[7 136 225],unit='ms'); plotEvent('Tone',0,[7 136 225],unit='ms',x=-1000); box off

% nexttile(3,[1 2]);
% plotSpikes(allToneSpikeRate,clusterList,timeRange,bluePurpleRed,average=false,unit='ms',smooth=true);
% plotEvent('Tone',0.25,[7 136 225],unit='ms'); plotEvent('Outcome',0,[7 136 225],unit='ms',x=1000); box off

% nexttile;
% plotSpikes(leftToneSpikeRate,clusterList,timeRange,bluePurpleRed,average=false,unit='ms',smooth=true);
% plotEvent('Left tone',0,[7 136 225],unit='ms'); box off

% nexttile;
% plotSpikes(rightToneSpikeRate,clusterList,timeRange,bluePurpleRed,average=false,unit='ms',smooth=true);
% plotEvent('Right tone',0,[7 136 225],unit='ms'); box off

saveFigures(gcf,'PSTH-outcome-unit55&162',session.path);

%% Plot average of all neurons

initializeFig(0.5,0.5);

unitList = 1:ap.nGoodClusters;

subplot(2,1,1)
plotSpikes(optoSpikeRateShort(:,1:400,:),unitList,timeRangeShort,bluePurpleRed,average=true,unit='ms',smooth=true);
plotEvent('Stim',0.02,color=[7 136 225],unit='ms');
xlabel('Time (s)'); ylabel('Spikes/s'); box off

subplot(2,1,2)
plotSpikes(optoSpikeRateLong(:,1:400,:),unitList,timeRangeLong,bluePurpleRed,average=true,unit='ms');
plotEvent('Stim',0.02,color=[7 136 225],unit='ms');
xlabel('Time (s)'); ylabel('Spikes/s'); box off

% save figure
% saveas(gcf,strcat(sessionpath,'\PSTH-allOpto_full\psth_stim_population_average.png'));

%% Find spikes for SALT

% Salt params
% Shun
salt_binSize = 0.001; % Should be 0.001
salt_window = 0.02; % 20ms
nBinsInHist = round(salt_window/salt_binSize);

% Find baseline segment
baselineIdx = round(rand([1 length(optoStimIdx)]) * length(airpuff(1:optoStimIdx(1))));
[spt_baseline,~,spike_baseline_params] = getSpikes([0,2],salt_binSize,baselineIdx,params);

% Find optostim segment
[spt_test,~,spt_test_params] = getSpikes([0,salt_window+0.005],salt_binSize,optoStimIdx,params,maxLatency=salt_window);

%% Run salt for each laser pulse for each cell

% Initialize result matrix
salt_p = nan(ap.nGoodClusters,1);
salt_i = nan(ap.nGoodClusters,1);
salt_hist = nan(ap.nGoodClusters,nBinsInHist+1,size(spt_baseline,3)/nBinsInHist + 1);

for neuron = 1:ap.nGoodClusters
    neuron_spt_baseline = logical(squeeze(spt_baseline(neuron,:,:)));
    neuron_spt_test = logical(squeeze(spt_test(neuron,:,:)));
    [salt_p(neuron),salt_i(neuron),salt_hist(neuron,:,:)] = salt(neuron_spt_baseline,neuron_spt_test,salt_binSize,salt_window);
    disp(['Finished: SALT for neuron ',num2str(neuron)]);
end

% initializeFig(0.5,0.5);
% subplot(1,2,1); histogram(salt_p,50); box off
% subplot(1,2,2); histogram(salt_i,50); box off

% Plot latency histogram
% i = 64;
% initializeFig(0.5,0.5);
% hist_bl = mean(squeeze(salt_hist(i,:,1:end-1)),2);
% x = categorical({'0','1','2','3','4','5','6','7','8','9','10',...
%     '11','12','13','14','15','16','17','18','19','20'});
% x = reordercats(x, {'0','1','2','3','4','5','6','7','8','9','10',...
%     '11','12','13','14','15','16','17','18','19','20'});
% bar(x, [hist_bl';squeeze(salt_hist(i,:,end))]);
% xlabel('Time (ms)'); ylabel('Number of trials'); legend('Baseline','Opto');
% set(gca, 'YScale', 'log');
% title(['Cluster: ',num2str(i),'; p-value: ',num2str(salt_p(i))]);
% autoArrangeFigures();
% 
% optotag_candidates = find(salt_p <= 0.01);
% box off

%% Save opto-triggered spikes

save(strcat(session.path,'/','optotag_',sessionName),'optoStimIdx',...
    'optoSpikeRateShort','optoSpikeShort','opto_short_params',...
    'optoSpikeRateLong','opto_long_params',...
    'spt_baseline','spike_baseline_params',...
    'spt_test','spt_test_params',...
    'salt_binSize','salt_window',...
    '-v7.3');

%% (With waveform) Plot every single neuron summary and save in a folder

optoStimDuration = 20; % Shun, in ms
% optoStimDuration = 0; % Ally

optoTrialRange = 1:600;

wvf_r = nan(1,ap.nGoodClusters);
for idx = 64%1:542%ap.nGoodClusters
    i = idx;
    initializeFig(1,1); textOn = true;
    tiledlayout(4,3);

    % plot raster
    nexttile(1,[1 2]);
    plotSpikeRaster(optoSpikeRateShort(:,optoTrialRange,:),i,[-0.005,salt_window],unit='ms',originalTimeRange=timeRangeShort,size=20);
    plotEvent('Stim',optoStimDuration,color=[7 136 225],unit='s');
    %yline([250],'--');
    %yline([600 1200 1800 2400],'-');
    
    nexttile(4,[1 2]);
    plotSpikeRaster(optoSpikeRateShort(:,optoTrialRange,:),i,timeRangeShort,unit='ms',size=20);
    plotEvent('Stim',optoStimDuration,color=[7 136 225],unit='s');
    %yline([250],'--');
    %yline([600 1200 1800 2400],'-');
    
    nexttile(7,[1 2]);
    plotSpikes(optoSpikeRateShort(:,optoTrialRange,:),i,timeRangeShort,bluePurpleRed,textOn=true,text_source=ap,unit='ms',smooth=true);
    plotEvent('Stim',optoStimDuration,color=[7 136 225],unit='s');
    
    nexttile(10,[1 2]);
    plotSpikes(optoSpikeRateLong(:,optoTrialRange,:),i,timeRangeLong,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
    plotEvent('Stim',optoStimDuration,color=[7 136 225],unit='s');


    % Get waveform
    timeRange = [-0.0020,0.0025];
    [neuron_optotag_wvf,neuron_control_wvf] = getSpikeWaveforms(mmf.Data.x,i,timeRange,params,...
                                                    event=optoStimIdx,maxLatency=salt_window,...
                                                    triggeredSpikeIdx=spt_test_params.triggeredSpikeIdx);
    t = linspace(timeRange(1),timeRange(2),size(neuron_control_wvf,3)) * 1000;

    nexttile(3,[2 1]);
    optotag_trace = rmmissing(squeeze(neuron_optotag_wvf(1,:,:)));
    baseline_trace = rmmissing(squeeze(neuron_control_wvf(1,:,:)));
    r = xcorr(mean(optotag_trace,1),mean(baseline_trace,1),0,'normalized');
    wvf_r(i) = r;
    % Plot avg traces
    plotSEM(t,baseline_trace,blueWhiteRed(1,:)); hold on
    plotSEM(t,optotag_trace,blueWhiteRed(500,:)); hold on
    xlabel('Time (ms)'); ylabel('Voltage (mV)'); 
    legend({['Baseline (n=',num2str(size(baseline_trace,1)),')'],...
        ['Opto (n=',num2str(size(optotag_trace,1)),')']});
    title(['r = ',num2str(r)]); box off

    nexttile(9,[2 1]);
    hist_bl = mean(squeeze(salt_hist(i,:,1:end-1)),2);
    x = categorical({'0','1','2','3','4','5','6','7','8','9','10',...
                '11','12','13','14','15','16','17','18','19','20'});
    x = reordercats(x, {'0','1','2','3','4','5','6','7','8','9','10',...
                '11','12','13','14','15','16','17','18','19','20'});
    bar(x, [hist_bl';squeeze(salt_hist(i,:,end))]);
    xlabel('First spike latency (ms)'); ylabel('Number of trials'); legend('Baseline','Opto');
    set(gca, 'YScale', 'log');
    title(['p-value: ',num2str(salt_p(i))]); box off
    
    % save figure
    % saveFigures(gcf,strcat('psth_stim_cluster',num2str(i)),strcat(sessionpath,'\PSTH-allOpto_1mW_pulse'),savePDF=true,saveFig=true);
    % close;
end

%% SALT p-value vs waveform correlation

% Plot scatter
initializeFig(.5,.5);

for i = 1:ap.nGoodClusters
    if salt_p(i) <= 0.01 && wvf_r(i) >= 0.5
        scatter(salt_p(i),wvf_r(i),'filled','MarkerFaceColor',blueWhiteRed(500,:)); hold on
        text(salt_p(i),wvf_r(i),num2str(i)); hold on
    else
        scatter(salt_p(i),wvf_r(i),'MarkerEdgeColor',blueWhiteRed(1,:)); hold on
        % text(salt_p(i),wvf_r(i),num2str(i)); hold on
    end
end

xlabel('p-value'); 
xlim([0 0.01]); 
%set(gca, 'XScale', 'log');
ylabel('Waveform correlation');
return

%% Optotag list Shun

optotagged = [55 57 64 65 67 81 84 85 86 92 97 98 99 103 104,...
                113 140 146 152 158 161 162 164 168,...
                194 195 196 201 213];
maybe_optotagged = [54 60 70 74 89 101 103 146 150 151 159 166 167 195 196 202 204 210 211 216];

%% Plot individual graph for a single neuron
clusterList = [64];
initializeFig(0.5,0.5);
plotSpikes(optoSpikeRateLong,clusterList,timeRangeLong,bluePurpleRed,average=false,textOn=true,text_source=ap,unit='ms');
plotEvent('Stim',0.02,[7 136 225]);

%% Plot different timeRange
initializeFig(0.5,0.5);
clusterList = 213; textOn = true;

% plot raster
subplot(4,1,1);
plotSpikeRaster(optoSpikeRateShort,clusterList,[-0.005,0.02],unit='ms',originalTimeRange=timeRangeShort);
plotEvent('Stim',0.02,[7 136 225],unit='ms');

subplot(4,1,2);
plotSpikeRaster(optoSpikeRateShort,clusterList,timeRangeShort,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');

subplot(4,1,3);
plotSpikes(optoSpikeRateShort,clusterList,timeRangeShort,bluePurpleRed,textOn=true,text_source=ap,unit='ms',smooth=true);
plotEvent('Stim',0.02,[7 136 225],unit='ms');

subplot(4,1,4);
plotSpikes(optoSpikeRateLong,clusterList,timeRangeLong,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');

%% Get waveform (calculate correlation between waveforms)

timeRange = [-0.001,0.001];  % Time (sec) before and after spiketime to include in waveform

[neuron_optotag_wvf,neuron_control_wvf] = getSpikeWaveforms(mmf.Data.x,optotag_candidates,timeRange,params,...
                                              event=optoStimIdx,maxLatency=0.02,...
                                              triggeredSpikeIdx=spt_test_params.triggeredSpikeIdx);

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

%% (Without waveform) Plot every single neuron summary and save in a folder

for idx = 1:ap.nGoodClusters
    i = idx;
    initializeFig(1,1); textOn = true;
    
    % plot raster
    subplot(4,1,1);
    plotSpikeRaster(optoSpikeRateShort,i,[-0.005,0.01],unit='ms',originalTimeRange=timeRangeShort);
    plotEvent('Stim',0.02,[7 136 225],unit='ms');
    
    subplot(4,1,2);
    plotSpikeRaster(optoSpikeRateShort,i,timeRangeShort,unit='ms');
    plotEvent('Stim',0.02,[7 136 225],unit='ms');
    
    subplot(4,1,3);
    plotSpikes(optoSpikeRateShort,i,timeRangeShort,bluePurpleRed,textOn=true,text_source=ap,unit='ms',smooth=true);
    plotEvent('Stim',0.02,[7 136 225],unit='ms');
    
    subplot(4,1,4);
    plotSpikes(optoSpikeRateLong,i,timeRangeLong,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
    plotEvent('Stim',0.02,[7 136 225],unit='ms');
    
    % save figure
    saveas(gcf,strcat(sessionpath,'\PSTH-Opto-afterSALT\psth_stim_cluster',num2str(i),'.png'));
    close;
end


%% Load LFP

% Load LFP waveforms
lfpwvf.fileName = strcat(session.name,'_t0.imec0.lf.bin'); 
fileName = fullfile(strcat(session.pathNidq,session.name,'_imec0\'),lfpwvf.fileName);  

% lfpwvf.fileName = strcat(session.name,'_tcat.imec0.lf.bin'); 
% fileName = fullfile(strcat(session.pathNidq,'catgt_',session.name),lfpwvf.fileName);  

lfpwvf.nCh = 385; lfpwvf.dataType = 'int16';
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, lfpwvf.dataType), 'uint8')); % Determine number of bytes per sample
nSamp = filenamestruct.bytes/(lfpwvf.nCh*dataTypeNBytes);  % Number of samples per channel
lfpmmf = memmapfile(fileName, 'Format', {lfpwvf.dataType, [lfpwvf.nCh nSamp], 'x'},'Writable',true);
disp(['Finished: ', lfpwvf.fileName,' for session ',sessionName,' loaded']);

%% get sample LFP waveform for quality inspection

timeRange = [0 0.5];
sampleLFP = squeeze(getLFPWaveforms(lfpmmf.Data.x,timeRange,params,event=2562367));
t = linspace(timeRange(1),timeRange(2),size(sampleLFP,2)) * 1000;

% Plot sample LFP trace
initializeFig(0.5,1);
max_value = max(max(sampleLFP,[],2));
for i = 1:100
    disp(['Plotting channel #',num2str(i)]);
    plot(t,i+((sampleLFP(i,:)-sampleLFP(i,1))/max_value),'LineWidth',1,'Color','k'); hold on
end
xlabel('Time (ms)'); ylabel('Channel');

%% Power spectrum of raw LFP

% LFPdata = squeeze(optoLFP_filtered(:,56,:));
LFPdata = sampleLFP; timeRange = [0, 0.5];

initializeFig(1,1);
t = (timeRange(1):1/params.sync.lfpFs:timeRange(2)-1/params.sync.lfpFs)';
pspectrum(LFPdata(1,:)',t);

%% Preprocess LFP

% Remove 60Hz noise
notch60 = designfilt('bandstopiir','FilterOrder',12, ...
               'HalfPowerFrequency1',57,'HalfPowerFrequency2',63, ...
               'DesignMethod','butter','SampleRate',params.sync.lfpFs);
%fvtool(notch60);
sampleLFP_filtered = filtfilt(notch60,sampleLFP')';

% Third-order Butterworth low-pass filter with 300 Hz cut-off frequency
% (based on Jun et al., Nature, 2017)
lpfilt = designfilt('lowpassiir','FilterOrder',8, 'PassbandFrequency',300,...
        'PassbandRipple',0.01, 'Samplerate',params.sync.lfpFs);
%fvtool(lpfilt);
sampleLFP_filtered = filtfilt(lpfilt,sampleLFP_filtered')';


% Plot sample LFP trace
initializeFig(0.5,1);
max_value = max(sampleLFP_filtered,[],'all','omitnan');
for i = 1:100
    disp(['Plotting channel #',num2str(i)]);
    % plot(t,i+((sampleLFP(i,:)-sampleLFP(i,1))/max_value),'LineWidth',1,'Color','b'); hold on
    plot(t,i+((sampleLFP_filtered(i,:)-sampleLFP_filtered(i,1))/max_value),'LineWidth',1,'Color','k'); hold on
end
xlabel('Time (ms)'); ylabel('Channel');

%% getLFPWaveforms

blueLaserON = find(blueLaser); % Find blue laser pulses
optoStimIdx = blueLaserON(1:end);

% For shun
% timeRange = [-0.01,0.03];
% optoLFP = getLFPWaveforms(lfpmmf.Data.x,timeRange,params,event=optoStimIdx);
% optoLFP(192,:,:) = nan; % Remove ref channel
% t = linspace(timeRange(1),timeRange(2),size(optoLFP,3)) * 1000;

% For ally 1/05
% First 100 pulse for each power: 1s on, 10s off
% Next 500 pulse for each power: 100ms on, 400ms off
% optoStim_025mW_1s = optoStimIdx(1:100); timeRange = [-1,3];
% optoStim_025mW_100ms = optoStimIdx(101:600); timeRange = [-0.05,0.15];
% optoStim_1mW_1s = optoStimIdx(1801:1900); timeRange = [-1,3];
% optoStim_1mW_100ms = optoStimIdx(1901:2400); timeRange = [-0.05,0.15];

% ally 01/07
% 0-100 stim: 0.25mW 1sec on 10sec off
% 100-200 stim: 0.5mW 1sec on 10sec off
% 201-300 stim: 0.75mW 1sec on 10sec off
% 301-400 stim: 1mW 1sec on 10sec off
% 401-650: 20ms pulse at 1mW, 1sec ITI
pulse1s1mW = optoStimIdx(201:300); timeRange = [-0.01,0.05];
% pulse20ms = optoStimIdx(401:650); timeRange = [-0.01,0.05];


% Retrieve LFP
optoLFP = getLFPWaveforms(lfpmmf.Data.x,timeRange,params,event=pulse1s1mW);
t = linspace(timeRange(1),timeRange(2),size(optoLFP,3)) * 1000;
optoStimID = 56; 

%% Preprocess optoLFP

% optoLFP_notch = filtfilt(notch60,optoLFP');
% optoLFP_filtered = filtfilt(lpfilt,optoLFP_notch)';

optoLFP_filtered = nan(size(optoLFP));
for i = 1:size(optoLFP,2)
    disp(['Filtering event #',num2str(i)]);
    stimTrace = squeeze(optoLFP(:,i,:));

    % CAR
    carTrace = stimTrace - median(stimTrace);
        
    % Filters
    %notchTrace = filtfilt(notch60,stimTrace);
    %lpTrace = filtfilt(lpfilt,notchTrace)';

    optoLFP_filtered(:,i,:) = carTrace;
end

optoLFP_filtered(192,:,:) = 0; % Remove ref channel

% Plot sample LFP trace
% initializeFig(0.5,1); 
% nCh = 385;
% max_value = max(optoLFP_filtered,[],'all','omitnan');
% for i = 1:nCh
%     disp(['Plotting channel #',num2str(i)]);
%     plot(t,i+(squeeze(optoLFP_filtered(i,optoStimID,:)-optoLFP_filtered(i,optoStimID,1))/max_value),'LineWidth',1,'Color','k'); hold on
% end
% xlabel('Time (ms)'); ylabel('Channel'); ylim([1 nCh]);

%% Plot LFP raw traces for all channels

initializeFig(.5,.5); tiledlayout(1,2);

% Individual optostim
nexttile; % optoStimID = 56; 
imagesc('XData',t,'YData',1:lfpwvf.nCh,'CData',squeeze(optoLFP_filtered(:,optoStimID,:)));
xlim([t(1) t(end)]); ylim([1 lfpwvf.nCh]);
colormap(blueWhiteRed); c = colorbar; c.Label.String = 'LFP';
xlabel('Time (ms)'); ylabel('Channel'); title(['Opto stim #',num2str(optoStimID)]);

% Average optostim
nexttile; 
imagesc('XData',t,'YData',1:lfpwvf.nCh,'CData',squeeze(mean(optoLFP_filtered,2)));
xlim([t(1) t(end)]); ylim([1 lfpwvf.nCh]);
colormap(blueWhiteRed); c = colorbar; c.Label.String = 'LFP';

xlabel('Time (ms)'); ylabel('Channel'); title('Averaged opto stim');

%% Plot variance of LFP for all channels

k = params.sync.lfpFs * 0.002; % movvar window around 2ms
optoLFP_movvar = movvar(optoLFP_filtered,k,0,3,'omitnan');

initializeFig(.5,.5); tiledlayout(1,2);

% Individual optostim
nexttile; % optoStimID = 56; 
imagesc('XData',t,'YData',1:lfpwvf.nCh,'CData',squeeze(optoLFP_movvar(:,optoStimID,:)));
xlim([t(1) t(end)]); ylim([1 lfpwvf.nCh]);
colormap(blueWhiteRed); c = colorbar; c.Label.String = 'Variance';
xlabel('Time (ms)'); ylabel('Channel'); title(['Opto stim #',num2str(optoStimID)]);

% Average optostim
nexttile; 
imagesc('XData',t,'YData',1:lfpwvf.nCh,'CData',squeeze(mean(optoLFP_movvar,2)));
xlim([t(1) t(end)]); ylim([1 lfpwvf.nCh]);
colormap(blueWhiteRed); c = colorbar; c.Label.String = 'Variance';

xlabel('Time (ms)'); ylabel('Channel'); title('Averaged opto stim var');
