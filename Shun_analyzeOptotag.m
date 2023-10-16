% Shun_analyzeOptotag
% Shun Li, 12/16/2022

%% Load data
clear; close all;
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

%% Load kilosort
[clusterLabel,spike_times,spike_clusters,cluster_info] = readNPYData(session.pathImec);

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
optoStimIdx = find(blueLaser); % Find blue laser pulses
% Define params
timeRangeLong = [-1, 5];
timeRangeShort = [-0.05, 0.2];
binSizeLong = 0.05; % in sec
binSizeShort = 0.0001;

% For shun
% blueLaserON = find(blueLaser); % Find blue laser pulses
% optoStimIdx = blueLaserON(1:end);
% redLaserON = find(redLaser);
% % Define params
% timeRangeLong = [-0.5, 1];
% timeRangeShort = [-0.05, 0.2];
% binSizeLong = 0.01; % in sec
% binSizeShort = 0.0001;

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

%% Histogram of modIdx

initializeFig();

histogram(water_modIdx,50); hold on
histogram(water_modIdx(clusterList),50); box off

% nexttile
% histogram(airpuff_modIdx,50); hold on
% histogram(airpuff_modIdx(clusterList),50); box off

%% Ally 20ms 

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0';

% For ally
optoStimIdx = find(blueLaser); % Find blue laser pulses
% Define params
timeRangeShort20 = [-0.02, 0.02];
binSizeShort = 0.0001;
[optoSpikeShort20,optoSpikeRateShort20,opto_short_params20] = getSpikes(timeRangeShort20,binSizeShort,optoStimIdx,params);

%separate the spikes out by power 0.25-1mW
optoSpikeRateShort20_025mw = optoSpikeRateShort20(:, 1:100, :);
optoSpikeRateShort20_050mw = optoSpikeRateShort20(:, 101:200, :);
optoSpikeRateShort20_075mw = optoSpikeRateShort20(:, 201:300, :);
optoSpikeRateShort20_1mw = optoSpikeRateShort20(:, 301:400, :);

%separate out pre and post laser stimulation 0.25mW
optoSpikeShort20_025mw_pre = optoSpikeRateShort20_025mw(:, :, 1:200);
optoSpikeShort20_025mw_post = optoSpikeRateShort20_025mw(:, :, 201:400);
[significance_025mW, pvalue_025mW] = ranksumSignificance(optoSpikeShort20_025mw_pre,optoSpikeShort20_025mw_post,1:ap.nGoodClusters);
modIdx_025mW = getModulationIdx(optoSpikeShort20_025mw_pre,optoSpikeShort20_025mw_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_025mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2);
xlabel('p values'); ylabel('Count'); box off
title('0.25mW P values');

ax2 = nexttile;
less_than_05_025mW = sum(pvalue_025mW < 0.05); 
num_non_p_less_05_025mW = length(ap.goodClusters) - less_than_05_025mW;
percent_mod_pvalue = [less_than_05_025mW, num_non_p_less_05_025mW];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod_pvalue)
% lgd = legend(labels);
% lgd.Layout.Tile = 'south';

ax3 = nexttile;
histogram(significance_025mW,10)
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off;
title('0.25mW Pre/Post Signficance');

ax4 = nexttile;
num_sign_025 = length(find(significance_025mW));
num_nonsign025 = length(ap.goodClusters) - num_sign_025;
percent_mod = [num_sign_025, num_nonsign025];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax5 = nexttile;
histogram(modIdx_025mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
title('0.25mW Modulation Index');

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0';
saveas(gcf,strcat(sessionpath,'\Pre_Post_20ms\0.25mW','.png'));
close; 

%separate out pre and post laser stimulation 0.5mW
optoSpikeShort20_050mw_pre = optoSpikeRateShort20_050mw(:, :, 1:200);
optoSpikeShort20_050mw_post = optoSpikeRateShort20_050mw(:, :, 201:400);
[significance_050mW, pvalue_050mW] = ranksumSignificance(optoSpikeShort20_050mw_pre,optoSpikeShort20_050mw_post,1:ap.nGoodClusters);
modIdx_050mW = getModulationIdx(optoSpikeShort20_050mw_pre,optoSpikeShort20_050mw_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_050mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2);
xlabel('p values'); ylabel('Count'); box off;
title('0.5mW P values');

ax2 = nexttile;
less_than_05 = sum(pvalue_050mW < 0.05); 
num_non_p_less_05 = length(ap.goodClusters) - less_than_05;
percent_mod_pvalue = [less_than_05, num_non_p_less_05];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod_pvalue)
% lgd = legend(labels);
% lgd.Layout.Tile = 'south';

ax3 = nexttile;
histogram(significance_050mW,10);
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off;
title('0.5mW Pre/Post Signficance');

ax4 = nexttile;
num_sign_050 = length(find(significance_050mW));
num_nonsign050 = length(ap.goodClusters) - num_sign_050;
percent_mod = [num_sign_050, num_nonsign050];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax5 = nexttile;
histogram(modIdx_050mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
title('0.50mW Modulation Index');

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0';
saveas(gcf,strcat(sessionpath,'\Pre_Post_20ms\0.50mW','.png'));
close; 


%separate out pre and post laser stimulation 0.75mW
optoSpikeShort15_075mw_pre = optoSpikeRateShort20_075mw(:, :, 1:200);
optoSpikeShort15_075mw_post = optoSpikeRateShort20_075mw(:, :, 201:400);
[significance_075mW, pvalue_075mW] = ranksumSignificance(optoSpikeShort15_075mw_pre,optoSpikeShort15_075mw_post,1:ap.nGoodClusters);
modIdx_075mW = getModulationIdx(optoSpikeShort15_075mw_pre,optoSpikeShort15_075mw_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_075mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2);
xlabel('p values'); ylabel('Count'); box off;
title('0.75mW P values');

ax2 = nexttile;
less_than_05 = sum(pvalue_075mW < 0.05); 
num_non_p_less_05 = length(ap.goodClusters) - less_than_05;
percent_mod_pvalue = [less_than_05, num_non_p_less_05];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod_pvalue)
% lgd = legend(labels);
% lgd.Layout.Tile = 'south';

ax3 = nexttile;
histogram(significance_075mW,10);
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off;
title('0.75mW Pre/Post Signficance');

ax4 = nexttile;
num_sign_075 = length(find(significance_075mW));
num_nonsign075 = length(ap.goodClusters) - num_sign_075;
percent_mod = [num_sign_075, num_nonsign075];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax5 = nexttile;
histogram(modIdx_075mW,20);
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
title('0.75mW Modulation Index');

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0';
saveas(gcf,strcat(sessionpath,'\Pre_Post_20ms\0.75mW','.png'));
close;

%separate out pre and post laser stimulation 1mW
optoSpikeShort20_1mW_pre = optoSpikeRateShort20_1mw(:, :, 1:200);
optoSpikeShort20_1mW_post = optoSpikeRateShort20_1mw(:, :, 201:400);
[significance_1mW, pvalue_1mW] = ranksumSignificance(optoSpikeShort20_1mW_pre,optoSpikeShort20_1mW_post,1:ap.nGoodClusters);
modIdx_1mW = getModulationIdx(optoSpikeShort20_1mW_pre,optoSpikeShort20_1mW_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_1mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2);
xlabel('p values'); ylabel('Count'); box off;
title('1mW P values')

ax2 = nexttile;
less_than_05_1mw = sum(pvalue_1mW < 0.05); 
num_non_p_less_05_1mW = length(ap.goodClusters) - less_than_05_1mw;
percent_mod_pvalue = [less_than_05_1mw, num_non_p_less_05_1mW];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod_pvalue)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax3 = nexttile;
histogram(significance_1mW,10);
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off;
title('1mW Pre/Post Signficance');

ax4 = nexttile;
num_sign_1 = length(find(significance_1mW));
num_nonsign1 = length(ap.goodClusters) - num_sign_1;
percent_mod = [num_sign_1, num_nonsign1];
labels = {'significantly modulated', 'not significantly modulated'};

pie(percent_mod)
%lgd = legend(labels);
%lgd.Layout.Tile = 'south';

ax5 = nexttile;
histogram(modIdx_1mW,20);
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
title('1mW Modulation Index');

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0';
saveas(gcf,strcat(sessionpath,'\Pre_Post_20ms\1mW','.png'));
close;


%Modulated Units
%using significance
sign_025_units = length(find(significance_025mW)); 
sign_1_units = length(find(significance_1mW)); 
same_mod_units = nnz((ismember(find(significance_025mW), find(significance_1mW))))
modulated_units = [sign_025_units, sign_1_units, same_mod_units]

initializeFig(1,1);
t = tiledlayout(2,2,'TileSpacing','compact');
ax1 = nexttile;
X = categorical({'# Sig Mod 0.25mW','# Sig Mod 1mW','# Overlap Units'});
X = reordercats(X,{'# Sig Mod 0.25mW','# Sig Mod 1mW','# Overlap Units'});
Y = modulated_units;
bar(X,Y)
ylabel("# of units")

ax2 = nexttile;
same_mod = [same_mod_units, (length(ap.goodClusters) - same_mod_units)];
pie(same_mod)
labels = {'% overlap units modulated at 0.25 and 1mW', '% units not modulated at 0.25 and 1mW'};
lgd = legend(labels);
lgd.Layout.Tile = 'east';
title('% overlap units (Bonfon Corrected)');

%using pvalue
sign_025_units_pvalue = length(find(pvalue_025mW < 0.05)); 
sign_1_units_pvalue = length(find(pvalue_1mW < 0.05)); 
%same_mod_units_p = nnz((ismember(find(pvalue_025mW < 0.05), find(pvalue_1mW < 0.05))))
same_mod_units_p = length(intersect(find(pvalue_025mW < 0.05), find(pvalue_1mW < 0.05)));
modulated_units = [sign_025_units_pvalue, sign_1_units_pvalue, same_mod_units_p];

%
intersect_units = intersect(find(pvalue_025mW < 0.05), find(pvalue_1mW < 0.05));

ax3 = nexttile;
X = categorical({'# Sig Mod p 0.25mW','# Sig Mod p 1mW','# Overlap Units p'});
X = reordercats(X,{'# Sig Mod p 0.25mW','# Sig Mod p 1mW','# Overlap Units p'});
Y = modulated_units;
bar(X,Y)
ylabel("# of units");

ax4 = nexttile;
same_mod = [same_mod_units_p, (length(ap.goodClusters) - same_mod_units_p)];
pie(same_mod)
labels = {'% overlap units modulated at 0.25 and 1mW (pvalue)', '% units not modulated at 0.25 and 1mW (pvalue)'};
lgd = legend(labels);
lgd.Layout.Tile = 'east';
title('% overlap units (pvalue)');

saveas(gcf,strcat(sessionpath,'\Pre_Post_20ms\# Same Modulated Units','.png'));
close;

%%

%Looking now at the same units that were modulated at 0.25mW and 1mW 
%plot a subset of the data around 20ms

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0';

optoStimDuration = 0; % Ally

clusterList = intersect_units;  

for idx = 1:length(clusterList);
    i = clusterList(idx);
    initializeFig(1,1); textOn = true;
    t = tiledlayout(2,2,'TileSpacing','compact');

    % plot raster
    ax1 = nexttile;
    plotSpikeRaster(optoSpikeRateShort20_025mw,i,[-0.02,0.02], unit='ms',originalTimeRange=timeRangeShort20,size=5);
    plotEvent('Stim',optoStimDuration,[7 136 225],unit='ms');
    yline([100],'--');
    title('0.25mW of Cluster', num2str(i));
    
    ax2 = nexttile;    
    modIdx_025mW = getModulationIdx(optoSpikeShort20_025mw_pre,optoSpikeShort20_025mw_post,i);
    histogram(modIdx_025mW,20)
    xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
    title('025mW Modulation Index');
    xlim([-1, 1]);

    ax3 = nexttile;      
    plotSpikeRaster(optoSpikeRateShort20_1mw,i,[-0.02,0.02],unit='ms',originalTimeRange=timeRangeShort20,size=5);
    plotEvent('Stim',optoStimDuration,[7 136 225],unit='ms');
    yline([100],'--');
    %yline([600 1200 1800 2400],'-');
    title('1mW of Cluster', num2str(i))
    
    ax4 = nexttile;   
    modIdx_1mW = getModulationIdx(optoSpikeShort20_1mW_pre,optoSpikeShort20_1mW_post,i);
    histogram(modIdx_1mW,20)
    xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
    title('1mW Modulation Index');
    xlim([-1, 1])

    % save figure
    saveas(gcf,strcat(sessionpath,'\Test\PSTH_Same_Mod_units',num2str(i),'.png'));
    close;
end 

modIdx_025mW = getModulationIdx(optoSpikeShort20_025mw_pre,optoSpikeShort20_025mw_post,clusterList);
histogram(modIdx_025mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
title('0.25mW Modulation Index Same Mod Units');
xlim([-1, 1])
saveas(gcf,strcat(sessionpath,'\Same_mod_units\0.25mW Modulation Index_sameModUnits','.png'));
close;

modIdx_1mW = getModulationIdx(optoSpikeShort20_1mW_pre,optoSpikeShort20_1mW_post,clusterList);
histogram(modIdx_1mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off;
title('1mW Modulation Index Same Mod Units');
xlim([-1, 1])
saveas(gcf,strcat(sessionpath,'\Same_mod_units\1mW Modulation Index_sameModUnits','.png'));
close;

%% Ally 15 ms 

% For ally
% optoStimIdx = find(blueLaser); % Find blue laser pulses
% % Define params
% timeRangeShort15 = [-0.015, 0.015];
% binSizeShort = 0.0001;
% [optoSpikeShort15,optoSpikeRateShort15,opto_short_params15] = getSpikes(timeRangeShort15,binSizeShort,optoStimIdx,params);

%separate the spikes out by power 0.25-1mW
optoSpikeRateShort15_025mw = optoSpikeRateShort15(:, 1:100, :);
optoSpikeRateShort15_050mw = optoSpikeRateShort15(:, 101:200, :);
optoSpikeRateShort15_075mw = optoSpikeRateShort15(:, 201:300, :);
optoSpikeRateShort15_1mw = optoSpikeRateShort15(:, 301:400, :);

%separate out pre and post laser stimulation 0.25mW
optoSpikeShort15_025mw_pre = optoSpikeRateShort15_025mw(:, :, 1:150);
optoSpikeShort15_025mw_post = optoSpikeRateShort15_025mw(:, :, 151:300);
[significance_025mW, pvalue_025mW] = ranksumSignificance(optoSpikeShort15_025mw_pre,optoSpikeShort15_025mw_post,1:ap.nGoodClusters);
modIdx_025mW = getModulationIdx(optoSpikeShort15_025mw_pre,optoSpikeShort15_025mw_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_025mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2)
xlabel('p values'); ylabel('Count'); box off
title('0.25mW P values')

ax2 = nexttile;
less_than_05_025mW = sum(pvalue_025mW < 0.05); 
num_non_p_less_05_025mW = length(ap.goodClusters) - less_than_05_025mW;
percent_mod_pvalue = [less_than_05_025mW, num_non_p_less_05_025mW];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod_pvalue)
% lgd = legend(labels);
% lgd.Layout.Tile = 'south';

ax3 = nexttile;
histogram(significance_025mW,10)
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off
title('0.25mW Pre/Post Signficance')

ax4 = nexttile;
num_sign_025 = length(find(significance_025mW));
num_nonsign025 = length(ap.goodClusters) - num_sign_025;
percent_mod = [num_sign_025, num_nonsign025];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax5 = nexttile;
histogram(modIdx_025mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off
title('0.25mW Modulation Index')

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0'
saveas(gcf,strcat(sessionpath,'\Pre_Post_15ms\0.25mW','.png'));
close; 

%separate out pre and post laser stimulation 0.5mW
optoSpikeShort15_050mw_pre = optoSpikeRateShort15_050mw(:, :, 1:150);
optoSpikeShort15_050mw_post = optoSpikeRateShort15_050mw(:, :, 151:300);
[significance_050mW, pvalue_050mW] = ranksumSignificance(optoSpikeShort15_050mw_pre,optoSpikeShort15_050mw_post,1:ap.nGoodClusters);
modIdx_050mW = getModulationIdx(optoSpikeShort15_050mw_pre,optoSpikeShort15_050mw_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_050mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2)
xlabel('p values'); ylabel('Count'); box off
title('0.5mW P values')

ax2 = nexttile;
less_than_05_050mW = sum(pvalue_050mW < 0.05); 
num_non_p_less_05_05mW = length(ap.goodClusters) - less_than_05_050mW;
percent_mod_pvalue = [less_than_05_050mW, num_non_p_less_05_05mW];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod_pvalue)
% lgd = legend(labels);
% lgd.Layout.Tile = 'south';

ax3 = nexttile;
histogram(significance_050mW,10)
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off
title('0.5mW Pre/Post Signficance')

ax4 = nexttile;
num_sign_050 = length(find(significance_050mW));
num_nonsign05 = length(ap.goodClusters) - num_sign_050;
percent_mod = [num_sign_050, num_nonsign05];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax5 = nexttile;
histogram(modIdx_050mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off
title('0.5mW Modulation Index')

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0'
saveas(gcf,strcat(sessionpath,'\Pre_Post_15ms\0.5mW','.png'));
close;

%separate out pre and post laser stimulation 0.75mW
optoSpikeShort15_075mw_pre = optoSpikeRateShort15_075mw(:, :, 1:150);
optoSpikeShort15_075mw_post = optoSpikeRateShort15_075mw(:, :, 151:300);
[significance_075mW, pvalue_075mW] = ranksumSignificance(optoSpikeShort15_075mw_pre,optoSpikeShort15_075mw_post,1:ap.nGoodClusters);
modIdx_075mW = getModulationIdx(optoSpikeShort15_075mw_pre,optoSpikeShort15_075mw_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_075mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2)
xlabel('p values'); ylabel('Count'); box off
title('0.75mW P values')

ax2 = nexttile;
less_than_05_075mW = sum(pvalue_075mW < 0.05); 
num_non_p_less_05_075mW = length(ap.goodClusters) - less_than_05_075mW;
percent_mod_pvalue = [less_than_05_075mW, num_non_p_less_05_075mW];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod_pvalue)
% lgd = legend(labels);
% lgd.Layout.Tile = 'south';

ax3 = nexttile;
histogram(significance_075mW,10)
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off
title('0.75mW Pre/Post Signficance')

ax4 = nexttile;
num_sign_075 = length(find(significance_075mW));
num_nonsign075 = length(ap.goodClusters) - num_sign_075;
percent_mod = [num_sign_075, num_nonsign075];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax5 = nexttile;
histogram(modIdx_075mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off
title('0.75mW Modulation Index')

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0'
saveas(gcf,strcat(sessionpath,'\Pre_Post_15ms\0.75mW','.png'));
close;

%separate out pre and post laser stimulation 1mW
optoSpikeShort15_1mW_pre = optoSpikeRateShort15_1mw(:, :, 1:150);
optoSpikeShort15_1mW_post = optoSpikeRateShort15_1mw(:, :, 151:300);
[significance_1mW, pvalue_1mW] = ranksumSignificance(optoSpikeShort15_1mW_pre,optoSpikeShort15_1mW_post,1:ap.nGoodClusters);
modIdx_1mW = getModulationIdx(optoSpikeShort15_1mW_pre,optoSpikeShort15_1mW_post,1:ap.nGoodClusters);

initializeFig(1,0.5);
t = tiledlayout(1,5,'TileSpacing','compact');

ax1 = nexttile;
histogram(pvalue_1mW, 50); hold on;
xline(0.05, 'r', '-', 'LineWidth', 2)
xlabel('p values'); ylabel('Count'); box off
title('1mW P values')

ax2 = nexttile;
less_than_05_1mW = sum(pvalue_1mW < 0.05); 
num_non_p_less_05_1mW = length(ap.goodClusters) - less_than_05_1mW;
percent_mod_pvalue = [less_than_05_1mW, num_non_p_less_05_1mW];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod_pvalue)
% lgd = legend(labels);
% lgd.Layout.Tile = 'south';

ax3 = nexttile;
histogram(significance_1mW,10)
xlabel('Significance (0 or 1), Bonferroni Correct'); ylabel('Count'); box off
title('1mW Pre/Post Signficance')

ax4 = nexttile;
num_sign_1 = length(find(significance_1mW));
num_nonsign1 = length(ap.goodClusters) - num_sign_1;
percent_mod = [num_sign_1, num_nonsign1];
labels = {'significantly modulated', 'not significantly modulated'}

pie(percent_mod)
lgd = legend(labels);
lgd.Layout.Tile = 'east';

ax5 = nexttile;
histogram(modIdx_1mW,20)
xlabel('Modulation Index (Fon - Foff / Fon + Foff)'); ylabel('Count'); box off
title('1mW Modulation Index')

sessionpath = '/Users/allygirasole/Desktop/Ephys/2023_01_07_MAA103/2023_01_07_MAA103_g0'
saveas(gcf,strcat(sessionpath,'\Pre_Post_15ms\1mW','.png'));
close;

%Modulated Units
sign_025_units = length(find(significance_025mW)); 
sign_1_units = length(find(significance_1mW)); 
same_mod_units = nnz((ismember(find(significance_025mW), find(significance_1mW))))
modulated_units = [sign_025_units, sign_1_units, same_mod_units]

initializeFig(1,1);
t = tiledlayout(2,2,'TileSpacing','compact');

ax1 = nexttile;
X = categorical({'# Sig Mod 0.25mW','# Sig Mod 1mW','# Overlap Units'});
X = reordercats(X,{'# Sig Mod 0.25mW','# Sig Mod 1mW','# Overlap Units'});
Y = modulated_units;
bar(X,Y)
ylabel("# of units")

ax2 = nexttile;
same_mod = [same_mod_units, (length(ap.goodClusters) - same_mod_units)]
pie(same_mod)
labels = {'% overlap units modulated at 0.25 and 1mW', '% units not modulated at 0.25 and 1mW'}
lgd = legend(labels);
lgd.Layout.Tile = 'east';
title('% overlap units (Bonfon Corrected)');

%
sign_025_units_pvalue = length(find(pvalue_025mW < 0.05)); 
sign_1_units_pvalue = length(find(pvalue_1mW < 0.05)); 
same_mod_units_p = nnz((ismember(find(pvalue_025mW < 0.05), find(pvalue_1mW < 0.05))))
modulated_units = [sign_025_units_pvalue, sign_1_units_pvalue, same_mod_units_p]

ax3 = nexttile;
X = categorical({'# Sig Mod p 0.25mW','# Sig Mod p 1mW','# Overlap Units p'});
X = reordercats(X,{'# Sig Mod p 0.25mW','# Sig Mod p 1mW','# Overlap Units p'});
Y = modulated_units;
bar(X,Y)
ylabel("# of units")

ax4 = nexttile;
same_mod = [same_mod_units_p, (length(ap.goodClusters) - same_mod_units_p)]
pie(same_mod)
labels = {'% overlap units modulated at 0.25 and 1mW (pvalue)', '% units not modulated at 0.25 and 1mW (pvalue)'}
lgd = legend(labels);
lgd.Layout.Tile = 'east';
title('% overlap units (pvalue)');

saveas(gcf,strcat(sessionpath,'\Pre_Post_15ms\# Same Modulated Units','.png'));
close;


%% Plot average of all neurons

initializeFig(0.5,0.5);

subplot(2,1,1)
plotSpikes(optoSpikeRateShort(:,1:400,:),1:ap.nGoodClusters,timeRangeShort,bluePurpleRed,average=true,unit='ms',smooth=true);
plotEvent('Stim',0.02,[7 136 225],unit='ms');
xlabel('Time (s)'); ylabel('Spikes/s'); box off

subplot(2,1,2)
plotSpikes(optoSpikeRateLong(:,1:400,:),1:ap.nGoodClusters,timeRangeLong,bluePurpleRed,average=true,unit='ms');
plotEvent('Stim',0.02,[7 136 225],unit='ms');
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

% optoStimDuration = 0.02; % Shun
optoStimDuration = 0; % Ally

wvf_r = nan(1,ap.nGoodClusters);
for idx = 1:542%ap.nGoodClusters
    i = idx;
    initializeFig(1,1); textOn = true;
    tiledlayout(4,3);

    % plot raster
    nexttile(1,[1 2]);
    plotSpikeRaster(optoSpikeRateShort(:,401:650,:),i,[-0.005,salt_window],unit='ms',originalTimeRange=timeRangeShort,size=5);
    plotEvent('Stim',optoStimDuration,[7 136 225],unit='ms');
    yline([250],'--');
    %yline([600 1200 1800 2400],'-');
    
    nexttile(4,[1 2]);
    plotSpikeRaster(optoSpikeRateShort(:,401:650,:),i,timeRangeShort,unit='ms',size=5);
    plotEvent('Stim',optoStimDuration,[7 136 225],unit='ms');
    yline([250],'--');
    %yline([600 1200 1800 2400],'-');
    
    nexttile(7,[1 2]);
    plotSpikes(optoSpikeRateShort(:,401:650,:),i,timeRangeShort,bluePurpleRed,textOn=true,text_source=ap,unit='ms',smooth=true);
    plotEvent('Stim',optoStimDuration,[7 136 225],unit='ms');
    
    nexttile(10,[1 2]);
    plotSpikes(optoSpikeRateLong(:,401:650,:),i,timeRangeLong,bluePurpleRed,textOn=true,text_source=ap,unit='ms');
    plotEvent('Stim',optoStimDuration,[7 136 225],unit='ms');


    % Get waveform
    timeRange = [-0.0025,0.0025];
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
    saveas(gcf,strcat(sessionpath,'\PSTH-allOpto_1mW_pulse\psth_stim_cluster',num2str(i),'.png'));
    close;
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
