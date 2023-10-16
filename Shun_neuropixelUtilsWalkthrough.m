%% Load data
clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
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

% Neuropixels util
session.fullImecPath = strcat(session.pathImec, session.apBin);
imec = Neuropixel.ImecDataset(session.fullImecPath);
fprintf('Duration of recording %s is %g minutes\n', imec.fileStem, imec.nSamplesAP / imec.fsAP / 60);
ks = Neuropixel.KilosortDataset(session.pathImec);
ks.load()

%% Mark bad channels
% imec.markBadChannels(channelIds); % manual marking
[channelInds, channelIds] = imec.lookup_channelIds(channelIds); % use channel ids not indices

imec.markBadChannelsByRMS('rmsRange', [3 100]); % [low high] range of RMS in uV

% Writing modified metadata back to disk
% imec.writeModifiedAPMeta();
% extraMeta.cleaned = true;
% extraMeta.cleaningAlgorithm = 'v1';
% imec.writeModifiedAPMeta('extraMeta', extraMeta);

%% Accessing raw AP data (use mine)

% methods(imec) % see full list of methods

% Read & plot specific time windows
idxWindow = [300000, 300500];
[data_partial, sampleIdx] = imec.readAP_idx(idxWindow);

% timeWindow = [20, 20.01]; % in seconds
% [data_partial, sampleIdx] = imec.readAP_timeWindow(timeWindow);

imec.inspectAP_idxWindow(idxWindow);
% imec.inspectAP_timeWindow(timeWindow);
%% Kilosort analysis
stats = ks.computeBasicStats();
ks.printBasicStats();

metrics = ks.computeMetrics();
% Drift maps
metrics.plotDriftmap();
% Cluster drift maps
cluster_ids = metrics.cluster_ids(1:5:metrics.nClusters);
metrics.plotDriftmap('tsi', tsi, 'exciseRegionsOutsideTrials', true, 'cluster_ids', cluster_ids);
metrics.plotDriftmap('tsi', tsi, 'exciseRegionsOutsideTrials', true, 'cluster_ids', cluster_ids, 'colorByAmp', true);
metrics.plotDriftmap('tsi', tsi, 'exciseRegionsOutsideTrials', true, 'cluster_ids', cluster_ids, ...
    'showSmooth', true, 'showIndividual', false, 'smoothWidthSeconds', 50);

% Plot cluster centroid
figure;
metrics.plotClusterWaveformAtCentroid();

% Plot cluster images
cluster_id = 23; figure;
metrics.plotClusterImage(cluster_id, 'best_n_channels', 20);

%% Extract raw waveforms

% Extracting datasnippets from ImecDataset
% times = 1e5:1e5:1e6;
% window = [-1000, 999]; % 1000 samples before plus 999 samples after
% snippetSet = imec.readAPSnippetSet(times, window);

cluster_id = 23;

% Extracting waveforms via KilosortDataset
snippetSet = ks.getWaveformsFromRawData('cluster_ids', cluster_id, ...
    'num_waveforms', 50, 'best_n_channels', 20, 'car', true);
snippetSet.plotAtProbeLocations('alpha', 0.8);

% snippetSet = ks.getWaveformsFromRawData('cluster_ids', cluster_id, ...
%     'num_waveforms', 50, 'best_n_channels', 20, 'car', true, ...
%     'subtractOtherClusters', true);
% snippetSet.plotAtProbeLocations('alpha', 0.8);