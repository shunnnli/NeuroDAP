%% Housekeeping stuff
% Code based on Bernardo's processSJData.m and SJ's getWaveForms_filtered.m
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% Contributed by C. Schoonover and A. Fink
% WaveForm extraction from temp_wh.data (filtered data from kilosort)

clear all;
close all;

% Load npy-matlab
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

%% Load NPY data to matlab

% Enter session name
sessionName = '20220330 SJ495-R_g0';
session.pathImec = strcat('/Volumes/Shun neuro data/Neuropixel/',sessionName, '/', sessionName, '_imec0/');
session.apBin = strcat(sessionName,'_t0.imec0.ap.bin');
session.lfpBin = strcat(sessionName,'_t0.imec0.lf.bin');
% myEventTimes = load('path'); % a vector of times in seconds of some event to align to

ap.meta = ReadMeta(session.apBin, session.pathImec);
ap.Fs=str2double(ap.meta.imSampRate); % LFP: 2500, AP: 30000
disp(['AP Sampling rate is ' num2str(ap.Fs)]);
lfp.meta = ReadMeta(session.lfpBin, session.pathImec);
lfp.Fs=str2double(lfp.meta.imSampRate); % LFP: 2500, AP: 30000
disp(['LFP Sampling rate is ' num2str(lfp.Fs)]);

[clusterLabel,spike_times,spike_clusters] = readNPYData(session.pathImec);

%% Ask to get .bin file
[binNameImec,pathImec] = uigetfile('*.bin', 'Select Binary File');
if isempty(binNameImec)
    error('pick a file')
end
disp(['Will read from ' binNameImec])

meta = ReadMeta(binNameImec, pathImec);
ap.Fs=str2double(meta.imSampRate); % LFP: 2500, AP: 30000
disp(['AP Sampling rate is ' num2str(ap.Fs)]);

lfpBinName = strcat(sessionName,'_t0.imec0.lf.bin');
metaLFP = ReadMeta(lfpBinName, pathImec);
lfp.Fs=str2double(metaLFP.imSampRate); % LFP: 2500, AP: 30000
disp(['LFP Sampling rate is ' num2str(lfp.Fs)]);

%% Setup data structure

gwfparams.fileName = 'temp_wh.dat';     % .dat file containing the raw
gwfparams.dataType = 'int16';           % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 383;                    % Number of channels that were streamed to disk in .dat file
gwfparams.nWf = 383;                    % Number of waveforms per unit to pull out

gwfparams.spikeTimes = spike_times;         % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = spike_clusters;   % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

gwfparams.keepTime = 0.0040;            % Time (sec) before and after spiketime to include in waveform
gwfparams.blockTime = 60;               % Time for analyzing

% AP specific data
ap.keepSamples=floor(gwfparams.keepTime*ap.Fs);   
ap.keepRange = -ap.keepSamples:ap.keepSamples;        % Range of samples around spiketime to include in waveform
ap.nSampPerBlock = floor(gwfparams.blockTime * ap.Fs);
ap.nSampInFile = str2double(ap.meta.fileSizeBytes) / (2*str2double(ap.meta.nSavedChans));
ap.nBlocks = 1; % Number of blocks used for analysis
ap.totalBlocks=floor(ap.nSampInFile/ap.nSampPerBlock);  % Total number of blocks

% LFP specific data
lfp.keepSamples=floor(gwfparams.keepTime*lfp.Fs);
lfp.keepRange = -lfp.keepSamples:lfp.keepSamples; % Range of samples around spiketime to include in waveform
lfp.nSampPerBlock = floor(gwfparams.blockTime * lfp.Fs);
lfp.nSampInFile = str2double(lfp.meta.fileSizeBytes) / (2*str2double(lfp.meta.nSavedChans));
lfp.nBlocks = 1;    % Number of blocks used for analysis
lfp.totalBlocks=floor(lfp.nSampInFile/lfp.nSampPerBlock);   % Total number of blocks

% Spike specific data
% KSLabel==good but group==unsorted means good but unstable units
ap.goodClusters=clusterLabel.cluster_id(find(clusterLabel.group=='good')); % Good units labeled by hand
ap.nGoodClusters=length(ap.goodClusters);
ap.clusterToGoodClusterIndex=zeros(max(ap.goodClusters), 1);  % Index column is cluster_id, first column is goodClusters id -> converst cluster_id to goodClusters id
for counter=1:length(ap.goodClusters)
    ap.clusterToGoodClusterIndex(ap.goodClusters(counter))=counter;
end

% Create separate array for cluster_id of good cluster spikes
goodClusterSpikeIndices=find(ismember(spike_clusters, ap.goodClusters));
% Spike times in samples (convert to seconds by dividing sample rate)
ap.goodSpikeTimes=spike_times(goodClusterSpikeIndices);
% Spike cluster_id in the order of spike occurance
ap.goodSpikeClusters=spike_clusters(goodClusterSpikeIndices);

%% Load filtered waveform from temp_wh.dat and KiloSort/Phy output
fileName = fullfile(session.pathImec,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh * dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});
chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1; % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

%% Read filtered AP spike time-centered waveforms
gwfparams.nWf = 383;            % Number of waveforms per unit to pull out
nBlocks = 3;
nSpikesPerCluster=zeros(ap.nGoodClusters,1);    % Number of spikes for each good cluster

% waveform = (#channels, activity around spiketime, #clusters_id)
% Put cluster_id at the back so we can draw image using waveform(:,:,cluster_id)
waveforms = zeros(gwfparams.nWf,length(ap.keepRange),ap.nGoodClusters);
commonSpikeTime = ap.goodSpikeTimes;
nSampPerBlock = ap.nSampPerBlock;
keepSamples = ap.keepSamples;
keepRange = ap.keepRange;
% goodSpikeClusters = ap.goodSpikeClusters;
% clusterToGoodClusterIndex = ap.clusterToGoodClusterIndex;
% nSpikesPerCluster = ap.nSpikesPerCluster;

figure
for counter=1:nBlocks
    disp(['Reading block ' num2str(counter) '/' num2str(nBlocks)]);

    startIndex=(counter-1)*nSampPerBlock; % First index of each block (used as offset)
    endIndex=startIndex+nSampPerBlock-1;  % Last index of each block
    
    blockWaveform = mmf.Data.x(1:gwfparams.nCh, startIndex+1:endIndex);
    blockSpikes=find((commonSpikeTime>(startIndex+keepSamples)) ...
                        & (commonSpikeTime<(endIndex-keepSamples)));
    
    if ~isempty(blockSpikes)
        disp(['   Found ' num2str(length(blockSpikes)) ' spikes']);
        
        for spikeCounter=1:length(blockSpikes)
            % Extract cluster_id for each spikes found
            spikeIndex=blockSpikes(spikeCounter);
            cluster=double(ap.goodSpikeClusters(spikeIndex));    % cluster_id
            clusterIndex=ap.clusterToGoodClusterIndex(cluster);  % goodcluster id
            
            % Relative spike time aligned to block start
            spikeTime=double(commonSpikeTime(spikeIndex));
            tempStart=spikeTime-startIndex;

            newWf=double(blockWaveform(:, tempStart+keepRange));
            newWf = newWf - mean(newWf,2); % raw V - mean of channel
            % [~,maxCh] = maxk(max(abs(newWf),[],2),gwfparams.nWf); % Select nWf channels with highest magnitude of Ve
            % newWfMax = newWf(maxCh,:);
            newWfMax = newWf;

            % Spike waveform for each cluster
            waveforms(:,:,clusterIndex)= ...
                (nSpikesPerCluster(clusterIndex)*waveforms(:,:,clusterIndex) + newWfMax)...
                / (nSpikesPerCluster(clusterIndex)+1);
            nSpikesPerCluster(clusterIndex) = nSpikesPerCluster(clusterIndex)+1;
        end
    end
    
    % Showing example cluster_id = 40
    % imagesc(waveforms(:,:,13));
    % colorbar;
    plot(keepRange,waveforms(3,:,1));
    xlabel('Time');
    ylabel('V'); 
    title(['Channel 3 V traces of a cluster']);
    drawnow
end

% Save waveforms and nBlocks
if AP_on 
    ap.waveforms = waveforms;
    ap.nSpikesPerCluster = nSpikesPerCluster;
    ap.nBlocks = nBlocks; lfp.nBlocks = nBlocks;
    disp(['Finished reading ', num2str(nBlocks), ' blocks']);
end

%% Read raw AP or LFP spike time-centered waveforms
gwfparams.nWf = 383;    % Number of waveforms per unit to pull out
nBlocks = 3;
AP_on = false;          % true->ap.bin; false->lfp.bin
nSpikesPerCluster=zeros(ap.nGoodClusters,1);    % Number of spikes for each good cluster

% waveform = (#channels, activity around spiketime, #clusters_id)
% Put cluster_id at the back so we can draw image using waveform(:,:,cluster_id)
if AP_on
    waveforms = zeros(gwfparams.nWf,length(ap.keepRange),ap.nGoodClusters);
    commonSpikeTime = ap.goodSpikeTimes;
    meta = ap.meta;
    binNameImec = session.apBin;
    nSampPerBlock = ap.nSampPerBlock;
    keepSamples = ap.keepSamples;
    keepRange = ap.keepRange;
else 
    waveforms = zeros(gwfparams.nWf,length(lfp.keepRange),ap.nGoodClusters);
    spikeToSample = lfp.Fs/30000;
    commonSpikeTime = floor(ap.goodSpikeTimes*spikeToSample);  % Sample step of each good spike
    meta = lfp.meta;
    binNameImec = session.lfpBin;
    nSampPerBlock = lfp.nSampPerBlock;
    keepSamples = lfp.keepSamples;
    keepRange = lfp.keepRange;
end

figure
for counter=1:nBlocks
    disp(['Reading block ' num2str(counter) '/' num2str(nBlocks)]);

    startIndex=(counter-1)*nSampPerBlock; % First index of each block (used as offset)
    endIndex=startIndex+nSampPerBlock-1;  % Last index of each block

    blockSpikes=find((commonSpikeTime>(startIndex+keepSamples)) ...
                        & (commonSpikeTime<(endIndex-keepSamples)));
                    
    raw = ReadBin(startIndex, nSampPerBlock, meta, binNameImec, session.pathImec);
    raw([192 193],:) = [];
    
    if ~isempty(blockSpikes)
        disp(['   Found ' num2str(length(blockSpikes)) ' spikes']);
        
        for spikeCounter=1:length(blockSpikes)
            % Extract cluster_id for each spikes found
            spikeIndex=blockSpikes(spikeCounter);
            cluster=double(ap.goodSpikeClusters(spikeIndex));    % cluster_id
            clusterIndex=ap.clusterToGoodClusterIndex(cluster);  % goodcluster id
            
            % Relative spike time aligned to block start
            spikeTime=double(commonSpikeTime(spikeIndex));
            tempStart=spikeTime-startIndex;

            % newWf=double(blockWaveform(:, tempStart+keepRange));
            newWf=double(raw(:, tempStart+keepRange));
            newWf = newWf - mean(newWf,2); % raw V - mean of channel
            % [~,maxCh] = maxk(max(abs(newWf),[],2),gwfparams.nWf); % Select nWf channels with highest magnitude of Ve
            % newWfMax = newWf(maxCh,:);
            newWfMax = newWf;

            % Spike waveform for each cluster
            waveforms(:,:,clusterIndex)= ...
                (nSpikesPerCluster(clusterIndex)*waveforms(:,:,clusterIndex) + newWfMax)...
                / (nSpikesPerCluster(clusterIndex)+1);
            nSpikesPerCluster(clusterIndex) = nSpikesPerCluster(clusterIndex)+1;
        end
    end
    
    % Showing example cluster_id = 40
    % imagesc(waveforms(:,:,13));
    % colorbar;
    plot(keepRange,waveforms(3,:,1));
    xlabel('Time');
    ylabel('V'); 
    title(['Channel 3 V traces of a cluster']);
    drawnow
end

% Save waveforms and nBlocks
if AP_on 
    ap.waveforms = waveforms;
    ap.nSpikesPerCluster = nSpikesPerCluster;
    ap.nBlocks = nBlocks; lfp.nBlocks = nBlocks;
    disp(['Finished reading ', num2str(nBlocks), ' blocks']);
else
    lfp.waveforms = waveforms;
    % lfp.nSpikesPerCluster = nSpikesPerCluster;
    ap.nBlocks = nBlocks; lfp.nBlocks = nBlocks;
    disp(['Finished reading ', num2str(nBlocks), ' blocks']);
end

%% Plot related waveform
figure
imagesc(newWf(:,:,1));
colorbar;
xlabel('Sample');
ylabel('Channel number'); 

figure
plot(ap.keepRange,waveforms(3,:,1));
xlabel('Sample');
ylabel('Ve'); 

%% Create spike summary and calculate spike properties

spikeTable = table;
spikeTable.cluster_id = ap.goodClusters;              % cluster_id
spikeTable.goodClusterIndex=(1:length(ap.nSpikesPerCluster))';   % goodcluster id
spikeTable.firingRate = ap.nSpikesPerCluster/(ap.nBlocks*60);     % firing rate

spikeTable.contactChannel = 0 * ap.nSpikesPerCluster;     % channel number that detects the smallest V during spike range
% spikeTable.spikeAmp=spikeTable.spikeContact;        % Amp = max(V) - min(V) during spike
% spikeTable.spikeMinTime=spikeTable.spikeContact;    % Sample step where potential reaches minimum
% spikeTable.spikeNegAmp=spikeTable.spikeContact;     % min(V) during spike


for counter=1:length(ap.nSpikesPerCluster)

    % Find spikeContact
    % mmm: min V detected across all channels at each sample step
    % mml: channel number that detects min V at each sample step
    [mmm,mml]=min(ap.waveforms(:,:,counter)); 
    % mmm2: min V detected across all channels of all sample step
    % mml2: sample step that mmm2 is at
    [mmm2, mml2]=min(mmm);
    spikeTable.contactChannel(counter)=mml(mml2);
    spikeTable.minVoltageTime(counter)=mml2;
    
    % Waveform -> voltage trace of the spike contact during spikes
    spikeWaveform=ap.waveforms(spikeTable.contactChannel(counter), :, counter);
    spikeWaveform=spikeWaveform-mean(spikeWaveform);
    % spikeTable.spikeWaveform(counter) = spikeWaveform;
    
    % t2pAmp -> trough-to-peak amplitude
    % t2pDuration -> trough-to-peak duration
    % Amp -> signed difference between either the peak or trough and the baseline
    [minpeak,minpeakloc]=min(spikeWaveform);
    [maxpeak,maxpeakloc]=max(spikeWaveform);
    spikeTable.troughAmp(counter) = minpeak;
    spikeTable.peakAmp(counter) = maxpeak;
    spikeTable.t2pAmp(counter) = max(spikeWaveform)-min(spikeWaveform); % pre: spikeAmp
    spikeTable.t2pDuration(counter) = maxpeakloc - minpeakloc;
    
    [firstpeak,firstpeakloc] = max(spikeWaveform(1:minpeakloc));
    [secondpeak,secondpeakloc] = max(spikeWaveform(minpeakloc:end));
    spikeTable.prePeak2TroughRatio(counter) = abs(firstpeak)/abs(minpeak); % pre: spikePrePeak
    spikeTable.postPeak2TroughRatio(counter) = abs(secondpeak)/abs(minpeak); % pre: spikePostPeak
    spikeTable.peak2peakTime(counter) = (secondpeakloc+minpeakloc-firstpeakloc)/ap.Fs;
    % spikeTable.spikePrePeak(counter)=-max(spikeWaveform(1:minpeakloc))/minpeak;
    % spikeTable.spikePostPeak(counter)=-max(spikeWaveform(minpeakloc:end))/minpeak;
    
    % Width, minToMaxWidth
    fff=find(spikeWaveform<minpeak/3);
    if isempty(fff)
        spikeTable.Width(counter) = NaN;
        spikeTable.minToMaxWidth(counter) = NaN;
    else
        spikeTable.Width(counter)=1000*(max(fff)-min(fff))/ap.Fs;
        [~,postMaxLoc]=max(spikeWaveform(minpeakloc:end));
        spikeTable.minToMaxWidth(counter)=1000*postMaxLoc/ap.Fs;
    end
end

%% Plot waveform along the probe (AP or LFP)

clusterList=1:59;
Fs = lfp.Fs;
waveforms = lfp.waveforms;

for indCounter=1:length(clusterList)
    cid=clusterList(indCounter);

    figure
    set(gcf, 'Position', [88          41         800        1296])

    hold on
    contact=spikeTable.contactChannel(cid);   % Contacting channel number
    nContacts=10;
    nSpots=2*nContacts+1;

    title(['cluster: ' num2str(spikeTable.cluster_id(cid)) ...
            ' Channel contact: ' num2str(contact) ...
            ' Firing rate: ' num2str(spikeTable.firingRate(cid)) ' Hz'])

    baseTraces=4*(-nContacts:nContacts)+contact-mod(contact,4);
    leftOuter=baseTraces+0;
    rightInner=baseTraces+1;
    leftInner=baseTraces+2;
    rightOuter=baseTraces+3;

    plotPlaces=zeros(2*nSpots,4);
    plotPlaces(:,1)=reshape([leftOuter' zeros(nSpots,1)]', 1, 2*nSpots)';
    plotPlaces(:,3)=reshape([rightInner' zeros(nSpots,1)]', 1, 2*nSpots)';
    plotPlaces(:,2)=reshape([zeros(nSpots,1) leftInner']', 1, 2*nSpots)';
    plotPlaces(:,4)=reshape([zeros(nSpots,1) rightOuter']', 1, 2*nSpots)';

    xStep_ap = max(1000*keepRange/Fs)-min(1000*keepRange/Fs)/2;
    yStep_ap = (max(max(squeeze(waveforms(:, :, cid))))- ...
                min(min(squeeze(waveforms(:, :, cid)))))/2;

    for xp=1:4
        for yp=1:2*nSpots
            if plotPlaces(yp,xp)>0 && plotPlaces(yp,xp)<gwfparams.nWf
                % Plot waveform traces
                x = (xp-1)*xStep_ap+1000*keepRange/Fs;
                y = (yp-1)*yStep_ap+waveforms(plotPlaces(yp,xp),:,cid);
                plot(x,y);
                % Plot channel number
                xtxtpos = min(x) + 0.5;
                ytxtpos = y(1) + 2;
                text(xtxtpos,ytxtpos,num2str(plotPlaces(yp,xp)));
            end
        end
    end
    drawnow
    pause
end
% autoArrangeFigures();

%% Plot waveform along the probe (AP and LFP)

clusterList=1;

for indCounter=1:length(clusterList)
    cid=clusterList(indCounter);

    figure
    set(gcf, 'Position', [88          41         800        1296])

    hold on
    contact=spikeTable.contactChannel(cid);   % Contacting channel number
    nContacts=10;
    nSpots=2*nContacts+1;

    title(['cluster: ' num2str(spikeTable.cluster_id(cid)) ...
            ' Channel contact: ' num2str(contact) ...
            ' Firing rate: ' num2str(spikeTable.firingRate(cid)) ' Hz'])

    baseTraces=4*(-nContacts:nContacts)+contact-mod(contact,4);
    leftOuter=baseTraces+0;
    rightInner=baseTraces+1;
    leftInner=baseTraces+2;
    rightOuter=baseTraces+3;

    plotPlaces=zeros(2*nSpots,4);
    plotPlaces(:,1)=reshape([leftOuter' zeros(nSpots,1)]', 1, 2*nSpots)';
    plotPlaces(:,3)=reshape([rightInner' zeros(nSpots,1)]', 1, 2*nSpots)';
    plotPlaces(:,2)=reshape([zeros(nSpots,1) leftInner']', 1, 2*nSpots)';
    plotPlaces(:,4)=reshape([zeros(nSpots,1) rightOuter']', 1, 2*nSpots)';

    xStep_ap = max(1000*ap.keepRange/ap.Fs)-min(1000*ap.keepRange/ap.Fs)/2;
    yStep_ap = (max(max(squeeze(ap.waveforms(:, :, cid))))- ...
                min(min(squeeze(ap.waveforms(:, :, cid)))))/2;
    xStep_lfp = max(1000*lfp.keepRange/lfp.Fs)-min(1000*lfp.keepRange/lfp.Fs)/2;
    yStep_lfp = (max(max(squeeze(lfp.waveforms(:, :, cid))))- ...
                min(min(squeeze(lfp.waveforms(:, :, cid)))))/2;

    for xp=1:4
        for yp=1:2*nSpots
            if plotPlaces(yp,xp)>0 && plotPlaces(yp,xp)<gwfparams.nWf
                disp(plotPlaces(yp,xp));
                % Plot waveform traces
                x_ap = (xp-1)*xStep_ap+1000*ap.keepRange/ap.Fs;
                y_ap = (yp-1)*yStep_ap+ap.waveforms(plotPlaces(yp,xp),:,cid);
                plot(x_ap,y_ap);
                % hold on
                x_lfp = (xp-1)*xStep_lfp+1000*lfp.keepRange/lfp.Fs;
                y_lfp = (yp-1)*yStep_lfp+lfp.waveforms(plotPlaces(yp,xp),:,cid);
                plot(x_lfp,y_lfp);
                % Plot channel number
                xtxtpos = min(x_ap) + 0.5;
                ytxtpos = y_ap(1) + 2;
                text(xtxtpos,ytxtpos,num2str(plotPlaces(yp,xp)));
            end
        end
    end
    drawnow
    % pause
end

%% SJ's code for reading individual spike time centered waveforms
%unitIDs = unique(gwfparams.spikeClusters);
unitIDs = 20;

numUnits = size(unitIDs,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
waveFormsMean = nan(numUnits,nChInMap,wfNSamples);

for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
        tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
    end
    waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
    disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.']);
end


% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;

% wf.cluster_id                     % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                 % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                      % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                  % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                   % nClu: number of different clusters in .spikeClusters
%                                   % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

%% SJ's code for ploting FFT of each channel for the snippet of time
L = size(wf.waveForms,4);
for i=1:10
    for j=1:gwfparams.nWf
        figure(i*2-1); plot(squeeze(wf.waveForms(1,j,i,:))); hold on;
    end
    figure(i*2-1); title(['post-processed data: ch',num2str(i)]); hold off;
    
    Y = squeeze(fft(wf.waveForms(1,1,i,:)));
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    %figure(i*2);
    figure(100);
    loglog(f,P1); hold on;
    title('Single-Sided Amplitude Spectrum of X(t)');
    xlabel('f (Hz)');
    ylabel('|P1(f)|');
end
%autoArrangeFigures();

%% General functions
% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [1 X nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
% IMPORTANT: samp0 and nSamp must be integers.
%
function dataArray = ReadBinByCh(samp0, nSamp, meta, binName, path, ch)

    nChan = str2double(meta.nSavedChans);
    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    samp0 = max(samp0, 0);
    nSamp = floor(min(nSamp, nFileSamp - samp0));
    sizeA = [1, nSamp];
    precision = 'int16=>double';
    skip = 2 * (nChan-1);

    fid = fopen(fullfile(path, binName), 'rb');
    fseek(fid, samp0 * 2 * nChan + 2 * (ch-1), 'bof'); %start reading from the specified channel
    dataArray = fread(fid, sizeA, precision, skip); %skip so that we are not reading other channels
    fclose(fid);
end % ReadBinByCh

% =========================================================
% Return an array [lines X timepoints] of uint8 values for
% a specified set of digital lines.
%
% - dwReq is the one-based index into the saved file of the
%    16-bit word that contains the digital lines of interest.
% - dLineList is a zero-based list of one or more lines/bits
%    to scan from word dwReq.
%
function digArray = ExtractDigital(dataArray, meta, dwReq, dLineList)
    % Get channel index of requested digital word dwReq
    if strcmp(meta.typeThis, 'imec')
        [AP, LF, SY] = ChannelCountsIM(meta);
        if SY == 0
            fprintf('No imec sync channel saved\n');
            digArray = [];
            return;
        else
            digCh = AP + LF + dwReq;
        end
    else
        [MN,MA,XA,DW] = ChannelCountsNI(meta);
        if dwReq > DW
            fprintf('Maximum digital word in file = %d\n', DW);
            digArray = [];
            return;
        else
            digCh = MN + MA + XA + dwReq;
        end
    end
    [~,nSamp] = size(dataArray);
    digArray = zeros(numel(dLineList), nSamp, 'uint8');
    for i = 1:numel(dLineList)
        digArray(i,:) = bitget(dataArray(digCh,:), dLineList(i)+1, 'int16');
    end
end % ExtractDigital


% =========================================================
% Return a multiplicative factor for converting 16-bit
% file data to voltage. This does not take gain into
% account. The full conversion with gain is:
%
%   dataVolts = dataInt * fI2V / gain.
%
% Note that each channel may have its own gain.
%
function fI2V = Int2Volts(meta)
    if strcmp(meta.typeThis, 'imec')
        if isfield(meta,'imMaxInt')
            maxInt = str2num(meta.imMaxInt);
        else
            maxInt = 512;
        end
        fI2V = str2double(meta.imAiRangeMax) / maxInt;
    else
        fI2V = str2double(meta.niAiRangeMax) / 32768;
    end
end % Int2Volts


% =========================================================
% Return array of original channel IDs. As an example,
% suppose we want the imec gain for the ith channel stored
% in the binary data. A gain array can be obtained using
% ChanGainsIM() but we need an original channel index to
% do the look-up. Because you can selectively save channels
% the ith channel in the file isn't necessarily the ith
% acquired channel, so use this function to convert from
% ith stored to original index.
%
% Note: In SpikeGLX channels are 0-based, but MATLAB uses
% 1-based indexing, so we add 1 to the original IDs here.
%
function chans = OriginalChans(meta)
    if strcmp(meta.snsSaveChanSubset, 'all')
        chans = (1:str2double(meta.nSavedChans));
    else
        chans = str2num(meta.snsSaveChanSubset);
        chans = chans + 1;
    end
end % OriginalChans


% =========================================================
% Return counts of each imec channel type that compose
% the timepoints stored in binary file.
%
function [AP,LF,SY] = ChannelCountsIM(meta)
    M = str2num(meta.snsApLfSy);
    AP = M(1);
    LF = M(2);
    SY = M(3);
end % ChannelCountsIM

% =========================================================
% Return counts of each nidq channel type that compose
% the timepoints stored in binary file.
%
function [MN,MA,XA,DW] = ChannelCountsNI(meta)
    M = str2num(meta.snsMnMaXaDw);
    MN = M(1);
    MA = M(2);
    XA = M(3);
    DW = M(4);
end % ChannelCountsNI


% =========================================================
% Return gain for ith channel stored in the nidq file.
%
% ichan is a saved channel index, rather than an original
% (acquired) index.
%
function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
    if ichan <= savedMN
        gain = str2double(meta.niMNGain);
    elseif ichan <= savedMN + savedMA
        gain = str2double(meta.niMAGain);
    else
        gain = 1;
    end
end % ChanGainNI


% =========================================================
% Return gain arrays for imec channels.
%
% Index into these with original (acquired) channel IDs.
%
function [APgain,LFgain] = ChanGainsIM(meta)

    if isfield(meta,'imDatPrb_type')
        probeType = str2num(meta.imDatPrb_type);
    else
        probeType = 0;
    end
    if (probeType == 21) || (probeType == 24)
        [AP,LF,~] = ChannelCountsIM(meta);
        % NP 2.0; APgain = 80 for all channels
        APgain = zeros(AP,1,'double');
        APgain = APgain + 80;
        % No LF channels, set gain = 0
        LFgain = zeros(LF,1,'double');
    else
        % 3A or 3B data?
        % 3A metadata has field "typeEnabled" which was replaced
        % with "typeImEnabled" and "typeNiEnabled" in 3B.
        % The 3B imro table has an additional field for the
        % high pass filter enabled/disabled
        if isfield(meta,'typeEnabled')
            % 3A data
            C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
        else
            % 3B data
            C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
        end
        APgain = double(cell2mat(C(1)));
        LFgain = double(cell2mat(C(2)));
    end
end % ChanGainsIM


% =========================================================
% Having acquired a block of raw nidq data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:MN]      % all MN chans (MN from ChannelCountsNI)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectNI(dataArray, chanList, meta)

    [MN,MA] = ChannelCountsNI(meta);
    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        conv = fI2V / ChanGainNI(j, MN, MA, meta);
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end


% =========================================================
% Having acquired a block of raw imec data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:AP]      % all AP chans (AP from ChannelCountsIM)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectIM(dataArray, chanList, meta)

    % Look up gain with acquired channel ID
    chans = OriginalChans(meta);
    [APgain,LFgain] = ChanGainsIM(meta);
    nAP = length(APgain);
    nNu = nAP * 2;

    % Common conversion factor
    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        k = chans(j);       % acquisition index
        if k <= nAP
            conv = fI2V / APgain(k);
        elseif k <= nNu
            conv = fI2V / LFgain(k - nAP);
        else
            continue;
        end
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end


