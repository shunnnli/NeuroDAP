%% (Required) Load sync data 
clear; close all;
% addpath(genpath('/Users/shunli/Downloads/Sabatini lab/Methods'))
addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));


% sessionName = '20220330-SJ495-R_g0';
% sessionName = '20220401-SJ495-R_g0';
% sessionName = '20220523-SJ518-L_g0';
sessionName = '20220525-SJ518-L_g0';
load(strcat('sync_',sessionName,'.mat'));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed] = loadColors;

% (Required) Load AP spike data
% session.pathImec = strcat('/Users/shunli/Downloads/Sabatini lab/',sessionName, '/', sessionName, '_imec0');
% session.pathNidq = strcat('/Volumes/T7/',sessionName, '/');
% session.pathImec = strcat('/Volumes/T7/',sessionName, '/', sessionName, '_imec0/');
% session.pathImec = strcat('D:\Shun\',sessionName, '\', sessionName, '_imec0\');
session.pathNidq = strcat('D:\Shun\Recordings\',sessionName, '\');
session.pathImec = strcat('D:\Shun\Recordings\',sessionName, '\catgt_', sessionName, '\');
[clusterLabel,spike_times,spike_clusters] = readNPYData(session.pathImec);

% Spike specific data
% ap = readSpikeData(clusterLabel,spike_times,spike_clusters,ap);
ap.goodClusters=clusterLabel.cluster_id(find(clusterLabel.KSLabel=='good')); % Good units
ap.nGoodClusters=length(ap.goodClusters);
ap.clusterToGoodClusterIndex=zeros(max(ap.goodClusters), 1);  % Index column is cluster_id, first column is goodClusters id -> converst cluster_id to goodClusters id
for counter=1:length(ap.goodClusters)
    ap.clusterToGoodClusterIndex(ap.goodClusters(counter))=counter;
end

% Create separate array for cluster_id of good cluster spikes
goodClusterSpikeIndices=find(ismember(spike_clusters, uint32(ap.goodClusters)));
% Spike times in samples
ap.goodSpikeTimes=spike_times(goodClusterSpikeIndices);
% Spike cluster_id in the order of spike occurance
ap.goodSpikeClusters=spike_clusters(goodClusterSpikeIndices);
clear goodClusterSpikeIndices

disp(['Session ', sessionName,' loaded']);

%% (Required) Generate cue-triggered firing rate for each good cluster

% ITI is random variable between 2-3.5s
timeRange = [-1,2]; binSize = 0.01; % in sec

% Generate left-cue-triggered PSTH
leftToneOnIdx = find(allTrials==1);
[leftCueSpike,leftCueSpikeRate] = getSpikes(timeRange,binSize,leftToneOnIdx,...
                                    ap,nidq,timeNI);
% leftCueSpikeRate = leftCueSpike / (binSize*length(leftToneOnIdx));

% Generate right-cue-triggered PSTH
rightToneOnIdx = find(allTrials==2);
[rightCueSpike,rightCueSpikeRate] = getSpikes(timeRange,binSize,rightToneOnIdx,...
                                    ap,nidq,timeNI);
% rightCueSpikeRate = rightCueSpike / (binSize*length(rightToneOnIdx));

%% Generate omission PSTH

% Left Omission
leftToneOnIdx = find(omissionTrials==1);
[leftOmissionSpike,leftOmissionSpikeRate] = getSpikes(timeRange,binSize,leftToneOnIdx,...
                                    ap,nidq,timeNI);
% leftOmissionSpikeRate = leftOmissionSpike / (binSize*length(leftToneOnIdx));

% Right-omission PSTH
rightToneOnIdx = find(omissionTrials==2);
[rightOmissionSpike,rightOmissionSpikeRate] = getSpikes(timeRange,binSize,rightToneOnIdx,...
                                    ap,nidq,timeNI);
% rightOmissionSpikeRate = rightOmissionSpike / (binSize*length(rightToneOnIdx));

% Remove neurons with high mean firing rate
threshold = 200;
leftRemoveIdx = find(mean(leftCueSpikeRate,2) > threshold);
leftCueSpikeRate(leftRemoveIdx,:) = [];
disp(['Removed ',num2str(length(leftRemoveIdx)),' left units']);
rightRemoveIdx = find(mean(rightCueSpikeRate,2) > threshold);
rightCueSpikeRate(rightRemoveIdx,:) = [];
disp(['Removed ',num2str(length(rightRemoveIdx)),' right units']);

clear leftToneOnIdx rightToneOnIdx threshold leftRemoveIdx rightRemoveIdx 

%% Select example neurons

% Neurons w/ mean firing rate at 0%, 25%, 50%, 75%, 100% 
% from the reference cue (left/right)
[~,sortedFRIdx] = sort(mean(leftCueSpikeRate,2),'ascend');
drawExampleIdx = floor(linspace(1,length(sortedFRIdx),10)); % Last input is the number of samples drawn
exampleClusters_fr = sortedFRIdx(drawExampleIdx);
clear drawExampleIdx

%% Plot separate trace per neuron (averaged cue-triggered)

clusterList = 1:ap.nGoodClusters;
% clusterList = exampleClusters;
% clusterList = both_lick_clusters;
textOn = true;

% Plot left cue
spikeRate = leftCueSpikeRate; figure; 
drawSpikeRate('Left cue',timeRange,spikeRate,clusterList,textOn,ap,colors)
subtitle(['Time bin: ',num2str(binSize*1000),' ms; '...
      'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
      
% Plot right cue
spikeRate = rightCueSpikeRate; figure; 
drawSpikeRate('Right cue',timeRange,spikeRate,clusterList,textOn,ap,colors)
subtitle(['Time bin: ',num2str(binSize*1000),' ms; '...
      'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
  
% Plot left omission cue
spikeRate = leftOmissionSpikeRate; figure; 
drawSpikeRate('Left omission',timeRange,spikeRate,clusterList,textOn,ap,colors)
subtitle(['Time bin: ',num2str(binSize*1000),' ms; '...
      'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
  
% Plot right omission cue
spikeRate = rightOmissionSpikeRate; figure; 
drawSpikeRate('Right omission',timeRange,spikeRate,clusterList,textOn,ap,colors)
subtitle(['Time bin: ',num2str(binSize*1000),' ms; '...
      'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
  
clear spikeRate

%% Plot averaged trace for all neuron (averaged cue-triggered)

% left = mean(leftCueSpikeRate);
% right = mean(rightCueSpikeRate);
% leftOmission = mean(leftOmissionSpikeRate);
% rightOmission = mean(rightOmissionSpikeRate);

% Plot left/right cue
figure;
t = linspace(timeRange(1),timeRange(2),size(leftCueSpikeRate,2));
% plot(t,smooth(left),'LineWidth',2,'Color',colors(1)); hold on
% plot(t,smooth(right),'LineWidth',2,'Color',colors(2));
plotCI(t,leftCueSpikeRate,colors(1)); hold on
plotCI(t,rightCueSpikeRate,colors(2)); hold on
legend({'Left','Right'});

xlabel('Time (sec)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Spikes/s'); % ylim([0,150]);
xline(0,'-','Tone','Color','r','LineWidth',2,'HandleVisibility','off');
title('Averged PSTH of tone cue for all neurons');
subtitle(['Time bin: ',num2str(binSize*1000),' ms; '...
          'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
      

%% Plot left/right omission cue
figure;
t = linspace(timeRange(1),timeRange(2),size(leftOmissionSpikeRate,2));
% plot(t,smooth(leftOmission),'LineWidth',2,'Color',colors(1)); hold on
% plot(t,smooth(rightOmission),'LineWidth',2,'Color',colors(2));
plotCI(t,leftOmissionSpikeRate,colors(1)); hold on
plotCI(t,rightOmissionSpikeRate,colors(2)); hold on
legend({'Left','Right'});

xlabel('Time (sec)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Spikes/s'); % ylim([0,150]);
xline(0,'-','Tone','Color','r','LineWidth',2,'HandleVisibility','off');
title('Averged PSTH of Omission trials for all neurons');
subtitle(['Time bin: ',num2str(binSize*1000),' ms; '...
          'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
      
% clear left right left_err right_err leftOmission rightOmission leftOmission_err rightOmission_err

%% Plot raster plot of licking

timeRange = [-1,2];

% Plot left-cue-triggered licks
eventIdx = find(allTrials==1);
figure; drawLicks('Left cue',timeRange,eventIdx,...
                        leftLick,rightLick,nidq,timeNI)

% Plot right-cue-triggered licks
eventIdx = find(allTrials==2);
figure; drawLicks('Right cue',timeRange,eventIdx,...
                        leftLick,rightLick,nidq,timeNI)
                    
% Plot left-omission-triggered licks
eventIdx = find(omissionTrials==1);
figure; drawLicks('Left omission',timeRange,eventIdx,...
                        leftLick,rightLick,nidq,timeNI) 
                    
% Plot right-omission-triggered licks
eventIdx = find(omissionTrials==2);
figure; drawLicks('Right omission',timeRange,eventIdx,...
                        leftLick,rightLick,nidq,timeNI)                    

%% Plot licking and spiking together
% Some neurons spikes rthymically with licking 

timeRange = [-1,2]; textOn = true;
% clusterList = 1:ap.nGoodClusters; % in goodcluster id
clusterList = exampleClusters_fr;

% Plot left cue trials
set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
spikeRate = leftCueSpikeRate; eventIdx = find(allTrials==1);
yyaxis left; drawLicks('Left cue',timeRange,eventIdx,...
                            leftLick,rightLick,nidq,timeNI);
yyaxis right; drawSpikeRate('Left cue',timeRange,spikeRate,clusterList,textOn,ap,colors);

% Plot right cue trials
set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
spikeRate = rightCueSpikeRate; eventIdx = find(allTrials==2);
yyaxis left; drawLicks('Right cue',timeRange,eventIdx,...
                            leftLick,rightLick,nidq,timeNI);
yyaxis right; drawSpikeRate('Right cue',timeRange,spikeRate,clusterList,textOn,ap,colors);

% Plot left omission trials
set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
spikeRate = leftCueSpikeRate; eventIdx = find(omissionTrials==1);
yyaxis left; drawLicks('Left omission',timeRange,eventIdx,...
                            leftLick,rightLick,nidq,timeNI);
yyaxis right; drawSpikeRate('Left omission',timeRange,spikeRate,clusterList,textOn,ap,colors);
title(['Averged PSTH of left omission trial per neuron']);

% Plot right omission trials
set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
spikeRate = rightCueSpikeRate; eventIdx = find(omissionTrials==2);
yyaxis left; drawLicks('Right omission',timeRange,eventIdx,...
                            leftLick,rightLick,nidq,timeNI);
yyaxis right; drawSpikeRate('Right omission',timeRange,spikeRate,clusterList,textOn,ap,colors);
title(['Averged PSTH of right omission trial per neuron']);

%% Plot single cue PSTH (optional)
leftOn = true;
if leftOn; eventIdx = find(allTrials=1);
else; eventIdx = find(allTrials==2);
end
timeRange = [-2,5]; i = 2;

% Find first & last NI index
niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidq.Fs);
niLastIdx = eventIdx(i) + floor(timeRange(2)*nidq.Fs);

% Find corresponding imec index
[~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
imecLastIdx = imecFirstIdx + floor(ap.Fs*(timeRange(2)-timeRange(1)));
time_axis = timeImec(imecFirstIdx:imecLastIdx);

% Find spike within timeRange
spikeTimesinRangeIdx = find...
    (ap.goodSpikeTimes>imecFirstIdx & ap.goodSpikeTimes<imecLastIdx);
spikeTimesinRange = ap.goodSpikeTimes(spikeTimesinRangeIdx);
spikeClustersinRange = ap.goodSpikeClusters(spikeTimesinRangeIdx);

% Generate spike raster array
figure
clusterlist = unique(spikeClustersinRange);
for c = 1:length(clusterlist)
    cluster_id = spikeClustersinRange(c);
    cSpikeIdx = find(spikeClustersinRange == cluster_id);
    cluster = ap.clusterToGoodClusterIndex(cluster_id);
    scatter(timeImec(spikeTimesinRange(cSpikeIdx))-timeNI(eventIdx(i)),...
        cluster,'filled'); hold on
end
xlabel('Time');
ylabel('Firing rate (spikes/sec)'); ylim([0,150]);
if leftOn
    xline(0,'-','Left tone','Color','r','LineWidth',2);
    title(['Averged PSTH of left cue']);
else 
    xline(0,'-','Right tone','Color','r','LineWidth',2);
    title(['Averged PSTH of right cue']);
end

%% ***************** Firirng rate PCA NOTES *******************
%{
    1. Calculate averaged firing rate of EACH neuron for EACH trial
    2. Concatenate average firing rate together (e.g. Left-Right-Omission)
    3. PCA
%}

%% PCA: PCA matrix (Neurons: row; timepoints: col)

% spikeAvg = left/rightCueSpikeRate + left/rightOmissionSpikeRate
spikeAvg = [leftCueSpikeRate rightCueSpikeRate ...
            leftOmissionSpikeRate rightOmissionSpikeRate];
[coeff,score,~,~,explained,~] = pca(zscore(spikeAvg,0,2));

% Plot top PC vs time
max_value = max(max(coeff,[],2)); nPC = 15;
figure;
for i = 1:nPC
    t = linspace(0,length(spikeAvg)*binSize,length(spikeAvg));
    plot(t,i+((coeff(:,i)-coeff(1,i))/max_value),'LineWidth',1.5); hold on
end
yticks([1:nPC]); xlabel('Time (s)');
set(gca,'yticklabel',{'PC1','PC2','PC3','PC4','PC5',...
                      'PC6','PC7','PC8','PC9','PC10',...
                      'PC11','PC12','PC13','PC14','PC15'});

% Plot PCA results: variance
figure;
plot(explained,'-x','LineWidth',2,'Color',colors(1)); hold on
plot(cumsum(explained),'-x','LineWidth',2,'Color',colors(2));
xlim([0,25]); legend({'Individual explained variance','Cumulative explained variance'});
xlabel('Principal component'); ylabel('% total variance explained');

% Plot PCA results: score
figure;
for i = 1:size(score,1)
    scatter3(score(i,1),score(i,2),score(i,3),'filled');
    text(score(i,1),score(i,2),score(i,3),num2str(ap.goodClusters(i)));
    hold on
end
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

%% PCA: plot top PC licking vs firing

nTopPC = 3; textOn = 1;
example_frpca = zeros(10,nTopPC);
for i = 1:nTopPC
    % Select top 5 and lowest 5
    [~,Itop] = maxk(score(:,i),5);
    example_frpca(1:5,i) = Itop;
    [~,Ilow] = mink(score(:,i),5);
    example_frpca(6:end,i) = Ilow;
    
    % Plot left cue trials
    % clusterList = exampleClusters_pca(:,i);
    set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
    spikeRate = leftCueSpikeRate; eventIdx = find(allTrials==1);
    yyaxis left; drawLicks('Left cue',timeRange,eventIdx,...
        leftLick,rightLick,nidq,timeNI);
    yyaxis right; drawSpikeRate('Left cue',timeRange,spikeRate,clusterList,textOn,ap,colors);

    % Plot right cue trials
    set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
    spikeRate = rightCueSpikeRate; eventIdx = find(allTrials==2);
    yyaxis left; drawLicks('Right cue',timeRange,eventIdx,...
        leftLick,rightLick,nidq,timeNI);
    yyaxis right; drawSpikeRate('Right cue',timeRange,spikeRate,clusterList,textOn,ap,colors);
end
% autoArrangeFigures()

%% dPCA:

%% ***************** Lick-modulated spiking analysis *******************
%{
    Recreating Figure 4 from Chen et al., Cell Reports, 2021 and Figure 1
    from Rossi et al, Nat Neuro, 2016

    1. Find all lick indexes
    2. Look -300 ms to 300 ms
    3. Find licks and spikes within this time range
    4. Plot time from spout contact vs lick/spike freq
    5. Compare lick rate vs spike rate

    Recreating Figure 2 from Rossi et al., Nat Neuro, 2016

%}

%% Lick vs spiking: plot

% Set params
timeRange = [-0.3,0.3]; binSize = 0.01; % in sec (should smaller than 100 ms)
nBinsPeriLick = (timeRange(2)-timeRange(1))/binSize;
totalTime = session.nBlocks*session.blockTime;
nBins = totalTime / binSize;
nSampPerBinNI = length(timeNI)/nBins;

% Find lick index by combining left/right licks
combine_lick = zeros(size(rightLick));
combine_lick(leftLick) = 1; combine_lick(rightLick) = 2;
% Down sample combine_lick
combineLickTime = find(combine_lick > 0);

% Create matrix to plot licks from spout contact
% Row: each lick instance, Col: time bin peri-lick
periLickImage = zeros(length(combineLickTime),nBinsPeriLick);

for i = 1:length(combineLickTime)
    % Find first and last index around lick
    
end

%% ***************** GLM NOTES *******************
%{
    Predictor name and respective column:
    1. l_cue    2. r_cue    3. l_solenoid   4. r_solenoid
    5. l_lick   6. r_lick   7. choice_lick  8. rewarded_lick 
    9. unrewarded_lick      10. velocity
    11. l_omission          12. r_omission
  
    Three ways to do things:
    1. Train GLM for each neuron with the same amount of inputs
    2. Train GLM for each neuron with 15 sample models
    3. Train GLM for each neuron and use backward selection to find the
    best model (lowest RMSE score)


%}

%% (Required) Create GLM: Set up predictor/response variables for GLM

% Set down sample params
binSize = 0.05; totalTime = session.nBlocks*session.blockTime;
nBins = totalTime / binSize;
nSampPerBinNI = length(timeNI)/nBins;
nSampPerBinImec = length(timeImec)/nBins;
% Create time shifted data
offset = -0.25:binSize:0.25; glm.offset = offset;

glm.binSize = binSize; glm.nBins = nBins;
glm.nSampPerBinNI = nSampPerBinNI;
glm.nSampPerBinImec = nSampPerBinImec;

% Set cue predictor
cueDurationinBin = 0.05/binSize;    % Auditory cue is 50ms
l_cue = zeros(nBins,1); cueOnBin = floor(find(allTrials==1)/nSampPerBinNI)+1;
l_cue(cueOnBin) = 1;
% l_cue(cueOnBin+cueDurationinBin) = 1;
r_cue = zeros(nBins,1); cueOnBin = floor(find(allTrials==2)/nSampPerBinNI)+1;
r_cue(cueOnBin) = 1;
% r_cue(cueOnBin+cueDurationinBin) = 1;
l_omission = zeros(nBins,1); cueOnBin = floor(find(omissionTrials==1)/nSampPerBinNI)+1;
l_omission(cueOnBin) = 1;
r_omission = zeros(nBins,1); cueOnBin = floor(find(omissionTrials==2)/nSampPerBinNI)+1;
r_omission(cueOnBin) = 1;

% Set water predictor
l_solenoid = zeros(nBins,1); solenoidOnBin = floor(find(leftSolenoid==1)/nSampPerBinNI)+1;
l_solenoid(solenoidOnBin) = 1;
r_solenoid = zeros(nBins,1); solenoidOnBin = floor(find(rightSolenoid==1)/nSampPerBinNI)+1;
r_solenoid(solenoidOnBin) = 1;

% Set lick predictor
l_lick = zeros(nBins,1); lickOnBin = floor(find(leftLick==1)/nSampPerBinNI)+1;
for i = 1:length(lickOnBin)
    l_lick(lickOnBin(i)) = l_lick(lickOnBin(i)) + 1;
end
r_lick = zeros(nBins,1); lickOnBin = floor(find(rightLick==1)/nSampPerBinNI)+1;
for i = 1:length(lickOnBin)
    r_lick(lickOnBin(i)) = r_lick(lickOnBin(i)) + 1;
end

% Splitting licking subtype: choice lick, rewarded lick, unrewarded lick,
% spontaneous lick
% Step 1: combine cue into one vector (1 is left cue, 2 is right cue)
combine_cue = zeros(nBins,1);
combine_cue(l_cue==1) = 1;
combine_cue(r_cue==1) = 2;
combine_cue(l_omission==1) = 3;
combine_cue(r_omission==1) = 4;
glm.combine_cue = combine_cue;
% Step 2: find all licks within each cue
choice_lick = zeros(nBins,1); % First lick after the cue
rewarded_lick = zeros(nBins,1); % Licks in solenoid on session
unrewarded_lick = zeros(nBins,1); % Licks in solenoid off session
% spontaneous_lick = zeros(nBins,1); % Licks during ITI (3s after cue before next cue)
% Current state: spontaneous licks durng iti or before first cue are ignored
for i = 1:size(glm.combine_cue,1)
    if combine_cue(i) == 1
        % Find next cue time
        next_cue_time = find(combine_cue(i+1:end),1)+i;
        if isempty(next_cue_time); next_cue_time = nBins; end
        % iti = (next_cue_time - i)*binSize;
        
        % Find choice lick
        first_l_lick = find(l_lick(i+1:end),1)+i;
        first_r_lick = find(r_lick(i+1:end),1)+i;
        if isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = first_r_lick;
        elseif ~isempty(first_l_lick) && isempty(first_r_lick)
            choice_lick_time = first_l_lick;
        elseif ~isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = min(first_l_lick,first_r_lick);
        else; continue; 
        end
        if choice_lick_time < next_cue_time
            choice_lick(choice_lick_time) = 1;
            % disp(['cl=',num2str(choice_lick_time),' i=',num2str(i)]);
        end
        
        % Find rewarded/unrewarded lick
        rewarded_lick(find(l_lick(choice_lick_time+1:next_cue_time))+i) = 1;
        unrewarded_lick(find(r_lick(choice_lick_time+1:next_cue_time))+i) = 1;
        
    elseif combine_cue(i) == 2
        % Find next cue time
        next_cue_time = find(combine_cue(i+1:end),1)+i;
        if isempty(next_cue_time); next_cue_time = nBins; end
        % iti = (next_cue_time - i)*binSize;
        
        % Find choice lick
        first_l_lick = find(l_lick(i+1:end),1)+i;
        first_r_lick = find(r_lick(i+1:end),1)+i;
        if isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = first_r_lick;
        elseif ~isempty(first_l_lick) && isempty(first_r_lick)
            choice_lick_time = first_l_lick;
        elseif ~isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = min(first_l_lick,first_r_lick);
        else; continue; 
        end
        if choice_lick_time < next_cue_time
            choice_lick(choice_lick_time) = 1;
            % disp(num2str(choice_lick_time));
        end
        
        % Find rewarded/unrewarded lick
        rewarded_lick(find(r_lick(choice_lick_time+1:next_cue_time))+i) = 1;
        unrewarded_lick(find(l_lick(choice_lick_time+1:next_cue_time))+i) = 1;
    
    % Left omission
    elseif combine_cue(i) == 3
        % Find next cue time
        next_cue_time = find(combine_cue(i+1:end),1)+i;
        if isempty(next_cue_time); next_cue_time = nBins; end
        % iti = (next_cue_time - i)*binSize;
        
        % Find choice lick
        first_l_lick = find(l_lick(i+1:end),1)+i;
        first_r_lick = find(r_lick(i+1:end),1)+i;
        if isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = first_r_lick;
        elseif ~isempty(first_l_lick) && isempty(first_r_lick)
            choice_lick_time = first_l_lick;
        elseif ~isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = min(first_l_lick,first_r_lick);
        else; continue; 
        end
        if choice_lick_time < next_cue_time
            choice_lick(choice_lick_time) = 1;
            % disp(['cl=',num2str(choice_lick_time),' i=',num2str(i)]);
        end
        
        % Find unrewarded lick (all are unrewarded due to omission)
        unrewarded_lick(find(l_lick(choice_lick_time+1:next_cue_time))+i) = 1;
        unrewarded_lick(find(r_lick(choice_lick_time+1:next_cue_time))+i) = 1;
        
    
    % Right omission trial
    elseif combine_cue(i) == 4
        % Find next cue time
        next_cue_time = find(combine_cue(i+1:end),1)+i;
        if isempty(next_cue_time); next_cue_time = nBins; end
        % iti = (next_cue_time - i)*binSize;
        
        % Find choice lick
        first_l_lick = find(l_lick(i+1:end),1)+i;
        first_r_lick = find(r_lick(i+1:end),1)+i;
        if isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = first_r_lick;
        elseif ~isempty(first_l_lick) && isempty(first_r_lick)
            choice_lick_time = first_l_lick;
        elseif ~isempty(first_l_lick) && ~isempty(first_r_lick)
            choice_lick_time = min(first_l_lick,first_r_lick);
        else; continue; 
        end
        if choice_lick_time < next_cue_time
            choice_lick(choice_lick_time) = 1;
            % disp(num2str(choice_lick_time));
        end
        
        % Find rewarded/unrewarded lick
        unrewarded_lick(find(r_lick(choice_lick_time+1:next_cue_time))+i) = 1;
        unrewarded_lick(find(l_lick(choice_lick_time+1:next_cue_time))+i) = 1;
        
    else
        continue
    end
end

% Calculate velocity and acceleration
xgyro = downsample(gyroX,floor(nSampPerBinNI))'; 
xgyro(length(xgyro)-(length(xgyro)-nBins)+1:end) = [];
xgyro = xgyro - mean(xgyro);
ygyro = downsample(gyroY,floor(nSampPerBinNI))';
ygyro(length(ygyro)-(length(ygyro)-nBins)+1:end) = [];
ygyro = ygyro - mean(ygyro);
zgyro = downsample(gyroZ,floor(nSampPerBinNI))';
zgyro(length(zgyro)-(length(zgyro)-nBins)+1:end) = [];
zgyro = zgyro - mean(zgyro);
% Movement velocity (a.u.)
velocity = sqrt(xgyro.^2 + ygyro.^2 + zgyro.^2)*1000; 
% Movement acceleration (a.u.)
acceleration = [0;diff(velocity)];

% Set up spiking matrix
% Create response matrix Y
clusterList = 1:ap.nGoodClusters; % Set spike response vector
Y = zeros(length(clusterList),nBins); % store y_fr
Y_mean = zeros(length(clusterList),1); % mean firing rate

for i = 1:length(clusterList)
    cluster_id = ap.goodClusters(clusterList(i));
    
    % Set up spike (response vector)
    y_spike = zeros(nBins,1);
    spikeOnBin = floor(double(ap.goodSpikeTimes(ap.goodSpikeClusters==cluster_id))/nSampPerBinImec)+1;
    for j = 1:length(spikeOnBin)
        if spikeOnBin(j) <= nBins
            y_spike(spikeOnBin(j)) = 1 + y_spike(spikeOnBin(j)); % Spike count
        end
    end
    
    % Calculate firing rate
    y_fr = y_spike/binSize;
    Y(i,:) = y_fr'; 
    Y_mean(i) = mean(y_fr);
end
clear lickOnBin combine_cue xgyro ygyro zgyro i j

%% Create single GLM: Combine into one predictor matrix
X_baseline = [l_cue r_cue l_solenoid r_solenoid l_lick r_lick];
% X_baseline = [choice_lick rewarded_lick unrewarded_lick];
% X_baseline = [l_cue r_cue l_lick r_lick];
% X_baseline = [l_cue r_cue l_solenoid r_solenoid l_lick r_lick choice_lick rewarded_lick unrewarded_lick velocity];
nPredictors = size(X_baseline,2)+1; % include spike history

% Create response matrix Y
clusterList = 1:ap.nGoodClusters; % Set spike response vector
Y = zeros(length(clusterList),nBins); % store y_fr
Y_mean = zeros(length(clusterList),1); % mean firing rate

for i = 1:length(clusterList)
    cluster_id = ap.goodClusters(clusterList(i));
    
    % Set up spike (response vector)
    y_spike = zeros(nBins,1);
    spikeOnBin = floor(double(ap.goodSpikeTimes(ap.goodSpikeClusters==cluster_id))/nSampPerBinImec)+1;
    for j = 1:length(spikeOnBin)
        if spikeOnBin(j) <= nBins
            y_spike(spikeOnBin(j)) = 1 + y_spike(spikeOnBin(j)); % Spike count
        end
    end
    
    % Calculate firing rate
    y_fr = y_spike/binSize;
    Y(i,:) = y_fr'; 
    Y_mean(i) = mean(y_fr);
end

% Combine all predictor data into one matrix
X = [];

for offsetIdx = 1:length(offset)
    
    X_shifted = X_baseline;
    offsetInBin = uint32(abs(offset(offsetIdx))/binSize); 
    padding = zeros(offsetInBin,nPredictors);
    
    % Shifts behavioral data based on offset time
    if offset(offsetIdx) > 0
        X_shifted(nBins-offsetInBin+1:end,:) = [];
        X_shifted = [padding;X_shifted];
    elseif offset(offsetIdx) < 0
        X_shifted(1:offsetInBin,:) = [];
        X_shifted = [X_shifted;padding];
    end
    
    % Combine time shifted data with baseline
    disp(['Combining behavioral data shifted by ',num2str(offset(offsetIdx)*1000),' ms']);
    X = [X, X_shifted];
end

%% Create single GLM: partition into training vs testing data

% 80% training, 20% testing
totalTestBins = 0.2*nBins; nTestChunk = 5;
nBinsPerChunk = totalTestBins/nTestChunk; % 960 Bins/chunk

% Generate test chunk index
cueBins = find(glm.combine_cue>0);
resplit = true; nSplitAttempt = 0;
while resplit
    % testing chunk ends 1 second before the next tone
    chunkEndBin = sort(cueBins(randperm(length(cueBins),nTestChunk))-20);
    nSplitAttempt = nSplitAttempt + 1;
    flag = 0;
    for i = 1:length(chunkEndBin)
        if chunkEndBin(i) - nBinsPerChunk < 0; flag=1; break
        elseif chunkEndBin(i) - nBinsPerChunk > 0
            if i == 1; continue
            elseif chunkEndBin(i) - nBinsPerChunk < chunkEndBin(i-1); flag=2; break
            end
        end
    end
    if flag == 0
        resplit = false;
        disp(['Found testing partitions in ',num2str(nSplitAttempt),' attempts']);
    end
end

% Genearte test data index
chunkStartBin = chunkEndBin - nBinsPerChunk + 1;
binIdx = linspace(1,nBins,nBins)'; chunkIdx = [];
for idx = 1:nTestChunk
    chunkIdx = [chunkIdx;binIdx(chunkStartBin(idx):chunkEndBin(idx))];
end

% Separate test and training data
X_test = X(chunkIdx,:);
X_training = X(setdiff(binIdx,chunkIdx),:);
Y_test = Y(:,chunkIdx);
Y_training = Y(:,setdiff(binIdx,chunkIdx));

%% Create single GLM: run GLM model for each neuron

% Set spike response vector
clusterList = 1:ap.nGoodClusters;
M = cell(length(clusterList),1); % for storing GLM model for each neuron

for i = 1:length(clusterList)
    cluster_id = ap.goodClusters(clusterList(i));
    y_fr_training = Y_training(i,:)';
    
    % Run GLM
    disp(['Generating GLM for neuron ',num2str(clusterList(i))]);
    model = fitglm(X_training,y_fr_training,'linear','Distribution','normal');
    M{i} = model;
end

% Set up coefficients and p-value table
modelSummary = zeros(length(clusterList),size(X,2)*2+2); % Last column is intercept and it's p-value
r2 = zeros(length(clusterList),2); % Column 1: ordinary; Column 2: adjusted

for i = 1:length(M) % iterate through each neuron
    
    r2(i,1) = M{i}.Rsquared.Ordinary;
    r2(i,2) = M{i}.Rsquared.Adjusted;
    
    coefficient = table2array(M{i}.Coefficients(2:end,:));
    for j = 1:size(coefficient,1) % iterate through each variable
        modelSummary(i,j) = coefficient(j,1);
        modelSummary(i,j+nPredictors) = coefficient(j,4); 
    end
    modelSummary(i,end-1) = table2array(M{i}.Coefficients(1,1));
    modelSummary(i,end) = table2array(M{i}.Coefficients(1,4));
end

% % Remove non-significant coefficients
% for i = 1:size(modelSummary,1)
%     for j = 10:size(modelSummary,2)
%         if modelSummary(i,j) >= 0.05
%             modelSummary(i,j-9) = 0;
%         end
%     end
% end

%% Create single GLM: cross validation on training set (K = 10)

% Prediction function given training/testing instances
fcn = @(Xtr, Ytr, Xte) predict(...
    GeneralizedLinearModel.fit(Xtr,Ytr,'linear','distr','normal'), ...
    Xte);

% Store averaged RMSE for each neuron
RMSE = zeros(length(clusterList),1);
for i = 1:length(clusterList)
    % Perform cross-validation, and return average MSE across folds
    mse = crossval('mse', X_training, Y_training(i,:)', 'Predfun',fcn, 'kfold',10);
    
    % Compute root mean squared error
    RMSE(i) = sqrt(mse);
end

%% Create single GLM: use function version

% predictorName = {'Left cue','Right cue','Left solenoid','Right solenoid','Left licks','Right licks'};
% predictorName = {'Left cue','Right cue','Left Omission','Right Omission',...
%     'Left solenoid','Right solenoid','Left licks','Right licks',...
%     'Choice licks','Rewarded licks','Unrewarded licks'};
predictorName = {'Left cue','Right cue','Left Omission','Right Omission',...
    'Left solenoid','Right solenoid','Left licks','Right licks',...
    'Choice licks','Rewarded licks','Unrewarded licks',...
    'Velocity'};


clusterList = 1:ap.nGoodClusters; % Set spike response vector
M = cell(length(clusterList),1); % for storing GLM model for each neuron
X_baseline = [l_cue r_cue l_omission r_omission l_solenoid r_solenoid ...
              l_lick r_lick choice_lick rewarded_lick unrewarded_lick velocity];
nPredictors = size(X_baseline,2);

for i = 1:length(clusterList)
    disp(['Generating GLM for neuron ',num2str(clusterList(i))]);
%     X_baseline = [X_baseline Y(clusterList(i),:)']; % include spike history
    [~,model] = runGLM(X_baseline,clusterList(i),glm,ap,true);
    M{i} = model;
end

% Set up coefficients and p-value table
modelSummary = zeros(length(clusterList),nPredictors*length(offset)*2+2); % Last column is intercept and it's p-value
r2 = zeros(length(clusterList),2); % Column 1: ordinary; Column 2: adjusted

for i = 1:length(M) % iterate through each neuron
    
    r2(i,1) = M{i}.Rsquared.Ordinary;
    r2(i,2) = M{i}.Rsquared.Adjusted;
    
    coefficient = table2array(M{i}.Coefficients(2:end,:));
    for j = 1:size(coefficient,1) % iterate through each variable
        modelSummary(i,j) = coefficient(j,1);
        modelSummary(i,j+nPredictors) = coefficient(j,4); 
    end
    modelSummary(i,end-1) = table2array(M{i}.Coefficients(1,1));
    modelSummary(i,end) = table2array(M{i}.Coefficients(1,4));
end

%% Iterate multiple GLM: select overall best model

nModels = 15;
clusterList = 1:ap.nGoodClusters; % Set spike response vector
R2 = zeros(length(clusterList),nModels); % Store r2 of each neuron for all model tested
RMSE = zeros(length(clusterList),nModels); % Store averaged rmse of each neuron for all model tested

for m = 1:nModels
    
    disp(['Start analyzing GLM model ',num2str(m)]);
    
    % STEP 1: Set up X_baseline based on current testing model
    if m == 1;      X_baseline = [l_cue];
    elseif m == 2;  X_baseline = [r_cue];
    elseif m == 3;  X_baseline = [l_cue r_cue];
    elseif m == 4;  X_baseline = [l_solenoid r_solenoid];
    elseif m == 5;  X_baseline = [l_lick r_lick];
    elseif m == 6;  X_baseline = [choice_lick rewarded_lick unrewarded_lick];
    elseif m == 7;  X_baseline = [l_cue r_cue l_solenoid r_solenoid];
    elseif m == 8;  X_baseline = [l_cue r_cue l_lick r_lick];
    elseif m == 9;  X_baseline = [l_cue r_cue choice_lick rewarded_lick unrewarded_lick];
    elseif m == 10; X_baseline = [l_cue r_cue l_solenoid r_solenoid l_lick r_lick];
    elseif m == 11; X_baseline = [l_cue r_cue l_solenoid r_solenoid choice_lick rewarded_lick unrewarded_lick];
    elseif m == 12; X_baseline = [l_cue r_cue l_solenoid r_solenoid l_lick r_lick choice_lick rewarded_lick unrewarded_lick];
    elseif m == 13; X_baseline = [l_cue r_cue l_solenoid r_solenoid l_lick r_lick xgyro ygyro zgyro];
    elseif m == 14; X_baseline = [l_cue r_cue l_solenoid r_solenoid choice_lick rewarded_lick unrewarded_lick xgyro ygyro zgyro];
    elseif m == 15; X_baseline = [l_cue r_cue l_solenoid r_solenoid l_lick r_lick choice_lick rewarded_lick unrewarded_lick xgyro ygyro zgyro];
    end
    nPredictors = size(X_baseline,2);
    
    % STEP 2: Add offset padding
    offset = -0.25:binSize:0.25;
    X = []; % Combine all predictor data into one matrix
    for offsetIdx = 1:length(offset)

        X_shifted = X_baseline;
        offsetInBin = uint32(abs(offset(offsetIdx))/binSize); 
        padding = zeros(offsetInBin,nPredictors);

        % Shifts behavioral data based on offset time
        if offset(offsetIdx) > 0
            X_shifted(nBins-offsetInBin+1:end,:) = [];
            X_shifted = [padding;X_shifted];
        elseif offset(offsetIdx) < 0
            X_shifted(1:offsetInBin,:) = [];
            X_shifted = [X_shifted;padding];
        end

        % Combine time shifted data with baseline
        disp(['     STEP 2: Combining behavioral data shifted by ',num2str(offset(offsetIdx)*1000),' ms']);
        X = [X, X_shifted];
    end
    
    % STEP 3: Create response matrix Y
    Y = zeros(length(clusterList),nBins); % store y_fr
    Y_mean = zeros(length(clusterList),1); % mean firing rate
    for i = 1:length(clusterList)
        cluster_id = ap.goodClusters(clusterList(i));

        % Set up spike (response vector)
        y_spike = zeros(nBins,1);
        spikeOnBin = floor(double(ap.goodSpikeTimes(ap.goodSpikeClusters==cluster_id))/nSampPerBinImec)+1;
        for j = 1:length(spikeOnBin)
            if spikeOnBin(j) <= nBins
                y_spike(spikeOnBin(j)) = 1 + y_spike(spikeOnBin(j)); % Spike count
            end
        end

        % Calculate firing rate
        y_fr = y_spike/binSize;
        Y(i,:) = y_fr'; 
        Y_mean(i) = mean(y_fr);
    end
    
    % STEP 4: Separate data into testing and training set
    % 80% training, 20% testing
    totalTestBins = 0.2*nBins; nTestChunk = 5;
    nBinsPerChunk = totalTestBins/nTestChunk; % 960 Bins/chunk

    % Generate test chunk index
    cueBins = find(glm.combine_cue>0);
    resplit = true; nSplitAttempt = 0;
    while resplit
        % testing chunk ends 1 second before the next tone
        chunkEndBin = sort(cueBins(randperm(length(cueBins),nTestChunk))-20);
        nSplitAttempt = nSplitAttempt + 1;
        flag = 0;
        for i = 1:length(chunkEndBin)
            if chunkEndBin(i) - nBinsPerChunk < 0; flag=1; break
            elseif chunkEndBin(i) - nBinsPerChunk > 0
                if i == 1; continue
                elseif chunkEndBin(i) - nBinsPerChunk < chunkEndBin(i-1); flag=2; break
                end
            end
        end
        if flag == 0
            resplit = false;
            disp(['     STEP 4: Found testing partitions in ',num2str(nSplitAttempt),' attempts']);
        end
    end

    % Genearte test data index
    chunkStartBin = chunkEndBin - nBinsPerChunk + 1;
    binIdx = linspace(1,nBins,nBins)'; chunkIdx = [];
    for idx = 1:nTestChunk
        chunkIdx = [chunkIdx;binIdx(chunkStartBin(idx):chunkEndBin(idx))];
    end

    % Separate test and training data
    X_test = X(chunkIdx,:);
    X_training = X(setdiff(binIdx,chunkIdx),:);
    Y_test = Y(:,chunkIdx);
    Y_training = Y(:,setdiff(binIdx,chunkIdx));
    
    % STEP 5: Run GLM for a given neuron and model
    M = cell(length(clusterList),1); % for storing GLM model for each neuron

    for i = 1:length(clusterList)
        cluster_id = ap.goodClusters(clusterList(i));
        y_fr_training = Y_training(i,:)';
        
        % Run GLM
        disp(['     STEP 5: Generating GLM model ',num2str(m),...
            ' for neuron ',num2str(clusterList(i))]);
        model = fitglm(X_training,y_fr_training,'linear','Distribution','normal');
        M{i} = model;
    end

    % STEP 6: record r2 for each neuron for the current model
    % modelSummary = zeros(length(clusterList),size(X,2)*2+2); % Last column is intercept and it's p-value
    % r2 = zeros(length(clusterList),1); % Column 1: ordinary; Column 2: adjusted
    for i = 1:length(M) % iterate through each neuron
        % r2(i,1) = M{i}.Rsquared.Ordinary;
        r2(i,1) = M{i}.Rsquared.Adjusted;
    end
    R2(:,m) = r2;
    
    % STEP 7: 10-fold CV to calculate averaged RMSE for each neuron
    % Prediction function given training/testing instances
    fcn = @(Xtr, Ytr, Xte) predict(...
        GeneralizedLinearModel.fit(Xtr,Ytr,'linear','distr','normal'), ...
        Xte);
    % Store averaged RMSE for each neuron
    for i = 1:length(clusterList)
        % Perform cross-validation, and return average MSE across folds
        disp(['     STEP 7: Running CV of GLM model ',num2str(m),...
            ' for neuron ',num2str(clusterList(i))]);
        mse = crossval('mse', X_training, Y_training(i,:)', 'Predfun',fcn, 'kfold',10);
        % Compute root mean squared error
        RMSE(i,m) = sqrt(mse);
    end
    
    % STEP 8: Finished!!
    disp(['Finished GLM model ',num2str(m)]);
end

%% Iterate multiple GLM: backward select best model for each neuron

X_baseline = [l_cue r_cue l_omission r_omission l_solenoid r_solenoid ...
        l_lick r_lick choice_lick rewarded_lick unrewarded_lick velocity];
glm.X_baseline = X_baseline;
clusterList = 1:ap.nGoodClusters; % Set spike response vector
offset = -0.25:binSize:0.25; glm.offset = offset;
full_model = [1 2 3 4 5 6 7 8 9 10 11 12]; % Including both velocity and acceleration causes error (rank deficient)

% Store them into array (selected predictors will have the min_rmse value)
bestModel = zeros(length(clusterList),length(full_model));

for neuron = 1:length(clusterList)
    % Set up X_baseline to include spike history
    % y_fr = Y(i,:)';
    % X_baseline = [X_baseline y_fr];
    
    % Run backward selection
    min_rmse = 999;
    [best_neuron_model,min_rmse] = backSelect(full_model,min_rmse,neuron,...
                                                X_baseline,glm,ap);
    if best_neuron_model ~= 0
        bestModel(neuron,best_neuron_model) = min_rmse;
    end
end

%% Analyze GLM: calculate proportion of each model

% Number of responsive predictors
nPredCount = zeros(size(bestModel,2),1);
% row is neuron, column is predictor, element is the number of total
% responsive predictors
nPredIdentity = zeros(size(bestModel));
for i = 1:size(bestModel,1)
    nPred = length(find(bestModel(i,:) > 0));
    nPredCount(nPred) = nPredCount(nPred) + 1;
    nPredIdentity(i,bestModel(i,:)>0) = nPred;
end

% Pie chart
nPredRatio = nPredCount/size(bestModel,1);
figure; labels = {'1','2','3','4','5','6','7','8','9','10','11','12'};
pie(nPredRatio); legend(labels);
title('Proportion of predictor responsiveness');

% Bar plot
figure; bar(1:length(nPredCount),nPredCount);
xlabel('Number of predictors'); ylabel('Number of neurons');
title('Number of predictors per neuron');
subtitle(['Total neuron analyzed: ',num2str(size(bestModel,1)),...
          '; Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
      

% Calculate predictor heatmap
% pred_heatmap = nPredictors x nPredictors
% Row is predictor, column means for this predictor, how many neurons are
% also cue/solenoid/... responsive
pred_heatmap = zeros(size(bestModel,2));
for i = 1:size(bestModel,2)
    predCol = bestModel(:,i);
    for j = 1:length(predCol) % for each neuron
        if predCol(j) == 0; continue; end
        neuron = bestModel(j,:);
        for k = 1:size(neuron,2)
            if neuron(k) == 0; continue; end
            pred_heatmap(i,k) = pred_heatmap(i,k) + 1;
        end   
    end
end

figure; imagesc(pred_heatmap); colorbar;
xticks([1:size(pred_heatmap,1)]);
set(gca,'xticklabel',{'Left cue','Right cue','Left omission','Right omission',...
                      'Left solenoid','Right solenoid',...
                      'Left licks','Right licks',...
                      'Choice licks','Rewarded licks','Unrewarded licks',...
                      'Velocity'});
yticks([1:size(pred_heatmap,1)]);
set(gca,'yticklabel',{'Left cue','Right cue','Left omission','Right omission',...
                      'Left solenoid','Right solenoid',...
                      'Left licks','Right licks',...
                      'Choice licks','Rewarded licks','Unrewarded licks',...
                      'Velocity'});
title('Heatmap of co-responsive predictors');
subtitle(['Total neuron analyzed: ',num2str(size(bestModel,1)),...
          '; Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);

      
% % Calculate proportion of multi-responsive neurons for each predictors
% % Row is predictor, column means for this predictor, how many neurons are
% % nPred (e.g. 1,2,3) responsive
stackedPred = zeros(size(bestModel,2));
for i = 1:size(nPredIdentity,2)
    predCol = nPredIdentity(:,i);
    for j = 1:length(predCol)
        if predCol(j) == 0; continue; end
        stackedPred(i,predCol(j)) = stackedPred(i,predCol(j)) + 1;    
    end

end
figure; bar(stackedPred,'stacked');
xticks([1:size(nPredIdentity,1)]);
set(gca,'xticklabel',{'Left cue','Right cue','Left Omission','Right Omission',...
                      'Left solenoid','Right solenoid',...
                      'Left licks','Right licks',...
                      'Choice licks','Rewarded licks','Unrewarded licks',...
                      'Velocity'});
ylabel('Number of neurons');legend({'1','2','3','4','5','6','7','8','9','10','11','12'});
title('Proportion of multi-responsive neurons for each predictor');
subtitle(['Total neuron analyzed: ',num2str(size(bestModel,1)),...
          '; Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);

%% Analyze GLM: Test model performance on testing data

% Single neuron
neuron = 1; % good: 43,98,105; bad: 1
disp(['Predicting cluster ',num2str(ap.goodClusters(i)),...
    ' with R^2 = ',num2str(r2(neuron,2))]);

predicted = predict(M{neuron},X_test);
residuals = Y_test(neuron,:)' - predicted;

%% Plot GLM results: plot coefficient vs offsets (beta)

% Set x axis
x = linspace(glm.offset(1),glm.offset(end),length(glm.offset));

% Store example neurons for each coefficient (> 1 std)
exampleClusters_glm = cell(nPredictors,1);

% Construct beta function for each neuron
beta = cell(size(modelSummary,1),nPredictors);
for n = 1:size(modelSummary,1)
    for p = 1:nPredictors
        beta_t = zeros(1,length(offset));
        for i = 1:length(offset)
            beta_t(i) = modelSummary(n,(i-1)*nPredictors+p);
        end
        
        % Store beta function for this neuron's coefficient
        beta{n,p} = beta_t;
    end
end

% Plot averaged coefficient vs offset time (averaged across neurons)
figure;
for p = 1:size(predictorName,2)
    x_mean = zeros(length(offset),1);
    CI = zeros(length(offset),1);
    for i = 1:length(offset)
        % Averaging coefficient across neurons
        x_mean(i) = mean(modelSummary(:,(i-1)*nPredictors+p),"omitnan");
        CI(i) = getCI(modelSummary(:,(i-1)*nPredictors+p),0.95);
    end
    subplot(4,3,p);
%     subplot(3,2,p);
    errorbar(x,x_mean,CI,'-','LineWidth',2);    % Plot CI as error bar

    xlabel('Offset time (s)'); ylabel('Coefficient estimates');
    title([predictorName{p},' coefficient']);
end

% Plot individual coefficient vs offset time
figure;
for p = 1:nPredictors
    
    subplot(4,3,p);
%     subplot(3,3,p);
    outlier = [];
    
    for n = 1:size(modelSummary,1)
        x_n = zeros(length(offset),1);
        for i = 1:length(offset)
            % Averaging coefficient across neurons
            x_n(i) = modelSummary(n,(i-1)*nPredictors+p);
        end
        
        plot(x,x_n,'-','LineWidth',2); hold on
        
        % Label outlier cluster
        [y_max,x_max] = max(abs(x_n));
        err = std(modelSummary(:,(i-1)*nPredictors+p),"omitnan");
        if abs(y_max) > 2.5 * err
            xtxtpos = offset(x_max);
            ytxtpos = x_n(x_max);
            cluster_id = ap.goodClusters(n);
            outlier = [outlier;cluster_id];
            text(xtxtpos,ytxtpos,num2str(cluster_id)); hold on
        end
        
        exampleClusters_glm{p} = outlier;
        
    end
    
    xlabel('Offset time (s)'); ylabel('Coefficient estimates');
    title([predictorName{p},' coefficient']);
end

%% Plot GLM results: Plot beta for one example neuron of choice
% neurons = [86;52;50;92;102;147;154]; % in good cluster id
neurons = [43];
% neurons = ap.clusterToGoodClusterIndex(exampleClusters_glm{1}(1:5));
% Set x axis
x = linspace(offset(1),offset(end),length(offset));

figure;
for p = 1:nPredictors
    
    subplot(3,2,p);
    % subplot(2,2,p);
    outlier = [];
    
    for n = 1:size(neurons,1)
        x_n = zeros(length(offset),1);
        gcid = neurons(n);
        for i = 1:length(offset)
            % Averaging coefficient across neurons
            x_n(i) = modelSummary(gcid,(i-1)*nPredictors+p);
        end

        plot(x,x_n,'-','LineWidth',2); hold on
    end
    
    xlabel('Offset time (s)'); ylabel('Coefficient estimates');
    title([predictorName{p},' coefficient']);
end

%% Plot GLM results: plot example clusters licking & spiking
% Some neurons spikes rthymically with licking 

timeRange = [-1,2]; textOn = true;
clusterList = exampleClusters_glm{2};
% clusterList = exampleClusters_pca(:,1);

% Plot left cue trials
clusterList = ap.clusterToGoodClusterIndex(clusterList);
set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
spikeRate = rightCueSpikeRate; eventIdx = find(allTrials==2);
yyaxis left; drawLicks('Left cue',timeRange,eventIdx,...
                            leftLick,rightLick,nidq,timeNI);
yyaxis right; drawSpikeRate('Left cue',timeRange,spikeRate,clusterList,textOn,ap,colors);

% Plot right cue trials
clusterList = ap.clusterToGoodClusterIndex(clusterList);
set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
spikeRate = rightCueSpikeRate; eventIdx = find(allTrials==2);
yyaxis left; drawLicks('Right cue',timeRange,eventIdx,...
                            leftLick,rightLick,nidq,timeNI);
yyaxis right; drawSpikeRate('Right cue',timeRange,spikeRate,clusterList,textOn,ap,colors);

%% Plot GLM results: plot mean normal model diagnostics: r^2, RMSE

figure;
for i = 1:size(R2,2)
    x = i*ones(1,length(clusterList));
    bar(i,mean(R2(:,i))); hold on
    swarmchart(x,R2(:,i),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5,'XJitterWidth',0.7); hold on
    errorbar(i,mean(R2(:,i)),getCI(R2(:,i),0.95),'k','LineWidth',1.5); hold on
end
xticks([1:15]);
set(gca,'xticklabel',{'Left cue','Right cue',...
    'Left/right cue','Left/right solenoid',...
    'Left/right licks','Choice/rewarded/unrewarded licks',...
    'Cues + solenoid','Cues + licks','Cues + lick subtypes',...
    'Cues + solenoid + licks','Cues + solenoid + lick subtypes',...
    'Cues + solenoid + licks + lick subtypes',...
    'Cues + solenoid + licks + velocity',...
    'Cues + solenoid + lick subtypes + velocity',...
    'Cues + solenoid + licks + lick subtypes + velocity'});
% xticklabels({'Left/right cue','Left/right cue/solenoid/lick'});
xlabel('Model'); ylabel('R^2'); ylim([0 0.7]);
title(['R^2 for GLMs fitted to firing rate of session 0331R']);

figure;
for i = 1:size(RMSE,2)
    x = i*ones(1,size(RMSE,1));
    bar(i,mean(RMSE(:,i))); hold on
    swarmchart(x,RMSE(:,i),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5,'XJitterWidth',0.7); hold on
    errorbar(i,mean(RMSE(:,i)),getCI(RMSE(:,i),0.95),'k','LineWidth',1.5); hold on
end
xticks([1:15]);
set(gca,'xticklabel',{'Left cue','Right cue',...
    'Left/right cue','Left/right solenoid',...
    'Left/right licks','Choice/rewarded/unrewarded licks',...
    'Cues + solenoid','Cues + licks','Cues + lick subtypes',...
    'Cues + solenoid + licks','Cues + solenoid + lick subtypes',...
    'Cues + solenoid + licks + lick subtypes',...
    'Cues + solenoid + licks + velocity',...
    'Cues + solenoid + lick subtypes + velocity',...
    'Cues + solenoid + licks + lick subtypes + velocity'});
% xticklabels({'Left/right cue','Left/right cue/solenoid/lick'});
xlabel('Model'); ylabel('RMSE'); % ylim([0 0.7]);
title(['RMSE for GLMs fitted to firing rate of session 0331R']);

%% Plot GLM results: plot individual RMSE across all models

% Plot r2
figure;
for i = 1:size(R2,1)
    plot(R2(i,:),'-x','LineWidth',1); hold on
end
xticks([1:15]);
set(gca,'xticklabel',{'Left cue','Right cue',...
    'Left/right cue','Left/right solenoid',...
    'Left/right licks','Choice/rewarded/unrewarded licks',...
    'Cues + solenoid','Cues + licks','Cues + lick subtypes',...
    'Cues + solenoid + licks','Cues + solenoid + lick subtypes',...
    'Cues + solenoid + licks + lick subtypes',...
    'Cues + solenoid + licks + velocity',...
    'Cues + solenoid + lick subtypes + velocity',...
    'Cues + solenoid + licks + lick subtypes + velocity'});
% xticklabels({'Left/right cue','Left/right cue/solenoid/lick'});
% xlabel('Model'); ylabel('R^2'); % ylim([0 0.7]);
title(['R^2 for GLMs fitted to firing rate of session 0331R']);

% % Plot RMSE
% figure;
% for i = 1:size(highR2Clusters,1)
%     plot(RMSE(i,:),'-x','LineWidth',1); hold on
%     text(0.5,RMSE(i,1),num2str(ap.goodClusters(highR2Clusters(i))));
% end
% xticks([1:15]);
% set(gca,'xticklabel',{'Left cue','Right cue',...
%     'Left/right cue','Left/right solenoid',...
%     'Left/right licks','Choice/rewarded/unrewarded licks',...
%     'Cues + solenoid','Cues + licks','Cues + lick subtypes',...
%     'Cues + solenoid + licks','Cues + solenoid + lick subtypes',...
%     'Cues + solenoid + licks + lick subtypes',...
%     'Cues + solenoid + licks + gyroXYZ',...
%     'Cues + solenoid + lick subtypes + gyroXYZ',...
%     'Cues + solenoid + licks + lick subtypes + gyroXYZ'});
% xticklabels({'Left/right cue','Left/right cue/solenoid/lick'});
% xlabel('Model'); ylabel('RMSE'); % ylim([0 0.7]);
% title(['RMSE for GLMs fitted to firing rate of session 0401R']);

%% Plot GLM results: plot predicted firing vs actual firing for first 50 s

% Single neuron
neuron = 1; % good: 43,98,105; bad: 1

disp(['Predicting cluster ',num2str(ap.goodClusters(i)),...
    ' with R^2 = ',num2str(r2(neuron,2))]);

[predicted,~] = predict(M{neuron},X_test,'Alpha',0.05,'Simultaneous',true);
residuals = Y_test(neuron,:)' - predicted;
% rmse

figure;
% subplot(1,2,1);
plot(Y_test(neuron,1:1000)); hold on; plot(predicted(1:1000));
% errorbar(1:1000,predicted(1:1000),ci(1:1000));
legend({'Actual spikes','Predicted spikes'});
xline(find(l_cue(1:1000)),'-','Left tone','Color','r','LineWidth',2,'HandleVisibility','off');
xline(find(r_cue(1:1000)),'-','Right tone','Color','r','LineWidth',2,'HandleVisibility','off');
xlabel('Time (50 ms)'); ylabel('Spikes/s');

% subplot(1,2,2);
% plot(residuals(1:1000));
% xline(find(l_cue(1:1000)),'-','Left tone','Color','r','LineWidth',2,'HandleVisibility','off');
% xline(find(r_cue(1:1000)),'-','Right tone','Color','r','LineWidth',2,'HandleVisibility','off');
% xlabel('Time (50 ms)'); ylabel('Residual');

%% Plot GLM results: cue-triggered averaged firing for predicted/observed data

timeRange = [-1,2]; neuron = 1;
nBinsPerTrial = round((timeRange(2)-timeRange(1))/glm.binSize);

% Generate predicted firing rate for testing data
Y_predicted = zeros(size(Y_test));
for n = 1:size(Y_test,1)
    [Y_predicted(n,:),~] = predict(M{n},X_test,'Alpha',0.05,'Simultaneous',true);
end

% Find cue on time
testLeftCueOnBin = find(X_test(:,(6-1)*nPredictors+1)==1);
testRightCueOnBin = find(X_test(:,(6-1)*nPredictors+2)==1);

psth_avg_fr_predicted_l = zeros(size(Y_test,1),nBinsPerTrial);
psth_avg_fr_observed_l = zeros(size(Y_test,1),nBinsPerTrial);
totalCueAnalyzed = 0;
for i = 1:length(testLeftCueOnBin)
    
    % Find start & end bin
    cueStartBin = testLeftCueOnBin(i) + timeRange(1)/glm.binSize + 1;
    cueEndBin = testLeftCueOnBin(i) + timeRange(2)/glm.binSize;
    if cueStartBin < 0; continue; end
    if cueEndBin > size(Y_test,2); continue; end

    % Sum PSTH for each neuron
    for n = 1:size(Y_test,1)
        psth_avg_fr_observed_l(n,:) = psth_avg_fr_observed_l(n,:) + ...
                                    Y_test(n,cueStartBin:cueEndBin);
        psth_avg_fr_predicted_l(n,:) = psth_avg_fr_predicted_l(n,:) + ...
                                    Y_predicted(n,cueStartBin:cueEndBin);
    end
    totalCueAnalyzed = totalCueAnalyzed + 1;
end
% Averaged PSTH for each neuron
psth_avg_fr_observed_l = psth_avg_fr_observed_l / totalCueAnalyzed;
psth_avg_fr_predicted_l = psth_avg_fr_predicted_l / totalCueAnalyzed;

% Plot left cue trials
% clusterList = ap.clusterToGoodClusterIndex(clusterList);
% set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
smoothed_fr_observed = smooth(psth_avg_fr_observed_l(neuron,:));
smoothed_fr_predicted = smooth(psth_avg_fr_predicted_l(neuron,:));
% [~,x_max] = max(abs(smoothed_fr));
t = linspace(timeRange(1),timeRange(2),size(psth_avg_fr_observed_l,2));
plot(t,smoothed_fr_observed,'-','LineWidth',1.5,'Color',colors(1)); hold on
plot(t,smoothed_fr_predicted,'-','LineWidth',1.5,'Color',colors(2)); hold on
drawnow

% Set plot related params
legend({'Observed','Predicted'});
xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Spikes/s'); % ylim([0,100]);
xline(0,'-','Left tone','Color','r','LineWidth',2,'HandleVisibility','off');
title(['Averged PSTH of left cue for cluster ',num2str(ap.goodClusters(neuron))]);
subtitle(['R^2: ',num2str(r2(neuron,2)),'  Total trial analyzed: ',num2str(length(testLeftCueOnBin))]);


% Analyze right trials
psth_avg_fr_predicted_r = zeros(size(Y_test,1),nBinsPerTrial);
psth_avg_fr_observed_r = zeros(size(Y_test,1),nBinsPerTrial);
totalCueAnalyzed = 0;
for i = 1:length(testRightCueOnBin)
    
    % Find start & end bin
    cueStartBin = testRightCueOnBin(i) + timeRange(1)/glm.binSize + 1;
    cueEndBin = testRightCueOnBin(i) + timeRange(2)/glm.binSize;
    if cueStartBin < 0; continue; end
    if cueEndBin > size(Y_test,2); continue; end
    
    % Sum PSTH for each neuron
    for n = 1:size(Y_test,1)
        psth_avg_fr_observed_r(n,:) = psth_avg_fr_observed_r(n,:) + ...
                                    Y_test(n,cueStartBin:cueEndBin);
        psth_avg_fr_predicted_r(n,:) = psth_avg_fr_predicted_r(n,:) + ...
                                    Y_predicted(n,cueStartBin:cueEndBin);
    end
    totalCueAnalyzed = totalCueAnalyzed + 1;
end
% Averaged PSTH for each neuron
psth_avg_fr_observed_r = psth_avg_fr_observed_r / totalCueAnalyzed;
psth_avg_fr_predicted_r = psth_avg_fr_predicted_r / totalCueAnalyzed;

% Plot right cue trials
figure;
% set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
smoothed_fr_observed = smooth(psth_avg_fr_observed_r(neuron,:));
smoothed_fr_predicted = smooth(psth_avg_fr_predicted_r(neuron,:));
% [~,x_max] = max(abs(smoothed_fr));
t = linspace(timeRange(1),timeRange(2),size(psth_avg_fr_observed_r,2));
plot(t,smoothed_fr_observed,'-','LineWidth',1.5,'Color',colors(1)); hold on
plot(t,smoothed_fr_predicted,'-','LineWidth',1.5,'Color',colors(2)); hold on
drawnow

% Set plot related params
legend({'Observed','Predicted'});
xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Spikes/s'); % ylim([0,100]);
xline(0,'-','Right tone','Color','r','LineWidth',2,'HandleVisibility','off');
title(['Averged PSTH of right cue for cluster ',num2str(ap.goodClusters(neuron))]);
subtitle(['R^2: ',num2str(r2(neuron,2)),'  Total trial analyzed: ',num2str(length(testRightCueOnBin))]);

%% Analyze GLM: coefficient PCA

% Run PCA for coefficient space
coefficient = modelSummary(:,1:nPredictors*length(offset));
% coefficient = zscore(modelSummary(:,1:9));
[~,score,~,~,explained,~] = pca(coefficient);

% Plot PCA results: variance
figure;
plot(explained,'-x','LineWidth',2,'Color',colors(1)); hold on
plot(cumsum(explained),'-x','LineWidth',2,'Color',colors(2));
xlim([0,25]); legend({'Individual explained variance','Cumulative explained variance'});
xlabel('Principal component'); ylabel('% total variance explained');

% Plot PCA results: score
figure;
for i = 1:size(modelSummary,1)
    scatter3(score(i,1),score(i,2),score(i,3),'filled');
    text(score(i,1),score(i,2),score(i,3),num2str(ap.goodClusters(i)));
    hold on
end
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

%% Analyze GLM: select example clusters using PCA score

% Find top 10 pca score 
exampleClusters_pca = zeros(10,6);
for i = 1:6
    % Select top 5 and lowest 5
    [~,Itop] = maxk(score(:,i),5);
    exampleClusters_pca(1:5,i) = Itop;
    [~,Ilow] = mink(score(:,i),5);
    exampleClusters_pca(6:end,i) = Ilow;
    
%     % Plot left cue trials
%     % clusterList = exampleClusters_pca(:,i);
%     set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
%     spikeRate = rightCueSpikeRate; eventIdx = find(allTrials==2);
%     yyaxis left; drawLicks('Left cue',timeRange,eventIdx,...
%     leftLick,rightLick,nidq,timeNI);
%     yyaxis right; drawSpikeRate('Left cue',timeRange,spikeRate,clusterList,textOn,ap,colors);
%
%     % Plot right cue trials
%     set(figure,'defaultAxesColorOrder',[[0.9290 0.6940 0.1250]; [0 0.4470 0.7410]]);
%     spikeRate = rightCueSpikeRate; eventIdx = find(allTrials==2);
%     yyaxis left; drawLicks('Right cue',timeRange,eventIdx,...
%     leftLick,rightLick,nidq,timeNI);
%     yyaxis right; drawSpikeRate('Right cue',timeRange,spikeRate,clusterList,textOn,ap,colors);
end

%% Analyze GLM: select clusters for R^2 > 0.4

r2_threshold = 0.4;

% Find above threshold clusters
highR2Clusters = find(r2(:,2) > r2_threshold);

% Set x axis
neurons = highR2Clusters;
x = linspace(offset(1),offset(end),length(offset));

figure;
for p = 1:nPredictors
    
    subplot(3,2,p);
    % subplot(2,2,p);
    outlier = [];
    
    for n = 1:size(neurons,1)
        x_n = zeros(length(offset),1);
        gcid = neurons(n);
        for i = 1:length(offset)
            % Averaging coefficient across neurons
            x_n(i) = modelSummary(gcid,(i-1)*nPredictors+p);
        end

        plot(x,x_n,'-','LineWidth',2); hold on
        [y_max,x_max] = max(abs(x_n));
        xtxtpos = offset(x_max);
        ytxtpos = x_n(x_max);
        cluster_id = ap.goodClusters(n);
        text(xtxtpos,ytxtpos,num2str(cluster_id)); hold on
    end
    
    xlabel('Offset time (s)'); ylabel('Coefficient estimates');
    title([predictorName{p},' coefficient']);
end

%% Analyze GLM: plot coefficient, deviance histogram

% Plot coefficient histogram
% figure;
% edges = [min(min(modelSummary)):5:max(max(modelSummary))];
% % histogram(modelSummary(:,1:nPredictors*length(offset)),edges);
% histogram(mean(modelSummary(:,1:nPredictors*length(offset)),2),edges);
% % xlim([-5000,5000]);

% Plot deviance histogram
figure;
% edges = [min(min(modelSummary)):5:max(max(modelSummary))];
% histogram(modelSummary(:,1:nPredictors*length(offset)),edges);
histogram(mean(modelSummary(:,1:nPredictors*length(offset)),2),edges);
% xlim([-5000,5000]);

%% Analyze GLM: identify unresonable clusters from modelSummary

for i = 1:size(modelSummary,1)
    coefficient_mean = mean(modelSummary(i,1:nPredictors*length(offset)));
    if abs(coefficient_mean) > 4000
        % modelSummary(i,:) = nan;
        disp(['Found cluster ',num2str(ap.goodClusters(i)),...
            ' which coefficient mean = ',num2str(coefficient_mean)]);
    end
end

% for i = 1:size(Y,1)
%     spikeCount = Y(i,:);
%     if ~isempty(find(spikeCount < 0, 1))
%         % disp(['Found negative spikes for neuron ',num2str(ap.goodClusters(i))]);
%         disp(['Found negative spikes for neuron ',num2str(i)]);
%     end
% end

% for i = 1:size(X,2)
%     design = X(:,i);
%     if ~isempty(design(rem(design,1)~=0))
%         % disp(['Found negative spikes for neuron ',num2str(ap.goodClusters(i))]);
%         disp(['Found non-integer values for predictor ',num2str(i)]);
%     end
% end

%% Save GLM data (for analysis in python)

save(strcat('glm_',sessionName),'bestModel');
% save(strcat('glm_',sessionName),'R2','RMSE','-append');

%% Read waveform data from temp_wh.dat

% Setup data structure
% session.pathImec = strcat('/Volumes/Shun neuro data/Neuropixel/',sessionName, '/', sessionName, '_imec0/');
% session.pathImec = strcat('/Volumes/T7/',sessionName, '/', sessionName, '_imec0/');
session.pathImec = strcat('D:\Shun\',sessionName, '\', sessionName, '_imec0\');
gwfparams.fileName = 'temp_wh.dat';     % .dat file containing the raw
% apD = dir(fullfile(session.pathImec, '*ap*.bin')); % AP band file from spikeGLX specifically
% gwfparams.fileName = apD(1).name;
gwfparams.dataType = 'int16';           % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 385;                    % Number of channels that were streamed to disk in .dat file
% nCh = 383 for dat file
gwfparams.nWf = 385;                    % Number of waveforms per unit to pull out

gwfparams.spikeTimes = spike_times;         % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = spike_clusters;   % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

% Load filtered waveform from temp_wh.dat and KiloSort/Phy output
fileName = fullfile(session.pathImec,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % Determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(fileName,'Format',{gwfparams.dataType, [gwfparams.nCh nSamp], 'x'},'Writable',true);
chMap = readNPY(fullfile(session.pathImec, 'channel_map.npy'))+1; % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

%% Find peri-lick AP waveforms

% Set up time variable
timeRangePerLick = [-0.05,0.05]; lickOnIdx = find(leftLick==1);
binSizePerLick = 0.001; nSampPerBin = floor(ap.Fs*binSizePerLick);
nBinsPerLick = (timeRangePerLick(2)-timeRangePerLick(1))/binSizePerLick;

lickWaveform = zeros(gwfparams.nWf,nBinsPerLick,length(lickOnIdx));

for i = 1:length(lickOnIdx)
    % Find first NI index
    niFirstIdx = lickOnIdx(i) + floor(timeRangePerLick(1)*nidq.Fs);
    % Find corresponding imec index
    [~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
    imecLastIdx = imecFirstIdx + nBinsPerLick*nSampPerBin;
    
    % Read waveform within range
    wf = mmf.Data.x(1:gwfparams.nCh, imecFirstIdx:imecLastIdx);
    wf_binned = downsample(wf',nSampPerBin)';
    % wf = wf - mean(wf,2); % V - mean of channel
    
    % Spike waveform for each lick
    wf_binned(:,length(wf_binned)-(length(wf_binned)-nBinsPerLick)+1:end) = [];
    lickWaveform(:,:,i)= wf_binned;
    
    disp(['Lick ',num2str(i),'/',num2str(length(lickOnIdx)),...
        ' (', num2str(i/length(lickOnIdx)*100),'%) is analyzed']);
end

%% Plot peri-lick AP waveforms averaged across licks

x = linspace(timeRangePerLick(1),timeRangePerLick(2),nBinsPerLick);

% Plot mean waveform per channel (averaged over licks)
lickWaveform_mean = mean(lickWaveform,3);
figure;
for i = 1:size(lickWaveform_mean,1)
    plot(x,lickWaveform_mean(i,:),'LineWidth',1); hold on
    % Plot channel number
    [y_max,x_max] = max(abs(lickWaveform_mean(i,:)));
    if abs(y_max) > 0.5
        xtxtpos = timeRangePerLick(1)+binSizePerLick*x_max;
        ytxtpos = lickWaveform_mean(i,x_max);
        text(xtxtpos,ytxtpos,num2str(i)); hold on
    end
end
xlabel('Time (s)'); ylabel('Voltage');
xline(0,'-','Lick','Color','r','LineWidth',1,'HandleVisibility','off');
title('Averaged peri-lick waveform per channel');
subtitle(['Time bin: ',num2str(binSizePerLick*1000),' ms; '...
          'Total licks analyzed: ',num2str(length(lickOnIdx)),'; '...
          'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
box off

%% Plot peri-lick waveform: all licks for 1 channel

ch = 20;
x = linspace(timeRangePerLick(1),timeRangePerLick(2),nBinsPerLick);

% Plot waveform for each individual lick for 1 channel
figure;
for i = 1:size(lickWaveform,3)
    plot(x,lickWaveform(ch,:,i),'LineWidth',1); hold on
end
plot(x,mean(lickWaveform(ch,:,:),3),'LineWidth',2,'Color',colors(1));
xlabel('Time (s)'); ylabel('Voltage'); %ylim([-5 5]);
xline(0,'-','Lick','Color','r','LineWidth',2,'HandleVisibility','off');
title(['Individual peri-lick waveform for channel ',num2str(ch)]);
subtitle(['Time bin: ',num2str(binSizePerLick*1000),' ms; '...
          'Total licks analyzed: ',num2str(length(lickOnIdx)),'; '...
          'Total time analyzed: ',num2str(session.nBlocks*session.blockTime),'s']);
box off

%% Find peri-stim AP waveforms


% Set up time variable
timeRange = [-0.01,0.01]; stimOnIdx = blueLaserON(701:800);
binSize = 1/ap.Fs; nSampPerBin = floor(ap.Fs*binSize);
nBins = (timeRange(2)-timeRange(1))/binSize;

stimWaveform = zeros(gwfparams.nCh,nBins,length(stimOnIdx));

for i = 1:length(stimOnIdx)
    % Find first NI index
    niFirstIdx = stimOnIdx(i) + floor(timeRange(1)*nidq.Fs);
    % Find corresponding imec index
    [~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
    imecLastIdx = imecFirstIdx + nBins*nSampPerBin;
    
    % Read waveform within range
    wf = mmf.Data.x(1:gwfparams.nCh, imecFirstIdx:imecLastIdx);
    % wf_binned = downsample(wf',nSampPerBin)';
    wf_binned = wf;
    % wf = wf - mean(wf,2); % V - mean of channel
    
    % Spike waveform for each stim
    wf_binned(:,length(wf_binned)-(length(wf_binned)-nBins)+1:end) = [];
    stimWaveform(:,:,i)= wf_binned;
    
    disp(['Stim ',num2str(i),'/',num2str(length(stimOnIdx)),...
        ' (', num2str(i/length(stimOnIdx)*100),'%) is analyzed']);
end

%% Plot peri-stim PSTH for each channel

x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);

% Plot mean waveform per channel
stimWaveform_mean = mean(stimWaveform,3);
% figure(1); max_value = max(max(stimWaveform_mean,[],1));
% for i = 1:size(stimWaveform_mean,1)
%     plot(x,i+(stimWaveform_mean(i,:)-stimWaveform_mean(1,:))/max_value,'LineWidth',1); hold on
% end
% xlabel('Time (s)'); ylabel('V_e');
% xline(0,'-','Stim','Color','r','LineWidth',1,'HandleVisibility','off');
% box off

figure(2); imagesc('XData',x,'YData',1:size(stimWaveform,1),'CData',stimWaveform_mean); 
colormap(colMap); colorbar; %caxis([-200 200]);
xlabel('Time (ms)'); ylabel('Channel'); 
xlim([x(1) x(end)]); ylim([1 size(stimWaveform,1)]);
xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
% patch('XData',[0 20*0.001],'YData',[size(stimWaveform,1),size(stimWaveform,1)],...
%       'FaceColor',colors(1),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
title('1mW 20ms pulse @ 1Hz');
box off

%% ************************ Optotag ******************************
%{
1. Use laser shuffle: type laser_shuffle in command window

%}

%% Generate input for laser_shuffle

table_name = strcat('Optotag/optotag_',sessionName,'_5mW_10ms.csv');

% Generate header row
header = cell(1,ap.nGoodClusters+1); header{1,1} = 'Laser_1';
for i = 1:ap.nGoodClusters
    cluster_header = strcat('sig',num2str(ap.goodClusters(i)));
    header{1,i+1} = cluster_header;
end
writecell(header,table_name);

% Generate event times column
spikeTable = cell(1,ap.nGoodClusters+1);
laserTime = timeNI(blueLaserON)';
spikeTable{1,1} = laserTime(1:100);

% Create spikeTable (spike timeImec of individual good cluster at each column)
for i = 1:ap.nGoodClusters
    cluster_id = ap.goodClusters(i);
    [row,~] = find(ap.goodSpikeClusters==cluster_id);
    spikeTimes = timeImec(ap.goodSpikeTimes(row))';
    spikeTable{1,i+1} = spikeTimes;
end

% Make other NaN
l = cellfun(@length,spikeTable);          % return length of each cell in cell array c
L = max(l);                               % and find the longest
n = arrayfun(@(l) nan(L-l,1),l,'uni',0);  % create vector of NaN to make each equal length 
for i = 1:size(spikeTable,2)              % and make the array uniform
  spikeTable(i)={[spikeTable{i};n{i}]};
end
writematrix(cell2mat(spikeTable),table_name,'WriteMode','append');
disp('laser_shuffle input generated');

%% Generate laser-triggered firing rate for each good cluster

% 20220523
% laser_pulse_duration(1:300) = 0.01;
% laser_pulse_duration(301:600) = 0.05;
% laser_pulse_duration(601:900) = 0.02;
% stim_per_pattern = 100; nPatterns = 9; % npulse_per_stim = 1; 

% 20220525
stim_per_pattern = 200; nPatterns = 5; % npulse_per_stim = 1; 

% Set up parameters
% timeRange = [-0.005,0.03]; 
timeRange = [-1, 10];
binSize = 0.005; % in sec
blueLaserON = find(blueLaser);
nStim = nPatterns * stim_per_pattern;
blueLaserOnIdx = blueLaserON(201:400);
laser_pulse_duration = zeros(1,nStim);

% Analysis window based on high variance
high_var_range = [0.005,0.025]; % in sec

% Generate PSTH
[optoSpikes,optoSpikeRate] = getSpikes(timeRange,binSize,blueLaserOnIdx,...
                                    ap,nidq,timeNI);

% Plot PSTH for each neuron across trials
clusterList = 1:ap.nGoodClusters; textOn = true;
figure; drawSpikeRate('Stim',timeRange,optoSpikeRate,clusterList,textOn,ap,colors);

% Example artifact unit
% 20220523: artifact: 56, 318
% 20220525: artifact: 741, 64; others: 225, 23 
clusterList = ap.clusterToGoodClusterIndex(64); textOn = true;
figure; drawSpikeRate('Stim',timeRange,optoSpikeRate,clusterList,textOn,ap,colors);

%% Plot distribution of first spike latency after stim 

timeRange = [0 0.2]; binSize = 0.005;
npulse_per_stim = 1; stim_per_pattern = 100; nPatterns = 9;
nStim = npulse_per_stim * stim_per_pattern * nPatterns;
laserOnset = find(redLaser);

% Intialize first_spike_latency: row is neurons, column is stim num
first_spike_latency = nan(ap.nGoodClusters,nStim);
nBins = (timeRange(2)-timeRange(1)) / (1/binSize);

% Looping over each stim
for i = 1:nStim
    spikes = getSpikes(timeRange,binSize,laserOnset(i),ap,nidq,timeNI);
    for neuron = 1:ap.nGoodClusters
        first_spike_bin = find(spikes(neuron,:),1);
        if isempty(first_spike_bin); continue; end
        first_spike_latency(neuron,i) = binSize * first_spike_bin;
    end
end

% Convert to ms
first_spike_latency = first_spike_latency .* 1000;

% Plot histogram distribution of first_spike_latency
ans = mean(first_spike_latency,2,'omitnan');
figure; histogram(mean(first_spike_latency,2,'omitnan'),50,'FaceColor',colors(1));
xlabel('First spike latency (ms)'); ylabel('Number of neurons');
box off

%{
% For analyzing each pattern
% for i = 701 % test one pattern only
for i = 1:stim_per_pattern:nStim % 1, 101, 201, ..., 801
    stimOnIdx = laserOnset(i:i+stim_per_pattern-1);
    stimWaveform = zeros(gwfparams.nCh,nBins,length(stimOnIdx));
    
    % Generate stimWaveform
    for j = 1:length(stimOnIdx)
        % Find first NI index
        niFirstIdx = stimOnIdx(j) + floor(timeRange(1)*nidq.Fs);
        [~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
        imecLastIdx = imecFirstIdx + nBins-1;
        % disp([imecFirstIdx,stimOnIdx(j),imecLastIdx]);

        % Read waveform within range
        wf = mmf.Data.x(1:gwfparams.nCh, imecFirstIdx:imecLastIdx);

        % Spike waveform for each stim
        stimWaveform(:,:,j)= wf;

        disp(['Stim ',num2str(i-1+j),'/',num2str(nStim),...
            ' (', num2str(((i-1+j)/nStim)*100),'%) is analyzed']);
    end
    
    % Store average wvf to artifact_wvf
    cur_pattern = ((i-1)/stim_per_pattern) + 1;
    artifact_wvf{cur_pattern} = mean(stimWaveform,3);
    
end
%}

%% Helper function

% ======================================= %
% function [spike,spikeRate] = getSpikes(timeRange,binSize,eventIdx,...
%                 ap,nidq,timeNI)
% 
% % The number of spikes for all good units at each time bin (binSize) within
% % time before/after the event (timeRange)
% 
% spike = zeros(ap.nGoodClusters, round((timeRange(2)-timeRange(1))/binSize));
% for i = 1:length(eventIdx)
%     % Find first & last NI index
%     niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidq.Fs);
%     % Find corresponding imec index
%     imecFirstIdx = round(timeNI(niFirstIdx)*ap.Fs); % make sure timeImec and timeNI is aligned to 0
%     % [~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
%     imecLastIdx = imecFirstIdx + floor(ap.Fs*(timeRange(2)-timeRange(1)));
% 
%     % Find spike within timeRange
%     spikeTimesinRangeIdx = find...
%         (ap.goodSpikeTimes>imecFirstIdx & ap.goodSpikeTimes<imecLastIdx);
%     spikeTimesinRange = ap.goodSpikeTimes(spikeTimesinRangeIdx);
%     spikeClustersinRange = ap.goodSpikeClusters(spikeTimesinRangeIdx);
%     disp(['Found ',num2str(length(spikeTimesinRange)),' spikes for event ', num2str(i)]);
% 
%     % Generate spike raster array
%     for s = 1:length(spikeTimesinRange)
%         % Find which neuron fires the spike
%         cluster_id = spikeClustersinRange(s);
%         goodCluster_id = ap.clusterToGoodClusterIndex(cluster_id);
%         
%         % Find the timebin of the spike
%         relativeSpikeSamp = spikeTimesinRange(s)-imecFirstIdx;
%         relativeSpikeTime = double(relativeSpikeSamp)/ap.Fs;
%         relativeSpikeBin = floor(relativeSpikeTime/binSize)+1;
%         
%         % Add spike to array
%         spike(goodCluster_id,relativeSpikeBin) = 1 + spike(goodCluster_id,relativeSpikeBin);
%     end
% end
% spikeRate = spike / (binSize*length(eventIdx));
% 
% end % getSpikes

% ======================================= %
% function [] = drawLicks(plotLeft,timeRange,eventIdx,...
%                 leftLick,rightLick,nidq,timeNI)
% 
% % Find licking events
% leftLickOnIdx = find(leftLick==1);
% rightLickOnIdx = find(rightLick==1);
% 
% for i = 1:length(eventIdx)
%     % Find first & last NI index
%     niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidq.Fs);
%     niLastIdx = eventIdx(i) + floor(timeRange(2)*nidq.Fs);
% 
%     % Find licks within timeRange
%     leftLickTimesinRange = leftLickOnIdx(leftLickOnIdx>niFirstIdx & leftLickOnIdx<niLastIdx);
%     rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
%     relativeLeftLickTime = timeNI(leftLickTimesinRange)-timeNI(eventIdx(i));
%     relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i));
%     
%     % Plot lick raster plot
%     scatter(relativeLeftLickTime,i,'filled','MarkerFaceColor','#6DBAA1'); hold on
%     scatter(relativeRightLickTime,i,'filled','MarkerFaceColor','#FFC25C'); hold on
%     drawnow
% end
% 
% xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
% ylabel('Trial'); ylim([0,length(eventIdx)]);
% if plotLeft
%     xline(0,'-','Left tone','Color','r','LineWidth',2,'HandleVisibility','off');
%     title(['Left cue-triggered licking']);
%     subtitle(['Total cue analyzed: ',num2str(length(eventIdx))]);
% else 
%     xline(0,'-','Right tone','Color','r','LineWidth',2,'HandleVisibility','off');
%     title(['Right cue-triggered licking']);
%     subtitle(['Total cue analyzed: ',num2str(length(eventIdx))]);
% end
% 
% end % drawLicks

% ======================================= %
% function [] = drawSpikeRate(plotLeft,timeRange,...
%                 spikeRate,clusterList,textOn,ap,color)
% 
% % Set colormap
% cmap = colormap(flipud(jet(length(clusterList))));
% 
% % Plot individual traces
% for i = 1:length(clusterList)
%     neuron = clusterList(i);
%     if sum(spikeRate(neuron,:)) > 10
%         smoothed_fr = smooth(spikeRate(neuron,:));
%         [~,x_max] = max(abs(smoothed_fr));
%         t = linspace(timeRange(1),timeRange(2),size(spikeRate,2));
%         if strcmp(color,'map')
%             plot(t,smoothed_fr,'-','LineWidth',1.5); hold on
%         else
%             plot(t,smoothed_fr,'-','LineWidth',1.5,'Color',color); hold on
%         end
%         drawnow
%         
%         if textOn
%             xtxtpos = t(x_max);
%             ytxtpos = smoothed_fr(x_max);
%             text(xtxtpos,ytxtpos,num2str(ap.goodClusters(neuron))); hold on
%         end
%     end
% end
% 
% % Set plot related params
% set(gca,'ColorOrder',cmap);
% xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
% ylabel('Spikes/s'); % ylim([0,100]);
% if plotLeft
%     xline(0,'-','Left tone','Color','r','LineWidth',2);
%     title(['Averged PSTH of left cue per neuron']);
% else 
%     xline(0,'-','Right tone','Color','r','LineWidth',2);
%     title(['Averged PSTH of right cue per neuron']);
% end
% 
% end % drawSpikeRate

% ======================================= %
function [rmse,model] = runGLM(X_bl,neuron,glm,ap,returnModel)
    
    % STEP 1: Set up X_baseline based on current testing model
    nPredictors = size(X_bl,2);
    % disp(['Start analyzing GLM model']);
    
    % STEP 2: Create response matrix Y
    cluster_id = ap.goodClusters(neuron);
    y_spike = zeros(glm.nBins,1); % Set up spike (response vector)
    spikeOnBin = floor(double(ap.goodSpikeTimes(ap.goodSpikeClusters==cluster_id))/glm.nSampPerBinImec)+1;
    for j = 1:length(spikeOnBin)
        if spikeOnBin(j) <= glm.nBins
            y_spike(spikeOnBin(j)) = 1 + y_spike(spikeOnBin(j)); % Spike count
        end
    end
    y_fr = y_spike/glm.binSize; % Calculate firing rate
    
    % STEP 3: Add offset padding for behavior variables
    offset = glm.offset;
    X = []; % Combine all predictor data into one matrix
    for offsetIdx = 1:length(offset)

        X_shifted = X_bl;
        offsetInBin = uint32(abs(offset(offsetIdx))/glm.binSize); 
        padding = zeros(offsetInBin,nPredictors);

        % Shifts behavioral data based on offset time
        if offset(offsetIdx) > 0
            X_shifted(glm.nBins-offsetInBin+1:end,:) = [];
            X_shifted = [padding;X_shifted];
        elseif offset(offsetIdx) < 0
            X_shifted(1:offsetInBin,:) = [];
            X_shifted = [X_shifted;padding];
        end

        % Combine time shifted data with baseline
        % disp(['     STEP 3: Combining behavioral data shifted by ',num2str(offset(offsetIdx)*1000),' ms']);
        X = [X, X_shifted];
    end
    
    % STEP 3.5: include spiking history to GLM
    % Equal to the number of offset, but only concerns history
    % i.e. if offset is 5 bins forward/backwards, include spike history 10
    % bins backwards
%     spikeHistory = [];
%     for i = 1:length(offset)
%         Y_shifted = y_fr; 
%         padding = zeros(i,1);
%         Y_shifted(1:i,1) = [];
%         Y_shifted = [Y_shifted;padding];
%         spikeHistory = [Y_shifted spikeHistory];
%     end
%     X = [X spikeHistory];
    
    
    % STEP 4: Separate data into testing and training set
    % 80% training, 20% testing
    totalTestBins = 0.2*glm.nBins; nTestChunk = 5;
    nBinsPerChunk = totalTestBins/nTestChunk; % 960 Bins/chunk

    % Generate test chunk index
    cueBins = find(glm.combine_cue > 0);
    resplit = true; nSplitAttempt = 0;
    while resplit
        % testing chunk ends 1 second before the next tone
        chunkEndBin = sort(cueBins(randperm(length(cueBins),nTestChunk))-20);
        nSplitAttempt = nSplitAttempt + 1;
        flag = 0;
        for i = 1:length(chunkEndBin)
            if chunkEndBin(i) - nBinsPerChunk < 0; flag=1; break
            elseif chunkEndBin(i) - nBinsPerChunk > 0
                if i == 1; continue
                elseif chunkEndBin(i) - nBinsPerChunk < chunkEndBin(i-1); flag=2; break
                end
            end
        end
        if flag == 0
            resplit = false;
            % disp(['     STEP 4: Found testing partitions in ',num2str(nSplitAttempt),' attempts']);
        end
    end

    % Genearte test data index
    chunkStartBin = chunkEndBin - nBinsPerChunk + 1;
    binIdx = linspace(1,glm.nBins,glm.nBins)'; chunkIdx = [];
    for idx = 1:nTestChunk
        chunkIdx = [chunkIdx;binIdx(chunkStartBin(idx):chunkEndBin(idx))];
    end

    % Separate test and training data
    X_training = X(setdiff(binIdx,chunkIdx),:);
    Y_training = y_fr(setdiff(binIdx,chunkIdx));
    
    % STEP 5: 10-fold CV to calculate averaged RMSE for each neuron
    % Prediction function given training/testing instances
    fcn = @(Xtr, Ytr, Xte) predict(...
        GeneralizedLinearModel.fit(Xtr,Ytr,'linear','distr','normal'), ...
        Xte);
    % Perform cross-validation, and return average MSE across folds
    mse = crossval('mse', X_training, Y_training, 'Predfun',fcn, 'kfold',10);
    % Compute root mean squared error
    rmse = sqrt(mse);
    disp(['     STEP 5: RMSE of current GLM model for neuron ',num2str(neuron), ...
          ' is ',num2str(rmse)]);
      
    % STEP 6: return model if returnModel is true
    if returnModel
        model = fitglm(X_training,Y_training,'linear','Distribution','normal');
    end
    
end % runGLM

% ======================================= %
function [best_model,min_rmse] = backSelect(best_model,min_rmse,neuron,...
                                            X_baseline,glm,ap)
% INPUT/OUTPUT:
    % best_model: an array storing predictorList for the best model 
    % (e.g. [1 2 4 5 6])
    % min_rmse: rmse score of the best model
    % X_baseline: baseline X including all variable
% DESCRIPTION: This function implements backward selection process for multiple GLM.
    % 1. The best model inputed (the full model is the best model for the
    % first iteration).
    % 2. backSelect calculate RMSE for the input best model, and also
    % models removing one predictors from best model.
    % 3. If the best model (with lowest RMSE) is the input best model, this
    % means that all predictors are important, therefore the selection
    % process ends.
    % 4. If the best model is the nPredictor-1 model, then the process runs
    % iteratively by calling another backSelect using the nPredictor-1
    % model as the best model.
    % 5. If two models are equal in rmse, then the same process is rerun
    % again to separete between two models.
    
nPredictors = length(best_model);
input_model = best_model; % Store input model separately (might changed when doing CV for nPredictor-1 models)
input_rmse = min_rmse;

if nPredictors == 1
    % Run with no predictor values (zero array)
    X_intercept = zeros(glm.nBins,1);
    rmse = runGLM(X_intercept,neuron,glm,ap,false);
    if rmse < min_rmse
        min_rmse = rmse;
        best_model = 0;
        disp(['Found best model for neuron ',num2str(neuron),...
              ' which RMSE = ',num2str(min_rmse),...
              ' using predictorList = [',num2str(best_model),']']);
        return
    else
        disp(['Found best model for neuron ',num2str(neuron),...
              ' which RMSE = ',num2str(min_rmse),...
              ' using predictorList = [',num2str(best_model),']']);
        return
    end
else
    % Run nPredictor-1 model
    for i = 1:nPredictors
        % Remove selected variable
        predictorList = input_model;
        disp(['Removing predictor ',num2str(predictorList(i)),...
              ' from predictorList = [',num2str(predictorList),']']);
        predictorList(i) = [];
        disp(['predictorList = [',num2str(predictorList),']']);
        
        % Run CV to obtain rmse
        X_remove = X_baseline(:,predictorList);
        rmse = runGLM(X_remove,neuron,glm,ap,false);
        if rmse < min_rmse
            min_rmse = rmse;
            best_model = predictorList;
            disp(['Current best model updated: min_rmse = ',num2str(min_rmse),...
                  '; predictorList = [',num2str(predictorList),']']);
        end
    end
    
    % Run the input model
    disp(['Testing input model with predictorList = [',num2str(input_model),']']);
    X_input = X_baseline(:,input_model);
    rmse = runGLM(X_input,neuron,glm,ap,false);
    if rmse < min_rmse
        min_rmse = rmse;
        best_model = input_model;
        disp(['Found best model for neuron ',num2str(neuron),...
              ' which RMSE = ',num2str(min_rmse),...
              ' using predictorList = [',num2str(best_model),']']);
        return
    elseif min_rmse == input_rmse
        best_model = input_model;
        disp(['Found best model for neuron ',num2str(neuron),...
              ' which RMSE = ',num2str(min_rmse),...
              ' using predictorList = [',num2str(best_model),']']);
        return
    else
        [best_model,min_rmse] = backSelect(best_model,min_rmse,neuron,X_baseline,glm,ap);
    end
end

end % backSelect
