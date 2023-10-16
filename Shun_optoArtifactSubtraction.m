%% For mean subtraction of each channel before kilosort
%{
1. Read raw data file, calculate mean waveform for each channel at laser
time pulse

2. Subtract raw signal with the mean waveform

3. Run kilosort
%}

%% Read waveform data from ap.bin

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed] = loadColors;

% Setup data structure
% sessionName = '20220523-SJ518-L_g0';
% sessionName = '20220525-SJ518-L_g0';
% sessionName = '20220612-SJ518-R_g0';
sessionName = '20220613-SJ518-R_g0';
load(strcat('sync_',sessionName,'.mat'),'blueLaser','redLaser','timeImec','timeNI'); % red laser for 0523 sessions

% session.pathImec = strcat('/Volumes/Shun neuro data/Neuropixel/',sessionName, '/', sessionName, '_imec0/');
session.pathNidq = strcat('D:\Shun\Recordings\',sessionName, '\');
session.pathImec = strcat('D:\Shun\Recordings\',sessionName, '\catgt_', sessionName, '\');
session.apBin = strcat(sessionName,'_tcat.imec0.ap.bin');
session.nidqBin = strcat(sessionName,'_t0.nidq.bin');
ap.meta = ReadMeta(session.apBin, session.pathImec);
ap.Fs = str2double(ap.meta.imSampRate);
nidq.meta = ReadMeta(session.nidqBin, session.pathNidq);
nidq.Fs = str2double(nidq.meta.niSampRate);

gwfparams.fileName = session.apBin; gwfparams.nCh = 385;
% gwfparams.fileName = 'temp_wh.dat'; gwfparams.nCh = 383;
gwfparams.dataType = 'int16';           % Data type of .dat file (this should be BP filtered)
gwfparams.nWf = 385;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes = spike_times;         % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = spike_clusters;   % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

% Load waveform
fileName = fullfile(session.pathImec,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % Determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'},'Writable',true);
chMap = readNPY(fullfile(session.pathImec, 'channel_map.npy'))+1; % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

disp([gwfparams.fileName,' for session ',sessionName,' loaded']);

%% Set up session stimulation pattern

% 20220523
% laser_pulse_duration = zeros(1,length(laserOnset));
% laser_pulse_duration(1:300) = 0.01;
% laser_pulse_duration(301:600) = 0.05;
% laser_pulse_duration(601:900) = 0.02;
% npulse_per_stim = 1; stim_per_pattern = 100; nPatterns = 9;
% Plot titles for 20220523
% titles = {'5mW 10ms pulse @ 1Hz','1mW 10ms pulse @ 1Hz','100uW 10ms pulse @ 1Hz',...
%           '5mW 50ms pulse @ 1Hz','1mW 50ms pulse @ 1Hz','100uW 50ms pulse @ 1Hz',...
%           '5mW 20ms pulse @ 1Hz','1mW 20ms pulse @ 1Hz','100uW 20ms pulse @ 1Hz'};

% 20220525
% npulse_per_stim = 1; stim_per_pattern = 200; nPatterns = 5;
% Plot titles for 20220525
% titles = {'1mW 20ms pulse @ 1Hz','300uW 20ms pulse @ 1Hz',...
%           '100uW 20ms pulse @ 1Hz','600uW 20ms pulse @ 1Hz',...
%           '0uW 20ms pulse @ 1Hz'};

% 20220612
% npulse_per_stim = 1; stim_per_pattern = 50; nPatterns = 12;
% titles = {'25uW','50uW','75uW','100uW',...
%           '125uW','150uW','175uW','200uW',...
%           '225uW','250uW','275uW','300uW'};

% 20220613
npulse_per_stim = 1; stim_per_pattern = 50; nPatterns = 17;
titles = {'25uW','50uW','75uW','100uW',...
          '125uW','150uW','175uW','200uW',...
          '225uW','250uW','275uW','300uW',...
          '0uW','50uW','100uW','150uW','200uW'};
% titles= {'0uW','50uW','100uW','150uW','200uW'};

%% (1) Get raw waveform for each pattern

timeRange = [-0.01 0.030]; % time range of mean subtraction
nStim = npulse_per_stim * stim_per_pattern * nPatterns;
% laserOnset = find(redLaser);
laserOnset = find(blueLaser);

% Generate average stim waveform
artifact_wvf = cell(3,nPatterns);
% nBins = double((timeRange(2)-timeRange(1))*ap.Fs);
nBins = (timeRange(2)-timeRange(1)) / (1/ap.Fs);

% for i = 201:stim_per_pattern:400 % test one pattern only
for i = 1:stim_per_pattern:nStim % 1, 101, 201, ..., 801
    stimOnIdx = laserOnset(i:i+stim_per_pattern-1);
    stimWaveform = zeros(gwfparams.nCh,nBins,length(stimOnIdx));
    
    % Generate stimWaveform
    for j = 1:length(stimOnIdx)

        %{
        Finds the wrong imec idx (not sure why)
        Can test this by type following in command window:
            format long g
            timeNI(stimOnIdx(j))
            timeImec(imecFirstIdx+300)

        % Find first Imec index
        [~, imecStimIdx] = min(abs(timeImec-timeNI(stimOnIdx(j))));
        % Find corresponding imec index
        imecFirstIdx = imecStimIdx - floor(timeRange(1)*ap.Fs);
        %}
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
    
    % Calculate variance of waveform for each channel at each time point
    % across all pulses
    

    % Store average wvf to artifact_wvf
    cur_pattern = ((i-1)/stim_per_pattern) + 1;
    artifact_wvf{1,cur_pattern} = stimWaveform;
    artifact_wvf{2,cur_pattern} = mean(stimWaveform,3);
    artifact_wvf{3,cur_pattern} = var(stimWaveform,0,3);
    
end

%% (1) Alternative: load preanalyzed artifact_wvf

sessionName = '20220525-SJ518-L_g0';
timeRange = [-0.01 0.030]; % time range of mean subtraction
nBins = (timeRange(2)-timeRange(1)) / (1/ap.Fs);
load(strcat('artifact_wvf_',sessionName));
npulse_per_stim = 1; stim_per_pattern = 200; nPatterns = 5;
% Plot titles for 20220525
titles = {'1mW 20ms pulse @ 1Hz','300uW 20ms pulse @ 1Hz',...
          '100uW 20ms pulse @ 1Hz','600uW 20ms pulse @ 1Hz',...
          '0uW 20ms pulse @ 1Hz'};

%% (1) Plot mean waveform
x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);

for i = 1:nPatterns
    % subplot(3,3,i);
    % subplot(2,3,i); 
    subplot(4,5,i);
    temp_wvf = artifact_wvf{2,i}; 
    % cmax = max(max(abs(temp_wvf([1:end-1],:))));
    cmax = max(max(abs(temp_wvf([1:191,193:end-1],:))));
    imagesc('XData',x,'YData',1:size(artifact_wvf{2,i},1),'CData',temp_wvf);
    colormap(blueWhiteRed); colorbar; caxis([-cmax cmax]);
    xlabel('Time (ms)'); ylabel('Channel');
    xlim([x(1) x(end)]); ylim([1 size(stimWaveform,1)]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    % patch('XData',[0 20*0.001],'YData',[size(stimWaveform,1),size(stimWaveform,1)],...
    %       'FaceColor',colors(1),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    title(titles(i));
    box off
end
autoArrangeFigures();

%% (1) Plot mean variance of each channel across stimulation

x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);
for i = 1:nPatterns
    %subplot(3,3,i);
    % subplot(2,3,i); 
    subplot(4,5,i);
    temp_wvf = artifact_wvf{3,i};
    cmax = max(max(abs(temp_wvf([1:191,193:end-1],:))));
    imagesc('XData',x,'YData',1:size(artifact_wvf{2,i},1),'CData',temp_wvf);
    colormap("parula"); colorbar; caxis([0 cmax]);
    xlabel('Time (ms)'); ylabel('Channel');
    xlim([x(1) x(end)]); ylim([1 size(stimWaveform,1)]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    % patch('XData',[0 20*0.001],'YData',[size(stimWaveform,1),size(stimWaveform,1)],...
    %       'FaceColor',colors(1),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    title(titles(i));
    box off
end

%% (1) Plot single channel waveform for all stim
channel = 102;

x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);
figure;
for i = 1:nPatterns
    %subplot(3,3,i);
    %subplot(2,3,i);
    subplot(4,5,i);
    temp_wvf = artifact_wvf{1,i}(channel,:,:); cmax = max(max(abs(temp_wvf)));
    imagesc('XData',x,'YData',1:size(artifact_wvf{1,i},3),'CData',squeeze(temp_wvf)');
    colormap(blueWhiteRed); colorbar; caxis([-cmax cmax]);
    xlabel('Time (ms)'); ylabel('Stimulation pulse');
    xlim([x(1) x(end)]); ylim([1 size(temp_wvf,3)]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    title(titles(i));
    box off
end

figure;
for i = 1:nPatterns
    % subplot(3,3,i);
    % subplot(2,3,i);
    figure(i+1);
    temp_wvf = artifact_wvf{1,i}(channel,:,:); cmax = max(max(abs(temp_wvf)));
    for j = 1:size(temp_wvf,3)
        y = (temp_wvf(1,:,j) / cmax) + j;
        plot(x,y,'Color','k'); hold on
    end
    xlabel('Time (ms)'); ylabel('Stimulation pulse');
    xlim([x(1) x(end)]); ylim([1 size(temp_wvf,3)+1]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    title(titles(i));
    box off
end
autoArrangeFigures();

%% (1) Plot same channel diagnostics at different power level

channel = 100;
x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);
subplot(1,3,1);
for i = 1:nPatterns
    plot(x,artifact_wvf{3,i}(channel,:),"LineWidth",1,"Color",blueRedYellow(i)); hold on
end
xlabel('Time (ms)'); ylabel('Variance');
legend(titles); box off

subplot(1,3,2);
for i = 1:nPatterns
    % plot(x,artifact_wvf{3,i}(channel,:)./artifact_wvf{2,i}(channel,:),"LineWidth",1,"Color",colors(i)); hold on
    % std_wvf = std(squeeze(artifact_wvf{1,i}(channel,:,:)),0,2);
    % plot(x,std_wvf,"LineWidth",1,"Color",colors(i)); hold on
    plot(x,artifact_wvf{2,i}(channel,:),"LineWidth",1,"Color",blueRedYellow(i)); hold on
end
xlabel('Time (ms)'); ylabel('Mean');
legend(titles); box off

subplot(1,3,3);
for i = 1:nPatterns
    % plot(abs(artifact_wvf{2,i}(channel,:)),artifact_wvf{3,i}(channel,:),"LineWidth",1,"Color",colors(i)); hold on
    scatter(abs(artifact_wvf{2,i}(channel,:)),artifact_wvf{3,i}(channel,:),'.',"Color",blueRedYellow(i)); hold on
end
xlabel('|Mean|'); ylabel('Variance');
legend(titles); box off

%% (2) Subtract average waveform from original data and blank out 300us data

corrected_wvf = cell(3,nPatterns);
eventSamp = abs(timeRange(1))*ap.Fs;
blankOutTime = 0.0003; blankOutSamp = blankOutTime*ap.Fs; % 300us
% openingTimeSamp = 0.0017*ap.Fs;

temp_wf = zeros(gwfparams.nCh,nBins);
% for i = 201
% for i = 201:stim_per_pattern:400
for i = 1:stim_per_pattern:nStim
    stimOnIdx = laserOnset(i:i+stim_per_pattern-1);
    avg_wvf = artifact_wvf{2,((i-1)/stim_per_pattern) + 1};
    correctedWaveform = zeros(gwfparams.nCh,nBins,length(stimOnIdx));
    disp(['Correcting pattern ',num2str(((i-1)/stim_per_pattern) + 1)]);
    
    for j = 1:length(stimOnIdx)

        niFirstIdx = stimOnIdx(j) + floor(timeRange(1)*nidq.Fs);
        [~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
        imecLastIdx = imecFirstIdx + nBins-1;

        % Read waveform within range
        wf = double(mmf.Data.x(1:gwfparams.nCh,imecFirstIdx:imecLastIdx));
        original_wf = wf;
        % figure(1); imagesc(wf); colorbar; title('Raw');
        % xline(abs(timeRange(1)*ap.Fs),'-','Color','r','LineWidth',2);
        % figure(2); imagesc(avg_wvf); colorbar; title('Avg');
        % xline(abs(timeRange(1)*ap.Fs),'-','Color','r','LineWidth',2);

        % Raw waveform - avg wvf for each channel
        wf = bsxfun(@minus,wf,avg_wvf);
        % figure(3); imagesc(wf); colorbar; title('Raw-avg');
        % xline(abs(timeRange(1)*ap.Fs),'-','Color','r','LineWidth',2);

        % P1: Blank out 300us: ignore shutter opening
        wf(:,eventSamp:eventSamp+blankOutSamp) = 0;
        % figure(4); imagesc(wf); colorbar; title('Blank out');
        % xline(abs(timeRange(1)*ap.Fs),'-','Color','r','LineWidth',2);

        % Add baseline (avg before shutter signal)
        wf = bsxfun(@plus,wf,mean(avg_wvf(:,1:abs(floor(timeRange(1)*ap.Fs))),2));
        % figure(5); imagesc(wf); colorbar; title('Added baselne');
        % xline(abs(timeRange(1)*ap.Fs),'-','Color','r','LineWidth',2);
        temp_wf = (temp_wf + wf)/j;
        % figure(6); imagesc(temp_wf); colorbar; title('Avg corrected');
        % xline(abs(timeRange(1)*ap.Fs),'-','Color','r','LineWidth',2);
        % autoArrangeFigures();

        % !!!Dangerous territory!!!: modify ap.bin
        %mmf.Data.x(1:gwfparams.nCh,imecFirstIdx:imecLastIdx) = wf;

        correctedWaveform(:,:,j)= wf;
        disp(['Stim ',num2str(i-1+j),'/',num2str(nStim),...
            ' (', num2str(((i-1+j)/nStim)*100),'%) is corrected']);

    end
    
    % Store corrected wvf to corrected_wvf
    cur_pattern = ((i-1)/stim_per_pattern) + 1;
    corrected_wvf{1,cur_pattern} = correctedWaveform;
    corrected_wvf{2,cur_pattern} = mean(correctedWaveform,3);
    corrected_wvf{3,cur_pattern} = var(correctedWaveform,0,3);
    
end

%% (2) Plot averaged corrected waveform
x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);

for i = 1:nPatterns
    subplot(2,3,i); temp_wvf = corrected_wvf{2,1};
    cmax = max(max(abs(temp_wvf([1:191,193:end-1],:))));
    imagesc('XData',x,'YData',1:size(corrected_wvf{2,i},1),'CData',temp_wvf);
    colormap(blueWhiteRed); colorbar; caxis([-cmax cmax]);
    xlabel('Time (ms)'); ylabel('Channel');
    xlim([x(1) x(end)]); ylim([1 size(correctedWaveform,1)]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    % patch('XData',[0 20*0.001],'YData',[size(stimWaveform,1),size(stimWaveform,1)],...
    %       'FaceColor',colors(1),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    title(titles(i));
    box off
end
autoArrangeFigures();

%% (2) Plot single channel waveform for all stim
channel = 123;

x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);
figure;
for i = 1:nPatterns
    %subplot(3,3,i);
    subplot(2,3,i);
    temp_wvf = corrected_wvf{1,i}(channel,:,:); cmax = max(max(abs(temp_wvf)));
    imagesc('XData',x,'YData',1:size(temp_wvf,3),'CData',squeeze(temp_wvf)');
    colormap(blueWhiteRed); colorbar; caxis([-cmax cmax]);
    xlabel('Time (ms)'); ylabel('Stimulation pulse');
    xlim([x(1) x(end)]); ylim([1 size(temp_wvf,3)]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    title(titles(i));
    box off
end

figure;
for i = 1:nPatterns
    % subplot(3,3,i);
    % subplot(2,3,i);
    figure(i+1);
    temp_wvf = corrected_wvf{1,i}(channel,:,:); cmax = max(max(abs(temp_wvf)));
    for j = 1:size(temp_wvf,3)
        y = (temp_wvf(1,:,j) / cmax) + j;
        plot(x,y,'Color','k'); hold on
    end
    xlabel('Time (ms)'); ylabel('Stimulation pulse');
    xlim([x(1) x(end)]); ylim([1 size(temp_wvf,3)+1]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    title(titles(i));
    box off
end
autoArrangeFigures();

%% (2) Plot variance of each channel across stimulation

x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);
for i = 1:nPatterns
    %subplot(3,3,i);
    subplot(2,3,i); temp_wvf = corrected_wvf{3,i};
    cmax = max(max(abs(temp_wvf([1:191,193:end-1],:))));
    imagesc('XData',x,'YData',1:size(corrected_wvf{2,i},1),'CData',temp_wvf);
    colormap(parula); colorbar; caxis([0 cmax]);
    xlabel('Time (ms)'); ylabel('Channel');
    xlim([x(1) x(end)]); ylim([1 size(correctedWaveform,1)]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    % patch('XData',[0 20*0.001],'YData',[size(stimWaveform,1),size(stimWaveform,1)],...
    %       'FaceColor',colors(1),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    title(titles(i));
    box off
end
autoArrangeFigures();

%% (2) Plot same channel variance at different power level, mean vs var

channel = 123;
x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);
subplot(1,3,1);
for i = 1:nPatterns
    plot(x,artifact_wvf{3,i}(channel,:),"LineWidth",1,"Color",colors(i)); hold on
end
xlabel('Time (ms)'); ylabel('Variance');
legend(titles); box off

subplot(1,3,2);
for i = 1:nPatterns
    plot(x,artifact_wvf{3,i}(channel,:)./artifact_wvf{2,i}(channel,:),"LineWidth",1,"Color",colors(i)); hold on
end
xlabel('Time (ms)'); ylabel('Variance/mean');
legend(titles); box off

subplot(1,3,3);
for i = 1:nPatterns
    plot(artifact_wvf{2,i}(channel,:),artifact_wvf{3,i}(channel,:),"LineWidth",1,"Color",colors(i)); hold on
end
xlabel('Mean'); ylabel('Variance');
legend(titles); box off

%% (2) Test writing directly to ap.bin

% Read in a test section
raw = [20,19,13,11,14,4,6,6,9,10];
mmf.Data.x(1,1:10) = 1;
modified = double(mmf.Data.x(1,1:10));
mmf.Data.x(1,1:10) = raw;
restored = double(mmf.Data.x(1,1:10));

%% (2) Re-plot original file with bigger time window

timeRange = [-0.003 0.01]; % time range of mean subtraction
npulse_per_stim = 1; stim_per_pattern = 100; nPatterns = 9;
nStim = npulse_per_stim * stim_per_pattern * nPatterns;
laserOnset = find(redLaser);
laser_pulse_duration = zeros(1,length(laserOnset));

% 20220523
laser_pulse_duration(1:300) = 0.01;
laser_pulse_duration(301:600) = 0.05;
laser_pulse_duration(601:900) = 0.02;

% Generate average stim waveform
artifact_wvf = cell(1,nPatterns);
% nBins = double((timeRange(2)-timeRange(1))*ap.Fs);
nBins = (timeRange(2)-timeRange(1)) / (1/ap.Fs);

% for i = 701 % test one pattern only
for i = 1:stim_per_pattern:nStim % 1, 101, 201, ..., 801
    stimOnIdx = laserOnset(i:i+stim_per_pattern-1);
    stimWaveform = zeros(gwfparams.nCh,nBins,length(stimOnIdx));
    
    % Generate stimWaveform
    for j = 1:length(stimOnIdx)

        %{
        Finds the wrong imec idx (not sure why)
        Can test this by type following in command window:
            format long g
            timeNI(stimOnIdx(j))
            timeImec(imecFirstIdx+300)

        % Find first Imec index
        [~, imecStimIdx] = min(abs(timeImec-timeNI(stimOnIdx(j))));
        % Find corresponding imec index
        imecFirstIdx = imecStimIdx - floor(timeRange(1)*ap.Fs);
        %}
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

% Plot mean waveform
x = linspace(timeRange(1)*1000,timeRange(2)*1000,nBins);

% Plot titles for 20220523
titles = {'5mW 10ms pulse @ 1Hz','1mW 10ms pulse @ 1Hz','100uW 10ms pulse @ 1Hz',...
          '5mW 50ms pulse @ 1Hz','1mW 50ms pulse @ 1Hz','100uW 50ms pulse @ 1Hz',...
          '5mW 20ms pulse @ 1Hz','1mW 20ms pulse @ 1Hz','100uW 20ms pulse @ 1Hz'};

for i = 1:nPatterns
    figure(i); 
    imagesc('XData',x,'YData',1:size(artifact_wvf{2,i},1),'CData',artifact_wvf{2,i});
    colormap(colormap_BlueWhiteRed); colorbar; %clim([-200 200]);
    xlabel('Time (ms)'); ylabel('Channel');
    xlim([x(1) x(end)]); ylim([1 size(stimWaveform,1)]);
    xline(0,'-','Stimulation','Color','r','LineWidth',2,'HandleVisibility','off');
    % patch('XData',[0 20*0.001],'YData',[size(stimWaveform,1),size(stimWaveform,1)],...
    %       'FaceColor',colors(1),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    title(titles(i));
    box off
end
autoArrangeFigures();