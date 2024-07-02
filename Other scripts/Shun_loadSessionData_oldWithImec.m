% Code based on SJ's ReadSGLXData_script.m
%% Load file
clear; close all;
addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
% sessionName = '20220330-SJ495-R_g0';
% sessionName = '20220401-SJ495-R_g0';
% sessionName = '20220518-SJ512-R-photometry only_g0';
% sessionName = '20220523-SJ518-L_g0';
% sessionName = '20220527-SJ512-L_g0';
% sessionName = '20220525-SJ518-L_g0';
% sessionName = '20220530-SJ512-L_g0';
% sessionName = '20220602-SJ513-L_g0';
% sessionName = '20220603-SJ513-L_g0';
% sessionName = '20220604-SJ513-L_g0';
% sessionName = '20220606-SJ513-R_g0';
% sessionName = "20220612-SJ518-R_g0";
% sessionName = "20220613-SJ518-R_g0";
%sessionName = "20220615-SJ518-R_g0";
sessionName = "20220908-SJ519-D4-VLS-L_g0";

% serverpath = '\\research.files.med.Harvard.edu\neurobio\MICROSCOPE\Kevin\3-Experiments\5-InVivo\3-Neuropixel_Phys\1-Raw_Data_KM\2022\8-Aug\';
localpath = 'D:\Shun\Recordings\'; 
session.path = localpath;

% session.pathImec = strcat('/Volumes/Shun neuro data/Neuropixel/',sessionName, '/', sessionName, '_imec0/');
session.pathImec = strcat(session.path,sessionName, '\catgt_', sessionName,'\');
session.pathNidq = strcat(session.path,sessionName, '\');
session.apBin = strcat(sessionName,'_tcat.imec0.ap.bin');
session.nidqBin = strcat(sessionName,'_t0.nidq.bin');
% session.apBin = strcat(sessionName,'_t0.imec0.ap.bin');

ap.meta = ReadMeta(session.apBin, session.pathImec);
ap.Fs = str2double(ap.meta.imSampRate);
% lfp.meta = ReadMeta(session.lfpBin, session.pathImec);
% lfp.Fs = str2double(lfp.meta.imSampRate);
nidq.meta = ReadMeta(session.nidqBin, session.pathNidq);
nidq.Fs = str2double(nidq.meta.niSampRate);

%% Set time blocks for NIDAQ data
% This is needed cause data is too big to load as a whole

% blockTime = 60; %s
% nBlocks = 20; % total 1200 sec
% ap.nSampPerBlock = floor(blockTime * ap.Fs);
% ap.totalSampIncluded = nBlocks * ap.nSampPerBlock;
% % lfp.nSampPerBlock = floor(blockTime * lfp.Fs);
% % lfp.totalSampIncluded = nBlocks * lfp.nSampPerBlock;
% nidq.nSampPerBlock = floor(blockTime * nidq.Fs);
% nidq.totalSampIncluded = nBlocks * nidq.nSampPerBlock; 

% Use the whole session
apTotalSecs = str2double(ap.meta.fileSizeBytes) / ap.Fs / str2double(ap.meta.nSavedChans) / 2;
ap.totalSampIncluded = floor(apTotalSecs * ap.Fs);
nidqTotalSecs = str2double(nidq.meta.fileSizeBytes) / nidq.Fs / str2double(nidq.meta.nSavedChans) / 2;
nidq.totalSampIncluded = floor(nidqTotalSecs * nidq.Fs);
blockTime = nidqTotalSecs;
nBlocks = 1;

session.blockTime = blockTime;
session.nBlocks = nBlocks;

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

return;

%% Read Sync Channel in IMEC data (#385)

% Normally takes ~10 seconds for 60s data; 230s for 1200s data
tic
syncImec = ReadBinByCh(0, ap.totalSampIncluded, ap.meta, session.apBin, session.pathImec, 385);
% syncImecLFP = ReadBinByCh(0, lfp.totalSampIncluded, lfp.meta, session.lfpBin, session.pathImec, 385);
display('time for reading sync signal from IMEC data');
toc

%% Read NIDAQ data

nidqData = ReadBin(0, nidq.totalSampIncluded, nidq.meta, session.nidqBin, session.pathNidq);

% Read Analog Channels in NIDAQ data
% For an analog channel: gain correct saved channel ch (1-based for MATLAB)
a_ch_start = str2num(nidq.meta.niXAChans1(1));
a_ch_end = str2num(nidq.meta.niXAChans1(3));
a_ch_list = [a_ch_start:a_ch_end];

analogNI = GainCorrectNI(nidqData, a_ch_list+1, nidq.meta); % remove +1 if analyzing session before 0512

for i = 1:numel(a_ch_list)
    figure(10+i);
    plot(analogNI(a_ch_list(i)+1,:)); % remove +1 if analyzing session before 0512
    title(['NIDAQ analog channel ',num2str(a_ch_list(i))]);
end

% Read Digital Channels in NIDAQ data
d_ch_start = str2num(nidq.meta.niXDChans1(1));
d_ch_end = str2num(nidq.meta.niXDChans1(3));
d_ch_list = [d_ch_start:d_ch_end];

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;
digitalNI = ExtractDigital(nidqData, nidq.meta, dw, d_ch_list);

syncNI = digitalNI(2,:); % Save sync channel from NIDAQ data separately

%% Plot sync channel data
figure 
plot(syncImecAP); hold on
plot(syncNI);
legend
title('sync pulse from IMEC data');
% autoArrangeFigures();

%% Save retreived data

save(strcat('sync_',sessionName),'sessionName','ap','nidq','session',...
                'syncImec','syncNI');
% save(strcat('sync_',sessionName),'sessionName','nidq','session',...
%                         'syncNI','analogNI','digitalNI');

%% Extract behavior from NIDAQ (Jay's task)

% Convert other analog signals to rising edge
leftLick = (analogNI(1,:)> (max(analogNI(1,:))/2));
temp = [false, diff(leftLick)];
leftLick = (temp==1);

rightLick = (analogNI(2,:)> (max(analogNI(2,:))/2));
temp = [false, diff(rightLick)];
rightLick = (temp==1);

% gyroX = analogNI(4,:);
% gyroY = analogNI(5,:);
% gyroZ = analogNI(6,:);
basket = single(analogNI(4:6,:));
photometry = analogNI(8,:);

% Convert other digital signals to rising edge.
temp = false(size(digitalNI));
for i=1:size(digitalNI,1)
    temp2 = [false, diff(digitalNI(i,:))];
    temp(i,:) = (temp2==1);
end
start_session = find(digitalNI(1,:)==1,1);
leftSolenoid = temp(3,:);
rightSolenoid = temp(4,:);
leftTone = temp(5,:);
rightTone = temp(6,:);
redLaser = temp(7,:);
blueLaser = temp(8,:);

% Set up allTrials: not including omissionTrials
% 1 is left trial, 2 is right trial
allTrials = zeros(size(leftTone));
allTrials(leftTone == 1) = 1;
allTrials(rightTone == 1) = 2;

figure(1);
plot(syncImec_diff); hold on
plot(syncNI_diff); title('sync pulse from Imec/NI');

figure(2);
plot(leftLick); title('left lick'); hold on
plot(leftSolenoid,'LineWidth',2); title('left reward'); hold on
plot(leftTone,'LineWidth',2); title('left tone');

figure(3);
plot(rightLick); title('right lick'); hold on
plot(rightSolenoid,'LineWidth',2); title('right reward'); hold on
plot(rightTone,'LineWidth',2); title('right tone');

%% Read Omission trials (Jay's task)
taskOutput = readtable('20220427 lick task.xlsx','Sheet','SJ498');
omissionRows = find(contains(taskOutput.Var3,"Omission"));
omissionTrialsNum = table2array(taskOutput(omissionRows-2,1));

% Set up ommissionTrials
% 1 is left omission trials, 2 is right omission trials
allTrialsIdx = find(allTrials > 0);
omissionTrialsNum = omissionTrialsNum(omissionTrialsNum <= length(allTrialsIdx));
omissionTrialsTime = allTrialsIdx(omissionTrialsNum);
omissionTrials = zeros(size(leftTone));
omissionTrials(omissionTrialsTime) = allTrials(omissionTrialsTime);

% Remove omission trials from allTrials
allTrials(omissionTrialsTime) = 0;

% Check omisison trial data
figure; 
startTime = omissionTrialsTime(1) - 50000;
endTime = startTime + 2000000;
plot(allTrials(startTime:endTime)); hold on
plot(leftSolenoid(startTime:endTime)); hold on
plot(rightSolenoid(startTime:endTime)); hold on
plot(omissionTrials(startTime:endTime),'LineWidth',1);
legend({'All trials','Left Solenoid','Right solenoid','Omission trials'});




%% Extract behavior from NIDAQ (valence switching task)

% Convert other analog signals to rising edge
leftLick = (analogNI(1,:)> (max(analogNI(1,:))/2));
temp = [false, diff(leftLick)];
leftLick = (temp==1);

rightLick = (analogNI(2,:)> (max(analogNI(2,:))/2));
temp = [false, diff(rightLick)];
rightLick = (temp==1);

basket = single(analogNI(4:6,:));
photometry = analogNI(8,:);

% Convert other digital signals to rising edge.
temp = false(size(digitalNI));
for i=1:size(digitalNI,1)
    temp2 = [false, diff(digitalNI(i,:))];
    temp(i,:) = (temp2==1);
end
airpuff = temp(1,:);
leftSolenoid = temp(3,:);
rightSolenoid = temp(4,:);
leftTone = temp(5,:);
rightTone = temp(6,:);
redLaser = temp(7,:);
blueLaser = temp(8,:);

% Set up allTrials: not including omissionTrials
% 1 is left trial, 2 is right trial
allTrials = zeros(size(leftTone));
allTrials(leftTone == 1) = 1;
allTrials(rightTone == 1) = 2;

figure(1);
plot(leftLick); title('left lick'); hold on
plot(leftSolenoid,'LineWidth',2); title('left reward'); hold on
plot(leftTone,'LineWidth',2); title('left tone');

figure(2);
plot(rightLick); title('right lick'); hold on
plot(rightSolenoid,'LineWidth',2); title('right reward'); hold on
plot(rightTone,'LineWidth',2); title('right tone');


%% Save behavioral data
save(strcat('sync_',sessionName),'leftLick','rightLick','basket', ...
    'leftSolenoid','rightSolenoid','allTrials',...
    'photometry','blueLaser','redLaser','-append');

% save(strcat('sync_',sessionName),'omissionTrials','-append');

%% (With camera) Extract rising edge and match sample freq

% Load data
camerapath = dir(strcat(session.path,sessionName,'\','times_cam1-*.csv'));
camera = load([camerapath.folder,'\',camerapath.name]);
cam_Fs = 143.93;

% Check skipped frames
skipped_frame = find(diff(camera(:,2)) > 1);
if ~isempty(skipped_frame)
    disp(['Found ', num2str(length(skipped_frame)), ' skipped frames!']);
end
% What should I do next???????????

% Convert sync signals to digital signal (boolean);
syncImec = (syncImec > (max(syncImec)/2));
syncNI = (syncNI==1);
syncCam = (camera(:,1) > max(camera(:,1)/2))';
% plot(camera(:,1));

% Extract location of rising edge
temp = [false,diff(syncImec)]; % diff(syncImec) records time of rising edge (1) and falling edge (0)
syncImec_diff = (temp==1);
temp = [false,diff(syncNI)];
syncNI_diff = (temp==1);
temp = [false,diff(syncCam)];
syncCam_diff = (temp==1);

% **Extract start and match sample freq**

% Loop over first n pulse of Imec/NI
x_idx = find(syncImec_diff>0,20);
y_idx = find(syncNI_diff>0,20);
z_idx = find(syncCam_diff>0,20);

% First 60s of sync channel signal
if 60*ap.Fs > ap.totalSampIncluded
    x = syncImec; 
    y = syncNI;
    z = syncCam;
else
    x = syncImec(1:floor(ap.Fs*60)); 
    y = syncNI(1:floor(nidq.Fs*60));
    z = syncCam(1:floor(cam_Fs*60));
end

% Upsample y,z to match sampling rate of x
c1 = floor(ap.Fs/nidq.Fs);
c2 = floor(ap.Fs/cam_Fs);
y_upsample = upsample(y,c1);
z_upsample = upsample(z,c2);
% Fill upsampled steps with 1 if during sync pulse
for i = 1:length(y_upsample)-c1
    if y_upsample(i) == 1 && y_upsample(i+c1) == 1
        for j = 1:c1-1
            y_upsample(i+j) = 1;
        end
    end
end
for i = 1:length(z_upsample)-c2
    if z_upsample(i) == 1 && z_upsample(i+c2) == 1
        for j = 1:c2-1
            z_upsample(i+j) = 1;
        end
    end
end

plot(x); hold on
plot(y_upsample); hold on
plot(z_upsample);
legend;

%% (With camera) Find common first pulse of imec and ni
maxdot = 0;
for i = 1:length(x_idx)
    for j = 1:length(y_idx)
        l = min(length(x)-x_idx(i),length(y_upsample)-c1*y_idx(j));
        dotprod = dot(x(x_idx(i):x_idx(i)+l),y_upsample(c1*y_idx(j):c1*y_idx(j)+l));
        disp(['dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot),...
            ', Imec pulse: ',num2str(i), ', NI pulse: ',num2str(j)]);

        if dotprod > maxdot
            x_first = x_idx(i);
            y_first = y_idx(j);
            maxdot = dotprod;

            figure;
            plot(x(x_idx(i):x_idx(i)+l)); hold on
            plot(2*y_upsample(c1*y_idx(j):c1*y_idx(j)+l)); legend({'Imec','NI'});
            title(['Imec pulse=',num2str(i), ', NI pulse=',num2str(j),...
                ', dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot)]);
            pause
        end
    end
end

syncImec_first = x_first;
syncNI_first = y_first;
close;

%% (With camera) Find common first pulse of imec and camera
maxdot = 0;
for i = 1:length(x_idx)
    for j = 1:length(z_idx)
        l = min(length(x)-x_idx(i),length(z_upsample)-c2*z_idx(j));
        dotprod = dot(x(x_idx(i):x_idx(i)+l),z_upsample(c2*z_idx(j):c2*z_idx(j)+l));
        disp(['dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot),...
            ', Imec pulse: ',num2str(i), ', Camera pulse: ',num2str(j)]);

        if dotprod > maxdot
            x_first = x_idx(i);
            z_first = z_idx(j);
            maxdot = dotprod;

            figure;
            plot(x(x_idx(i):x_idx(i)+l)); hold on
            plot(2*z_upsample(c2*z_idx(j):c2*z_idx(j)+l));
            title(['Imec pulse=',num2str(i), ', Camera pulse=',num2str(j),...
                ', dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot)]);
            pause
        end
    end
end

if syncImec_first ~= x_first
    error('First common sync pulse between Imec&NI and Imec&Cam is different');
end
syncCam_first = z_first;
close;

%% (With camera) Assign each timestep of IMEC, NI, camera with common time stamp

% Filling in the time for IMEC timebin using the sampling rate
timeImec=zeros(1,length(syncImec));
for i=1:length(syncImec)
    % Time of each timebin of IMEC in sec
    timeImec(i)=(i-syncImec_first)*(1/ap.Fs);
end
idx_Imec = find(syncImec_diff(syncImec_first:end)==1)+syncImec_first-1;

% Filling in the time for NI using the sampling rate and interpolation
% reference = IMEC
timeNI=zeros(1,length(syncNI));
% Index of all sync pulses starting from the first pulse
idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;

% Filling in the time for camera using the sampling rate and interpolation
% reference = IMEC
timeCamera=zeros(1,length(syncCam));
% Index of all sync pulses starting from the first pulse
idx_Cam = find(syncCam_diff(syncCam_first:end)==1)+syncCam_first-1;


% Fill in time
idx_previous_NI = 0;
idx_previous_cam = 0;
l = min([length(idx_Imec),length(idx_NI),length(idx_Cam)]);
for i=1:l
    % Sync the time using the rising edge of the current sync pulse
    timeNI(idx_NI(i)) = timeImec(idx_Imec(i));
    timeCamera(idx_Cam(i)) = timeImec(idx_Imec(i));
    % Filled in time in between two pulses
    if idx_previous_NI < idx_NI(i)
        timeNI(idx_previous_NI+1 : idx_NI(i)) = ...
            timeNI(idx_NI(i)) + (1/nidq.Fs)*[idx_previous_NI-idx_NI(i)+1 :1: 0];
    end
    if idx_previous_cam < idx_Cam(i)
        timeCamera(idx_previous_cam+1 : idx_Cam(i)) = ...
            timeCamera(idx_Cam(i)) + (1/cam_Fs)*[idx_previous_cam-idx_Cam(i)+1 :1: 0];
    end
    % Filled in time after the last syncNI pulse
    if i == l
        timeNI(idx_NI(i):end) = ...
            timeNI(idx_NI(i)) + (1/nidq.Fs)*[0 :1: length(timeNI)-idx_NI(i)];
        timeCamera(idx_Cam(i):end) = ...
            timeCamera(idx_Cam(i)) + (1/cam_Fs)*[0 :1: length(timeCamera)-idx_Cam(i)];
    end

    idx_previous_NI = idx_NI(i);
    idx_previous_cam = idx_Cam(i);
end

% Align timeNI and timeImec to start at 0
t0 = timeImec(1);
timeImec = timeImec - t0; % make timeImec to start at 0s
timeNI = timeNI - t0; % align timeNI to timeImec
timeCamera = timeCamera - t0; % align timeCamera to timeImec

subplot(3,3,1); plot(timeImec); title('Imec time in sec');
subplot(3,3,2); plot(timeNI); title('NI time in sec');
subplot(3,3,3); plot(timeCamera); title('Camera time in sec');
subplot(3,3,4); plot(diff(timeImec)); title('d(timeImec)');
subplot(3,3,5); plot(diff(timeNI)); title('d(timeNI)');
subplot(3,3,6); plot(diff(timeCamera)); title('d(timeCamera)');
subplot(3,3,7); plot(timeImec(idx_Imec(1:l))-timeNI(idx_NI(1:l))); title('Imec-NI');
subplot(3,3,8); plot(timeImec(idx_Imec(1:l))-timeCamera(idx_Cam(1:l))); title('Imec-cam');
subplot(3,3,9); plot(timeNI(idx_NI(1:l))-timeCamera(idx_Cam(1:l))); title('NI-cam');

%% (With camera) Save sync time

% sync.syncImec_first = syncImec_first;
% sync.syncNI_first = syncNI_first;
% sync.syncCam_first = syncCam_first;
% sync.syncImec_diff = syncImec_diff;
% sync.syncNI_diff = syncNI_diff;
% sync.syncCam_diff = syncCam_diff;

sync.imecFs = ap.Fs;
sync.nidqFs = nidq.Fs;
sync.camFs = cam_Fs;

save(strcat('sync_',sessionName),'sync','timeImec','timeNI','timeCamera','-append');

%% (No camera) Find the first common pulse between Imec and NI

% For NIDAQ only session
% syncImec = syncNI;
% syncImec_diff = syncNI_diff;
% ap.Fs = nidq.Fs;

% Find indices of first two pulses
% x_idx=find(syncImec_diff>0,2);
% y_idx=find(syncNI_diff>0,2);
% 
% pair_idx=[x_idx(1) y_idx(1);x_idx(1) y_idx(2);x_idx(2) y_idx(1)];
% maxdot = 0;
% x=syncImec(1:floor(ap.Fs*60)); % First 60s of sync channel signal
% y=syncNI(1:floor(nidq.Fs*60));
% for i=1:size(pair_idx,1)
%     l = min(length(x)-pair_idx(i,1),length(y)-pair_idx(i,2));
%     dotprod = dot(x(pair_idx(i,1):pair_idx(i,1)+l),y(pair_idx(i,2):pair_idx(i,2)+l));
%     disp(['dotprod: ',num2str(dotprod),' maxdot: ',num2str(maxdot),' i: ',num2str(i)]);
%     if dotprod > maxdot
%         x_first=pair_idx(i,1);
%         y_first=pair_idx(i,2);
%         maxdot = dotprod;
%     end
% end

% Convert sync signals to digital signal (boolean);
syncImec = (syncImec > (max(syncImec)/2));
syncNI = (syncNI==1);

% Extract location of rising edge
temp = [false,diff(syncImec)]; % diff(syncImec) records time of rising edge (1) and falling edge (0)
syncImec_diff = (temp==1);
temp = [false,diff(syncNI)];
syncNI_diff = (temp==1);

% Loop over first 10 pulse of Imec/NI
x_idx=find(syncImec_diff>0,20);
y_idx=find(syncNI_diff>0,20);
x = syncImec(1:floor(ap.Fs*60)); % First 60s of sync channel signal
y = syncNI(1:floor(nidq.Fs*60));

% Upsample y to match sampling rate of x
c = floor(ap.Fs/nidq.Fs);
y_upsample = upsample(y,c);
for i = 1:length(y_upsample)-3
    if y_upsample(i) == 1 && y_upsample(i+3) == 1
        y_upsample(i+1) = 1;
        y_upsample(i+2) = 1;
    end
end

maxdot = 0;
for i = 1:length(x_idx)
    for j = 1:length(y_idx)
        l = min(length(x)-x_idx(i),length(y_upsample)-c*y_idx(j));
        dotprod = dot(x(x_idx(i):x_idx(i)+l),y_upsample(c*y_idx(j):c*y_idx(j)+l));
        disp(['dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot),...
            ', Imec pulse: ',num2str(i), ', NI pulse: ',num2str(j)]);

        if dotprod > maxdot
            x_first = x_idx(i);
            y_first = y_idx(j);
            maxdot = dotprod;

            figure;
            plot(x(x_idx(i):x_idx(i)+l)); hold on
            plot(2*y_upsample(y_idx(j):y_idx(j)+l));
            title(['Imec pulse=',num2str(i), ', NI pulse=',num2str(j),...
                ', dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot)]);
            pause
        end
    end
end

% figure(9); plot(x(x_first:end),'b'); hold on;
% plot(y(y_first:end),'r'); hold on;

syncImec_first = x_first;
syncNI_first = y_first;

%% (No camera) Assign each timestep of IMEC and NI with common time stamp

% Filling in the time for IMEC timebin using the sampling rate
timeImec=zeros(1,length(syncImec));
for i=1:length(syncImec)
    % Time of each timebin of IMEC in sec
    timeImec(i)=(i-syncImec_first)*(1/ap.Fs);
end

% Filling in the time for NI using the sampling rate and interpolation
% reference = IMEC
timeNI=zeros(1,length(syncNI));
% Index of all sync pulses starting from the first pulse
idx_Imec = find(syncImec_diff(syncImec_first:end)==1)+syncImec_first-1;
idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;

% idx_previous = 1;
idx_previous = 0;
for i=1:min(length(idx_NI),length(idx_Imec))
    % Sync the time using the rising edge of the current sync pulse
    timeNI(idx_NI(i))=timeImec(idx_Imec(i));
    % Filled in time in between two pulses
    if idx_previous < idx_NI(i)
        timeNI(idx_previous+1 : idx_NI(i)) = ...
            timeNI(idx_NI(i)) + (1/nidq.Fs)*[idx_previous-idx_NI(i)+1 :1: 0];
    end
    % Filled in time after the last syncNI pulse
    if i == length(idx_NI)
        timeNI(idx_NI(i):end) = ...
            timeNI(idx_NI(i)) + (1/nidq.Fs)*[0 :1: length(timeNI)-idx_NI(i)];
    end
    idx_previous = idx_NI(i);
end

% % Align timeNI and timeImec to start at 0
t0 = timeImec(1);
timeImec = timeImec - t0; % make timeImec to start at 0s
timeNI = timeNI - t0; % align timeNI to timeImec

figure(10); plot(timeImec); title('Imec Time in sec'); hold on
figure(11); plot(timeNI); title('NI Time in sec');

l=min(length(idx_Imec),length(idx_NI));
figure(12); plot(timeImec(idx_Imec(1:l))-timeNI(idx_NI(1:l)));
% figure(12); plot(timeImec(1:l) - timeNI(1:l));
% autoArrangeFigures();

%% (No camera) Save sync time

sync.syncImec_first = syncImec_first;
sync.syncNI_first = syncNI_first;
sync.syncImec_diff = syncImec_diff;
sync.syncNI_diff = syncNI_diff;

save(strcat('sync_',sessionName),'sync','timeImec','timeNI','-append');



