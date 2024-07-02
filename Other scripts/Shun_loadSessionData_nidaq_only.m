%% Load file
clear; close all;
%addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
% addpath(genpath('/Volumes/neurobio/MICROSCOPE/Shun/Neuropixel analysis/Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Neuropixel analysis\Methods'));

% Enter session variables
sessionName = "20221106-SL037-D4_g0";
% sessionName = "20221004-TP174-D9_g0";
videoName = 'times_cam1-2022-09-18T17_26_25.csv';

% session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Project risk\Recordings\';
session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Project valence\Recordings\';
% session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Neuropixel VLS\Recordings\';
% session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Test photometry\';
session.pathNidq = strcat(session.path,sessionName, '\');
session.nidqBin = strcat(sessionName,'_t0.nidq.bin');
% session.pathNidq = strcat(session.path,sessionName, '\catgt_',sessionName);
% session.nidqBin = strcat(sessionName,'_tcat.nidq.bin');

nidq.meta = ReadMeta(session.nidqBin, session.pathNidq);
nidq.Fs = str2double(nidq.meta.niSampRate);

disp(sessionName);

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
nidqTotalSecs = str2double(nidq.meta.fileSizeBytes) / nidq.Fs / str2double(nidq.meta.nSavedChans) / 2;
nidq.totalSampIncluded = floor(nidqTotalSecs * nidq.Fs);
blockTime = nidqTotalSecs;
nBlocks = 1;

session.blockTime = blockTime;
session.nBlocks = nBlocks;

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

autoArrangeFigures();

%% Save retrieved data
save(strcat(session.path,sessionName,'/','sync_',sessionName),...
    'sessionName','nidq','session','syncNI');
disp("Session data saved");

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

basket = single(analogNI(4:6,:));
photometry = analogNI(8,:);

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
redLaser = temp(7,:); % or ENL
blueLaser = temp(8,:);

% Set up allTrials: not including omissionTrials
% 1 is left trial, 2 is right trial
allTrials = zeros(size(leftTone));
allTrials(leftTone == 1) = 1;
allTrials(rightTone == 1) = 2;

% figure(1);
% plot(syncNI_diff); title('sync pulse from Imec/NI');
% 
% figure(2);
% plot(leftLick); title('left lick'); hold on
% plot(leftSolenoid,'LineWidth',2); hold on
% plot(leftTone,'LineWidth',2);

% figure(3);
% plot(rightLick); title('right lick'); hold on
% plot(rightSolenoid,'LineWidth',2); title('right reward'); hold on
% plot(rightTone,'LineWidth',2); title('right tone');

figure;
plot(rightSolenoid);

%% Save behavioral data (valence switching task)
save(strcat(session.path,sessionName,'\','sync_',sessionName),...
    'airpuff','leftLick','rightLick',...
    'leftSolenoid','rightSolenoid','allTrials','temperature',...
    'photometry','blueLaser','redLaser','-append');

% save(strcat('sync_',sessionName),'omissionTrials','-append');
disp("NIDAQ data saved");

%% (Without camera) Save sync time

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

save(strcat(session.path,sessionName,'\','sync_',sessionName),'timeNI','-append');
disp("Time data saved");

%% (With camera) Extract rising edge and match sample freq

% Load data
camerapath = dir(strcat(session.path,sessionName,'\',videoName));
camera = load([camerapath.folder,'\',camerapath.name]);
cam_Fs = 143.93;

% Check skipped frames
skipped_frame = find(diff(camera(:,2)) > 1);
if ~isempty(skipped_frame)
    disp(['Found ', num2str(length(skipped_frame)), ' skipped frames!']);
end
% What should I do next???????????

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
    z = syncCam(1:floor(cam_Fs*60));
end

% Upsample camera to match sampling rate of nidq
c2 = floor(nidq.Fs/cam_Fs);
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

%% (With camera) Find common first pulse of NI and camera
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
            pause
        end
    end
end

syncNI_first = y_first;
syncCam_first = z_first;
close;

%% (With camera) Assign each timestep of NI and camera with common time stamp

% Filling in the time for IMEC timebin using the sampling rate
timeNI = zeros(1,length(syncNI));
for i=1:length(syncNI)
    % Time of each timebin of IMEC in sec
    timeNI(i)=(i-syncNI_first)*(1/nidq.Fs);
end
idx_NI = find(syncNI_diff(syncNI_first:end)==1)+syncNI_first-1;

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
            timeCamera(idx_Cam(i)) + (1/cam_Fs)*[idx_previous_cam-idx_Cam(i)+1 :1: 0];
    end
    % Filled in time after the last syncNI pulse
    if i == l
        timeCamera(idx_Cam(i):end) = ...
            timeCamera(idx_Cam(i)) + (1/cam_Fs)*[0 :1: length(timeCamera)-idx_Cam(i)];
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

%% (With camera) Save sync time

sync.nidqFs = nidq.Fs;
sync.camFs = cam_Fs;

save(strcat(session.path,sessionName,'\','sync_',sessionName),...
    'sync','timeNI','timeCamera','-append');


