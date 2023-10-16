%% Load data
clear; close all;
% addpath(genpath('/Users/shunli/Downloads/Sabatini lab/Methods'))
addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));

session_path = 'D:\Shun\Recordings\';
sessionName = '20220810-sync-test_g0';
camera = load(strcat(session_path,sessionName,'\','times_cam1-2022-08-10T09_51_06.csv'));
cam_Fs = 119.13;

%% Check skipped frames
skipped_frame = find(diff(camera(:,2)) > 1);
if ~isempty(skipped_frame)
    disp(['Found ', num2str(length(skipped_frame)), ' skipped frames!']);
end

% What should I do next???????????

%% Load imec and ni data
session.pathNidq = strcat(session_path,sessionName, '\');
session.pathImec = strcat(session_path,sessionName, '\', sessionName, '_imec0\');
session.apBin = strcat(sessionName,'_t0.imec0.ap.bin');
session.nidqBin = strcat(sessionName,'_t0.nidq.bin');

ap.meta = ReadMeta(session.apBin, session.pathImec);
ap.Fs = str2double(ap.meta.imSampRate);
% lfp.meta = ReadMeta(session.lfpBin, session.pathImec);
% lfp.Fs = str2double(lfp.meta.imSampRate);
nidq.meta = ReadMeta(session.nidqBin, session.pathNidq);
nidq.Fs = str2double(nidq.meta.niSampRate);

%% Load sync data from imec
% Set time blocks for NIDAQ data
% This is needed cause data is too big to load as a whole

% Use the whole session
apTotalSecs = str2double(ap.meta.fileSizeBytes) / ap.Fs / str2double(ap.meta.nSavedChans) / 2;
ap.totalSampIncluded = floor(apTotalSecs * ap.Fs);
nidqTotalSecs = str2double(nidq.meta.fileSizeBytes) / nidq.Fs / str2double(nidq.meta.nSavedChans) / 2;
nidq.totalSampIncluded = floor(nidqTotalSecs * nidq.Fs);
blockTime = nidqTotalSecs;
nBlocks = 1;

% Read Sync Channel in IMEC data (#385)
tic
syncImec = ReadBinByCh(0, ap.totalSampIncluded, ap.meta, session.apBin, session.pathImec, 385);
% syncImecLFP = ReadBinByCh(0, lfp.totalSampIncluded, lfp.meta, session.lfpBin, session.pathImec, 385);
disp('time for reading sync signal from IMEC data');
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


%% Extract rising edge
% Convert sync signals to digital signal (boolean);
syncImec = (syncImec > (max(syncImec)/2));
syncNI = (syncNI==1);
syncCam = (camera(:,1) > max(camera(:,1)/2))';
% camera(camera(:,1) < middle) = 0;
% camera(camera(:,1) >= middle) = 1;
% plot(camera(:,1));

% Extract location of rising edge
temp = [false,diff(syncImec)]; % diff(syncImec) records time of rising edge (1) and falling edge (0)
syncImec_diff = (temp==1);
temp = [false,diff(syncNI)];
syncNI_diff = (temp==1);
temp = [false,diff(syncCam)];
syncCam_diff = (temp==1);

%% Plot sync channel data
figure 
plot(syncImec); hold on
plot(syncNI); hold on
plot(syncCam);
legend
title('Sync pulses');

%% Find the first common pulse between Imec and NI

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

%% Find common first pulse of imec and ni
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
            plot(2*y_upsample(c1*y_idx(j):c1*y_idx(j)+l));
            title(['Imec pulse=',num2str(i), ', NI pulse=',num2str(j),...
                ', dotprod: ',num2str(dotprod),', maxdot: ',num2str(maxdot)]);
            pause
        end
    end
end

syncImec_first = x_first;
syncNI_first = y_first;

%% Find common first pulse of imec and camera
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

%% Assign each timestep of IMEC and NI with common time stamp

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

subplot(2,3,1); plot(timeImec); title('Imec time in sec');
subplot(2,3,2); plot(timeNI); title('NI time in sec');
subplot(2,3,3); plot(timeCamera); title('Camera time in sec');
subplot(2,3,4); plot(diff(timeImec)); title('Diff timeImec');
subplot(2,3,5); plot(diff(timeNI)); title('Diff timeNI');
subplot(2,3,6); plot(diff(timeCamera)); title('Diff timeCamera');

figure(13); plot(timeImec(idx_Imec(1:l))-timeNI(idx_NI(1:l))); hold on
plot(timeImec(idx_Imec(1:l))-timeCamera(idx_Cam(1:l))); hold on
plot(timeNI(idx_NI(1:l))-timeCamera(idx_Cam(1:l)));

