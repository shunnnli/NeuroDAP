% Shun_analyzeRTPP

% 2025/01/31

% Plots trajectory of RTPP centroid and performs other analysis

%% Load the CSV file

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
[twoColors,~,~,~,~,~,bluePurpleRed] = loadColors;

filename = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project misc/Recordings'),...
                        'Prompt','Select an animal');

sessionList = dir(filename{1});
sessionList = sessionList(~ismember({sessionList.name},{'.','..'}));
nSessions = length(sessionList);

%% Process data

session_cutoffs = [0, 0, 0]; % in min;
Fs = 20; % frame rate

sessions = struct([]);
for s = 1:nSessions
    cur_session = dir(fullfile(sessionList(s).folder,sessionList(s).name));
    tableName = cur_session(contains({cur_session.name},'times-')).name;
    tablePath = cur_session(contains({cur_session.name},'times-')).folder;
    cur_data = readtable(fullfile(tablePath,tableName));

    dirsplit = split(tablePath,filesep);
    sessions(s).name = dirsplit{end};
    sessions(s).path = tablePath;
    sessions(s).data = cur_data;
    sessions(s).cutoff = session_cutoffs(s);

    % Add time column
    sessions(s).data.Item1 = datetime(sessions(s).data.Item1, ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ','TimeZone', 'UTC');
    startTime = sessions(s).data.Item1(1);
    minuteSinceStart = minutes(sessions(s).data.Item1 - startTime);

    % Align to cutoff
    if session_cutoffs(s) >= 0; offset = session_cutoffs(s);
    else; offset = max(minuteSinceStart) + session_cutoffs(s);
    end
    sessions(s).data.time = minuteSinceStart - offset;
end
disp('Finished: animal data loaded');

%% Plot summary figure

initializeFig(0.5,1); tiledlayout(1,nSessions+1);
stim_side = 'left';
noMovementThreshold = 20;

leftColor = [156, 219, 17]./255;
rightColor = [144, 126, 171]./255;
stimColor = [7, 162, 222]./255;
ctrlColor = [.3 .3 .3];

for s = 1:nSessions
    cur_data = sessions(s).data(sessions(s).data.time >= 0,:);
    cur_name = sessions(s).name;

    X_raw = cur_data.Item2_X;
    Y_raw = cur_data.Item2_Y;
    Y_midpoint = 350; %(min(Y_raw)+max(Y_raw))/2;

    % Drop potential sleep time
    staticWindow = getStaticPeriod(X_raw, Y_raw,noMovementThreshold=3,windowDuration=30);
    cur_data_clean = cur_data(~staticWindow,:);
    X = X_raw(~staticWindow);
    Y = Y_raw(~staticWindow);
    % X = X_raw; Y = Y_raw;
    
    % Plot the trajectory
    nexttile;
    if ~contains(cur_name,'RTPP')
        plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), Color=rightColor, LineWidth=2); hold on;
        plot(X(Y<Y_midpoint), Y(Y<Y_midpoint), Color=leftColor, LineWidth=2);
        legend({'Right', 'Left'}, 'Location', 'northeast');
    elseif strcmpi(stim_side,'right')
        plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), Color=stimColor, LineWidth=2); hold on;
        plot(X(Y<Y_midpoint), Y(Y<Y_midpoint), Color=ctrlColor, LineWidth=2);
        legend({'Stim OFF', 'Stim ON'}, 'Location', 'northeast');
    elseif strcmpi(stim_side,'left')
        plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), Color=ctrlColor, LineWidth=2); hold on;
        plot(X(Y<Y_midpoint), Y(Y<Y_midpoint), Color=stimColor, LineWidth=2);
        legend({'Stim OFF', 'Stim ON'}, 'Location', 'northeast');
    end
    
    xlim([300,650]); xlabel('X Position');
    ylim([0,700]);ylabel('Y Position');
    title(cur_name);

    % Plot time in each chamber
    % stim = sum(strcmp(cur_data.Item3,'True'));
    if strcmpi(stim_side,'right')
        stim = find(Y>=Y_midpoint);
    elseif strcmpi(stim_side,'left')
        stim = find(Y<=Y_midpoint);
    end
    stim_pct(s) = length(stim)/length(Y) * 100;
end

% Plot bar plot
nexttile;
for s = 1:nSessions
    if ~contains(sessions(s).name,'RTPP'); color = ctrlColor;
    else; color = stimColor; end
    plotScatterBar(s,stim_pct(s),Color=color,style='bar'); hold on
end
xticks(1:nSessions); xticklabels({sessions.name});
ylabel('Time spent in stimulated side (%)');

%% Heatmap (ongoing)

% Define heatmap grid
squareSize = 5; % in pixels
edgesX = linspace(300, 640, round(range(X)/squareSize));
edgesY = linspace(0, 700, round(range(Y)/squareSize));

% Count occurences of points in each grid bin
[counts, edgesX, edgesY] = histcounts2(X, Y, edgesX, edgesY);
timeInBins = counts * 1/Fs;
% log_counts = log(counts);
% log_times = log(timeInBins);
gamma_counts = counts .^ 0.5;

%4. plot
figure;
imagesc(edgesX(1:end-1), edgesY(1:end-1), gamma_counts');
% % set(gca, 'YDir', 'normal'); 
% set(gca,'ColorScale','log');
colormap(sky); colorbar;
xlim([300,640]); xlabel('X Position');
ylim([0,700]);ylabel('Y Position');
title ('Heatmap');

%% Velocity

% velocity calculations
dx = diff(X);
dy = diff(Y);
dt = 1/Fs;
velocity = sqrt(dx.^2 + dy.^2) ./ dt;

% normalize velocity
velocity_norm = (velocity - min(velocity)) / (max(velocity) - min(velocity));

% colormap
cmap = sky;
colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, velocity_norm);

% make plot
figure;
hold on;

% color gradient
for i = 1:length(velocity)
    plot(X(i:i+1), Y(i:i+1), 'Color', colors(i, :), 'LineWidth', 2)
end

% 7. axis label and titles
xlabel('X-Position')
ylabel('Y-Position')
title('RTPP Practice Analysis')
clim([min(velocity) max(velocity)]); % Set the color axis based on velocity range
colormap(sky);
h = colorbar;  % Create a colorbar
ylabel(h, 'Velocity (units per second)');

hold off;

%% (Test) Show removed static periods

windowDuration = 30;
noMovementThrehsold = 5;
staticWindow = getStaticPeriod(X_raw, Y_raw,...
    noMovementThreshold=noMovementThrehsold,windowDuration=windowDuration);

figure; tiledlayout(1,2);
nexttile;
% Plot the full trajectory in blue
plot(X_raw, Y_raw, 'b.-');  
hold on;
% Overlay the points marked for removal in red circles
plot(X_raw(staticWindow), Y_raw(staticWindow), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('X Position');
ylabel('Y Position');
legend('Recorded Trajectory','Removed Data');
title('Animal Trajectory with Removed Data Highlighted');

nexttile;
movingWindowInSamples = windowDuration * 20;
diffXinWindow = movsum(diff(X_raw),movingWindowInSamples);
diffYinWindow = movsum(diff(Y_raw),movingWindowInSamples);
distanceInWindow = diffXinWindow + diffYinWindow;
time = (1:length(distanceInWindow))./20 / 60;
plot(time, distanceInWindow); hold on;
yline(50,'r',LineWidth=2); hold on 
yline(-50,'r',LineWidth=2); hold on; 
plot(time(staticWindow),distanceInWindow(staticWindow),'r');
xlabel('Time (min)');
ylabel('Distane travel in 30s moving window');


%% (Test) Show static period criteria on result

cur_data = sessions(2).data(sessions(2).data.time >= 0,:);
for s = 1:100
    % Drop potential sleep time
    staticWindow = getStaticPeriod(X_raw,Y_raw,noMovementThreshold=s,windowDuration=30);

    Y_raw = cur_data.Item2_Y;
    Y = Y_raw(~staticWindow);
    stim = find(Y<=Y_midpoint);
    stim_pct(s) = length(stim)/length(Y) * 100;
end

plot(stim_pct)