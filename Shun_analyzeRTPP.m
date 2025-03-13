% Shun_analyzeRTPP

% 2025/01/31

% Plots trajectory of RTPP centroid and performs other analysis

%% Load the CSV file

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

filename = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'));

sessionList = dir(filename{1});
sessionList = sessionList(~ismember({sessionList.name},{'.','..'}));
nSessions = length(sessionList);

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
end
Fs = 30; % frame rate

%% Plot summary figure

initializeFig(0.5,1); tiledlayout(1,nSessions);

for s = 1:nSessions
    cur_data = sessions(s).data;
    cur_name = sessions(s).name;

    X = cur_data.Item3_X;
    Y = cur_data.Item3_Y;
    Y_midpoint = (min(Y)+max(Y))/2;
    
    % Plot the trajectory
    nexttile;
    if contains(cur_name,'Habit')
        plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), 'k', 'LineWidth', 2); hold on; % Black for OFF
        plot(X(Y<Y_midpoint), Y(Y<Y_midpoint), 'b', 'LineWidth', 2); % Red for ON
        legend({'Right', 'Left'}, 'Location', 'northeast');
    else
        plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), 'b', 'LineWidth', 2); hold on; % Black for OFF
        plot(X(Y<Y_midpoint), Y(Y<Y_midpoint), 'k', 'LineWidth', 2); % Red for ON
        legend({'Stim OFF', 'Stim ON'}, 'Location', 'northeast');
    end
    
    xlim([300,650]); xlabel('X Position');
    ylim([0,700]);ylabel('Y Position');
    title(cur_name);
end


%% Heatmap

% Define heatmap grid
squareSize = 10; % in pixels
edgesX = linspace(300, 640, round(range(X)/squareSize));
edgesY = linspace(0, 700, round(range(Y)/squareSize));

% Count occurences of points in each grid bin
[counts, edgesX, edgesY] = histcounts2(X, Y, edgesX, edgesY);
timeInBins = counts * 1/Fs;
% log_counts = log(counts);
% log_times = log(timeInBins);
gamma_counts = counts .^ 0.5;

% Convert count to time/

%4. plot
figure;
imagesc(edgesX(1:end-1), edgesY(1:end-1), gamma_counts');
set(gca, 'YDir', 'normal'); 
% set(gca,'ColorScale','log');
colormap(sky); colorbar;
xlim([300,640]); xlabel('X Position');
ylim([0,700]);ylabel('Y Position');
title ('Heatmap');


%% Bar plot of %time at opto side

stim_pct = zeros(1,nSessions);
for s = 1:nSessions
    cur_data = sessions(s).data;
    cur_name = sessions(s).name;

    X = cur_data.Item3_X;
    Y = cur_data.Item3_Y;
    Y_midpoint = (min(Y)+max(Y))/2;

    % stim = strcmp(cur_data.Item2,'True'); % Assuming Stim is a binary column (0 = Off, 1 = On)
    stim = find(Y>=Y_midpoint);
    stim_pct(s) = length(stim)/length(Y) * 100;
end

% Plot bar plot
initializeFig(0.3,0.5);
bar(1,stim_pct(1),'k'); hold on
bar(2,stim_pct(2),'b'); hold on
xticks([1 2]);
xticklabels({'Habituation','RTPP1'});
ylabel('Time spent in stimulated side (%)');


