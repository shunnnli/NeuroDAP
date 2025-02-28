% Shun_analyzeRTPP

% 2025/01/31

% Plots trajectory of RTPP centroid and performs other analysis

%% Load the CSV file

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

filename = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'));
data = readtable(filename{1});

Fs = 30; % frame rate

%% Extract columns
X = data.Item3_X;
Y = data.Item3_Y;
stim = strcmp(data.Item2,'True'); % Assuming Stim is a binary column (0 = Off, 1 = On)

% Remove artifact (to remove wall edges)
% x_artifact = X < 5.5 & X > 5;
% y_artifact = Y < 650 & Y > 600;
% artifact_idx = intersect(find(x_artifact),find(y_artifact));
% X(artifact_idx) = nan;
% Y(artifact_idx) = nan;

% Plot the trajectory
figure; hold on;
plot(X(stim == 0), Y(stim == 0), 'k', 'LineWidth', 1); % Black for OFF
plot(X(stim == 1), Y(stim == 1), 'b', 'LineWidth', 1); % Red for ON

% Labels and title
xlim([300,650]); xlabel('X Position');
ylim([0,700]);ylabel('Y Position');
title('Animal Trajectory');
legend({'Stim OFF', 'Stim ON'}, 'Location', 'northeast');


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

stim_pct = sum(stim)/length(stim) * 100;

% % Plot bar plot
% figure; hold on;
