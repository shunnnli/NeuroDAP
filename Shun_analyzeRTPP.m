% Shun_analyzeRTPP

% 2025/01/31

% Plots trajectory of RTPP centroid and performs other analysis

%% Load the CSV file
filename = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project valence/Recordings'));
data = readtable(filename{1});

%% Extract columns
X = data.Item3_X;
Y = data.Item3_Y;
stim = strcmp(data.Item2,'True'); % Assuming Stim is a binary column (0 = Off, 1 = On)

% Plot the trajectory
figure;
hold on;

% Plot segments where stimulation is OFF
plot(X(stim == 0), Y(stim == 0), 'k', 'LineWidth', 1); % Black for OFF

% Plot segments where stimulation is ON
plot(X(stim == 1), Y(stim == 1), 'b', 'LineWidth', 1); % Red for ON

% Labels and title
xlabel('X Position');
ylabel('Y Position');
title('Animal Trajectory with Optogenetic Stimulation');

% Legend
legend({'Stimulation OFF', 'Stimulation ON'}, 'Location', 'best');

hold off;