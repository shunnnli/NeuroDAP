% Shun_PokeToUnclamp

%% Load FED3 data

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Load CSV (unfinished)
sessionPath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'));

% Concatenate data
data = cell(size(sessionPath));
for i = 1:length(sessionPath)
    data{i} = readtable(sessionPath{i});
end

%% Count number of left & right poke

for i = 1:length(data)
    cur_data = data{i};
    nLeftPoke = sum(strcmpi(cur_data.Event,'Left'));
    nRightPoke = sum(strcmpi(cur_data.Event,'Right'));
    disp(['Session ',num2str(i),...
          ': nLeftPoke = ',num2str(nLeftPoke),...
          '; nRightPoke = ',num2str(nRightPoke)]);
end

%% Calculate poke rate

binWindow = 5; % in min

% Plot data
initializeFig(0.5,0.5); tiledlayout(length(data),2);

for i = 1:length(data)
    poke_rate = getPokeRate(data{i},method='bin',binWindow=binWindow);
    poke_rate_smoothed = getPokeRate(data{i},binWindow=binWindow);
    x = binWindow * (1:length(poke_rate));

    nexttile;
    plot(x,poke_rate);
    xlabel('Time in Minutes');
    ylabel('Poke Rate')
    title ('Poke Rate Over Time Using Bins');
    
    nexttile;
    plot(x,poke_rate_smoothed);
    xlabel('Time in Minutes');
    ylabel('Poke Rate');
    title('Poke Rate Over Time Using Moving Window');
end