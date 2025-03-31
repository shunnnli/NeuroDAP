% Shun_PokeToUnclamp

%% Load FED3 data

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Load CSV (unfinished)
sessionPath = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings/202503-PokeToUnclamp/SL341'));

%%
% Concatenate data
data = cell(size(sessionPath));
for i = 1:length(sessionPath)
    cur_data = readtable(sessionPath{i});
    cur_data.Time = minutes(cur_data.MM_DD_YYYYHh_mm_ss - cur_data.MM_DD_YYYYHh_mm_ss(1));
    data{i} = cur_data;
end

%% Histogram of number of left & right poke across sessions

close all;

leftColor = [233 34 216]./255;
rightColor = [23 134 216]./255;

initializeFig(.5, .5);
for i = 1:length(data)
    cur_data = data{i};
    nLeftPoke = sum(strcmpi(cur_data.Event,'Left'));
    nRightPoke = sum(strcmpi(cur_data.Event,'Right'));
    disp(['Session ',num2str(i),...
          ': nLeftPoke = ',num2str(nLeftPoke),...
          '; nRightPoke = ',num2str(nRightPoke)]);
    plotScatterBar(2*i-0.45,nLeftPoke,style='bar',color=leftColor);
    plotScatterBar(2*i+0.45,nRightPoke,style='bar',color=rightColor);
end

xlabel('Session Number');
ylabel('Number of Pokes')
title ('Number of Left/Right Pokes by Session');
xticklabels(1:length(data));
legend({'Left Poke', 'Right Poke'}, 'Location', 'northeast');

%% Count number of left & right poke, first two hours only
    
% set color here:
leftColor = [233 34 216]./255;
rightColor = [23 134 216]./255;
timeCutoff = 120; % in minute

% Get number of nosepokes before cutoff
nNosePoke = zeros(nSessions,2); % first col: left poke; second col: right poke

for i = 1:length(data)

    
    % Filter data: Keep only rows within the first 120 minutes
    filteredData = cur_data(cur_data.Time <= timeCutoff, :);
    nLeftPoke = sum(strcmpi(filteredData.Event,'Left'));
    nRightPoke = sum(strcmpi(filteredData.Event,'Right'));
    disp(['Session ',num2str(i),...
          ': nLeftPoke = ',num2str(nLeftPoke),...
          '; nRightPoke = ',num2str(nRightPoke)]);
end 

% Plot data
initializeFig(.5, .5);
for i = 1:length(data)
    nLeftPoke = sum(strcmpi(cur_data.Event,'Left'));
    nRightPoke = sum(strcmpi(cur_data.Event,'Right'));
    disp(['Session ',num2str(i),...
          ': nLeftPoke = ',num2str(nLeftPoke),...
          '; nRightPoke = ',num2str(nRightPoke)]);
    plotScatterBar(2*i-0.45,nLeftPoke,style='bar',color=leftColor);
    plotScatterBar(2*i+0.45,nRightPoke,style='bar',color=rightColor);
end

xlabel('Session Number');
ylabel('Number of Pokes')
title ('Number of Left/Right Pokes by Session');
xticklabels({'1','2','3','4', '5', '6', '7'});
legend({'Left Poke', 'Right Poke'}, 'Location', 'northeast');

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