% Shun_PokeToUnclamp

%% Load FED3 data

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Load CSV (unfinished)
animalList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings/'),...
                         'Prompt','Select animal folder')';

%%
% Concatenate data
animals = cell(length(animalList),1);

% all animal folders
baseDir = '/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings/202503-PokeToUnclamp';

for a = 1:length(animalList)
    % Get all the path to csv files inside the folder
    sessionList = dir(fullfile(animalList{a}));
    sessionList = sessionList(~ismember({sessionList.name},{'.','..'}));
    % sessionList = ...
    data = cell(length(sessionList));
    for i = 1:length(sessionList)
        % Concatenate with sessionPath to get the full path
        sessionpath = strcat(animalList{a},filesep,sessionList(i).name);
        cur_data = readtable(sessionpath);
        cur_data.Time = minutes(cur_data.MM_DD_YYYYHh_mm_ss - cur_data.MM_DD_YYYYHh_mm_ss(1));
        data{i} = cur_data;
    end
    % Add all session data to the animal
    animals{a} = data;
end

%% Histogram of number of left & right poke across sessions

close all;

baseDir = '/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings/202503-PokeToUnclamp';
Average = true;

totalLeft = 0;
totalRight = 0;
totalSessions = 0; % why do these need to be initialized?animalDirs = dir(parentDir);

sessionList = dir(fullfile(animalList{a}));
sessionList = sessionList(~ismember({sessionList.name},{'.','..'})); % why would does this need to be stated again if in previous block

leftColor = [233 34 216]./255;
rightColor = [23 134 216]./255;

figure; hold on;

%initializeFig(.5, .5);
for a = 1: length(animalList)
    for a = 1:length(sessionList)
        sessionList = dir(fullfile(animalList{a}));
        sessionpath = strcat(animalList{a},filesep,sessionList(i).name); 
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
    end
      totalLeft = totalLeft + nLeftPoke;
        totalRight = totalRight + nRightPoke;
        totalSessions = totalSessions + 1;
end
 
xlabel('Session Number');
ylabel('Number of Pokes')
title ('Number of Left/Right Pokes by Session');
xticks(2*(1:length(data)));             
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