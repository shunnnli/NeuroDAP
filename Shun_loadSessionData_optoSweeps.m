% Shun_loadSessionData
% Shun Li, 4/2/2022

% 2023/07/25
% Modified from original code for analyzing optoPair sessions

%% Load files for optoSweeps

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun'));

% Set session params
taskList = cell(size(sessionList));
stimPatternList = cell(size(sessionList));
for s = 1:length(sessionList)
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    
    % Determine stim pattern
    stim_prompt = {'Enter tested frequencies: ','Enter repetitions per pattern: '};
    definput = {'[30,35,40,45,50]','30'};
    stimPatternList{s} = inputdlg(stim_prompt,'Opto sweep params',[1,35],definput);
end


% Set edge cages
ni_photometry = true;
reload = false;

% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList taskList stimPatternList ni_photometry reload
    
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{s},ni_photometry=ni_photometry,reload=reload);
        analyzeSessions_optoSweeps(sessionList{s},taskList{s},stimPatternList{s});
    catch
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end
end

return