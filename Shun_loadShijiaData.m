%% Sample analysis code for Shijia
% Shun Li, 2023/11/20

%% Load sessions

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun'));


%% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList plotPhotometry reloadAll plotLicks
    
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{s},reloadAll=true,...
            labjackSetup="Shijia",...
            invertStim=false,...
            withPhotometryNI=false,photometryNI_mod=false,...
            recordLJ=[1 0 0],...
            rollingWindowTime=180);
        analyzeSessions_Shijia(sessionList{s},...
            redo=false,...
            analyzeTraces=true,...
            plotPhotometry=true,...
            plotLicks=true);
    catch ME
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end