%% Sample analysis code for Shijia
% Shun Li, 2023/11/20

%% Load sessions

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
% addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\NeuroDAP\Methods'));
% addpath(genpath(osPathSwitch('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shijia\A_Code\Methods')));

% Select sessions via uipickfiles
% sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun'));
% sessionList = uipickfiles('FilterSpec',osPathSwitch('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shijia\'));
sessionList = uipickfiles('FilterSpec',osPathSwitch('C:\Shijia\Recordings'));

%% Ben - please uncomment the line below, change the folder if needed
% sessionList = uipickfiles('FilterSpec',osPathSwitch('C:\Shijia\Recordings'));

%% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList plotPhotometry reloadAll plotLicks
    
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{s},reloadAll=true,...
            labjackSetup="Shijia",...
            recordLJ=[1 0 0],...
            rollingWindowTime=180); % it has the concatLabjack_shijia which cleans the cue signal 
        analyzeSessions_Shijia(sessionList{s},...
            redo=false,...
            analyzeTraces=true,...
            plotPhotometry=true,...
            plotLicks=true); % it uses the cleaned labjack.cue and then clean the manual water, then feeds the events to getTrialTable_shijiaCatch where I clean the licks in each trial. also independently clean licks over the whole session to update labjack.lick in analyzeSessions_Shijia, before analyzeTraces 
    catch ME
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end