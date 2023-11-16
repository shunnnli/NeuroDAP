% Shun_loadSessionData
% Shun Li, 4/2/2022


%% Load files

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun'));

% Select anlaysis params
[analysisParams,canceled] = inputAnalysisParams(sessionList,...
                                reloadAll=false,...
                                recordLJ='[1 0 0]',...
                                plotPhotometry=true,...
                                plotLicks=true,...
                                withPhotometryNI=false);
if canceled; return; end
for s = 1:length(sessionList)
    analysisParams(s).rollingWindowTime = str2double(analysisParams(s).rollingWindowTime);
    analysisParams(s).recordLJ = eval(analysisParams(s).recordLJ);
end 

% Select session params
[sessionParams,canceled] = inputSessionParams(sessionList,...
                                paradigm=3,...
                                reactionTime=2);
if canceled; return; end
taskList = cell(size(sessionList));
taskOptions = {"random","reward pairing","punish pairing"};
stimPatternList = cell(size(sessionList));
for s = 1:length(sessionList)
    taskList{s} = taskOptions{sessionParams(s).Paradigm};
    stimPatternList{s} = {sessionParams(s).OptoPulseFreq,sessionParams(s).OptoPulseDuration,sessionParams(s).OptoStimDuration};
    sessionParams(s).ReactionTime = str2double(sessionParams(s).ReactionTime);
    sessionParams(s).minLicks = str2double(sessionParams(s).minLicks);
end


% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList analysisParams sessionParams taskList stimPatternList withPhotometryNI plotPhotometry reloadAll plotLicks
    
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{s},reloadAll=analysisParams(s).reloadAll,...
            invertStim=sessionParams(s).OptoInverted,...
            withPhotometryNI=analysisParams(s).withPhotometryNI,photometryNI_mod=false,...
            recordLJ=analysisParams(s).recordLJ,...
            rollingWindowTime=analysisParams(s).rollingWindowTime);
        analyzeSessions_optoPair(sessionList{s},taskList{s},stimPatternList{s},...
            redo=false,round=false,performing=false,...
            analyzeTraces=true,...
            plotPhotometry=analysisParams(s).plotPhotometry,...
            plotLicks=analysisParams(s).plotLicks,...
            pavlovian=sessionParams(s).Pavlovian,reactionTime=sessionParams(s).ReactionTime);
    catch ME
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end

return
