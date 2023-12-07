% Shun_loadSessionData
% Shun Li, 4/2/2022


%% Single session analysis

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'));
errorSessionIdx = [];

% Select anlaysis params
[analysisParams,canceled] = inputAnalysisParams(sessionList,...
                                reloadAll=false,...
                                recordLJ='[1 0 0]',...
                                plotPhotometry=true,...
                                plotBehavior=true,...
                                withPhotometryNI=false);
if canceled; return; end
for s = 1:length(sessionList)
    analysisParams(s).rollingWindowTime = str2double(analysisParams(s).rollingWindowTime);
    analysisParams(s).recordLJ = eval(analysisParams(s).recordLJ);
end 

% Select session params
[sessionParams,canceled] = inputSessionParams(sessionList,...
                                paradigm=2,...
                                reactionTime=2);
if canceled; return; end
taskList = cell(size(sessionList));
taskOptions = {'random','reward pairing','punish pairing'};
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
    clearvars -except s sessionList errorSessionIdx retryAttempt analysisParams sessionParams taskList stimPatternList withPhotometryNI plotPhotometry reloadAll plotLicks
    
    idx = s; % if loop through all
    % idx = errorSessionIdx(s); % if loop through error sessions

    dirsplit = strsplit(sessionList{idx},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{idx},reloadAll=analysisParams(idx).reloadAll,...
            invertStim=sessionParams(idx).OptoInverted,...
            withPhotometryNI=analysisParams(idx).withPhotometryNI,photometryNI_mod=false,...
            recordLJ=analysisParams(idx).recordLJ,...
            rollingWindowTime=analysisParams(idx).rollingWindowTime);
        analyzeSessions_optoPair(sessionList{idx},...
            task=taskList{idx},...
            stimPattern=stimPatternList{idx},...
            redo=false,round=false,performing=false,...
            analyzeTraces=true,...
            plotPhotometry=analysisParams(idx).plotPhotometry,...
            plotBehavior=analysisParams(idx).plotBehavior,...
            pavlovian=sessionParams(idx).Pavlovian,...
            reactionTime=sessionParams(idx).ReactionTime);
    catch ME
        errorSessionIdx = [errorSessionIdx;idx];
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end

return

%% Bulk analysis

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun/Project valence/Recordings'));
errorSessionIdx = [];
error()
% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList errorSessionIdx

    idx = s; % if loop through all
    % idx = errorSessionIdx(s); % if loop through error sessions
    
    dirsplit = strsplit(sessionList{idx},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{idx},reloadAll=true,...
            invertStim=false,...
            withPhotometryNI=false,photometryNI_mod=false,...
            recordLJ=[1 1 0],...
            rollingWindowTime=180);
        analyzeSessions_optoPair(sessionList{idx},...
            redo=false,round=false,performing=false,...
            analyzeTraces=true,...
            plotPhotometry=false,...
            plotBehavior=false,...
            pavlovian=false,...
            reactionTime=2);
    catch ME
        errorSessionIdx = [errorSessionIdx;idx];
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end
