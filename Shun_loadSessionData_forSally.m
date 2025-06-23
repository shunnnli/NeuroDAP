% Shun_loadSessionData_forSally
% Shun Li, 4/2/2022

%% Single session analysis

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project misc/Recordings'))';
errorSessionIdx = [];

% Select anlaysis params
[analysisParams,canceled] = inputAnalysisParams(sessionList,...
                                reloadAll=false,...
                                recordLJ='[1 1 0]',...
                                plotPhotometry=true,...
                                plotBehavior=true,...
                                withPhotometryNI=false);
if canceled; return; end
for s = 1:length(sessionList)
    analysisParams(s).rollingWindowTime = str2double(analysisParams(s).rollingWindowTime);
    analysisParams(s).recordLJ = eval(analysisParams(s).recordLJ);
end


% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList errorSessionIdx analysisParams sessionParams taskList redStimPatternList blueStimPatternList withPhotometryNI plotPhotometry reloadAll
    
    idx = s; % if loop through all
    % idx = errorSessionIdx(s); % if loop through error sessions

    dirsplit = strsplit(sessionList{idx},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{idx},reloadAll=analysisParams(idx).reloadAll,...
            invertStim=true,...
            withPhotometryNI=analysisParams(idx).withPhotometryNI,photometryNI_mod=false,...
            recordLJ=analysisParams(idx).recordLJ,...
            rollingWindowTime=analysisParams(idx).rollingWindowTime,...
            followOriginal=false);
        analyzeSessions_forSally(sessionList{idx},...
            task='',...
            redo=true,round=false,...
            analyzeTraces=true,...
            plotPhotometry=analysisParams(idx).plotPhotometry,...
            plotBehavior=analysisParams(idx).plotBehavior);
    catch ME
        errorSessionIdx = [errorSessionIdx;idx];
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end
close all
return
