% Shun_loadOptoSweepData.m
% 2025/09/18


%% Single session analysis

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project clamping/Recordings'))';
errorSessionIdx = [];

% Select anlaysis params
[analysisParams,canceled] = inputAnalysisParams(sessionList,...
                                reloadAll=false,...
                                recordLJ='[1 0 0]',...
                                plotPhotometry=true,...
                                withPhotometryNI=false);
if canceled; return; end
for s = 1:length(sessionList)
    analysisParams(s).rollingWindowTime = str2double(analysisParams(s).rollingWindowTime);
    analysisParams(s).recordLJ = eval(analysisParams(s).recordLJ);
end


% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList errorSessionIdx analysisParams taskList redStimPatternList blueStimPatternList withPhotometryNI plotPhotometry reloadAll
    
    idx = s; % if loop through all
    % idx = errorSessionIdx(s); % if loop through error sessions

    dirsplit = strsplit(sessionList{idx},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{idx},reloadAll=analysisParams(idx).reloadAll,...
            NISetup='clamp',...
            invertStim=false,...
            withPhotometryNI=analysisParams(idx).withPhotometryNI,photometryNI_mod=false,...
            recordLJ=analysisParams(idx).recordLJ,...
            rollingWindowTime=analysisParams(idx).rollingWindowTime,...
            followOriginal=false);
        analyzeSessions_optoSweeps(sessionList{idx},...
            laserSource='clamp',...
            durationList=[0.1, 0.5, 1],...
            bluePower=[20,30,40],...
            redPower=[45, 55],...
            redo=true,round=false,...
            plotPhotometry=analysisParams(idx).plotPhotometry,...
            shortTimeRange=[-1,5],longTimeRange=[-5,10]);
    catch ME
        errorSessionIdx = [errorSessionIdx;idx];
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end
close all
return