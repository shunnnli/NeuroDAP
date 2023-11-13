% Shun_loadSessionData
% Shun Li, 4/2/2022


%% Load files

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun'));

% Select anlaysis params
[analysisParams,canceled] = inputAnalysisParams(sessionList);
if canceled; return; end
recordLJList = nan(length(sessionList),3);
for s = 1:length(sessionList)
    analysisParams(s).rollingWindowTime = str2double(analysisParams(s).rollingWindowTime);
    analysisParams(s).recordLJ = eval(analysisParams(s).recordLJ);
end


% prompt = {'NI photometry:';'Reload All:';'Plot photometry:';'Plot licks:'};
% format = struct('type', {}, 'style', {}, 'items', {}, 'size', {});
% format(1,1).type   = 'list'; format(1,1).style = 'popupmenu';
% format(1,1).items  = {'true','false'}; format(1,1).size  = [150 20];
% format(2,1).type   = 'list'; format(2,1).style = 'popupmenu';
% format(2,1).items  = {'true','false'}; format(2,1).size  = [150 20];
% format(3,1).type   = 'list'; format(3,1).style = 'popupmenu';
% format(3,1).items  = {'true','false'}; format(3,1).size  = [150 20];
% format(4,1).type   = 'list'; format(4,1).style = 'popupmenu';
% format(4,1).items  = {'true','false'}; format(4,1).size  = [150 20];
% defaultanswers = {2, 2, 1, 1}; %index from the items in the list
% [answer, ~] = inputsdlg(prompt, 'Set analysis params', format, defaultanswers);
% Set analysis params
% withPhotometryNI = (answer{1} == 1);
% reloadAll = (answer{2} == 1);
% plotPhotometry = (answer{3} == 1); 
% plotLicks = (answer{4} == 1); 

% Select session params
[sessionParams,canceled] = inputSessionParams(sessionList);
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
            redo=true,round=false,performing=false,...
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
