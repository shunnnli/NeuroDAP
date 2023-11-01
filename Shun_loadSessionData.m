% Shun_loadSessionData
% Shun Li, 4/2/2022


%% Load files

clear; close all;
addpath(genpath(osPathSwitch('/Volumes/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));

% Select sessions via uipickfiles
sessionList = uipickfiles('FilterSpec',osPathSwitch('/Volumes/MICROSCOPE/Shun'));

% Select anlaysis params
prompt = {'NI photometry:';'Reload All:';'Plot photometry:';'Plot licks:'};
format = struct('type', {}, 'style', {}, 'items', {}, 'size', {});
format(1,1).type   = 'list'; format(1,1).style = 'popupmenu';
format(1,1).items  = {'true','false'}; format(1,1).size  = [150 20];
format(2,1).type   = 'list'; format(2,1).style = 'popupmenu';
format(2,1).items  = {'true','false'}; format(2,1).size  = [150 20];
format(3,1).type   = 'list'; format(3,1).style = 'popupmenu';
format(3,1).items  = {'true','false'}; format(3,1).size  = [150 20];
format(4,1).type   = 'list'; format(4,1).style = 'popupmenu';
format(4,1).items  = {'true','false'}; format(4,1).size  = [150 20];
defaultanswers = {2, 2, 1, 1}; %index from the items in the list
[answer, ~] = inputsdlg(prompt, 'Set analysis params', format, defaultanswers);
% Set analysis params
withPhotometryNI = (answer{1} == 1);
reloadAll = (answer{2} == 1);
plotPhotometry = (answer{3} == 1); 
plotLicks = (answer{4} == 1); 

% Select session params
taskList = cell(size(sessionList));
stimPatternList = cell(size(sessionList));
invertStim = nan(size(sessionList));
% Build dlg box
prompt_task = cell(length(sessionList),1);
prompt_opto = cell(length(sessionList),1);
format_task = struct('type', {}, 'style', {}, 'items', {}, 'size', {});
format_opto = struct('type', {}, 'style', {}, 'items', {}, 'size', {});
default_task = num2cell(ones(size(sessionList)));
default_opto = num2cell(ones(size(sessionList)));
for s = 1:length(sessionList)
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit

    prompt_task(s) = {sessionName(1:end-3)};
    prompt_opto(s) = {sessionName(1:end-3)};

    format_task(s,1).type   = 'list';
    format_task(s,1).style = 'popupmenu';
    format_task(s,1).items  = {'Random','Reward pairing','Punish pairing'};
    format_task(s,1).size  = [150 20];

    format_opto(s,1).type   = 'list';
    format_opto(s,1).style = 'popupmenu';
    format_opto(s,1).items  = {'Pattern inverted','Pattern normal','Trigger'};
    format_opto(s,1).size  = [150 20];
end
[answer_task, canceled1] = inputsdlg(prompt_task, 'Session task', format_task, default_task);
[answer_opto, canceled2] = inputsdlg(prompt_opto, 'Opto signal type', format_opto, default_opto);
if any([canceled1,canceled2]); return; end

% Set session params
for s = 1:length(sessionList)
    switch answer_task{s}
        case 1; taskList{s} = "random";
        case 2; taskList{s} = "reward pairing";
        case 3; taskList{s} = "punish pairing";
    end
    switch answer_opto{s}
        case 1; stimPatternList{s} = {'20','5','500'}; invertStim(s) = true;
        case 2; stimPatternList{s} = {'20','5','500'}; invertStim(s) = false;
        case 3; stimPatternList{s} = {'2','500','500'}; invertStim(s) = false;
    end
end

% Run each session
for s = 1:length(sessionList)
    close all;
    clearvars -except s sessionList taskList stimPatternList withPhotometryNI plotPhotometry reloadAll plotLicks invertStim
    
    dirsplit = strsplit(sessionList{s},filesep); 
    sessionName = dirsplit{end}; clear dirsplit
    try
        loadSessions(sessionList{s},reloadAll=reloadAll,invertStim=invertStim(s),...
            withPhotometryNI=withPhotometryNI,photometryNI_mod=false,...
            rollingWindowTime=180);
        analyzeSessions_optoPair(sessionList{s},taskList{s},stimPatternList{s},...
            redo=false,round=false,performing=false,...
            plotPhotometry=plotPhotometry,plotLicks=plotLicks,...
            pavlovian=true,reactionTime=1);
    catch ME
        disp(getReport(ME));
        warning(['Session ', sessionName, ' have an error, skipped for now!!!!']);
        continue
    end 
end

return
