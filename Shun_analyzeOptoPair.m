% Shun_analyzeOptoPair
% 2023/04/25

% Outputs a .mat file with aligned traces for each animal and behavior

%% Setup

clear; close all;
% addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Analysis\Methods'));
[~,~,~,blueGreenYellow,blueWhiteRed,~,bluePurpleRed,purpleWhiteRed] = loadColors;
r2p_cmap = getColormap([255, 50, 58],[0 0 0],500,'midcol',[255 255 255]);
p2r_cmap = getColormap([241 160 255],[0 0 0],500,'midcol',[255 255 255]);

% Define result directory
resultspath = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\project valence\Results\';
% resultspath = uigetdir('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\','Select Results folder directory');

% For pairing, sessions needs arrange in days then in animal number

%% Notes for each animal

%{
Common stuff:
    1. need to remove not performaning trials
    2. try group by days or group by performance
%}

%{
EP-LHb ChrimsonR, NAc dLight
- SL043
    - Baseline: 20221205(water, stim), 20221206 (water, stim only), 20230102,
    20230103 (water, airpuff)
    - Reward: 20221210, 20221211, 20221212 (remove surprise trial for long)
    - Punish: 20221219, 20221220, 20221221
    - Reward: 20230104 (not learn), 0105 (learned), 0106

- SL044
    - Baseline: 20221205(water, stim), 20221206 (water, stim only), 20230102,
    20230103 (water, weak airpuff)
    - Reward: 20221210 (not learn), 20221211 (ramp), 20221212 (ramp)
    - Baseline-Punish: 20230104 (only pair learn), 0105, 0106

- SL045
    - Baseline: 20220106 (water, stim), 20221219 (water)
    - Reward: water decrase DA during reward, exclude from analysis (without enough conditioning)
    - Punish: exclude from analysis due to the same reasoning

- SL046
    - Baseline: 20220106 (stim), 20221216 (water, 1-60), 20230102, 03 (water, airpuff)
    - Reward: 20221210 (not learn), 20221211 (not learn), 20221212 (ramp)
    - Punish: 20221215 (not learn), 20221216 (flat), 1217 (slight up)
    - Reward: 20230104 (ramp), 0105 (ramp), 0106 (ramp)
%}

%{
EP-LHb ChrimsonR, LHb GCaMP, NAc dLight
- SL060
    - Baseline: 20230326, 20230327 (no DA recording), 20230328 (no DA dip)
    - Reward: 20230415, 20230416, 20230417
    - Punish: 20230418, 20230419, 20230420
    - Reward: 20230425, 20230426, 20230427

- SL061
    - Baseline: 
    - Reward:
    - Punish:
    - Reward:

- SL062
    - Baseline: 20230326, 20230327, 20230328 (0328 with DA, other just LHb)
    - Reward: 20230415, 20230416, 20230417
    - Punish: 20230418, 20230419, 20230420
    - Reward:

- SL063
    - Baseline: 20230326, 20230327, 20230328 (just LHb)
    - Reward: 20230415, 20230416, 20230417
    - Punish: 20230418, 20230419, 20230420
    - Reward: 20230426 (weird LHb), 20230427, 20230428

- SL064
    - Baseline: 20230326, 20230327, 20230328
    - Reward: 20230415, 20230416, 20230417
    - Punish: 20230418, 20230419, 20230420
    - Reward: 20230424, 20230425, 20240426 (only LHb for 0426)

- SL065
    - Baseline: 
    - Reward:
    - Punish:
    - Reward:

- SL066
    - Baseline: 20230327, 20230328, 20230326 (optional: werid DA signal)
    - Reward: 20230415, 20230416, 20230417
    - Punish: 20230418, 20230419, 20230420
    - Reward: 20230423, 20230424, 20230425

- SL067
    - Baseline: 20230326, 20230327, 20230328 (LHb water signal weird)
    - Reward: 20230415 (ramp)
    - exclude from analysis

- SL068
    - Baseline: 20230326, 20230328
    - Reward: 20230329 (not learn), 0330 (only learned pair), 0331, 0401
    - Punish: 20230402, 0403, 0404 (not learned)
    - Reward: 20230405 (not learn), 0406, 0407

- SL069
    - No LHb signal, increase DA for airpuff, exclude from analysis

%}


%% ALL: Load processed list and struct 
load(strcat(resultspath,'\','sessions.mat'),'b2rList','r2pList','p2rList','baselineList');

% load(strcat(resultspath,'\','baseline.mat'),'baseline');
% load(strcat(resultspath,'\','pairing.mat'),'pairing');
% load(strcat(resultspath,'\','animals.mat'),'animals');

%% ALL: add sessions

addList = uipickfiles('FilterSpec','\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');

% Decide which task it is
tasklist = {'baseline','baseline->reward','reward->punish','punish->reward'};
[idx,selected] = listdlg('PromptString',{'Select a task for these sessions',...
    'Only one task can be selected at a time',''},...
    'SelectionMode','single',...
    'ListString',tasklist);

if selected
    if strcmp(tasklist{idx},'baseline')
        % Update session list
        load(strcat(resultspath,'\','sessions.mat'),'baselineList');
        baselineList = horzcat(baselineList,addList);
        save(strcat(resultspath,'\','sessions'),'baselineList','-append');
        
        % Update related struct
        if ~exist('baseline','var')
            load(strcat(resultspath,'\','baseline.mat'),'baseline');
            disp('Finished: baseline struct loaded');
        end
        baseline = addSessionsToStruct(addList,baseline,'baseline');
        
    elseif strcmp(tasklist{idx},'baseline->reward')
        % Update session list
        load(strcat(resultspath,'\','sessions.mat'),'b2rList');
        b2rList = horzcat(b2rList,addList);
        save(strcat(resultspath,'\','sessions'),'b2rList','-append');
        
        % Update related struct
        if ~exist('pairing','var')
            load(strcat(resultspath,'\','pairing.mat'),'pairing');
            disp('Finished: pairing struct loaded');
        end
        pairing = addSessionsToStruct(addList,pairing,'baseline->reward');
        
    elseif strcmp(tasklist{idx},'reward->punish')
        % Update session list
        load(strcat(resultspath,'\','sessions.mat'),'r2pList');
        r2pList = horzcat(r2pList,addList);
        save(strcat(resultspath,'\','sessions'),'r2pList','-append');
        
        % Update related struct
        if ~exist('pairing','var')
            load(strcat(resultspath,'\','pairing.mat'),'pairing');
            disp('Finished: pairing struct loaded');
        end
        pairing = addSessionsToStruct(addList,pairing,'reward->punish');
        
    elseif strcmp(tasklist{idx},'punish->reward')
        % Update session list
        load(strcat(resultspath,'\','sessions.mat'),'p2rList');
        p2rList = horzcat(p2rList,addList);
        save(strcat(resultspath,'\','sessions'),'p2rList','-append');
        
        % Update related struct
        if ~exist('pairing','var')
            load(strcat(resultspath,'\','pairing.mat'),'pairing');
            disp('Finished: pairing struct loaded');
        end
        pairing = addSessionsToStruct(addList,pairing,'punish->reward');
    end
end

%% Baseline: Select baseline files
%{
Analyze baseline, reward, punish separately, meaning one sessionList for
each condition
%}
baselineList = uipickfiles('FilterSpec','\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
save(strcat(resultspath,'\','baseline'),'baselineList','-v7.3');

%% Baseline: Load or create baseline (random reward/punish/stim) struct
%{
Similar to the cells struct in analyzeSlice
each row is a session
%}

baseline = struct();

% Loading all data into cell
for i = 1:length(baselineList)
    
    % Get session information
    dirsplit = strsplit(baselineList{i},filesep); 
    sessionName = dirsplit{end};
    namesplit = strsplit(sessionName,'-');
    sessionAnimal = namesplit{2};
    
    % Store session information in struct
    baseline(i).mouse = sessionAnimal;
    baseline(i).session = sessionName;

    % Load recording data
    % Some session might not have LJ or NI photometry, initialize related
    % variable to 0 to avoid error
    rollingGreenLP = 0; timePhotometry = 0; photometryNI = 0;
    trials = 0; blueLaser = 0; airpuff_rounded = 0; firstPulse = 0;
    load(strcat(baselineList{i},'\','sync_',sessionName,'.mat'),...
        'trials','params','rightLick','blueLaser',...
        'timeNI','timePhotometry','rollingGreenLP','photometryNI',...
        'rightSolenoid','leftTone','airpuff_rounded','airpuff',...
        'firstPulse','redLaser');
    if airpuff_rounded == 0; airpuff_rounded = airpuff; end
    if firstPulse == 0; firstPulse = find(redLaser); end
    
    % For old sessions
    if ~isfield(params.sync,'timeNI'); params.sync.timeNI = timeNI; end
    if ~isfield(params.sync,'timePhotometry'); params.sync.timePhotometry = timePhotometry; end
    if ~isfield(params.sync,'behaviorFs'); params.sync.behaviorFs = 10000; end
    
    % Store session recordings in struct
    baseline(i).params = params;
    baseline(i).trials = trials;
    baseline(i).licks = find(rightLick);
    baseline(i).water = find(rightSolenoid);
    baseline(i).airpuff = find(airpuff_rounded);
    baseline(i).tone = find(leftTone);
    baseline(i).stim = firstPulse;
    baseline(i).randomShutter = round(rand([500,1])*length(blueLaser));
    baseline(i).photometryLJ = rollingGreenLP;
    baseline(i).photometryNI = photometryNI;
    
    % Initialize analysis params (window to analyze)
    % Behavior analysis window
    targetStruct(i).analysis.behavior.water = 1:length(targetStruct(i).water);
    targetStruct(i).analysis.behavior.airpuff = 1:length(targetStruct(i).airpuff);
    targetStruct(i).analysis.behavior.stim = 1:length(targetStruct(i).stim);
    targetStruct(i).analysis.behavior.tone = 1:length(targetStruct(i).tone);
    % DA analysis window
    targetStruct(i).analysis.DA.water = 1:length(targetStruct(i).water);
    targetStruct(i).analysis.DA.airpuff = 1:length(targetStruct(i).airpuff);
    targetStruct(i).analysis.DA.stim = 1:length(targetStruct(i).stim);
    targetStruct(i).analysis.DA.tone = 1:length(targetStruct(i).tone);
    % LHb analysis window
    targetStruct(i).analysis.LHb.water = 1:length(targetStruct(i).water);
    targetStruct(i).analysis.LHb.airpuff = 1:length(targetStruct(i).airpuff);
    targetStruct(i).analysis.LHb.stim = 1:length(targetStruct(i).stim);
    targetStruct(i).analysis.LHb.tone = 1:length(targetStruct(i).tone);

    clear dirsplit namesplit params rollingGreenLP photometryNI ...
        trials rightLick airpuff_rounded timeNI timePhotometry leftTone ...
        rightSolenoid firstPulse airpuff redLaser blueLaser
    disp(['Session ',sessionName,' loaded']);
end

%% Baseline: loop through struct

for i = 1:length(baselineList)
    % Get session information
    dirsplit = strsplit(baselineList{i},filesep); 
    sessionName = dirsplit{end};
    load(strcat(baselineList{i},'\','sync_',sessionName,'.mat'),'rightLick');
    baseline(i).licks = find(rightLick==1);
    disp(i);
end

%% Baseline: 1.0 Behavior: find lick traces for baseline

timeRange = [-1,3]; lick_binSize = 0.1;

for i = 1:length(baseline)  
    trials = baseline(i).trials; 
    params = baseline(i).params; 
    rightLick = baseline(i).licks;

    % getLicks by trial type
    [waterLickRate,~,~] = getLicks(timeRange,baseline(i).water,lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
    [airpuffLickRate,~,~] = getLicks(timeRange,baseline(i).airpuff,lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
    [stimLickRate,~,~] = getLicks(timeRange,baseline(i).stim,lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
    [toneLickRate,~,~] = getLicks(timeRange,baseline(i).tone,lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
                            
    % Save results for this session
    baseline(i).lickRate.timeRange = timeRange;
    %baseline(i).lickRaster.timeRange = timeRange;
    baseline(i).lickRate.lick_binSize = lick_binSize;
    %baseline(i).lickRaster.lick_binSize = lick_binSize;
    baseline(i).lickRate.water = waterLickRate;
    %baseline(i).lickRaster.water = waterLicks;
    baseline(i).lickRate.airpuff = airpuffLickRate;
    %baseline(i).lickRaster.airpuff = airpuffLicks;
    baseline(i).lickRate.stim = stimLickRate;
    %baseline(i).lickRaster.stim = stimLicks;
    baseline(i).lickRate.tone = toneLickRate;
    %baseline(i).lickRaster.tone = toneLicks;

    disp(['Session ',baseline(i).session,' calculated']);
    clear stimLickRate toneLickRate waterLickRate airpuffLickRate
end
disp('Lick rate traces genearted');

%% Baseline: 1.1 Behavior: concat lick rate vs time

waterLickRate_all = [];
airpuffLickRate_all = [];
stimLickRate_all = [];
toneLickRate_all = [];

for i = 1:length(baseline)
    % Water
    analysisRange = baseline(i).analysis.behavior.water;
    waterLickRate_all = [waterLickRate_all; baseline(i).lickRate.water(analysisRange,:)];
    % Airpuff
    analysisRange = baseline(i).analysis.behavior.airpuff;
    airpuffLickRate_all = [airpuffLickRate_all; baseline(i).lickRate.airpuff(analysisRange,:)];
    % Stim
    analysisRange = baseline(i).analysis.behavior.stim;
    stimLickRate_all = [stimLickRate_all; baseline(i).lickRate.stim(analysisRange,:)];
    % Tone
    analysisRange = baseline(i).analysis.behavior.tone;
    toneLickRate_all = [toneLickRate_all; baseline(i).lickRate.tone(analysisRange,:)];
end

% Remove nan rows
waterLickRate_all = rmmissing(waterLickRate_all);
airpuffLickRate_all = rmmissing(airpuffLickRate_all);
stimLickRate_all = rmmissing(stimLickRate_all);
toneLickRate_all = rmmissing(toneLickRate_all);

disp('Lick rate traces concatenated');

%% Baseline: 2.1.0 Photometry: get baseline DA/LHb traces

timeRange = [-1,3];

for i = 1:length(baseline)  
    trials = baseline(i).trials; 
    params = baseline(i).params; 

    [DA_water,~] = plotTraces(baseline(i).water,timeRange,baseline(i).photometryLJ,bluePurpleRed(1,:),params,plot=false);
    [DA_airpuff,~] = plotTraces(baseline(i).airpuff,timeRange,baseline(i).photometryLJ,[0.2, 0.2, 0.2],params,plot=false);
    [DA_tone,~] = plotTraces(baseline(i).tone,timeRange,baseline(i).photometryLJ,bluePurpleRed(350,:),params,plot=false);
    [DA_stim,~] = plotTraces(baseline(i).stim,timeRange,baseline(i).photometryLJ,bluePurpleRed(end,:),params,plot=false);
    [DA_randomShutter,~] = plotTraces(baseline(i).randomShutter,timeRange,baseline(i).photometryLJ,bluePurpleRed(end,:),params,plot=false);
    
    [LHb_water,~] = plotTraces(baseline(i).water,timeRange,baseline(i).photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni',plot=false);
    [LHb_airpuff,~] = plotTraces(baseline(i).airpuff,timeRange,baseline(i).photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni',plot=false);
    [LHb_tone,~] = plotTraces(baseline(i).tone,timeRange,baseline(i).photometryNI,bluePurpleRed(350,:),params,photometrySystem='ni',plot=false);
    [LHb_stim,~] = plotTraces(baseline(i).stim,timeRange,baseline(i).photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni',plot=false);
    [LHb_randomShutter,~] = plotTraces(baseline(i).randomShutter,timeRange,baseline(i).photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni',plot=false);
    
    % Store results
    baseline(i).traces.DA_water = DA_water;
    baseline(i).traces.DA_airpuff = DA_airpuff;
    baseline(i).traces.DA_tone = DA_tone;
    baseline(i).traces.DA_stim = DA_stim;
    baseline(i).traces.DA_randomShutter = DA_randomShutter;
    baseline(i).traces.LHb_water = LHb_water;
    baseline(i).traces.LHb_airpuff = LHb_airpuff;
    baseline(i).traces.LHb_tone = LHb_tone;
    baseline(i).traces.LHb_stim = LHb_stim;
    baseline(i).traces.LHb_randomShutter = LHb_randomShutter;
    baseline(i).traces.timeRange = timeRange;

    disp(['Session ',baseline(i).session,' calculated']);
    clear DA_water DA_airpuff DA_tone DA_stim LHb_water LHb_airpuff LHb_tone LHb_stim
end

disp('DA/LHb traces generated');

%% Baseline: 2.1.1 Photometry: concat baseline DA/LHb traces

DA_water_all = []; LHb_water_all = [];
DA_airpuff_all = []; LHb_airpuff_all = [];
DA_tone_all = []; LHb_tone_all = [];
DA_stim_all = []; LHb_stim_all = [];
DA_randomShutter_all = []; LHb_randomShutter_all = [];

for i = 1:length(baseline)
    % For DA
    % Water
    analysisRange = baseline(i).analysis.DA.water;
    DA_water_all = [DA_water_all; baseline(i).traces.DA_water(analysisRange,:)];
    % Airpuff
    analysisRange = baseline(i).analysis.DA.airpuff;
    DA_airpuff_all = [DA_airpuff_all; baseline(i).traces.DA_airpuff(analysisRange,:)];
    % Tone
    analysisRange = baseline(i).analysis.DA.tone;
    DA_tone_all = [DA_tone_all; baseline(i).traces.DA_tone(analysisRange,:)];
    % Stim
    analysisRange = baseline(i).analysis.DA.stim;
    DA_stim_all = [DA_stim_all; baseline(i).traces.DA_stim(analysisRange,:)];
    % Random shutter
    DA_randomShutter_all = [DA_randomShutter_all; baseline(i).traces.DA_randomShutter];
    
    % For LHb
    % Water
    analysisRange = baseline(i).analysis.LHb.water;
    LHb_water_all = [LHb_water_all; baseline(i).traces.LHb_water(analysisRange,:)];
    % Airpuff
    analysisRange = baseline(i).analysis.LHb.airpuff;
    LHb_airpuff_all = [LHb_airpuff_all; baseline(i).traces.LHb_airpuff(analysisRange,:)];
    % Tone
    analysisRange = baseline(i).analysis.LHb.tone;
    LHb_tone_all = [LHb_tone_all; baseline(i).traces.LHb_tone(analysisRange,:)];
    % Stim
    analysisRange = baseline(i).analysis.LHb.stim;
    LHb_stim_all = [LHb_stim_all; baseline(i).traces.LHb_stim(analysisRange,:)];
    % Random shutter
    LHb_randomShutter_all = [LHb_randomShutter_all; baseline(i).traces.LHb_randomShutter];
end

% Remove nan entries
DA_water_all = rmmissing(DA_water_all);
DA_airpuff_all = rmmissing(DA_airpuff_all);
DA_tone_all = rmmissing(DA_tone_all);
DA_stim_all = rmmissing(DA_stim_all);
DA_randomShutter_all = rmmissing(DA_randomShutter_all);
LHb_water_all = rmmissing(LHb_water_all);
LHb_airpuff_all = rmmissing(LHb_airpuff_all);
LHb_tone_all = rmmissing(LHb_tone_all);
LHb_stim_all = rmmissing(LHb_stim_all);
LHb_randomShutter_all = rmmissing(LHb_randomShutter_all);

disp('DA/LHb traces concatenated');

%% Baseline: 3.0 Plot all traces (3 subplots of DA, LHb, Licks)
timeRange = baseline(1).traces.timeRange;
mouseNum = numel(unique({baseline.mouse}));

initializeFig(1,0.33);
tiledlayout(1,3);
nexttile;
binNumber = size(baseline(1).traces.DA_water,2);
t = linspace(timeRange(1),timeRange(2),binNumber);
plotSEM(t,DA_water_all,bluePurpleRed(1,:));
plotSEM(t,DA_airpuff_all,[0.2, 0.2, 0.2]);
plotSEM(t,DA_tone_all,bluePurpleRed(350,:));
plotSEM(t,DA_stim_all,bluePurpleRed(end,:));
plotEvent('',0,'r');
taskLegend = {['Water (n=',num2str(size(DA_water_all,1)),', ',num2str(mouseNum),'mice)'],...
            ['Airpuff (n=',num2str(size(DA_airpuff_all,1)),', ',num2str(mouseNum),'mice)'],...
            ['Tone (n=',num2str(size(DA_tone_all,1)),', ',num2str(mouseNum),'mice)'],...
            ['Stim (n=',num2str(size(DA_stim_all,1)),', ',num2str(mouseNum),'mice)']};
xlabel('Time (s)'); ylabel('NAc dLight1.3b signal (z-score)'); 
legend(taskLegend,'Location','Northeast');
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');

nexttile;
binNumber = size(baseline(1).traces.LHb_water,2);
t = linspace(timeRange(1),timeRange(2),binNumber);
plotSEM(t,LHb_water_all,bluePurpleRed(1,:));
plotSEM(t,LHb_airpuff_all,[0.2, 0.2, 0.2]);
plotSEM(t,LHb_tone_all,bluePurpleRed(350,:));
plotSEM(t,LHb_stim_all,bluePurpleRed(end,:));
plotEvent('',0,'r');
taskLegend = {['Water (n=',num2str(size(LHb_water_all,1)),', ',num2str(mouseNum),'mice)'],...
            ['Airpuff (n=',num2str(size(LHb_airpuff_all,1)),', ',num2str(mouseNum),'mice)'],...
            ['Tone (n=',num2str(size(LHb_tone_all,1)),', ',num2str(mouseNum),'mice)'],...
            ['Stim (n=',num2str(size(LHb_stim_all,1)),', ',num2str(mouseNum),'mice)']};
xlabel('Time (s)'); ylabel('LHb GCaMP8m signal (z-score)'); 
legend(taskLegend,'Location','Northeast');
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');

% Plot lick rate
nexttile
timeRange = baseline(1).lickRate.timeRange;
binNumber = size(baseline(1).lickRate.water,2);
t = linspace(timeRange(1),timeRange(2),binNumber);
plotSEM(t,waterLickRate_all,bluePurpleRed(1,:));
plotSEM(t,airpuffLickRate_all,[0.2, 0.2, 0.2]);
plotSEM(t,stimLickRate_all,bluePurpleRed(350,:));
plotSEM(t,toneLickRate_all,bluePurpleRed(end,:));
plotEvent('',0,'r');
taskLegend = {['Water (n=',num2str(length(waterLickRate_all)),', ',num2str(mouseNum),'mice)'],...
            ['Airpuff (n=',num2str(length(airpuffLickRate_all)),', ',num2str(mouseNum),'mice)'],...
            ['Tone (n=',num2str(length(toneLickRate_all)),', ',num2str(mouseNum),'mice)'],...
            ['Stim (n=',num2str(length(stimLickRate_all)),', ',num2str(mouseNum),'mice)']};
xlabel('Time (s)'); ylabel('Licks/s'); 
legend(taskLegend,'Location','Northeast');
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
saveFigures(gcf,'baseline_all',resultspath);

%% Baseline: 3.1 Plot all traces (table form)
initializeFig(1,1);
tiledlayout(3,4);
mouseNum = numel(unique({baseline.mouse}));

% DA
timeRange = baseline(1).traces.timeRange;
binNumber = size(baseline(1).traces.DA_water,2);
t = linspace(timeRange(1),timeRange(2),binNumber);
ylimit = [-1,2]; 
% textloc_randomShutter = [timeRange-0.1, -0.05];
% DA, water
nexttile;
plotSEM(t,DA_randomShutter_all,[.7 .7 .7]); 
plotSEM(t,DA_water_all,bluePurpleRed(1,:));
plotEvent('',0,bluePurpleRed(1,:)); 
taskLegend = {['Water (n=',num2str(size(DA_water_all,1)),', ',num2str(mouseNum),' mice)']};
% xlabel('Time (s)'); ylabel('NAc dLight1.3b signal (z-score)'); 
ylim(ylimit); axis off
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% DA, airpuff
nexttile;
plotSEM(t,DA_randomShutter_all,[.7 .7 .7]);
plotSEM(t,DA_airpuff_all,[0.2, 0.2, 0.2]); 
plotEvent('',0,[0.2, 0.2, 0.2]);
taskLegend = ['Airpuff (n=',num2str(size(DA_airpuff_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('NAc dLight1.3b signal (z-score)'); 
ylim(ylimit); axis off
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% DA, stim
nexttile;
plotSEM(t,DA_randomShutter_all,[.7 .7 .7]); 
plotSEM(t,DA_stim_all,bluePurpleRed(end,:));
taskLegend = ['Stim (n=',num2str(size(DA_stim_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('NAc dLight1.3b signal (z-score)'); 
ylim(ylimit); axis off
plotEvent('',0.5,'r');
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% DA, tone
nexttile;
plotSEM(t,DA_randomShutter_all,[.7 .7 .7]);
plotSEM(t,DA_tone_all,bluePurpleRed(300,:)); 
taskLegend = ['Tone (n=',num2str(size(DA_tone_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('NAc dLight1.3b signal (z-score)'); 
ylim(ylimit); axis off
plotEvent('',0.5,bluePurpleRed(300,:));
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');


% LHb
timeRange = baseline(1).traces.timeRange;
binNumber = size(baseline(1).traces.LHb_water,2);
t = linspace(timeRange(1),timeRange(2),binNumber);
ylimit = [-1,2]; textloc_randomShutter = [timeRange-0.1, -0.05];
% LHb, water
nexttile;
plotSEM(t,LHb_randomShutter_all,[.7 .7 .7]); 
plotSEM(t,LHb_water_all,bluePurpleRed(1,:));
taskLegend = ['Water (n=',num2str(size(LHb_water_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('LHb GCaMP8m signal (z-score)'); 
ylim(ylimit); axis off
plotEvent('',0,bluePurpleRed(1,:));
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% LHb, airpuff
nexttile;
plotSEM(t,LHb_randomShutter_all,[.7 .7 .7]); 
plotSEM(t,LHb_airpuff_all,[0.2, 0.2, 0.2]);
taskLegend = ['Airpuff (n=',num2str(size(LHb_airpuff_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('LHb GCaMP8m signal (z-score)'); 
ylim(ylimit); axis off
plotEvent('',0,[0.2, 0.2, 0.2]);
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% LHb, stim
nexttile;
plotSEM(t,LHb_randomShutter_all,[.7 .7 .7]); 
plotSEM(t,LHb_stim_all,bluePurpleRed(end,:));
taskLegend = ['Stim (n=',num2str(size(LHb_stim_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('LHb GCaMP8m signal (z-score)'); 
ylim(ylimit); axis off
plotEvent('',0.5,'r');
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% LHb, tone
nexttile;
plotSEM(t,LHb_randomShutter_all,[.7 .7 .7]); 
plotSEM(t,LHb_tone_all,bluePurpleRed(300,:));
taskLegend = ['Tone (n=',num2str(size(LHb_tone_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('LHb GCaMP8m signal (z-score)'); 
ylim(ylimit); axis off
plotEvent('',0.5,bluePurpleRed(300,:));
%text(textloc_randomShutter,textloc_randomShutter,'\fontname{Arial} Random shutter','Color',[.7 .7 .7],'HorizontalAlignment','right');
legend(taskLegend,'Location','Northeast'); 
scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');


% Licks
timeRange = baseline(1).lickRate.timeRange;
binNumber = size(baseline(1).lickRate.water,2);
t = linspace(timeRange(1),timeRange(2),binNumber);
ylimit = [0,10]; 
% lick, water
nexttile;
plotSEM(t,waterLickRate_all,bluePurpleRed(1,:));
yline(0,'-','Color',[.7 .7 .7],'LineWidth',2);
plotEvent('',0,bluePurpleRed(1,:));
taskLegend = ['Water (n=',num2str(size(waterLickRate_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('Licks/s'); 
ylim(ylimit); axis off
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% lick, airpuff
nexttile;
plotSEM(t,airpuffLickRate_all,[0.2, 0.2, 0.2]);
yline(0,'-','Color',[.7 .7 .7],'LineWidth',2);
plotEvent('',0,[0.2, 0.2, 0.2]);
taskLegend = ['Airpuff (n=',num2str(size(airpuffLickRate_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('Licks/s'); 
ylim(ylimit); axis off
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% lick, stim
nexttile;
plotSEM(t,stimLickRate_all,bluePurpleRed(end,:));
yline(0,'-','Color',[.7 .7 .7],'LineWidth',2);
taskLegend = ['Stim (n=',num2str(size(stimLickRate_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('Licks/s'); 
ylim(ylimit); axis off
plotEvent('',0.5,'r');
legend(taskLegend,'Location','Northeast'); 
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% lick, tone
nexttile;
plotSEM(t,toneLickRate_all,bluePurpleRed(300,:));
yline(0,'-','Color',[.7 .7 .7],'LineWidth',2);
taskLegend = ['Tone (n=',num2str(size(toneLickRate_all,1)),', ',num2str(mouseNum),' mice)'];
% xlabel('Time (s)'); ylabel('Licks/s'); 
ylim(ylimit); axis off
plotEvent('',0.5,bluePurpleRed(300,:));
legend(taskLegend,'Location','Northeast'); 
scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','Hz');

% Save figure
% set(gcf,'PaperType','A4');
saveFigures(gcf,'baseline_all_table',resultspath);

%% Baseline: 4.0 Num of licks vs DA, LHb traces

%% Baseline: 5. save baseline struct
save(strcat(resultspath,'\','baseline'),'baseline','-append','-v7.3');
disp('baseline struct saved');

%% Pairing: Select files for baseline-reward or punish-reward

b2rList = uipickfiles('FilterSpec','\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
p2rList = uipickfiles('FilterSpec','\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
r2pList = uipickfiles('FilterSpec','\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun');
save(strcat(resultspath,'\','sessions'),'b2rList','p2rList','r2pList','-append');

%% Pairing: Load or create pairing struct
%{
Similar to the cells struct in analyzeSlice
each row is a session
%}

pairing = struct();
% max_misses_allowed = 3; windowlength = 10;

% Baseline -> reward
pairing = addSessionsToStruct(b2rList,pairing,'baseline->reward',newStruct=true);
disp('Baseline->reward finished');
pairing = addSessionsToStruct(r2pList,pairing,'reward->punish');
disp('Reward->punish finished');
pairing = addSessionsToStruct(p2rList,pairing,'punish->reward');
disp('Punish->reward finished');

%% Pairing: Calculate cutoff of trials based on performance for reward
%{
Several criteria:
    1. total trial number have to be more than 100
    2. >=n misses within 10 trials (moving window)
%}

% max_misses_allowed = 3; windowlength = 10;

for i = 1:length(pairing)
    trials = pairing(i).trials;
    disp(['Session calculating: ',pairing(i).session]);
    [trials,cutoff_sample] = getSessionCutoff(trials,pairing(i).type);
    pairing(i).params.analysis.cutoff_sample = cutoff_sample;
    pairing(i).trials = trials;
end

%% Pairing: save struct
save(strcat(resultspath,'\','sessions'),'b2rList','p2rList','r2pList','baselineList','-append');
save(strcat(resultspath,'\','pairing'),'pairing','-v7.3');
save(strcat(resultspath,'\','animals'),'animals','-v7.3');

%% Pairing: sort rows
pairing = sortrows_struct(pairing,{'mouse','session'});

%% Pairing: loop through struct

for i = 1:length(pairing)
    trials = pairing(i).trials;
    pairing(i).stim = trials{trials.isTone == 0 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","performing"]};
    pairing(i).tone = trials{trials.isTone == 1 & trials.isStim == 0,["TrialNumber","CueTime","OutcomeTime","performing"]};
    pairing(i).pair = trials{trials.isTone == 1 & trials.isStim == 1,["TrialNumber","CueTime","OutcomeTime","performing"]};
end

%% Pairing: 1.1 Behavior: create animal struct and get performance params
%{
Performance params
    - reward latency to cue offset
    - first lick latency
    - nAnticipatoryLicks (licks before outcome)

Should return 3 matrix (stim, pair, tone)
    - row is each trial 
    - column is trial number in session + each params
%}

% Initialize value
animals = struct();
mouselist = unique({pairing.mouse});
timeRange = [-1,3];
targetSize = 8209; %size(cur_DA_allTrials,2);
lick_binSize = 0.1;

prev_mouse_row = 0;
for i = 1:length(mouselist)
    cur_mouse = mouselist{i};
    mouse_row = find(strcmp({pairing.mouse},cur_mouse)); % corresponding row of that animal
    
    % Store in animal struct
    animals(i).mouse = cur_mouse;
    animals(i).sessions = mouse_row; % row in pairing struct
    
    %% Initialize stuffs
    % Record task transitions
    b2rEnd_allTrials = 0; r2pEnd_allTrials = 0; p2rEnd_allTrials = 0; 
    b2rEnd_stim = 0; r2pEnd_stim = 0; p2rEnd_stim = 0; 
    b2rEnd_tone = 0; r2pEnd_tone = 0; p2rEnd_tone = 0; 
    b2rEnd_pair = 0; r2pEnd_pair = 0; p2rEnd_pair = 0; 
    
    % Initialize performance matrix
    rl_allTrials = []; fll_allTrials = []; nal_allTrials = [];
    rl_stim = []; fll_stim = []; nal_stim = [];
    rl_tone = []; fll_tone = []; nal_tone = [];
    rl_pair = []; fll_pair = []; nal_pair = [];
    
    % Initialize lick traces
    allTrialsLickRate = []; waterLickRate = []; airpuffLickRate = [];
    stimLickRate = []; toneLickRate = []; pairLickRate = []; randomShutterLickRate = [];
    
    % Initialize DA traces
    DA_allTrials = []; DA_stim = []; 
    DA_tone = []; DA_pair = []; 
    DA_water = []; DA_airpuff = []; DA_randomShutter = [];
    
    % Initialize LHb traces
    LHb_allTrials = []; LHb_stim = []; 
    LHb_tone = []; LHb_pair = []; 
    LHb_water = []; LHb_airpuff = []; LHb_randomShutter = [];
    
    % Initialize DA stats
    % time window: -1-0s, 0-0.5s, 0.5-1s, 1-2s, 2-3s
    
    %% For each session of that mouse
    for j = 1:length(mouse_row)
        %% Current session info
        sessionID = prev_mouse_row + j;
        params = pairing(sessionID).params;
        trials = pairing(sessionID).trials; 
        Fs = pairing(sessionID).params.sync.behaviorFs;
        type = pairing(sessionID).type;
        rightLick = pairing(sessionID).licks;

        % Get analysisRange
        cutoff_sample = pairing(sessionID).params.analysis.cutoff_sample;
        analysisRange_allTrials = find(trials.performing);
        analysisRange_stim = find(pairing(sessionID).stim(:,4));
        analysisRange_tone = find(pairing(sessionID).tone(:,4));
        analysisRange_pair = find(pairing(sessionID).pair(:,4));
        analysisRange_water = find(pairing(sessionID).water < cutoff_sample);
        analysisRange_airpuff = find(pairing(sessionID).airpuff < cutoff_sample);
        
        % Get event index
        allTrialsIdx = trials.CueTime(analysisRange_allTrials);
        waterIdx = pairing(sessionID).water(analysisRange_water);
        airpuffIdx = pairing(sessionID).airpuff(analysisRange_airpuff);
        stimIdx = pairing(sessionID).stim(analysisRange_stim,2);
        toneIdx = pairing(sessionID).tone(analysisRange_tone,2);
        pairIdx = pairing(sessionID).pair(analysisRange_pair,2);
        
        %% Get behavior performance
        
        % Calculate reward latency to cue offset
        cur_rl_allTrials = (trials.OutcomeTime - 5000)./Fs*1000;
        rl_allTrials = [rl_allTrials;cur_rl_allTrials(analysisRange_allTrials)];
        rl_stim = [rl_stim;cur_rl_allTrials(pairing(sessionID).stim(analysisRange_stim,1))];
        rl_tone = [rl_tone;cur_rl_allTrials(pairing(sessionID).tone(analysisRange_tone,1))];
        rl_pair = [rl_pair;cur_rl_allTrials(pairing(sessionID).pair(analysisRange_pair,1))];
        
        % Calculate first lick latency
        cur_fll_allTrials = trials.ReactionTime./Fs*1000;
        fll_allTrials = [fll_allTrials;cur_fll_allTrials(analysisRange_allTrials)];
        fll_stim = [fll_stim;cur_fll_allTrials(pairing(sessionID).stim(analysisRange_stim,1))];
        fll_tone = [fll_tone;cur_fll_allTrials(pairing(sessionID).tone(analysisRange_tone,1))];
        fll_pair = [fll_pair;cur_fll_allTrials(pairing(sessionID).pair(analysisRange_pair,1))];
        
        % Calculate num of anticipatory licks
        cur_nal_allTrials = trials.nAnticipatoryLicks;
        nal_allTrials = [nal_allTrials;cur_nal_allTrials(analysisRange_allTrials)];
        nal_stim = [nal_stim;cur_nal_allTrials(pairing(sessionID).stim(analysisRange_stim,1))];
        nal_tone = [nal_tone;cur_nal_allTrials(pairing(sessionID).tone(analysisRange_tone,1))];
        nal_pair = [nal_pair;cur_nal_allTrials(pairing(sessionID).pair(analysisRange_pair,1))];
        
        % Update trials of each type
        if strcmp(type,'baseline->reward')
            b2rEnd_allTrials = b2rEnd_allTrials + length(analysisRange_allTrials);
            b2rEnd_stim = b2rEnd_stim + length(analysisRange_stim);
            b2rEnd_tone = b2rEnd_tone + length(analysisRange_tone);
            b2rEnd_pair = b2rEnd_pair + length(analysisRange_pair);
        elseif strcmp(type,'reward->punish')
            r2pEnd_allTrials = r2pEnd_allTrials + length(analysisRange_allTrials);
            r2pEnd_stim = r2pEnd_stim + length(analysisRange_stim);
            r2pEnd_tone = r2pEnd_tone + length(analysisRange_tone);
            r2pEnd_pair = r2pEnd_pair + length(analysisRange_pair);
        elseif strcmp(type,'punish->reward')
            p2rEnd_allTrials = p2rEnd_allTrials + length(analysisRange_allTrials);
            p2rEnd_stim = p2rEnd_stim + length(analysisRange_stim);
            p2rEnd_tone = p2rEnd_tone + length(analysisRange_tone);
            p2rEnd_pair = p2rEnd_pair + length(analysisRange_pair);
        else; error('Unknown trial type!');
        end
        
        %% Get event-triggered lick traces
        [cur_allTrialsLickRate,~,~] = getLicks(timeRange,allTrialsIdx,lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
        [cur_waterLickRate,~,~] = getLicks(timeRange,waterIdx,lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
        [cur_airpuffLickRate,~,~] = getLicks(timeRange,airpuffIdx,lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
        [cur_stimLickRate,~,~] = getLicks(timeRange,stimIdx,lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
        [cur_toneLickRate,~,~] = getLicks(timeRange,toneIdx,lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
        [cur_pairLickRate,~,~] = getLicks(timeRange,pairIdx,lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
        [cur_randomShutterLickRate,~,~] = getLicks(timeRange,pairing(sessionID).randomShutter,lick_binSize,[],rightLick,...
                                    params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
                                
        % Combine
        allTrialsLickRate = [allTrialsLickRate; cur_allTrialsLickRate];
        waterLickRate = [waterLickRate; cur_waterLickRate];
        airpuffLickRate = [airpuffLickRate; cur_airpuffLickRate];
        stimLickRate = [stimLickRate; cur_stimLickRate];
        toneLickRate = [toneLickRate; cur_toneLickRate];
        pairLickRate = [pairLickRate; cur_pairLickRate];
        randomShutterLickRate = [randomShutterLickRate; cur_randomShutterLickRate];
                            
        %% Get DA/LHb traces
        [cur_DA_allTrials,DA_timestamp] = plotTraces(allTrialsIdx,timeRange,pairing(sessionID).photometryLJ,bluePurpleRed(1,:),params,plot=false);
        [cur_DA_water,~] = plotTraces(waterIdx,timeRange,pairing(sessionID).photometryLJ,bluePurpleRed(1,:),params,plot=false);
        [cur_DA_airpuff,~] = plotTraces(airpuffIdx,timeRange,pairing(sessionID).photometryLJ,[0.2, 0.2, 0.2],params,plot=false);
        [cur_DA_tone,~] = plotTraces(toneIdx,timeRange,pairing(sessionID).photometryLJ,bluePurpleRed(350,:),params,plot=false);
        [cur_DA_stim,~] = plotTraces(stimIdx,timeRange,pairing(sessionID).photometryLJ,bluePurpleRed(end,:),params,plot=false);
        [cur_DA_pair,~] = plotTraces(pairIdx,timeRange,pairing(sessionID).photometryLJ,bluePurpleRed(end,:),params,plot=false);
        [cur_DA_randomShutter,~] = plotTraces(pairing(sessionID).randomShutter,timeRange,pairing(sessionID).photometryLJ,bluePurpleRed(end,:),params,plot=false);

        [cur_LHb_allTrials,LHb_timestamp] = plotTraces(allTrialsIdx,timeRange,pairing(sessionID).photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni',plot=false);
        [cur_LHb_water,~] = plotTraces(waterIdx,timeRange,pairing(sessionID).photometryNI,bluePurpleRed(1,:),params,photometrySystem='ni',plot=false);
        [cur_LHb_airpuff,~] = plotTraces(airpuffIdx,timeRange,pairing(sessionID).photometryNI,[0.2, 0.2, 0.2],params,photometrySystem='ni',plot=false);
        [cur_LHb_tone,~] = plotTraces(toneIdx,timeRange,pairing(sessionID).photometryNI,bluePurpleRed(350,:),params,photometrySystem='ni',plot=false);
        [cur_LHb_stim,~] = plotTraces(stimIdx,timeRange,pairing(sessionID).photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni',plot=false);
        [cur_LHb_pair,~] = plotTraces(pairIdx,timeRange,pairing(sessionID).photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni',plot=false);
        [cur_LHb_randomShutter,~] = plotTraces(pairing(sessionID).randomShutter,timeRange,pairing(sessionID).photometryNI,bluePurpleRed(end,:),params,photometrySystem='ni',plot=false);
        
        % For sessions before I changed to unmod laser in the middle, 
        % add columns in front to match the unmod number of samples 
        % (SpectualWindowNew = 1, so unmod have one more samples than mod)
        
        % Pad current trace
        pad_num = targetSize - size(cur_DA_allTrials,2);
        if pad_num > 0
            disp(['Number of columns padded to current trace: ', num2str(pad_num)]);
            cur_DA_allTrials = [ones(size(cur_DA_allTrials,1),pad_num).*cur_DA_allTrials(:,1), cur_DA_allTrials];
            cur_DA_water = [ones(size(cur_DA_water,1),pad_num).*cur_DA_water(:,1), cur_DA_water];
            cur_DA_airpuff = [ones(size(cur_DA_airpuff,1),pad_num).*cur_DA_airpuff(:,1), cur_DA_airpuff];
            cur_DA_tone = [ones(size(cur_DA_tone,1),pad_num).*cur_DA_tone(:,1), cur_DA_tone];
            cur_DA_stim = [ones(size(cur_DA_stim,1),pad_num).*cur_DA_stim(:,1), cur_DA_stim];
            cur_DA_pair = [ones(size(cur_DA_pair,1),pad_num).*cur_DA_pair(:,1), cur_DA_pair];
            cur_DA_randomShutter = [ones(size(cur_DA_randomShutter,1),pad_num).*cur_DA_randomShutter(:,1), cur_DA_randomShutter];
        end
        
        % Combine
        DA_allTrials = [DA_allTrials; cur_DA_allTrials];
        DA_water = [DA_water; cur_DA_water];
        DA_airpuff = [DA_airpuff; cur_DA_airpuff];
        DA_tone = [DA_tone; cur_DA_tone];
        DA_stim = [DA_stim; cur_DA_stim];
        DA_pair = [DA_pair; cur_DA_pair];
        DA_randomShutter = [DA_randomShutter; cur_DA_randomShutter];
        
        LHb_allTrials = [LHb_allTrials; cur_LHb_allTrials];
        LHb_water = [LHb_water; cur_LHb_water];
        LHb_airpuff = [LHb_airpuff; cur_LHb_airpuff];
        LHb_tone = [LHb_tone; cur_LHb_tone];
        LHb_stim = [LHb_stim; cur_LHb_stim];
        LHb_pair = [LHb_pair; cur_LHb_pair];
        LHb_randomShutter = [LHb_randomShutter; cur_LHb_randomShutter];
        
    end
    prev_mouse_row = prev_mouse_row + j;
    
    %% Store values to animals
    % Store task trials
    animals(i).taskRange.baseline2reward.allTrials = 1:b2rEnd_allTrials;
    animals(i).taskRange.baseline2reward.stim = 1:b2rEnd_stim;
    animals(i).taskRange.baseline2reward.tone = 1:b2rEnd_tone;
    animals(i).taskRange.baseline2reward.pair = 1:b2rEnd_pair;
    
    animals(i).taskRange.reward2punish.allTrials = (1:r2pEnd_allTrials) + b2rEnd_allTrials;
    animals(i).taskRange.reward2punish.stim = (1:r2pEnd_stim) + b2rEnd_stim;
    animals(i).taskRange.reward2punish.tone = (1:r2pEnd_tone) + b2rEnd_tone;
    animals(i).taskRange.reward2punish.pair = (1:r2pEnd_pair) + b2rEnd_pair;
    
    animals(i).taskRange.punish2reward.allTrials = (1:p2rEnd_allTrials) + (b2rEnd_allTrials + r2pEnd_allTrials);
    animals(i).taskRange.punish2reward.stim = (1:p2rEnd_stim) + (b2rEnd_stim + r2pEnd_stim);
    animals(i).taskRange.punish2reward.tone = (1:p2rEnd_tone) + (b2rEnd_tone + r2pEnd_tone);
    animals(i).taskRange.punish2reward.pair = (1:p2rEnd_pair) + (b2rEnd_pair + r2pEnd_pair);
    
    % Store calculated latency
    animals(i).performance.allTrials.reward_latency = rl_allTrials;
    animals(i).performance.stim.reward_latency = rl_stim;
    animals(i).performance.tone.reward_latency = rl_tone;
    animals(i).performance.pair.reward_latency = rl_pair;
    
    % Store first lick latency
    animals(i).performance.allTrials.first_lick_latency = fll_allTrials;
    animals(i).performance.stim.first_lick_latency = fll_stim;
    animals(i).performance.tone.first_lick_latency = fll_tone;
    animals(i).performance.pair.first_lick_latency = fll_pair;
    
    % Store num of anticipatory licks
    animals(i).performance.allTrials.n_anticipatory_licks = nal_allTrials;
    animals(i).performance.stim.n_anticipatory_licks = nal_stim;
    animals(i).performance.tone.n_anticipatory_licks = nal_tone;
    animals(i).performance.pair.n_anticipatory_licks = nal_pair;
    
    % Store photometry traces
    animals(i).traces.DA_allTrials = DA_allTrials;
    animals(i).traces.DA_water = DA_water;
    animals(i).traces.DA_airpuff = DA_airpuff;
    animals(i).traces.DA_tone = DA_tone;
    animals(i).traces.DA_stim = DA_stim;
    animals(i).traces.DA_pair = DA_pair;
    animals(i).traces.DA_randomShutter = DA_randomShutter;
    animals(i).traces.LHb_allTrials = LHb_allTrials;
    animals(i).traces.LHb_water = LHb_water;
    animals(i).traces.LHb_airpuff = LHb_airpuff;
    animals(i).traces.LHb_tone = LHb_tone;
    animals(i).traces.LHb_stim = LHb_stim;
    animals(i).traces.LHb_pair = LHb_pair;
    animals(i).traces.LHb_randomShutter = LHb_randomShutter;
    animals(i).traces.timeRange = timeRange;
    animals(i).traces.DA_timestamp = linspace(timeRange(1),timeRange(2),size(DA_allTrials,2));
    animals(i).traces.LHb_timestamp = LHb_timestamp;
    
    % Store lick traces
    animals(i).lickRate.timeRange = timeRange;
    animals(i).lickRate.lick_binSize = lick_binSize;
    animals(i).lickRate.timestamp = timeRange(1):lick_binSize:timeRange(2);
    animals(i).lickRate.allTrials = allTrialsLickRate;
    animals(i).lickRate.water = waterLickRate;
    animals(i).lickRate.airpuff = airpuffLickRate;
    animals(i).lickRate.tone = toneLickRate;
    animals(i).lickRate.stim = stimLickRate;
    animals(i).lickRate.pair = pairLickRate;
    animals(i).lickRate.randomShutter = randomShutterLickRate;
    
    %% Analyze traces stat in windows
    timeWindows = [[-1,0];[0.25,0.5];[0.75,1];[1,2];[2,3]];
    animals(i).stats.timeWindows = timeWindows;
    
    % Get stats for DA
    animals(i).stats.DA_allTrials = getTracesStats(DA_allTrials,timeWindows,timestamp=DA_timestamp);
    animals(i).stats.DA_water = getTracesStats(DA_water,timeWindows,timestamp=DA_timestamp);
    animals(i).stats.DA_airpuff = getTracesStats(DA_airpuff,timeWindows,timestamp=DA_timestamp);
    animals(i).stats.DA_tone = getTracesStats(DA_tone,timeWindows,timestamp=DA_timestamp);
    animals(i).stats.DA_stim = getTracesStats(DA_stim,timeWindows,timestamp=DA_timestamp);
    animals(i).stats.DA_pair = getTracesStats(DA_pair,timeWindows,timestamp=DA_timestamp);
    animals(i).stats.DA_randomShutter = getTracesStats(DA_randomShutter,timeWindows,timestamp=DA_timestamp);
    
    % Get stats for LHb
    animals(i).stats.LHb_allTrials = getTracesStats(LHb_allTrials,timeWindows,timestamp=LHb_timestamp);
    animals(i).stats.LHb_water = getTracesStats(LHb_water,timeWindows,timestamp=LHb_timestamp);
    animals(i).stats.LHb_airpuff = getTracesStats(LHb_airpuff,timeWindows,timestamp=LHb_timestamp);
    animals(i).stats.LHb_tone = getTracesStats(LHb_tone,timeWindows,timestamp=LHb_timestamp);
    animals(i).stats.LHb_stim = getTracesStats(LHb_stim,timeWindows,timestamp=LHb_timestamp);
    animals(i).stats.LHb_pair = getTracesStats(LHb_pair,timeWindows,timestamp=LHb_timestamp);
    animals(i).stats.LHb_randomShutter = getTracesStats(LHb_randomShutter,timeWindows,timestamp=LHb_timestamp);
    
    disp(['Animal ',animals(i).mouse,' calculated']);
end

clear nal_* fll_* rl_* trials Fs mouse_row sessionID b2r* r2p* p2r* cur* ...
    DA* LHb* analysisRange* *Idx *LickRate cutoff_sample prev* ...
    rightLick type pad*
disp('Finished: animal struct calculated');

% Save struct
save(strcat(resultspath,'\','pairing'),'pairing','-v7.3');
save(strcat(resultspath,'\','animals'),'animals','-v7.3');
disp('Pairings and aniamls struct stored');

%% Pairing: 1.3.1 Behavior: plot performance params vs trial number for each animal

animalList = {animals([1,3,5:9,11]).mouse};

initializeFig(1,1);
tiledlayout(length(animalList),2);
movmean_window = 5;

% Plot baseline2reward -> reward2punish
for i = 1:length(animalList)
    % Get row of that animal
    row = contains({animals.mouse},animalList{i});
    % Get window
    b2rWindow = animals(row).taskRange.baseline2reward.stim;
    r2pWindow = animals(row).taskRange.reward2punish.stim;
    p2rWindow = animals(row).taskRange.punish2reward.stim;
    
    nexttile
    % Anticipatory lick
    trace = animals(row).performance.stim.n_anticipatory_licks';
%     plot(trialWindow,trace(trialWindow)); hold on
    plotSEM(b2rWindow,movmean(trace(b2rWindow),movmean_window),bluePurpleRed(1,:)); hold on
    plotSEM(r2pWindow,movmean(trace(r2pWindow),movmean_window),bluePurpleRed(end,:)); hold on
    plotSEM(p2rWindow,movmean(trace(p2rWindow),movmean_window),bluePurpleRed(200,:)); hold on
    xlabel('Trial number'); ylabel('Anticipatory licks'); ylim([0 8]);
    title(animals(row).mouse); 
%     axis off; 
%     if i == length(animalList)
%         scalebar('XLen',10,'YLen',1,'XUnit','Trial','YUnit','lick');
%     end
    
%     nexttile
%     % reward_latency
%     trace = animal(row).performance.stim.reward_latency';
% %     plot(trialWindow,trace(trialWindow)); hold on
%     plotSEM(b2rWindow,movmean(trace(b2rWindow),movmean_window),bluePurpleRed(1,:)); hold on
%     plotSEM(r2pWindow,movmean(trace(r2pWindow),movmean_window),bluePurpleRed(end,:)); hold on
%     plotSEM(p2rWindow,movmean(trace(p2rWindow),movmean_window),bluePurpleRed(1,:)); hold on
%     xlabel('Trial number'); ylabel('Reward latency (ms)'); % ylim([50 2000]);
    
    nexttile
    % first lick latency
    trace = animals(row).performance.stim.first_lick_latency';
%     plot(trialWindow,trace(trialWindow)); hold on
    plotSEM(b2rWindow,movmean(trace(b2rWindow),movmean_window),bluePurpleRed(1,:)); hold on
    plotSEM(r2pWindow,movmean(trace(r2pWindow),movmean_window),bluePurpleRed(end,:)); hold on
    plotSEM(p2rWindow,movmean(trace(p2rWindow),movmean_window),bluePurpleRed(200,:)); hold on
    xlabel('Trial number'); ylabel('First lick latency (ms)'); ylim([0 1500]);
    title(animals(row).mouse); 
%     axis off; 
%     if i == length(animalList)
%         scalebar('XLen',10,'YLen',250,'XUnit','Trial','YUnit','ms');
%     end
end
% saveFigures(gcf,'performance_vs_trial',resultspath);

%% Pairing: 1.3.2. Average across all animals
animalList = {animals.mouse};

movmean_window = 5;

% Initalize matrix
b2r_nal = nan(length(animalList),1000); b2r_fll = nan(length(animalList),1000);
r2p_nal = nan(length(animalList),1000); r2p_fll = nan(length(animalList),1000);
p2r_nal = nan(length(animalList),1000); p2r_fll = nan(length(animalList),1000);

% loop through all animals
for i = 1:length(animalList)
    % Get row of that animal
    row = contains({animals.mouse},animalList{i});
    % Get window
    b2rWindow = animals(row).taskRange.baseline2reward.stim;
    r2pWindow = animals(row).taskRange.reward2punish.stim;
    p2rWindow = animals(row).taskRange.punish2reward.stim;
    
    % Anticipatory lick
    trace = animals(row).performance.stim.n_anticipatory_licks';
    b2r_nal(i,1:size(b2rWindow,2)) = trace(b2rWindow); 
    r2p_nal(i,1:size(r2pWindow,2)) = trace(r2pWindow);
    p2r_nal(i,1:size(p2rWindow,2)) = trace(p2rWindow);
    
    % First lick latency
    trace = animals(row).performance.stim.first_lick_latency';
    b2r_fll(i,1:size(b2rWindow,2)) = trace(b2rWindow); 
    r2p_fll(i,1:size(r2pWindow,2)) = trace(r2pWindow);
    p2r_fll(i,1:size(p2rWindow,2)) = trace(p2rWindow);
end

initializeFig(0.5,0.67); tiledlayout(2,3);
% Plot num anticipatory licks average
nexttile
plotWindow = 1:150;
plotSEM(plotWindow,b2r_nal(:,plotWindow),bluePurpleRed(1,:));
xlabel('Trial number'); ylabel('Anticipatory licks'); ylim([0 6])
nexttile
plotWindow = 1:150;
plotSEM(plotWindow,r2p_nal(:,plotWindow),bluePurpleRed(end,:));
xlabel('Trial number'); ylabel('Anticipatory licks'); ylim([0 6])
nexttile
plotWindow = 1:150;
plotSEM(plotWindow,p2r_nal(:,plotWindow),bluePurpleRed(200,:));
xlabel('Trial number'); ylabel('Anticipatory licks'); ylim([0 6])

% Plot first lick latency average
nexttile
plotWindow = 1:150;
plotSEM(plotWindow,b2r_fll(:,plotWindow),bluePurpleRed(1,:));
xlabel('Trial number'); ylabel('First lick latency (ms)'); ylim([0 1500]);
nexttile
plotWindow = 1:150;
plotSEM(plotWindow,r2p_fll(:,plotWindow),bluePurpleRed(end,:));
xlabel('Trial number'); ylabel('First lick latency (ms)'); ylim([0 1500]);
nexttile
plotWindow = 1:150;
plotSEM(plotWindow,p2r_fll(:,plotWindow),bluePurpleRed(200,:));
xlabel('Trial number'); ylabel('First lick latency (ms)'); ylim([0 1500]);

saveFigures(gcf,'performance_vs_trial_avg',resultspath);

%% Pairing: 1.4. group traces by performance, define analysis list

% SL060: normal DA response, LHb response to CS-punish is always dip
% SL062: normal DA response, LHb response to CS-punish is decreasing but
% still a dip
% SL063: weird LHb traces (response profile changes suddenly during 
% CS-punsih) and DA traces(DA signal decreases very fast during reward pairing)
% SL064: B2R session 3 spout moved, DA signal increase during airpuff, LHb
% strongly inhibits during CS-punish; DA dip during p2r session 2
% SL066: normal response for all
% SL068: failed to learn punishment in turns of LHb/DA traces, slow
% learning for reward task

goodDA = {'SL043','SL044','SL046','SL060','SL062','SL063','SL066','SL068'};
goodLHb = {'SL056','SL060','SL062','SL066','SL068'};
goodAll = {'SL043','SL044','SL046','SL056','SL060','SL062','SL066','SL067','SL068'};
all = {'SL043','SL044','SL046','SL056','SL060','SL062','SL063','SL064','SL066','SL067','SL068'};
bothDALHb = {'SL060','SL062','SL063','SL064','SL066','SL067','SL068'};

%% Pairing: 1.5.0. Plot pair photometry signal (important!)

animalList = {'SL066'};
regionList = {'DA','LHb'};
eventName = 'stim';
taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'}; 

earlyWindow = 1:20; lateWindow = 100:120;
midWindow_reward = 50:70; midWindow_punish = 30:50;

initializeFig(1,1); tiledlayout(2,3);
for r = 1:length(regionList)
    for i = 1:length(taskNames)
        % Plot traces
        if strcmp(taskNames{i},'baseline2reward')
            nexttile;
            [early,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=earlyWindow,color=blueWhiteRed(150,:));
            [middle,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=midWindow_reward,color=blueWhiteRed(75,:));
            [late,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=lateWindow,color=blueWhiteRed(1,:));
            plotEvent('Stim',0.5,'r');
            xlabel('Time (s)'); ylabel('z-score'); 
            legend({['Trial ',num2str(earlyWindow(1)),'-',num2str(earlyWindow(end)),' (n=', num2str(size(early,1)),', ',num2str(length(animalList)),' mice)'],...
                ['Trial ',num2str(midWindow_reward(1)),'-',num2str(midWindow_reward(end)),' (n=', num2str(size(middle,1)),', ',num2str(length(animalList)),' mice)'],...
                ['Trial ',num2str(lateWindow(1)),'-',num2str(lateWindow(end)),' (n=', num2str(size(late,1)),', ',num2str(length(animalList)),' mice)']});
        elseif strcmp(taskNames{i},'reward2punish')
            nexttile;
            [early,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=earlyWindow,color=blueWhiteRed(350,:));
            [middle,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=midWindow_punish,color=blueWhiteRed(425,:));
            [late,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=lateWindow,color=blueWhiteRed(end,:));
            plotEvent('Stim',0.5,'r');
            xlabel('Time (s)'); ylabel('z-score'); 
            legend({['Trial ',num2str(earlyWindow(1)),'-',num2str(earlyWindow(end)),' (n=', num2str(size(early,1)),', ',num2str(length(animalList)),' mice)'],...
                ['Trial ',num2str(midWindow_punish(1)),'-',num2str(midWindow_punish(end)),' (n=', num2str(size(middle,1)),', ',num2str(length(animalList)),' mice)'],...
                ['Trial ',num2str(lateWindow(1)),'-',num2str(lateWindow(end)),' (n=', num2str(size(late,1)),', ',num2str(length(animalList)),' mice)']});
        elseif strcmp(taskNames{i},'punish2reward')
            nexttile;
            [early,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=earlyWindow,color=purpleWhiteRed(150,:));
            [middle,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=midWindow_reward,color=purpleWhiteRed(75,:));
            [late,~] = plotTracesFromAnimals(animals,animalList,regionList{r},eventName,task=taskNames{i},traceRange=lateWindow,color=purpleWhiteRed(1,:));
            plotEvent('Stim',0.5,'r');
            xlabel('Time (s)'); ylabel('z-score'); 
            legend({['Trial ',num2str(earlyWindow(1)),'-',num2str(earlyWindow(end)),' (n=', num2str(size(early,1)),', ',num2str(length(animalList)),' mice)'],...
                ['Trial ',num2str(midWindow_reward(1)),'-',num2str(midWindow_reward(end)),' (n=', num2str(size(middle,1)),', ',num2str(length(animalList)),' mice)'],...
                ['Trial ',num2str(lateWindow(1)),'-',num2str(lateWindow(end)),' (n=', num2str(size(late,1)),', ',num2str(length(animalList)),' mice)']});
        end
    end
end
% saveFigures(gcf,'pair_all',resultspath);
clear *Window animalList 

%% Pairing: 1.6. activity changes across trials (table summary)
%{
Plot in tables where row is one stat (e.g. mean/var/max...), column is the
task (baseline->reward, reward->punish ...)
%}
trialRange = 1:100;
animalList = goodAll;
eventName = 'stim';
statsParams_DA = {'avg','avg','avg'};
statsParams_LHb = {'avg','avg','avg'};
smoothWindow = 5;

taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'};
% legendList = {['Baseline'],['CS'],['US'],['Post-US'],['End']};
initializeFig(1,1);
tiledlayout(2,3);

% Initialize empty matrix to store data
cur_DA_traces = nan(length(animalList),length(trialRange));
cur_LHb_traces = nan(length(animalList),length(trialRange));

timeWindows = animals(1).stats.timeWindows;
DA_traces = cell(size(timeWindows,1),1);
LHb_traces = cell(size(timeWindows,1),1);

for i = 1:length(taskNames)
    for w = 1:size(timeWindows,1)
        % Retreive values and store them
        for m = 1:length(animalList)
            row = contains({animals.mouse},animalList{m});
            taskRange = animals(row).taskRange.(taskNames{i}).(eventName);
            if length(trialRange) >= length(taskRange); trials = taskRange;
            else; trials = taskRange(trialRange); end

            individual_trace_DA = animals(row).stats.(['DA_',eventName]).(statsParams_DA{i})(trials,w)';
            individual_trace_LHb = animals(row).stats.(['LHb_',eventName]).(statsParams_LHb{i})(trials,w)';
            cur_DA_traces(m,1:length(individual_trace_DA)) = individual_trace_DA;
            cur_LHb_traces(m,1:length(individual_trace_LHb)) = individual_trace_LHb;
        end
        DA_traces{w} = cur_DA_traces;
        LHb_traces{w} = cur_LHb_traces;
    end

    % Plot
    nexttile(i);
    plotSEM(trialRange,movmean(DA_traces{1},smoothWindow,2),[.7, .7, .7]); % baseline
    plotSEM(trialRange,movmean(DA_traces{2},smoothWindow,2),blueWhiteRed(end,:)); % CS+
    plotSEM(trialRange,movmean(DA_traces{3},smoothWindow,2),blueWhiteRed(1,:)); % US
    %plotSEM(trialRange,movmean(DA_traces{4},smoothWindow,2),blueWhiteRed(350,:)); % post-US
    %plotSEM(trialRange,movmean(DA_traces{5},smoothWindow,2),blueWhiteRed(300,:)); % end
    ylabel([statsParams_DA{i},' DA CS signal']);
    xlabel('Trial'); %ylim([-1, 2]);
    legend({['Baseline (',num2str(size(rmmissing(DA_traces{1},'MinNumMissing',trialRange(end)),1)),' mice)'],...
            ['CS (',num2str(size(rmmissing(DA_traces{2},'MinNumMissing',trialRange(end)),1)),' mice)'],...
            ['US (',num2str(size(rmmissing(DA_traces{3},'MinNumMissing',trialRange(end)),1)),' mice)']});
    
    nexttile(i+3);
    plotSEM(trialRange,movmean(LHb_traces{1},smoothWindow,2),[.7, .7, .7]); % baseline
    plotSEM(trialRange,movmean(LHb_traces{2},smoothWindow,2),blueWhiteRed(end,:)); % CS+
    plotSEM(trialRange,movmean(LHb_traces{3},smoothWindow,2),blueWhiteRed(1,:)); % US
    %plotSEM(trialRange,movmean(LHb_traces{4},smoothWindow,2),blueWhiteRed(350,:)); % post-US
    %plotSEM(trialRange,movmean(LHb_traces{5},smoothWindow,2),blueWhiteRed(300,:)); % end
    ylabel([statsParams_LHb{i},' LHb CS signal']);
    xlabel('Trial'); %ylim([-2.5, 1.5]);
    legend({['Baseline (',num2str(size(rmmissing(LHb_traces{1},'MinNumMissing',trialRange(end)),1)),' mice)'],...
            ['CS (',num2str(size(rmmissing(LHb_traces{2},'MinNumMissing',trialRange(end)),1)),' mice)'],...
            ['US (',num2str(size(rmmissing(LHb_traces{3},'MinNumMissing',trialRange(end)),1)),' mice)']});
end
saveFigures(gcf,'peak_trial_stim_goodavg_100',resultspath);

clear taskNames cur_* DA_traces LHb_traces taskNames statsParams* *Range ...
     trace trials

%% Pairing: 1.6.1 actvity vs trial for single animal scatter plot

animalList = {animals.mouse};
eventName = 'stim';
statsParams_DA = {'avg','avg','avg'};
statsParams_LHb = {'avg','avg','avg'};
maxTrial = 200;
colorList = [bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:);...
    bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:)];

% Loop through all animal
for m = 1:length(animalList)
    % Setup params
    row = contains({animals.mouse},animalList{m});
    taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'};
    initializeFig(1,1); tiledlayout(2,3);

    for i = 1:length(taskNames)
        % Retreive values and store them
        row = contains({animals.mouse},animalList{m});
        taskRange = animals(row).taskRange.(taskNames{i}).(eventName);
        individual_trace_DA = animals(row).stats.(['DA_',eventName]).(statsParams_DA{i})(taskRange,2)';
        individual_trace_LHb = animals(row).stats.(['LHb_',eventName]).(statsParams_LHb{i})(taskRange,2)';
        
        % Fit data
        fit_DA = false; fit_LHb = false;
        if ~any(isnan(individual_trace_DA)) && ~isempty(individual_trace_DA)
            x_DA = (1:size(individual_trace_DA,2))';
            f_DA = fit(x_DA,individual_trace_DA','exp1');
            fit_DA = true;
        end
        if ~any(isnan(individual_trace_LHb)) && ~isempty(individual_trace_LHb)
            x_LHb = (1:size(individual_trace_LHb,2))';
            f_LHb = fit(x_LHb,individual_trace_LHb','exp1');
            fit_LHb = true;
        end

        % Plot
        nexttile(i);
        scatter(1:size(individual_trace_DA,2),individual_trace_DA,'filled',...
            'MarkerFaceColor',colorList(i,:),'MarkerFaceAlpha',0.5); hold on
        if fit_DA
            h = plot(f_DA,x_DA,individual_trace_DA'); delete(h(1));
            h(2).ColorMode = 'manual'; h(2).Color = colorList(i,:);
            h(2).LineWidth = 2;
        end
        if ~isempty(individual_trace_DA); xlim([0,size(individual_trace_DA,2)]); end
        ylabel([statsParams_DA{i},' DA CS signal']);
        xlabel('Trial');

        nexttile(i+3);
        scatter(1:size(individual_trace_LHb,2),individual_trace_LHb,'filled',...
            'MarkerFaceColor',colorList(i,:),'MarkerFaceAlpha',0.5); hold on
        if fit_LHb
            h = plot(f_LHb,x_LHb,individual_trace_LHb'); delete(h(1));
            h(2).ColorMode = 'manual'; h(2).Color = colorList(i,:);
            h(2).LineWidth = 2;
        end
        if ~isempty(individual_trace_LHb); xlim([0,size(individual_trace_LHb,2)]); end
        ylabel([statsParams_LHb{i},' LHb CS signal']);
        xlabel('Trial');
    end
    
    % Save figure
    saveFigures(gcf,['peak_trial_scatter_',num2str(animals(row).mouse)],[resultspath,'\Individual_peak_trial_scatter_mean']);
    close;
end

clear individual* x_* w r i j k

%% Pairing: 1.6.2. Scatter plot of avg photometry response vs trial for all

animalList = goodAll;
eventName = 'stim';
statsParams_DA = {'avg','avg','avg'};
statsParams_LHb = {'avg','avg','avg'};
polyfitdegree = 1;
maxTrial = 90;
colorList = [bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:);...
    bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:)];

initializeFig(0.67,0.67); tiledlayout(2,3);
% Loop through all animal
for i = 1:length(taskNames)
    % Plot DA
    traces_DA = nan(length(animalList),maxTrial);
    nexttile(i);
    for m = 1:length(animalList)
        % Setup params
        row = contains({animals.mouse},animalList{m});
        taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'};
        
        % Retreive values and store them
        taskRange = animals(row).taskRange.(taskNames{i}).(eventName);
        individual_trace_DA = animals(row).stats.(['DA_',eventName]).(statsParams_DA{i})(taskRange,2)';
        selectRange = 1:min(maxTrial,size(individual_trace_DA,2));
        traces_DA(m,selectRange) = individual_trace_DA(selectRange);
        
        scatter(1:size(individual_trace_DA,2),individual_trace_DA,'filled',...
            'MarkerFaceColor',colorList(i,:),'MarkerFaceAlpha',0.3); hold on
    end
    if isempty(maxTrial); xlim([0,size(individual_trace_DA,2)]); 
    else; xlim([1,maxTrial]); end
    ylabel([statsParams_DA{i},' DA CS signal']);
    xlabel('Trial');
    
    % Analyze DA correlation
    mean_DA = mean(traces_DA,1,'omitnan');
    %plotSEM(selectRange, traces_DA, colorList(i,:));
    x_trials = ones(size(traces_DA)).*selectRange; % selectRange
    [r_DA,p_DA] = corr(x_trials',traces_DA','Type','Spearman');
    % Fit DA
    x = selectRange;
    [p,S] = polyfit(x,mean_DA,polyfitdegree);
    [yfit,delta] = polyconf(p,x,S,'alpha',0.05);
    plotSEM(x,yfit,colorList(i,:),delta=delta);
    title(['Fit: ',texlabel(polystr(p))]);

    % Plot LHb
    traces_LHb = nan(length(animalList),maxTrial);
    nexttile(i+3);
    for m = 1:length(animalList)
        % Setup params
        row = contains({animals.mouse},animalList{m});
        taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'};
        
        % Retreive values and store them
        taskRange = animals(row).taskRange.(taskNames{i}).(eventName);
        individual_trace_LHb = animals(row).stats.(['LHb_',eventName]).(statsParams_LHb{i})(taskRange,2)';
        selectRange = 1:min(maxTrial,size(individual_trace_LHb,2));
        traces_LHb(m,selectRange) = individual_trace_LHb(selectRange);
        
        scatter(1:size(individual_trace_LHb,2),individual_trace_LHb,'filled',...
            'MarkerFaceColor',colorList(i,:),'MarkerFaceAlpha',0.3); hold on
    end    
    if isempty(maxTrial); xlim([0,size(individual_trace_LHb,2)]);
    else; xlim([1,maxTrial]); end
    ylabel([statsParams_LHb{i},' LHb CS signal']);
    xlabel('Trial');
    
    % Analyze LHb correlation
    mean_LHb = mean(traces_LHb,1,'omitnan');
    %plotSEM(selectRange, traces_LHb, colorList(i,:));
    x_trials = ones(size(traces_LHb)).*selectRange; % selectRange
    [r_LHb,p_LHb] = corr(x_trials',traces_LHb','Type','Spearman');
    x = selectRange;
    [p,S] = polyfit(x,mean_LHb,polyfitdegree);
    [yfit,delta] = polyconf(p,x,S,'alpha',0.05);
    plotSEM(x,yfit,colorList(i,:),delta=delta);
    title(['Fit: ',texlabel(polystr(p))]);
end

% saveFigures(gcf,'avg_response_vs_trials',resultspath);
clear individual* x* w r i j k mean_* delta p S yfit

%% Pairing: 1.7.0 Get photometry stats vs performance params
animalList = {animals.mouse};
timeWindow = 2; % CS
eventName = 'stim';
performance = 'n_anticipatory_licks';
statsParams = {'max','min','avg'};

% Initialize params
DA_params = cell(length(animalList),length(statsParams));
LHb_params = cell(length(animalList),length(statsParams));

% loop through all selected animals
for i = 1:length(animalList)
    % Setup
    row = contains({animals.mouse},animalList{i});
    %initializeFig(1,1); tiledlayout(2,3);
    
    % Extract data
    perf = animals(row).performance.(eventName).(performance);
    for j = 1:length(statsParams)
        stats_DA = animals(row).stats.(['DA_',eventName]).(statsParams{j})(:,timeWindow);
        stats_LHb = animals(row).stats.(['LHb_',eventName]).(statsParams{j})(:,timeWindow);
        DA_params{i,j} = [stats_DA, perf];
        LHb_params{i,j} = [stats_LHb, perf];
    end
end
clear perf stats_* performance statsParams regionList i j

%% Pairing: 1.7.1. Plot scatter plot for each animal
animalList = {animals.mouse};
timeWindow = 2; % CS
eventName = 'stim';
performance = 'n_anticipatory_licks';
statsParams = {'max','min','avg'};
smoothWindow = 8;
colorList = [bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:);...
    bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:)];
taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'};

for i = 1:length(animalList)
    % Setup
    row = contains({animals.mouse},animalList{i});
    initializeFig(1,0.6); tiledlayout(3,8);
    
    % loop through params
    for j = 1:length(statsParams)
        % All DA
        nexttile
        scatter(DA_params{row,j}(:,1),DA_params{row,j}(:,2),'filled',...
            'MarkerFaceColor',[.3,.3,.3],'MarkerFaceAlpha',0.5);
        xlabel([statsParams{j},' DA signal']); ylabel('Anticipatory licks');
        % Analyze DA correlation
        [r_DA,p_DA] = corr(DA_params{row,j}(:,1),DA_params{row,j}(:,2),'Type','Spearman');
        title(['r=',num2str(round(r_DA,3)),', p=',num2str(round(p_DA,3))]);
        
        % All LHb
        nexttile
        scatter(LHb_params{row,j}(:,1),LHb_params{row,j}(:,2),'filled',...
            'MarkerFaceColor',[.3,.3,.3],'MarkerFaceAlpha',0.5);
        xlabel([statsParams{j},' LHb signal']); ylabel('Anticipatory licks');
        % Analyze LHb correlation
        [r_LHb,p_LHb] = corr(LHb_params{row,j}(:,1),LHb_params{row,j}(:,2),'Type','Spearman');
        title(['r=',num2str(round(r_LHb,3)),', p=',num2str(round(p_LHb,3))]);
        
        % task
        for k = 1:length(taskNames)
            taskRange = animals(row).taskRange.(taskNames{k}).(eventName);
            % DA
            nexttile
            scatter(DA_params{row,j}(taskRange,1),DA_params{row,j}(taskRange,2),'filled',...
                'MarkerFaceColor',colorList(k,:),'MarkerFaceAlpha',0.5);
            xlabel([statsParams{j},' DA signal']); ylabel('Anticipatory licks');
            if ~isempty(taskRange)
                % Analyze DA correlation
                [r_DA,p_DA] = corr(DA_params{row,j}(taskRange,1),DA_params{row,j}(taskRange,2),'Type','Spearman');
                title(['r=',num2str(round(r_DA,3)),', p=',num2str(round(p_DA,3))]);
            end
            % LHb
            nexttile
            scatter(LHb_params{row,j}(taskRange,1),LHb_params{row,j}(taskRange,2),'filled',...
                'MarkerFaceColor',colorList(k,:),'MarkerFaceAlpha',0.5);
            xlabel([statsParams{j},' LHb signal']); ylabel('Anticipatory licks');
            if ~isempty(taskRange)
                % Analyze LHb correlation
                [r_LHb,p_LHb] = corr(LHb_params{row,j}(taskRange,1),LHb_params{row,j}(taskRange,2),'Type','Spearman');
                title(['r=',num2str(round(r_LHb,3)),', p=',num2str(round(p_LHb,3))]);
            end
        end
    end
    
    % Save figure
    saveFigures(gcf,['photometry_nal_scatter_',num2str(animals(row).mouse)],[resultspath,'\Individual_photometry_nal_scatter'],pngOnly=true);
    close;
end
clear r_* p_* taskRange statsParams i j k idx row

%% Pairing: 1.7.2. Plot scatter plot for each animal in one plot
animalList = goodAll;
timeWindow = 2; % CS
eventName = 'stim';
performance = 'n_anticipatory_licks';
statsParams = 'avg';
smoothWindow = 8;
colorList = [bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:);...
    bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:)];
taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'};

initializeFig(1,1); tiledlayout(length(animalList),8);
for i = 1:length(animalList)
    % Setup
    row = contains({animals.mouse},animalList{i});
    
    j = 3;
    % All DA
    nexttile
    scatter(DA_params{row,j}(:,1),DA_params{row,j}(:,2),'filled',...
        'MarkerFaceColor',[.3,.3,.3],'MarkerFaceAlpha',0.5);
    xlabel([statsParams,' DA signal']); ylabel('Anticipatory licks');
    % Analyze DA correlation
    [r_DA,p_DA] = corr(DA_params{row,j}(:,1),DA_params{row,j}(:,2),'Type','Spearman');
    title(['r=',num2str(round(r_DA,3)),', p=',num2str(round(p_DA,3))]);

    % All LHb
    nexttile
    scatter(LHb_params{row,j}(:,1),LHb_params{row,j}(:,2),'filled',...
        'MarkerFaceColor',[.3,.3,.3],'MarkerFaceAlpha',0.5);
    xlabel([statsParams,' LHb signal']); ylabel('Anticipatory licks');
    % Analyze LHb correlation
    [r_LHb,p_LHb] = corr(LHb_params{row,j}(:,1),LHb_params{row,j}(:,2),'Type','Spearman');
    title(['r=',num2str(round(r_LHb,3)),', p=',num2str(round(p_LHb,3))]);

    % task
    for k = 1:length(taskNames)
        taskRange = animals(row).taskRange.(taskNames{k}).(eventName);
        % DA
        nexttile
        scatter(DA_params{row,j}(taskRange,1),DA_params{row,j}(taskRange,2),'filled',...
            'MarkerFaceColor',colorList(k,:),'MarkerFaceAlpha',0.5);
        xlabel([statsParams,' DA signal']); ylabel('Anticipatory licks');
        if ~isempty(taskRange)
            % Analyze DA correlation
            [r_DA,p_DA] = corr(DA_params{row,j}(taskRange,1),DA_params{row,j}(taskRange,2),'Type','Spearman');
            title(['r=',num2str(round(r_DA,3)),', p=',num2str(round(p_DA,3))]);
        end
        % LHb
        nexttile
        scatter(LHb_params{row,j}(taskRange,1),LHb_params{row,j}(taskRange,2),'filled',...
            'MarkerFaceColor',colorList(k,:),'MarkerFaceAlpha',0.5);
        xlabel([statsParams,' LHb signal']); ylabel('Anticipatory licks');
        if ~isempty(taskRange)
            % Analyze LHb correlation
            [r_LHb,p_LHb] = corr(LHb_params{row,j}(taskRange,1),LHb_params{row,j}(taskRange,2),'Type','Spearman');
            title(['r=',num2str(round(r_LHb,3)),', p=',num2str(round(p_LHb,3))]);
        end
    end
end

Save figure
saveFigures(gcf,'meanPhotometry_nal_scatter_all',resultspath,pngOnly=false);
close;

clear i j k r m r_* p_*

%% Pairing: 1.7.3. Plot scatter for all animals
animalList = goodAll;
timeWindow = 2; % CS
eventName = 'stim';
performance = 'n_anticipatory_licks';
statsParams = 'avg';
corrType = 'Spearman';
colorList = [bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:);...
    bluePurpleRed(1,:); bluePurpleRed(end,:); bluePurpleRed(200,:)];
taskNames = {'baseline2reward', 'reward2punish', 'punish2reward'};

% Get summary matrix
avgDA_nal_all = []; avgLHb_nal_all = [];
avgDA_nal_b2r = []; avgLHb_nal_b2r = [];
avgDA_nal_r2p = []; avgLHb_nal_r2p = [];
avgDA_nal_p2r = []; avgLHb_nal_p2r = [];

for i = 1:length(animalList)
    % Setup
    row = contains({animals.mouse},animalList{i});
    j = 3;
    
    % All DA/LHb
    avgDA_nal_all = [avgDA_nal_all; DA_params{row,j}(:,1:2)];
    avgLHb_nal_all = [avgLHb_nal_all; LHb_params{row,j}(:,1:2)];
    
    % B2R
    taskRange = animals(row).taskRange.baseline2reward.(eventName);
    avgDA_nal_b2r = [avgDA_nal_b2r; DA_params{row,j}(taskRange,1:2)];
    avgLHb_nal_b2r = [avgLHb_nal_b2r; LHb_params{row,j}(taskRange,1:2)];
    
    % R2P
    taskRange = animals(row).taskRange.reward2punish.(eventName);
    avgDA_nal_r2p = [avgDA_nal_r2p; DA_params{row,j}(taskRange,1:2)];
    avgLHb_nal_r2p = [avgLHb_nal_r2p; LHb_params{row,j}(taskRange,1:2)];
    
    % P2R
    taskRange = animals(row).taskRange.punish2reward.(eventName);
    avgDA_nal_p2r = [avgDA_nal_p2r; DA_params{row,j}(taskRange,1:2)];
    avgLHb_nal_p2r = [avgLHb_nal_p2r; LHb_params{row,j}(taskRange,1:2)];
end

% Plot
initializeFig(.67,.5); tiledlayout(2,4);
% Plot scatter
ax1 = nexttile; data = rmmissing(avgDA_nal_all);
dscatter(data(:,1),data(:,2)); colormap(ax1,getColormap([.3,.3,.3],[255,255,255],250));
xlabel([statsParams,' DA signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax2 = nexttile; data = rmmissing(avgDA_nal_b2r);
dscatter(data(:,1),data(:,2)); colormap(ax2, blueWhiteRed(1:250,:));
xlabel([statsParams,' DA signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax3 = nexttile; data = rmmissing(avgDA_nal_r2p);
dscatter(data(:,1),data(:,2)); colormap(ax3,r2p_cmap(1:250,:));
xlabel([statsParams,' DA signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax4 = nexttile; data = rmmissing(avgDA_nal_p2r);
dscatter(data(:,1),data(:,2)); colormap(ax4,p2r_cmap(1:250,:));
xlabel([statsParams,' DA signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax5 = nexttile; data = rmmissing(avgLHb_nal_all);
dscatter(data(:,1),data(:,2)); colormap(ax5,getColormap([.3,.3,.3],[255,255,255],250));
xlabel([statsParams,' LHb signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax6 = nexttile; data = rmmissing(avgLHb_nal_b2r);
dscatter(data(:,1),data(:,2)); colormap(ax6,blueWhiteRed(1:250,:));
xlabel([statsParams,' LHb signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax7 = nexttile; data = rmmissing(avgLHb_nal_r2p);
dscatter(data(:,1),data(:,2)); colormap(ax7,r2p_cmap(1:250,:));
xlabel([statsParams,' LHb signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax8 = nexttile; data = rmmissing(avgLHb_nal_p2r);
dscatter(data(:,1),data(:,2)); colormap(ax8,p2r_cmap(1:250,:));
xlabel([statsParams,' LHb signal']); ylabel('Anticipatory licks');
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

%Save figure
saveFigures(gcf,'meanPhotometry_nal_scatter_oneplot',resultspath,pngOnly=false);
close;

clear i j k r m r_* p_* taskRange ax* f avg*

%% Pairing: 1.7.4 scatter plot avg LHb vs DA for all

animalList = bothDALHb;
eventName = 'stim';
timeWindow = 2;
statsParams = 'avg';


DALHb_all = []; DALHb_b2r = []; DALHb_r2p = []; DALHb_p2r = [];
for i = 1:length(animalList)
    % Setup
    row = contains({animals.mouse},animalList{i});
    j = 3;
    
    % Extract
    stats_DA = animals(row).stats.(['DA_',eventName]).(statsParams)(:,timeWindow);
    stats_LHb = animals(row).stats.(['LHb_',eventName]).(statsParams)(:,timeWindow);
    cur = [stats_DA, stats_LHb];
    DALHb_all = [DALHb_all; cur];

    % Trial type
    taskRange = animals(row).taskRange.baseline2reward.(eventName);
    cur = [stats_DA(taskRange), stats_LHb(taskRange)];
    DALHb_b2r = [DALHb_b2r; cur];
    taskRange = animals(row).taskRange.reward2punish.(eventName);
    cur = [stats_DA(taskRange), stats_LHb(taskRange)];
    DALHb_r2p = [DALHb_r2p; cur];
    taskRange = animals(row).taskRange.punish2reward.(eventName);
    cur = [stats_DA(taskRange), stats_LHb(taskRange)];
    DALHb_p2r = [DALHb_p2r; cur];
end

% Plot 
initializeFig(0.8,0.3); tiledlayout(1,4);

ax1 = nexttile; data = DALHb_all;
dscatter(data(:,1),data(:,2)); colormap(ax1,getColormap([.3,.3,.3],[255,255,255],250));
xlabel([statsParams,' DA signal']); ylabel([statsParams,' LHb signal']);
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax2 = nexttile; data = DALHb_b2r;
dscatter(data(:,1),data(:,2)); colormap(ax2, blueWhiteRed(1:250,:));
xlabel([statsParams,' DA signal']); ylabel([statsParams,' LHb signal']);
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax3 = nexttile; data = DALHb_r2p;
dscatter(data(:,1),data(:,2)); colormap(ax3,r2p_cmap(1:250,:));
xlabel([statsParams,' DA signal']); ylabel([statsParams,' LHb signal']);
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

ax4 = nexttile; data = DALHb_p2r;
dscatter(data(:,1),data(:,2)); colormap(ax4,p2r_cmap(1:250,:));
xlabel([statsParams,' DA signal']); ylabel([statsParams,' LHb signal']);
[r,p] = corr(data(:,1),data(:,2),'Type',corrType);
title(['r=',num2str(round(r,3)),', p=',num2str(round(p,3))]);

%Save figure
saveFigures(gcf,'meanDA_LHb_scatter_all',resultspath,pngOnly=false);
close;

clear i j k r m r_* p_* taskRange ax* f avg* data DALHb*

%% Pairing: loop through animals struct

for i = 1:length(animals)

end

%% Photometry: 1.1 DA and LHb anticorrelation analysis

DA_LHb_sessions = {};
params = pairing(64).params;
DA = pairing(64).photometryLJ(params.sync.commonStartPhotometry:end); 
LHb = pairing(64).photometryNI(params.sync.commonStartNI:end);

% Option 1: downsample photometry
DA_resampled = resample(DA,params.sync.ni_photometryFs,params.sync.photometryFs);
DA_matched = DA_resampled(1:length(LHb));

% Cross corr
[xx, xl]=xcorr(normalize(DA_matched), normalize(LHb), 200, 'normalized');
plot(xl, xx);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Not used %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pairing: 1.0 Behavior: find lick traces for baseline

timeRange = [-1,3]; lick_binSize = 0.1;

for i = 1:length(pairing)  
    trials = pairing(i).trials; 
    params = pairing(i).params; 
    rightLick = pairing(i).licks;
    

    % getLicks by trial type
    [waterLickRate,~,~] = getLicks(timeRange,pairing(i).water,lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
    [airpuffLickRate,~,~] = getLicks(timeRange,pairing(i).airpuff,lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
    [stimLickRate,~,~] = getLicks(timeRange,pairing(i).stim(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
    [toneLickRate,~,~] = getLicks(timeRange,pairing(i).tone(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
    [pairLickRate,~,~] = getLicks(timeRange,pairing(i).pair(:,2),lick_binSize,[],rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,inputLickIdx=true);
                            
    % Save results for this session
    pairing(i).lickRate.timeRange = timeRange;
    %baseline(i).lickRaster.timeRange = timeRange;
    pairing(i).lickRate.lick_binSize = lick_binSize;
    %baseline(i).lickRaster.lick_binSize = lick_binSize;
    pairing(i).lickRate.water = waterLickRate;
    %baseline(i).lickRaster.water = waterLicks;
    pairing(i).lickRate.airpuff = airpuffLickRate;
    %baseline(i).lickRaster.airpuff = airpuffLicks;
    pairing(i).lickRate.stim = stimLickRate;
    %baseline(i).lickRaster.stim = stimLicks;
    pairing(i).lickRate.tone = toneLickRate;
    %baseline(i).lickRaster.tone = toneLicks;
    pairing(i).lickRate.pair = pairLickRate;
    %baseline(i).lickRaster.tone = toneLicks;

    disp(['Session ',pairing(i).session,' calculated']);
    clear stimLickRate toneLickRate waterLickRate airpuffLickRate pairLickRate
end
disp('Lick rate traces genearted');

%% Pairing: 1.1 Behavior: concat lick rate vs time
waterLickRate_pairing = [];
airpuffLickRate_pairing = [];
stimLickRate_pairing = [];
toneLickRate_pairing = [];
pairLickRate_pairing = [];

for i = 1:length(pairing)
    cutoff_sample = pairing(i).params.analysis.cutoff_sample;
    % Water
    analysisRange = find(pairing(i).water < cutoff_sample);
    waterLickRate_pairing = [waterLickRate_pairing; pairing(i).lickRate.water(analysisRange,:)];
    % Airpuff
    analysisRange = find(pairing(i).airpuff < cutoff_sample);
    airpuffLickRate_pairing = [airpuffLickRate_pairing; pairing(i).lickRate.airpuff(analysisRange,:)];
    % Stim
    analysisRange = find(pairing(i).stim(:,4)==1);
    if ~isempty(pairing(i).lickRate.stim)
        stimLickRate_pairing = [stimLickRate_pairing; pairing(i).lickRate.stim(analysisRange,:)];
    end
    % Tone
    analysisRange = find(pairing(i).tone(:,4)==1);
    if ~isempty(pairing(i).lickRate.tone)
        toneLickRate_pairing = [toneLickRate_pairing; pairing(i).lickRate.tone(analysisRange,:)];
    end
    % Pair
    analysisRange = find(pairing(i).pair(:,4)==1);
    if ~isempty(pairing(i).lickRate.pair)
        pairLickRate_pairing = [pairLickRate_pairing; pairing(i).lickRate.pair(analysisRange,:)];
    end
end

% Remove nan rows
waterLickRate_pairing = rmmissing(waterLickRate_pairing);
airpuffLickRate_pairing = rmmissing(airpuffLickRate_pairing);
stimLickRate_pairing = rmmissing(stimLickRate_pairing);
toneLickRate_pairing = rmmissing(toneLickRate_pairing);
pairLickRate_pairing = rmmissing(pairLickRate_pairing);

disp('Lick rate traces concatenated');

%% Pairing: 1.1 Behavior: example lick raster for B-R and R-P 

%% Pairing: 1.2.0 Behavior: (not used) lick rate vs time for water, airpuff, stim, tone, pair 
% (just to check my code is working)

initializeFig(.5,.5);
tiledlayout(1,2);
mouseNum = numel(unique({pairing.mouse}));

nexttile;
timeRange = pairing(1).lickRate.timeRange;
binNumber = size(pairing(1).lickRate.water,2);
t = linspace(timeRange(1),timeRange(2),binNumber);
plotSEM(t,waterLickRate_pairing,bluePurpleRed(1,:));
plotSEM(t,airpuffLickRate_pairing,[0.2, 0.2, 0.2]);
taskLegend = {['Water (n=',num2str(length(waterLickRate_pairing)),', ',num2str(mouseNum),'mice)'],...
            ['Airpuff (n=',num2str(length(airpuffLickRate_pairing)),', ',num2str(mouseNum),'mice)']};
xlabel('Time (s)'); ylabel('Licks/s'); 
legend(taskLegend,'Location','Northeast');

nexttile;
plotSEM(t,toneLickRate_pairing,bluePurpleRed(300,:));
plotSEM(t,pairLickRate_pairing,bluePurpleRed(150,:));
plotSEM(t,stimLickRate_pairing,bluePurpleRed(end,:));
plotEvent('',0.5,'r');
taskLegend = {['Tone (n=',num2str(length(toneLickRate_pairing)),', ',num2str(mouseNum),'mice)'],...
              ['Pair (n=',num2str(length(pairLickRate_pairing)),', ',num2str(mouseNum),'mice)'],...
              ['Stim (n=',num2str(length(stimLickRate_pairing)),', ',num2str(mouseNum),'mice)']};
xlabel('Time (s)'); ylabel('Licks/s'); 
legend(taskLegend,'Location','Northeast');
% scalebar('XLen',0.5,'YLen',0.5,'XUnit','sec','YUnit','z');
% saveFigures(gcf,'pairing_lickRate',resultspath);

%% 2.4 Photometry: reward pair LHb

%% 2.5 Photometry: punish pair DA

%% 2.6 Photometry: punish pair LHb

%% 2.7 Photometry: performance params vs DA during stim

%% 2.8 Photometry: performacne params vs LHb during stim

%% 2.9 Photometry: DA stim AUC vs trials during B-R

%% 2.10 Photometry: DA stim AUC vs trials during R-P

%% 2.11 Photometry: DA stim AUC vs trials during P-R

%% 2.12 Photometry: LHb stim AUC vs trials during B-R

%% 2.13 Photometry: LHb stim AUC vs trials during R-P

%% 2.14 Photometry: LHb stim AUC vs trials during P-R

%% 2.15 Photometry: bar plot DA AUC at baseline, reward end, punish end

%% 2.16 Photometry: bar plot LHb AUC at baseline, reward end, punish end

