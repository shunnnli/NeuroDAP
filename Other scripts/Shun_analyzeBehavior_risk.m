% Risk taking behavior analysis

clear; close all;
% addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));
addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Neuropixel analysis\Methods'));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed,blueGreenPurple] = loadColors;
             

sessionName = "20221014-SL030-D10_g2";
session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Project risk\Recordings\';
% session.path = '\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Test photometry\';
load(strcat(session.path,sessionName,'\','sync_',sessionName,'.mat'));
disp(['Session ',sessionName,' loaded']);

timeRange = [-0.5,2]; minLicks = 2; reactionTimeSamp = 1 * nidq.Fs;

%% Generate trial tables

% Find tones
allTrialsTime = find(allTrials);
leftSolenoidON = find(leftSolenoid);
leftLickON = find(leftLick);
rightSolenoidON = find(rightSolenoid);
rightLickON = find(rightLick);
airpuffON = find(airpuff);

% Initialize trial table
% Selection time: time of last choice lick (before outcome)
% Reaction time: time of first lick
varTypes = {'double','string','string',...
            'double','double','double','double','double'};
varNames = {'TrialNumber','Choice','Outcome',...
            'CueTime','ReactionTime','OutcomeTime','LastLickTime','NextCue'};

% Initialize trial subtypes
trials = table('Size',[length(allTrialsTime) length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);
goLeft = []; goRight = []; noChoice = [];
gotReward = []; gotPunish = []; gotNothing = []; 


% Loop over all trials
for i=1:length(allTrialsTime)
    cur_cue = allTrialsTime(i);
    if i == length(allTrialsTime); next_cue = length(allTrials);
    else; next_cue = allTrialsTime(i+1); end

    % Initialize trial table value
    reactionTime = 0; lastLickTime = 0;

    % Trial licks
    trial_l_licks = (leftLickON(leftLickON > cur_cue & leftLickON < next_cue) - cur_cue)';
    trial_l_licks = [trial_l_licks, ones(length(trial_l_licks),1)];
    trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
    trial_r_licks = [trial_r_licks, 2*ones(length(trial_r_licks),1)];
    % first col is lick time in sample, second col is side of lick (1->left, 2->right)
    trialLicks = sortrows([trial_l_licks;trial_r_licks]);
    if ~isempty(trialLicks)
        reactionTime = trialLicks(1,1); lastLickTime = trialLicks(end,1); 
    end

    % Trial outcome time
    leftReward = leftSolenoidON(leftSolenoidON > cur_cue & leftSolenoidON < next_cue) - cur_cue;
    rightReward = rightSolenoidON(rightSolenoidON > cur_cue & rightSolenoidON < next_cue) - cur_cue;
    punishment = airpuffON(airpuffON > cur_cue & airpuffON < next_cue)- cur_cue;
    isLeftReward = ~isempty(leftReward); 
    isRightReward = ~isempty(rightReward);
    isPunishment = ~isempty(punishment);
    outcomeTime = min([reactionTimeSamp, leftReward, rightReward, punishment]);

    % Trial choice
    if isempty(trialLicks); choiceLicks = [];
    else
        choiceLicks = trialLicks(trialLicks(:,1)<=outcomeTime,:); 
        % Get consecutive lick number
        d = [true; diff(choiceLicks(:,2))~=0; true];  % TRUE if values change
        n = diff(find(d));
        consecLicks = repelem(n, n);
    end

    if isempty(choiceLicks); noChoice = [noChoice;cur_cue]; choice = 'N';
    else
        if consecLicks(end) >= minLicks && choiceLicks(end,2) == 1
            goLeft = [goLeft;cur_cue]; choice = 'L';
        elseif consecLicks(end) >= minLicks && choiceLicks(end,2) == 2
            goRight = [goRight;cur_cue]; choice = 'R';
        else; noChoice = [noChoice;cur_cue]; choice = 'N';
        end
    end

    % Update outcome time for nothing trials (no outcome but have choice)
    if ~(isLeftReward || isRightReward || isPunishment) && ~strcmp(choice,'N')
        choiceLicks_idx = find(consecLicks>=minLicks,minLicks);
        choiceLicks = choiceLicks(1:choiceLicks_idx(end),:);
        outcomeTime = min(outcomeTime,choiceLicks(end,1));
    end

    % Trial outcome
    if ~(isLeftReward || isRightReward || isPunishment)
        if strcmp(choice,'N'); outcome = 'NC';
        else; outcome = 'NO'; gotNothing = [gotNothing;cur_cue]; end
    elseif isPunishment; outcome = 'P'; gotPunish = [gotPunish;cur_cue];
    elseif isLeftReward || isRightReward
        gotReward = [gotReward;cur_cue];
        if choice == 'N' % free reward
            if isLeftReward && ~isRightReward; outcome = 'LF';
            elseif ~isLeftReward && isRightReward; outcome = 'RF';
            elseif isLeftReward && isRightReward; outcome = 'BF';
            end
        else
            if isLeftReward && ~isRightReward; outcome = 'LR';
            elseif ~isLeftReward && isRightReward; outcome = 'RR';
            elseif isLeftReward && isRightReward; outcome = 'BR';
            end
        end
    end

    % Trial table
    trials(i,:) = {i,choice,outcome,cur_cue,...
                reactionTime,outcomeTime,lastLickTime,next_cue};
end

% Sanity check
disp(['Total reward = ',num2str(length(gotReward))]);
disp(['Total punishment = ',num2str(length(gotPunish))]);
disp(['Total left = ',num2str(length(goLeft))]);
disp(['Total right = ',num2str(length(goRight))]);
disp(['No choice = ',num2str(length(noChoice))]);
% disp(['Free reward = ',num2str(length(free_reward))]);

%% Trial history

lick_summary_fig = figure('Position', get(0,'Screensize'));
eventIdx = find(allTrials == 1); 
% drawLicks(timeRange,eventIdx,leftLick,rightLick,nidq,timeNI,leftSolenoid,rightSolenoid,airpuff);
drawTrials(timeRange,eventIdx,trials,leftLick,rightLick,nidq,timeNI,leftSolenoid,rightSolenoid,airpuff);
ylimit = ylim;
patch([0 .1 .1 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
        'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
xline(0,'-','Tone','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off'); 
box off

saveas(lick_summary_fig,strcat(session.path,sessionName,'\summary_lick_',sessionName,'.png'));
return

%% Reaction time analysis

safe = 'L'; risk = 'R';

% Reaction time after reward vs punish in prev trials
reward_reaction = []; punish_reaction = []; nothing_reaction = [];
safe_reward_reaction = []; risk_reward_reaction = [];
safe_punish_reaction = []; risk_punish_reaction = [];
safe_nothing_reaction = []; risk_nothing_reaction = [];

for i = 2:size(trials,1)
    prev_outcome = trials{i-1,"Outcome"};
    prev_choice = trials{i-1,"Choice"};
    if contains(prev_outcome,'R') && ~contains(prev_outcome,'F')
        reward_reaction = [reward_reaction; trials{i,"ReactionTime"}];
        if strcmp(prev_choice,safe)
            safe_reward_reaction = [safe_reward_reaction; trials{i,"ReactionTime"}];
        elseif strcmp(prev_choice,risk)
            risk_reward_reaction = [risk_reward_reaction; trials{i,"ReactionTime"}];
        end
    elseif contains(prev_outcome,'P')
        punish_reaction = [punish_reaction; trials{i,"ReactionTime"}];
        if strcmp(prev_choice,safe)
            safe_punish_reaction = [safe_punish_reaction; trials{i,"ReactionTime"}];
        elseif strcmp(prev_choice,risk)
            risk_punish_reaction = [risk_punish_reaction; trials{i,"ReactionTime"}];
        end
    elseif strcmp(prev_outcome,'NO')
        nothing_reaction = [nothing_reaction; trials{i,"ReactionTime"}];
        if strcmp(prev_choice,safe)
            safe_nothing_reaction = [safe_nothing_reaction; trials{i,"ReactionTime"}];
        elseif strcmp(prev_choice,risk)
            risk_nothing_reaction = [risk_nothing_reaction; trials{i,"ReactionTime"}];
        end
    end
end

bar(1,mean(reward_reaction)); hold on
errorbar(1,mean(reward_reaction),getCI(reward_reaction,0.95),'k','LineWidth',1.5); hold on
%swarmchart(1,reward_reaction,20,'filled'); hold on

bar(2,mean(punish_reaction)); hold on
errorbar(2,mean(punish_reaction),getCI(punish_reaction,0.95),'k','LineWidth',1.5); hold on
%swarmchart(2,punish_reaction,20,'filled'); hold on

bar(3,mean(nothing_reaction)); hold on
errorbar(3,mean(nothing_reaction),getCI(nothing_reaction,0.95),'k','LineWidth',1.5); hold on
%swarmchart(3,nothing_reaction,20,'filled'); hold on

bar(4,mean(safe_reward_reaction)); hold on
errorbar(4,mean(safe_reward_reaction),getCI(safe_reward_reaction,0.95),'k','LineWidth',1.5); hold on
bar(4.5,mean(risk_reward_reaction)); hold on
errorbar(4.5,mean(risk_reward_reaction),getCI(risk_reward_reaction,0.95),'k','LineWidth',1.5); hold on

bar(6,mean(safe_punish_reaction)); hold on
errorbar(6,mean(safe_punish_reaction),getCI(safe_punish_reaction,0.95),'k','LineWidth',1.5); hold on
bar(6.5,mean(risk_punish_reaction)); hold on
errorbar(6.5,mean(risk_punish_reaction),getCI(risk_punish_reaction,0.95),'k','LineWidth',1.5); hold on

bar(8,mean(safe_nothing_reaction)); hold on
errorbar(8,mean(safe_nothing_reaction),getCI(safe_nothing_reaction,0.95),'k','LineWidth',1.5); hold on
bar(8.5,mean(risk_nothing_reaction)); hold on
errorbar(8.5,mean(risk_nothing_reaction),getCI(risk_nothing_reaction,0.95),'k','LineWidth',1.5); hold on


% xticks([1:15]);
% set(gca,'xticklabel',{'Left cue','Right cue',...
%     'Left/right cue','Left/right solenoid',...
%     'Left/right licks','Choice/rewarded/unrewarded licks',...
%     'Cues + solenoid','Cues + licks','Cues + lick subtypes',...
%     'Cues + solenoid + licks','Cues + solenoid + lick subtypes',...
%     'Cues + solenoid + licks + lick subtypes',...
%     'Cues + solenoid + licks + velocity',...
%     'Cues + solenoid + lick subtypes + velocity',...
%     'Cues + solenoid + licks + lick subtypes + velocity'});