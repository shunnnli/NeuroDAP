%% (Required) Load sync data 
clear; close all;
% addpath(genpath('/Users/shunli/Downloads/Sabatini lab/Methods'));
addpath(genpath('D:\Shun\Neuropixel analysis\Methods'));

% sessionName = '20220612-SJ518-R_g0';
sessionName = '20220613-SJ518-R_g0';
% sessionName = '20220615-SJ518-R_g0';
load(strcat('sync_',sessionName,'.mat'));
[twoColors,colors,blueRedYellow,blueGreenYellow,blueWhiteRed] = loadColors;
% twoColors = {'#4B5CCC';'#C355BC'};

%% Extract event time (in samples)

blockLength = 50;
allTrialsTime = find(allTrials); % in sample
allTrialsLR = zeros(size(allTrialsTime)); % trial type (left/right) of each trial
allTrialsRP = zeros(size(allTrialsTime)); % trial type (reward/punish) of each trial
leftTrialsTime = find(allTrials == 1); % Left tone trials time in samples
rightTrialsTime = find(allTrials == 2); % Right tone trials time in sample

% Extract outcome delivery time (in samples)
rewardDelivery = find(leftSolenoid); % Reward delivery
punishDelivery = find(airpuff); % Punishment delivery

% Extract reward/punishment cue/trial time (in samples)
rewardTrials = []; punishTrials = [];
rewardOmissionTrials = []; punishOmissionTrials = [];
block = 1; 
for i = 1:length(allTrialsTime)
    disp(['Analyzing trial ',num2str(i),'/',num2str(length(allTrialsTime))]);
    cur_cue = allTrialsTime(i);
    allTrialsLR(i) = allTrials(cur_cue);
    
    if (allTrials(cur_cue) == 1 && mod(block,2) ~= 0)
        % Odd block: left -> reward
        allTrialsRP(i) = 1;
        rewardTrials = [rewardTrials; cur_cue];
    elseif (allTrials(cur_cue) == 2 && mod(block,2) ~= 0)
        % Odd block: right -> punishment
        allTrialsRP(i) = 2;
        punishTrials = [punishTrials; cur_cue];
    elseif (allTrials(cur_cue) == 1 && mod(block,2) == 0)
        % Even block: left -> punishment
        allTrialsRP(i) = 2;
        punishTrials = [punishTrials; cur_cue];
    elseif (allTrials(cur_cue) == 2 && mod(block,2) == 0)
        % Even block: right -> reward
        allTrialsRP(i) = 1;
        rewardTrials = [rewardTrials; cur_cue];
        
    elseif (allTrials(cur_cue) == 3 && mod(block,2) ~= 0)
        % Odd block: left -> reward
        allTrialsRP(i) = 3;
        rewardOmissionTrials = [rewardOmissionTrials; cur_cue];
    elseif (allTrials(cur_cue) == 4 && mod(block,2) ~= 0)
        % Odd block: right -> punishment
        allTrialsRP(i) = 4;
        punishOmissionTrials = [punishOmissionTrials; cur_cue];
    elseif (allTrials(cur_cue) == 3 && mod(block,2) == 0)
        % Even block: left -> punishment
        allTrialsRP(i) = 3;
        punishOmissionTrials = [punishOmissionTrials; cur_cue];
    elseif (allTrials(cur_cue) == 4 && mod(block,2) == 0)
        % Even block: right -> reward
        allTrialsRP(i) = 3;
        rewardOmissionTrials = [rewardOmissionTrials; cur_cue];
        
    else
        disp(num2str(allTrialsTime(i)));
    end

    % Update block if needed
    if mod(i,blockLength) == 0
        block = block + 1;
    end
end

leftTrials = find(allTrialsLR == 1); % Left tone trials
rightTrials = find(allTrialsLR == 2); % Right tone trials

% Sanity check
if length(punishTrials) ~= length(punishDelivery)
    disp("Error: punish trials != punish delivery; check blockLength");
end

%% Extract lick subtypes
% Generate a list of timestamps where corresponding events happened

thresholdSamp = 0.2 * nidq.Fs; % maiximum ILI to be in a bout
reactionTimeSamp = 5 * nidq.Fs;
delayTimeSamp = 1 * nidq.Fs;

allTrialsTime = find(allTrials);
leftLickTime = find(leftLick);

% Initialize arrays for individual lick events
anticipatory_lick = []; % Any lick between cue and delivery
reward_lick = []; % First bout of lick after reward delivery
punish_lick = []; % First bout of lick after punishment delivery
spontaneous_lick = []; % Spontaneous licks after first lick bout
reward_omitted_lick = [];
punish_omitted_lick = [];
miss_reward_trial = []; % Reward trials where animal do not have any licks

% Initialize arrays for lick diagnostic stats
nAnticipatoryLicks = zeros(1,length(allTrialsTime)); % #anticipatory licks for each trial
nTotalLicks = zeros(1,length(allTrialsTime)); % #total licks for each trial
hit = 0; % hit rate for reward trials
fa = 0; % false alarm rate for punish trials

total_analyzed = 0;
for i = 1:length(allTrialsTime)
    % Sanity check
    %total_analyzed = total_analyzed + length(correct_lick)+length(incorrect_lick);
    %disp('----------')
    %disp(i-1);
    %disp(num2str(total_analyzed));
    %disp(num2str(length(rewarded_lick)+length(unrewarded_lick)+...
    %    length(choice_lick)+length(spontaneous_lick)));
    %disp('----------');
    
    cur_cue = allTrialsTime(i);
    
    % Licks before first cue as spontaneous licks
    if i == 1
        spon = leftLickTime((1<leftLickTime)&(leftLickTime<cur_cue));
        spontaneous_lick = [spontaneous_lick,spon];
    end
    
    % Reward cue
    if ismember(cur_cue,rewardTrials)
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick & delivery time for the current trial
        currentLicks = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        currentDelivery = rewardDelivery((cur_cue<rewardDelivery)&(rewardDelivery<next_cue));
        nTotalLicks(i) = length(currentLicks);
        
        % Assign licks to subtypes
        if ~isempty(currentLicks)
            if isempty(currentDelivery)
                % No reward delivered, all licks are spontaneous
                spontaneous_lick = [spontaneous_lick,currentLicks];
            elseif length(currentDelivery) == 1
                hit = hit + 1;
                % Anticipatory lick
                anticipatory = currentLicks(currentLicks<currentDelivery);
                anticipatory_lick = [anticipatory_lick,anticipatory];
                nAnticipatoryLicks(i) = length(anticipatory);
                % Reward lick
                [~,col] = find(currentLicks>currentDelivery,1);
                if col<=1
                    last_anticipatory_lick = 0;
                else
                    last_anticipatory_lick = currentLicks(col-1);
                end
                [bout,spontaneous] = findLickBout(thresholdSamp,last_anticipatory_lick,currentLicks,[]);
                reward_lick = [reward_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                disp(['More than 1 reward delivery found: : trial ',num2str(i),...
                    '; #delivery', num2str(length(currentDelivery))]);
                pause;
            end
        else
            miss_reward_trial = [miss_reward_trial, cur_cue];
        end

    % Punishment cue
    elseif ismember(cur_cue,punishTrials)
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick & delivery time for the current trial
        currentLicks = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        currentDelivery = punishDelivery((cur_cue<punishDelivery)&(punishDelivery<next_cue));
        nTotalLicks(i) = length(currentLicks);
        
        % Assign licks to subtypes
        if ~isempty(currentLicks)
            if isempty(currentDelivery)
                % Display error since there should always be punishment
                disp(['No punishment delivery found: trial ',num2str(i)]);
                pause;
            elseif length(currentDelivery) == 1
                fa = fa + 1;
                % Anticipatory lick
                anticipatory = currentLicks(currentLicks<currentDelivery);
                anticipatory_lick = [anticipatory_lick,anticipatory];
                nAnticipatoryLicks(i) = length(anticipatory);
                % Reward lick
                [~,col] = find(currentLicks>currentDelivery,1);
                if col<=1
                    last_anticipatory_lick = 0;
                else
                    last_anticipatory_lick = currentLicks(col-1);
                end
                [bout,spontaneous] = findLickBout(thresholdSamp,last_anticipatory_lick,currentLicks,[]);
                punish_lick = [punish_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            else
                disp(['More than 1 punish delivery found: : trial ',num2str(i),...
                    '; #delivery', num2str(length(currentDelivery))]);
                pause;
            end
        else
            % miss_reward_trial = [miss_reward_trial, cur_cue];
        end
    
    % Reward omission
    elseif ismember(cur_cue,rewardOmissionTrials)
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick & delivery time for the current trial
        currentLicks = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        currentDelivery = rewardDelivery((cur_cue<rewardDelivery)&(rewardDelivery<next_cue));
        nTotalLicks(i) = length(currentLicks);
        
        % Assign licks to subtypes
        if ~isempty(currentLicks)
            if isempty(currentDelivery)
                % Should be the case since its reward omission
                % Anticipatory lick
                fakeDelivery = cur_cue + delayTimeSamp;
                anticipatory = currentLicks(currentLicks<fakeDelivery);
                anticipatory_lick = [anticipatory_lick,anticipatory];
                nAnticipatoryLicks(i) = length(anticipatory);
                % Reward omission lick
                [~,col] = find(currentLicks>fakeDelivery,1);
                if col<=1
                    last_anticipatory_lick = 0;
                else
                    last_anticipatory_lick = currentLicks(col-1);
                end
                [bout,spontaneous] = findLickBout(thresholdSamp,last_anticipatory_lick,currentLicks,[]);
                reward_omitted_lick = [reward_omitted_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            elseif length(currentDelivery) == 1
                disp('Reward delivery found for omission trials');
                pause;
            else
                disp(['More than 1 reward delivery found: ',...
                    num2str(length(currentDelivery))]);
                pause;
            end
        else
            % miss_reward_trial = [miss_reward_trial, cur_cue];
        end
    
    % Punishment omission trial
    elseif ismember(cur_cue,punishOmissionTrials)
        % Find next cue time
        if i == length(allTrialsTime)
            next_cue = length(allTrials); 
        else
            next_cue = allTrialsTime(i+1);
        end
        
        % Find lick & delivery time for the current trial
        currentLicks = leftLickTime((cur_cue<leftLickTime)&(leftLickTime<next_cue));
        currentDelivery = rewardDelivery((cur_cue<rewardDelivery)&(rewardDelivery<next_cue));
        nTotalLicks(i) = length(currentLicks);
        
        % Assign licks to subtypes
        if ~isempty(currentLicks)
            if isempty(currentDelivery)
                % Should be the case since its punish omission
                % Anticipatory lick
                fakeDelivery = cur_cue + delayTimeSamp;
                anticipatory = currentLicks(currentLicks<fakeDelivery);
                anticipatory_lick = [anticipatory_lick,anticipatory];
                nAnticipatoryLicks(i) = length(anticipatory);
                % Punish omission lick
                [~,col] = find(currentLicks>fakeDelivery,1);
                if col<=1
                    last_anticipatory_lick = 0;
                else
                    last_anticipatory_lick = currentLicks(col-1);
                end
                [bout,spontaneous] = findLickBout(thresholdSamp,last_anticipatory_lick,currentLicks,[]);
                punish_omitted_lick = [punish_omitted_lick,bout];
                spontaneous_lick = [spontaneous_lick,spontaneous];
            elseif length(currentDelivery) == 1
                disp('Punishment delivery found for omission trials');
                pause;
            else
                disp(['More than 1 punishment delivery found: ',...
                    num2str(length(currentDelivery))]);
                pause;
            end
        else
            % miss_reward_trial = [miss_reward_trial, cur_cue];
        end
        
    else
        continue
    end
    
end

% Sanity check
disp(num2str(length(leftLickTime)));
disp(num2str(length(reward_lick)+length(punish_lick)+...
    length(spontaneous_lick)+length(anticipatory_lick)));


%% Plot lick trend across block

% Plot left-trial licks
eventIdx = leftTrials; lastTrialsInBlock = [];
subplot(2,1,1);
for i = blockLength:blockLength:length(allTrialsTime)
  trialsInBlock = find(eventIdx>=i-blockLength & eventIdx<i);
  plotTrials = [lastTrialsInBlock,trialsInBlock];
  if mod(floor(i/blockLength),2) ~= 0
    plot(plotTrials,nAnticipatoryLicks(eventIdx(plotTrials)),...
      '-','LineWidth',1.5,'Color',twoColors(1)); hold on
    plot(plotTrials,nTotalLicks(eventIdx(plotTrials)),...
      '--','LineWidth',1.5,'Color',twoColors(1)); hold on
  else
    plot(plotTrials,nAnticipatoryLicks(eventIdx(plotTrials)),...
      '-','LineWidth',1.5,'Color',twoColors(2)); hold on
    plot(plotTrials,nTotalLicks(eventIdx(plotTrials)),...
      '--','LineWidth',1.5,'Color',twoColors(2)); hold on
  end
  lastTrialsInBlock = trialsInBlock(end);
end
legend('Anticipatory licks (reward)','Total licks (reward)',...
    'Anticipatory licks (airpuff)','Total licks (airpuff)');
ylabel('Licks'); xlim([1 length(eventIdx)-1]);

% Plot right-trial licks
eventIdx = rightTrials; lastTrialsInBlock = [];
subplot(2,1,2);
for i = blockLength:blockLength:length(allTrialsTime)
  trialsInBlock = find( eventIdx>=i-blockLength & eventIdx<i);
  plotTrials = [lastTrialsInBlock,trialsInBlock];
  if mod(floor(i/blockLength),2) == 0
    plot(plotTrials,nAnticipatoryLicks(eventIdx(plotTrials)),...
      '-','LineWidth',1.5,'Color',twoColors(1)); hold on
    plot(plotTrials,nTotalLicks(eventIdx(plotTrials)),...
      '--','LineWidth',1.5,'Color',twoColors(1)); hold on
  else
    plot(plotTrials,nAnticipatoryLicks(eventIdx(plotTrials)),...
      '-','LineWidth',1.5,'Color',twoColors(2)); hold on
    plot(plotTrials,nTotalLicks(eventIdx(plotTrials)),...
      '--','LineWidth',1.5,'Color',twoColors(2)); hold on
  end
  lastTrialsInBlock = trialsInBlock(end);
end
legend('Anticipatory licks (airpuff)','Total licks (airpuff)',...
    'Anticipatory licks (reward)','Total licks (reward)');
ylabel('Licks'); xlim([1 length(eventIdx)-1]);

%% Plot lick raster plot for reward/punish trial
timeRange = [-1,3];

% Plot left-cue-triggered licks
eventIdx = leftTrialsTime;
figure; drawLicks('Left cue',timeRange,eventIdx,...
                        leftLick,[],nidq,timeNI);
xline(1,'-','Outcome delivery','Color','b','LineWidth',2,'HandleVisibility','off');

% Plot right-cue-triggered licks
eventIdx = rightTrialsTime;
figure; drawLicks('Right cue',timeRange,eventIdx,...
                        leftLick,[],nidq,timeNI);
xline(1,'-','Outcome delivery','Color','b','LineWidth',2,'HandleVisibility','off');

% Plot reward trials licks
eventIdx = rewardTrials;
figure; drawLicks('Reward cue',timeRange,eventIdx,...
                        leftLick,[],nidq,timeNI);
xline(1,'-','Reward delivery','Color','b','LineWidth',2,'HandleVisibility','off');

% Plot punish trials licks
eventIdx = punishTrials;
figure; drawLicks('Punish cue',timeRange,eventIdx,...
                        leftLick,[],nidq,timeNI);
xline(1,'-','Airpuff delivery','Color','b','LineWidth',2,'HandleVisibility','off');

autoArrangeFigures();

%% Compare licks of different contexts

% #Anticipatory lick for reward vs punish trials
subplot(1,2,1);
anticipatory_reward = nAnticipatoryLicks(allTrialsRP==1);
anticipatory_punish = nAnticipatoryLicks(allTrialsRP==2);

i = 1; x = i*ones(1,length(anticipatory_reward));
swarmchart(x,anticipatory_reward,[],hex2rgb(twoColors(1)),'filled',...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5,'XJitterWidth',0.7); hold on
errorbar(i,mean(anticipatory_reward),getCI(anticipatory_reward,0.95),'k','LineWidth',1.5); hold on
scatter(i,mean(anticipatory_reward),100,'red','LineWidth',2); hold on

i = 2; x = i*ones(1,length(anticipatory_punish));
% bar(i,mean(anticipatory_punish)); hold on
swarmchart(x,anticipatory_punish,[],hex2rgb(twoColors(2)),'filled',...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5,'XJitterWidth',0.7); hold on
errorbar(i,mean(anticipatory_punish),getCI(anticipatory_punish,0.95),'k','LineWidth',1.5); hold on
scatter(i,mean(anticipatory_punish),100,'red','LineWidth',2); hold on

xticks([1:2]); ylabel('Anticipatory licks'); box off
set(gca,'xticklabel',{'Water','Airpuff'});

% #Total licks for reward vs punish trials
subplot(1,2,2);
totalLicks_reward = nTotalLicks(allTrialsRP==1);
totalLicks_punish = nTotalLicks(allTrialsRP==2);

i = 1; x = i*ones(1,length(totalLicks_reward));
% bar(i,mean(totalLicks_reward)); hold on
swarmchart(x,totalLicks_reward,[],hex2rgb(twoColors(1)),'filled',...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5,'XJitterWidth',0.7); hold on
errorbar(i,mean(totalLicks_reward),getCI(totalLicks_reward,0.95),'k','LineWidth',1.5); hold on
scatter(i,mean(totalLicks_reward),100,'red','LineWidth',2); hold on

i = 2; x = i*ones(1,length(totalLicks_punish));
% bar(i,mean(totalLicks_punish)); hold on
swarmchart(x,totalLicks_punish,[],hex2rgb(twoColors(2)),'filled',...
    'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5,'XJitterWidth',0.7); hold on
errorbar(i,mean(totalLicks_punish),getCI(totalLicks_punish,0.95),'k','LineWidth',1.5); hold on
scatter(i,mean(totalLicks_punish),100,'red','LineWidth',2); hold on

xticks([1:2]); ylabel('Total licks'); box off
set(gca,'xticklabel',{'Water','Airpuff'});
