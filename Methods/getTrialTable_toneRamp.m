function trials = getTrialTable_toneRamp(task,events,rightSolenoid_rounded,airpuff_rounded,options)

arguments
    task string
    events cell
    rightSolenoid_rounded double
    airpuff_rounded double

    options.behaviorFs double = 10000 % in Hz
    options.reactionTime double = 2 % in seconds
    options.minLicks double = 2 
    options.pavlovian logical = true

    options.toneNames = {'ramp_short','ramp_long','jump_short','jump_long','sine'}
end

% Load events as separte vectors
allTrials = events{1};
airpuffON = events{2};
rightSolenoidON = events{3};
rightLickON = events{4};
toneON = events{5};
stimON = events{6};
stim2ON = events{7};
toneTypes = events{8};

% Initialize trial table
% Selection time: time of last choice lick (before outcome)
% Reaction time: time of first lick
% BaselineLicks (in session time), others are in trial time
varTypes = {'double','double','string',...
            'logical','string','logical','logical','logical','logical',...
            'double','double','double','double','double','double','double','double',...
            'double','double','double','double','double','double',...
            'cell','cell','cell','cell'};
varNames = {'TrialNumber','Choice','Outcome',...
            'isTone','toneType','isStim','isReward','isPunishment','isStim2',...
            'CueTime','NextCue','ENL','ITI','ReactionTime','OutcomeTime','OutcomeReactionTime','LastLickTime',...
            'RewardSize','PunishSize','nLicks','nAnticipatoryLicks','nOutcomeLicks','nBaselineLicks',...
            'TrialLicks','AnticipatoryLicks','OutcomeLicks','BaselineLicks'};

% Initialize trial subtypes
trials = table('Size',[length(allTrials) length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

if ~options.pavlovian
    go = []; nogo = [];
    hit = []; miss = []; fa = []; cr = [];
end

% Define session related params
gracePeriod = floor(0.2 * options.behaviorFs);
reactionTimeSamp = options.reactionTime * options.behaviorFs;

% Loop over all trials
for i=1:length(allTrials)
    cur_cue = allTrials(i);
    if i == length(allTrials); next_cue = NaN;
    else; next_cue = allTrials(i+1); end
    ITI = (next_cue - cur_cue) / options.behaviorFs;

    % Trial related
    outcome = nan; choice = nan; outcomeReactionTime = nan;

    % Tone types
    toneTypeIdx = cellfun(@(v) ismember(cur_cue, v), toneTypes);
    toneType  = options.toneNames{find(toneTypeIdx,1)};

    % Trial cue/stim
    toneTime = toneON(toneON >= cur_cue-gracePeriod & toneON < next_cue-gracePeriod) - cur_cue;
    isTone = ~isempty(toneTime);
    stimTime = stimON(stimON >= cur_cue-gracePeriod & stimON < next_cue-gracePeriod) - cur_cue;
    isStim = ~isempty(stimTime);
    stim2Time = stim2ON(stim2ON >= cur_cue-gracePeriod & stim2ON < next_cue-gracePeriod) - cur_cue;
    isStim2 = ~isempty(stim2Time);

    % Trial outcome
    rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue-gracePeriod & rightSolenoidON<next_cue-gracePeriod) - cur_cue;
    isReward = ~isempty(rewardTime);
    if isReward; rewardSize = sum(rightSolenoid_rounded(rewardTime + cur_cue));
    else; rewardSize = 0; rewardTime = reactionTimeSamp; end
    punishTime = airpuffON(airpuffON>=cur_cue-gracePeriod & airpuffON<next_cue-gracePeriod) - cur_cue;
    isPunishment = ~isempty(punishTime);
    if isPunishment; punishSize = sum(airpuff_rounded(punishTime + cur_cue));
    else; punishSize = 0; punishTime = reactionTimeSamp; end
    outcomeTime = min([reactionTimeSamp, rewardTime, punishTime]);

    % Trial licks
    trial_r_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
    trialLicks = sortrows([trial_r_licks, 2*ones(length(trial_r_licks),1)]);
    % Initialize trial table value
    reactionTime = 0; lastLickTime = 0; 
    if ~isempty(trialLicks)
        reactionTime = min([trialLicks(1,1),outcomeTime]); 
        lastLickTime = trialLicks(end,1); 
    end

    % Separate licks (i.e. anticipatory licks: licks before outcome)
    if isempty(trialLicks)
        choiceLicks = [];
        outcomeLicks = []; 
        baselineLicks = [];
    else
        choiceLicks = trialLicks(trialLicks(:,1)<= outcomeTime,:); 
        consecLicks = getConsecutive(choiceLicks(:,2)); % Get consecutive lick number
        
        % Separate OutcomeLicks vs BaselineLicks
        % Opt 1: OutcomeLicks are licks within 5 seconds of outcome
        outcomeWindow = 5*options.behaviorFs;
        outcomeLicks = trialLicks(trialLicks(:,1)>=outcomeTime+outcomeWindow,:);
        baselineLicks = setdiff(trialLicks,[choiceLicks;outcomeLicks],'rows');
        baselineLicks(:,1) = baselineLicks(:,1) + cur_cue;

        % Opt 2: outcomelicks are first lickbout after outcome
        % remainingLickBout = getLickBout(remainingLicks(:,1));
        % outcomeLicks = remainingLicks(remainingLicks(:,1)>=remainingLickBout(1,1) & remainingLicks(:,1) < remainingLicksBout(2,1)); 
        % baselineLicks = getLickBout(remainingLicks,consummatoryLicks);
    end

    % For storing in trial table
    nLicks = size(trialLicks,1);
    nAnticipatoryLicks = size(choiceLicks,1);
    

    % Calculate ENL
    lastLickPrevTrial = rightLickON(rightLickON < cur_cue); 
    if isempty(lastLickPrevTrial)
        ENL = 0;
    else
        if i == 1; last_cue = 0; else; last_cue = allTrials(i-1); end
        ENL = cur_cue - max(lastLickPrevTrial(end),last_cue); 
    end

    % Trial choice (1: go; 0: nogo)
    if ~options.pavlovian
        if isempty(choiceLicks); nogo = [nogo;cur_cue]; choice = 0;
        else
            if consecLicks(end) >= options.minLicks
                go = [go;cur_cue]; 
                choice = 1;
            else; nogo = [nogo;cur_cue]; choice = 0;
            end
        end

        % Update outcomeTime for omission trials (no outcome but licked)
        if ~(isReward || isPunishment) && (choice ~= 0)
            choiceLicks_idx = find(consecLicks>=options.minLicks,options.minLicks);
            choiceLicks = choiceLicks(1:choiceLicks_idx(end),:);
            outcomeTime = min(outcomeTime,choiceLicks(end,1));
        end
    end


    % Trial outcome & outcome reaction time
    if contains(task,'reward')
        if ~options.pavlovian
            if choice == 1; hit = [hit; cur_cue]; outcome = 'H';
            elseif choice == 0; miss = [miss; cur_cue]; outcome = 'M'; end
        end
    elseif contains(task,'punish')
        if ~options.pavlovian
            if choice == 0; hit = [hit; cur_cue]; outcome = 'H';
            elseif choice == 1; miss = [miss; cur_cue]; outcome = 'M'; end
        end
    else
        if ~options.pavlovian; outcome = 'H'; end
    end
    outcomeReactionTime = rightLickON(find(rightLickON>outcomeTime+cur_cue,1)) - (cur_cue+outcomeTime);
    if isempty(outcomeReactionTime); outcomeReactionTime = nan; end

    % Trial table
    trials(i,:) = {i,choice,outcome,...
        isTone,toneType,isStim,isReward,isPunishment,isStim2,...
        cur_cue,next_cue,ENL,ITI,reactionTime,outcomeTime,outcomeReactionTime,lastLickTime,...
        rewardSize,punishSize,nLicks,nAnticipatoryLicks,size(outcomeLicks,1),size(baselineLicks,1),...
        num2cell(trialLicks,[1 2]),num2cell(choiceLicks+cur_cue,[1 2]),num2cell(outcomeLicks+cur_cue,[1 2]),num2cell(baselineLicks+cur_cue,[1 2])};
end

% Add a row for outside events
trials(i+1,:) = {0,0,0,0,0,0,0,0,0,...
                 0,0,0,0,0,0,0,0,...
                 0,0,0,0,0,0,...
                 {},{},{},{}};

% Sanity check
% disp(['     Total hit = ',num2str(length(hit))]);
% disp(['     Total miss = ',num2str(length(miss))]);
% disp(['     Total FA = ',num2str(length(fa))]);
% disp(['     Total CR = ',num2str(length(cr))]);

end