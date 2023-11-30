function trials = getTrialTable_shijia(events,options)

arguments
    events cell

    options.behaviorFs double = 2000
    options.reactionTime double = 1.5
end

% Load events as separte vectors
allTrials = events{1};
rightSolenoidON = events{2};
rightLickON = events{3};

% Initialize trial table
% Selection time: time of last choice lick (before outcome)
% Reaction time: time of first lick
% BaselineLicks (in session time), others are in trial time
varTypes = {'double','logical','string','logical',...
            'double','double','double','double','double',...
            'double','double','double','double',...
            'cell','cell','cell','cell',...
            'cell'};
varNames = {'TrialNumber','Outcome','TrialType','isReward',...
            'CueTime','NextCue','ReactionTime','FirstBoutLastLickTime','TrialLastLickTime',...
            'nLicks','nChoiceLicks','nOutcomeLicks','nBaselineLicks',...
            'AllTrialLicks','ChoiceLicks','OutcomeLicks','BaselineLicks',...
            'boutStartIdx'};


% Initialize trial subtypes
trials = table('Size',[length(allTrials) length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

% Define session related params
reactionTimeSamp = options.reactionTime * options.behaviorFs;

% Loop over all trials
for i=1:length(allTrials)
    cur_cue = allTrials(i);
    if i == length(allTrials); next_cue = NaN;
    else; next_cue = allTrials(i+1); end
    % ITI = (next_cue - cur_cue) / options.behaviorFs;

    % Initialize trial table params
    outcome = 0;
    reactionTime = 0; firstBoutLastLickTime = 0; trialLastLickTime = 0;
    choiceLicks = []; outcomeLicks = []; baselineLicks = [];

    % Trial outcome
    rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue & rightSolenoidON<next_cue) - cur_cue;
    isReward = ~isempty(rewardTime);
    outcomeTime = min([reactionTimeSamp, rewardTime]);

    % Trial licks
    trial_licks = (rightLickON(rightLickON > cur_cue & rightLickON < next_cue) - cur_cue)';
    trialLicks_all = sortrows([trial_licks, 2*ones(length(trial_licks),1)]);
    trialLicks_all = trialLicks_all(:,1);
    nLicks = size(trialLicks_all,1);

    % Separate licks (i.e. anticipatory licks: licks before outcome)
    if ~isempty(trialLicks_all)
        reactionTime = trialLicks_all(1,1);
        trialLastLickTime = trialLicks_all(end,1);

        % Find bout cutoff
        ILI = [100; diff(trialLicks_all)/options.behaviorFs];
        boutStartIdx = find(ILI>1);

        % Select lick within response window
        trialLicks_inWindow = trialLicks_all(trialLicks_all<=reactionTimeSamp,:);
        if size(trialLicks_inWindow,1) >= 4
            choiceLicks = trialLicks_inWindow(1:3);
            remainingLicks = setdiff(trialLicks_all,choiceLicks);

            % Remove spontaneous licks afterwards
            ILI = [10; diff(remainingLicks)/options.behaviorFs];
            boutStart = find(ILI>1,2);
            if length(boutStart) == 1; outcomeLicks = remainingLicks;
            else; outcomeLicks = remainingLicks(boutStart(1):boutStart(2)-1); end
            baselineLicks = setdiff(remainingLicks,outcomeLicks);
            
            if ~isempty(outcomeLicks); firstBoutLastLickTime = outcomeLicks(end);
            else; firstBoutLastLickTime = remainingLicks(end); end

            % define outcome
            outcome = 1;
        else
            choiceLicks = trialLicks_inWindow;
            if ~isempty(choiceLicks); firstBoutLastLickTime = choiceLicks(end); end
            remainingLicks = setdiff(trialLicks_all,choiceLicks);
            
            if ~isempty(remainingLicks)
                % Remove spontaneous licks afterwards
                ILI = [10; diff(remainingLicks)/options.behaviorFs];
                boutStart = find(ILI>1,2);
                if length(boutStart) == 1; outcomeLicks = remainingLicks;
                else; outcomeLicks = remainingLicks(boutStart(1):boutStart(2)-1); end
                baselineLicks = setdiff(remainingLicks,outcomeLicks);
            end

            % define outcome
            outcome = 0;
        end
    end


    % Calculate ENL
%     lastLickPrevTrial = rightLickON(rightLickON < cur_cue); 
%     if isempty(lastLickPrevTrial)
%         ENL = 0;
%     else
%         if i == 1; last_cue = 0; else; last_cue = allTrials(i-1); end
%         ENL = cur_cue - max(lastLickPrevTrial(end),last_cue); 
%     end

    % Get trial type
    if outcome && isReward; trialType = 'Rewarded';
    elseif outcome && ~isReward; trialType = 'Omission';
    else; trialType = 'Miss';
    end

    % Trial table
    trials(i,:) = {i,outcome,trialType,isReward,...
        cur_cue,next_cue,reactionTime,firstBoutLastLickTime,trialLastLickTime,...
        nLicks,size(choiceLicks,1),size(outcomeLicks,1),size(baselineLicks,1),...
        num2cell(trialLicks_all,[1 2]),num2cell(choiceLicks,[1 2]),num2cell(outcomeLicks,[1 2]),num2cell(baselineLicks,[1 2]),...
        num2cell(boutStartIdx,[1 2])};
end

% Add a row for outside events
trials(i+1,:) = {0,0,0,0,...
                0,0,0,0,0,...
                0,0,0,0,...
                {},{},{},{},{}};

end