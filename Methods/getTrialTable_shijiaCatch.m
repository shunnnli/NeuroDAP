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
varTypes = {'double','logical','string','logical','string',...
    'double','double','double','double','double','double',...
    'double','double','double','double',...
    'cell','cell','cell','cell',...
    'cell'};
varNames = {'TrialNumber','Outcome','TrialType','isReward','RewardType',...
    'CueTime','NextCue','RewardTime','ReactionTime','FirstBoutLastLickTime','TrialLastLickTime',...
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
    rewardTimeInTrial = 0; reactionTime = 0; firstBoutLastLickTime = 0; trialLastLickTime = 0; 
    rewardTimeUse = 0; reactionTimeUse = 0; firstBoutLastLickTimeUse = 0; trialLastLickTimeUse = 0; 

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

    %     % plot the licks to see if they are good
    %     figure;
    %     data = zeros(size(trialLicks_all));
    %     data(trialLicks_all)=1;
    %     scatter(1:length(data),data,'b');
    %     legend('original licks');
    %     hold on;
    %     data2 = zeros(size(trialLicks_all_cleaned));
    %     data2(trialLicks_all_cleaned)=1;
    %     scatter(1:length(data2),data2,'r');
    %     legend('cleaned licks');


    % clean the licks with 90ms lick ITI
    lickITICutoffArduino = 90; % in ms

    % Separate licks (i.e. anticipatory licks: licks before outcome)
    if ~isempty(trialLicks_all)

        % clean the licks with 90ms lick ITI - no need since we did it in concatLabJack_shijia
        previousLickOnset = trialLicks_all(1);
        trialLicks_all_cleaned = [previousLickOnset];
        for p = 1:length(trialLicks_all)
            if trialLicks_all(p)-previousLickOnset>(lickITICutoffArduino/1000*options.behaviorFs)
                previousLickOnset = trialLicks_all(p);
                trialLicks_all_cleaned = [trialLicks_all_cleaned previousLickOnset];
            end
        end
        trialLicks_all = trialLicks_all_cleaned';

%         if isReward

            %     %     % plot the licks and reward onset to see if they are good
            %     figure;
            %     data = zeros(size(trialLicks_all));
            %     data(trialLicks_all)=1;
            %     scatter(1:length(data),data,'b');
            %     legend('original licks');
            %     hold on;
            %     data2 = zeros(size(rewardTime));
            %     data2(rewardTime)=1;
            %     scatter(1:length(data2),data2,'r');
            %     legend('water reward');
%         end


        reactionTime = trialLicks_all(1,1);
        trialLastLickTime = trialLicks_all(end,1);
        reactionTimeUse = trialLicks_all(1,1) + cur_cue;
        trialLastLickTimeUse = trialLicks_all(end,1) + cur_cue;

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
            if ~isempty(choiceLicks); firstBoutLastLickTime = choiceLicks(end); firstBoutLastLickTimeUse = choiceLicks(end)+cur_cue; end
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


    % Get trial type
    if outcome && isReward; trialType = 'Rewarded';
        rewardTimeInTrial = rewardTime(1);
        [~,closestLickToReward] = min(abs(rewardTimeInTrial-trialLicks_all));
        if closestLickToReward == 2
            RewardType = 'lick2';
        elseif closestLickToReward == 6
            RewardType = 'lick6';
        else
            RewardType = 'lick4';
        end
        rewardTimeUse = rewardTimeInTrial + cur_cue;
    elseif outcome && ~isReward; trialType = 'Omission'; RewardType = 'NoReward';
    else; trialType = 'Miss'; RewardType = 'NoReward';
    end

    % Trial table
    trials(i,:) = {i,outcome,trialType,isReward,RewardType,...
        cur_cue,next_cue,rewardTimeUse,reactionTimeUse,firstBoutLastLickTimeUse,trialLastLickTimeUse,... % add 'Use': add the time of cur_cue
        nLicks,size(choiceLicks,1),size(outcomeLicks,1),size(baselineLicks,1),... % for analyzing only the first bout: last lick: outcomeLicks(end); non-last lick: [choiceLicks + outcomeLicks(1:end-1)]
        num2cell(trialLicks_all+cur_cue,[1 2]),num2cell(choiceLicks+cur_cue,[1 2]),num2cell(outcomeLicks+cur_cue,[1 2]),num2cell(baselineLicks+cur_cue,[1 2]),...
        num2cell(boutStartIdx,[1 2])};
end

% Add a row for outside events
trials(i+1,:) = {0,0,0,0,0,...
    0,0,0,0,0,0,...
    0,0,0,0,...
    {},{},{},{},{}};

end