function trials = getTrialTable_shijiaCatch(events,options)

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
% spontaneousLicks (in session time), others are in trial time
varTypes = {'double','logical','string','logical','string',...
    'double','double','double','double','double','double','double',...
    'double','double','double','double',...
    'cell','cell','cell','cell',...
    'cell'};
varNames = {'TrialNumber','Outcome','TrialType','isReward','RewardType',...
    'CueTime','NextCue','RewardTime','RewardLickTime','ReactionTime','FirstBoutLastLickTime','TrialLastLickTime',...
    'nLicks','nChoiceLicks','nOutcomeLicks','nSpontaneousLicks',...
    'AllTrialLicks','ChoiceLicks','OutcomeLicks','SpontaneousLicks',...
    'boutStartIdx'};


% Initialize trial subtypes
trials = table('Size',[length(allTrials) length(varNames)],...
    'VariableTypes',varTypes,'VariableNames',varNames);

% Define session related params
reactionTimeSamp = options.reactionTime * options.behaviorFs;
%%
initializeFig(1,1);
% Loop over all trials
for i=1:length(allTrials)
    cur_cue = allTrials(i);
    if i == length(allTrials); next_cue = NaN;
    else; next_cue = allTrials(i+1); end
    % ITI = (next_cue - cur_cue) / options.behaviorFs;

    % Initialize trial table params
    outcome = 0;
    rewardTimeInTrial = 0; reactionTime = 0; trialLastLickTime = 0;
    rewardTimeUse = 0; reactionTimeUse = 0; firstBoutLastLickTimeUse = 0; trialLastLickTimeUse = 0;

    choiceLicks = []; outcomeLicks = []; spontaneousLicks = []; rewardTime = []; rewardLickTime = 0;
    boutStartIdx = [];

    % Trial outcome
    rewardTime = rightSolenoidON(rightSolenoidON>=cur_cue & rightSolenoidON<next_cue) - cur_cue;
    isReward = ~isempty(rewardTime);
    outcomeTime = min([reactionTimeSamp, rewardTime]);

    % rewardTimeDebug = rightSolenoidON(rightSolenoidON>=cur_cue & rightSolenoidON<next_cue);
    % if ~isempty(rewardTimeDebug)
    %     rewardTimeDebugPlot = rewardTimeDebug(1);
    % end

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
            % 
            %     if isReward
            % 
            % %     % plot the licks and reward onset to see if they are good
            % figure;
            % data = zeros(size(trialLicks_all));
            % data(trialLicks_all)=1;
            % scatter(1:length(data),data,'b');
            % legend('original licks');
            % hold on;
            % data2 = zeros(size(rewardTime));
            % data2(rewardTime)=1;
            % scatter(1:length(data2),data2,'r');
            % legend('water reward');
            %     end


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
            spontaneousLicks = setdiff(remainingLicks,outcomeLicks);

            if ~isempty(outcomeLicks); firstBoutLastLickTimeUse = outcomeLicks(end)+cur_cue;
            else; firstBoutLastLickTimeUse = remainingLicks(end)+cur_cue; end

            % define outcome
            outcome = 1;
        else
            choiceLicks = trialLicks_inWindow;
            if ~isempty(choiceLicks); firstBoutLastLickTimeUse = choiceLicks(end)+cur_cue; end
            remainingLicks = setdiff(trialLicks_all,choiceLicks);

            if ~isempty(remainingLicks)
                % Remove spontaneous licks afterwards
                ILI = [10; diff(remainingLicks)/options.behaviorFs];
                boutStart = find(ILI>1,2);
                if length(boutStart) == 1; outcomeLicks = remainingLicks;
                else; outcomeLicks = remainingLicks(boutStart(1):boutStart(2)-1); end
                spontaneousLicks = setdiff(remainingLicks,outcomeLicks);
            end

            % define outcome
            outcome = 0;
        end
    end

    closestLickToReward = [];

    % Get trial type
    if outcome && isReward; trialType = 'Rewarded';
        rewardTimeInTrial = rewardTime(1);
        [~,closestLickToReward] = min(abs(rewardTimeInTrial-trialLicks_all));
        if closestLickToReward == 2
            RewardType = 'lick2';
            if length(trialLicks_all)>2
            rewardLickTime = trialLicks_all(closestLickToReward+1)  + cur_cue;
            else
            rewardLickTime = trialLicks_all(closestLickToReward)  + cur_cue;
            end
        elseif closestLickToReward == 6
            RewardType = 'lick6';
            if length(trialLicks_all)>6
            rewardLickTime = trialLicks_all(closestLickToReward+1) + cur_cue;
            else
            rewardLickTime = trialLicks_all(closestLickToReward) + cur_cue;
            end
        else
            RewardType = 'lick4';
            if length(trialLicks_all)>4
            rewardLickTime = trialLicks_all(closestLickToReward+1) + cur_cue;
            else
            rewardLickTime = trialLicks_all(closestLickToReward) + cur_cue;
            end
        end
        rewardTimeUse = rewardTimeInTrial + cur_cue;
    elseif outcome && ~isReward; trialType = 'Omission'; RewardType = 'NoReward';
    else; trialType = 'Miss'; RewardType = 'NoReward';
    end

    % Trial table
    trials(i,:) = {i,outcome,trialType,isReward,RewardType,...
        cur_cue,next_cue,rewardTimeUse,rewardLickTime,reactionTimeUse,firstBoutLastLickTimeUse,trialLastLickTimeUse,... % add 'Use': add the time of cur_cue
        nLicks,size(choiceLicks,1),size(outcomeLicks,1),size(spontaneousLicks,1),... % for analyzing only the first bout: last lick: outcomeLicks(end); non-last lick: [choiceLicks + outcomeLicks(1:end-1)]
        num2cell(trialLicks_all+cur_cue,[1 2]),num2cell(choiceLicks+cur_cue,[1 2]),num2cell(outcomeLicks+cur_cue,[1 2]),num2cell(spontaneousLicks+cur_cue,[1 2]),...
        num2cell(boutStartIdx,[1 2])};


    % For every trial, plot the licks and reward onset to see if they are good
    exampleTrialIdx = floor(linspace(1,length(allTrials)-1,20));

    if ismember(i,exampleTrialIdx)
        nexttile;
        if ~isempty(trialLicks_all)
            data = zeros(size(trialLicks_all));
            data(trialLicks_all)=1;
            scatter(1:length(data),data,'b');
            hold on;
        end

        if strcmp(trialType,'Rewarded')
            data2 = zeros(size(rewardTime));
            data2(rewardTime)=1;
            scatter(1:length(data2),data2,'r');
            hold on;

            data3 = zeros(size(rewardTimeInTrial));
            data3(rewardTimeInTrial)=1;
            scatter(1:length(data3),data3,'.y');
        end

        legend({'original licks','water reward','rewardtimeintrial'});
        xlim([0 options.behaviorFs*2]);
        if ~isempty(closestLickToReward)
            % if exist("closestLickToReward","var")
            title(['Trial' num2str(i) '_' num2str(closestLickToReward)]);
        else
            title(['Trial' num2str(i)]);
        end
    end

    % Add a row for outside events
    trials(i+1,:) = {0,0,0,0,0,...
        0,0,0,0,0,0,0,...
        0,0,0,0,...
        {},{},{},{},{}};


end

% cd('\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Shijia\BE28 BE31 Task_photometry\BE29 M164-167\20231205_M167L_4W40O40C10D10_g0');
% saveas(gcf, strcat(extractBefore(sessionName,'_g'),'_ExampleTrialEventScatterPlot.tif'));
