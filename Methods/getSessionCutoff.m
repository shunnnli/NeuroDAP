function [trials,cutoff_sample] = getSessionCutoff(trials,type,options)

%{
Several criteria:
    1. total trial number have to be more than 100
    2. >=n misses within 10 trials (moving window)
%}

% max_misses_allowed = 3; windowlength = 10;

arguments
    trials table
    type string
    options.max_misses_allowed double = 3
    options.windowlength double = 10
    options.minTrials double = 80 % start after 80 trials
end

max_misses_allowed = options.max_misses_allowed;
windowlength = options.windowlength;

if contains(type,'->reward') %strcmp(type,'baseline->reward') || strcmp(type,'punish->reward')
    rewardTrials = trials{trials.isReward == 1,["TrialNumber","Outcome"]}; % trials.isStim before 07/21/23
    outcome = strcmp(rewardTrials(:,2),'H');
    rewardTrials = [str2double(rewardTrials(:,1)) outcome];
    hitnum = movsum(rewardTrials(:,2),[0 windowlength-1],"Endpoints","discard");
    
    % Find candidate cutoffs
    candidate_cutoffs = find(hitnum < (windowlength - max_misses_allowed));
    candidate_cutoffs = candidate_cutoffs(candidate_cutoffs >= options.minTrials); 
    
    % Choose the first one as cutoff
    if ~isempty(candidate_cutoffs)
        cutoff_row = rewardTrials(candidate_cutoffs(1),1);
        % Create a new column in trials called performing
        performing = zeros(size(trials,1),1);
        performing(1:cutoff_row) = 1;
        trials.performing = performing;
        % Get sample of cutoff
        cutoff_sample = trials{cutoff_row,'NextCue'};
    else
        trials.performing = ones(size(trials,1),1);
        cutoff_sample = trials{end,"CueTime"} + 10*10000;
    end

    % Initialize a new column in trials called stage (i.e. learning stage)
    trials.stage = zeros(size(trials,1),1);
else
    % Create a new column in trials called performing
    trials.performing = ones(size(trials,1),1);
    cutoff_sample = trials{end,"CueTime"} + 10*10000;
    % Initialize a new column in trials called stage (i.e. learning stage)
    trials.stage = zeros(size(trials,1),1);
end

end