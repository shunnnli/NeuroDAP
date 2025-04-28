function trialNumber = findTrials(eventTime,trialTable)

% Given a event time, find which trial this event belongs to
arguments
    eventTime double
    trialTable table
end

trialNumber = zeros(size(eventTime));

for i = 1:length(eventTime)
    target = eventTime(i);
    trialNum = trialTable{trialTable.CueTime<=target & trialTable.NextCue>target,"TrialNumber"};

    if ~isempty(trialNum); trialNumber(i) = trialNum; 
    else; trialNumber(i) = height(trialTable); end
end

end