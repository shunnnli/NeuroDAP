function targetIdx = findCorrespondingTime(eventIdx,timeRef,timeTarget)

% Find sample index of the target system corresponding to reference system

targetIdx = zeros(size(eventIdx));

for i = 1:length(eventIdx)
    % if behaviorFs is not integer, not rounding will raise error
    [~, eventTime_target] = min(abs(timeTarget-timeRef(round(eventIdx(i))))); 
    targetIdx(i) = eventTime_target;
end

end