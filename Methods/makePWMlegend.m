function legendStr = makePWMlegend(events, color)
% events: table from extractPWMStim
% color : "blue" or "red"

if isempty(events)
    legendStr = '';
    return;
end

[grp, dutyVals] = findgroups(events.duty);
counts = splitapply(@numel, events.onset, grp);

legendStr = strings(numel(dutyVals),1);
for i = 1:numel(dutyVals)
    legendStr(i) = sprintf("%s PWM %d (n=%d)", color, dutyVals(i), counts(i));
end
end