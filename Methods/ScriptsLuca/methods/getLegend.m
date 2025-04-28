function legend = getLegend(events,labels,options)

arguments
    events cell
    labels cell
    options.ignore double = []% indexes to ignore
end

% Check events and labels are the same size
if length(events) ~= length(labels)
    error('Number of events and labels does not match!!');
end

% Find out empty or ignored labels/events
legend = {};
for i = 1:length(events)
    if isempty(events{i}) || ismember(i,options.ignore)
        continue
    else
        eventLegend = [labels{i},' (n=',num2str(length(events{i})),')'];
        legend(end+1) = {eventLegend};
    end
end

end