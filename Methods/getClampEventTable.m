function eventTable = getClampEventTable(behavior,timeNI)
%GETCLAMPEVENTTABLE Collect water / airpuff / tone deliveries with clamp state.
%
% Reuses the classification saved by analyzeSessions_clamp (waterClampIdx,
% waterUnclampIdx, ...), whose clampON_event marks clamp activity within 3 s
% after each delivery. Returns a table sorted by sample with:
%   Sample    native NI sample of the delivery
%   Type      "Water", "Airpuff" or "Tone"
%   IsClamp   true when the delivery was clamped
%   Time_sec  seconds since the first NI sample

arguments
    behavior (1,1) struct
    timeNI {mustBeNonempty}
end

eventFields = {'waterClampIdx','waterUnclampIdx','airpuffClampIdx', ...
               'airpuffUnclampIdx','toneClampIdx','toneUnclampIdx'};
if ~all(isfield(behavior,eventFields))
    error('getClampEventTable:MissingClassification', ...
        'Run analyzeSessions_clamp for this random session to save water/airpuff/tone clamp classifications.');
end

eventTypes = ["Water","Water","Airpuff","Airpuff","Tone","Tone"];
eventSamples = []; eventType = strings(0,1); eventIsClamp = false(0,1);
for i = 1:length(eventFields)
    idx = unique(double(behavior.(eventFields{i})(:)));
    idx = idx(isfinite(idx) & idx >= 1 & idx <= numel(timeNI) & idx == fix(idx));
    eventSamples = [eventSamples;idx]; %#ok<AGROW>
    eventType = [eventType;repmat(eventTypes(i),length(idx),1)]; %#ok<AGROW>
    eventIsClamp = [eventIsClamp;repmat(mod(i,2) == 1,length(idx),1)]; %#ok<AGROW>
end

eventTable = table(eventSamples,eventType,eventIsClamp, ...
    'VariableNames',{'Sample','Type','IsClamp'});
eventTable = sortrows(eventTable,'Sample');
if isempty(eventTable)
    eventTable.Time_sec = zeros(0,1);
else
    eventTable.Time_sec = reshape(timeNI(eventTable.Sample)-timeNI(1),[],1);
end
end
