function [eventTable,extractRange] = scoreClampEventRMS(eventTable,channel,params,options)
%SCORECLAMPEVENTRMS Size each event response on one photometry channel.
%
% Adds two columns to the table from getClampEventTable:
%   RMS    root mean square of the baseline-subtracted response window. RMS is
%          direction-independent, so peaks and dips both count as responses;
%          a flat (well clamped) response gives a small RMS.
%   Valid  true when the whole extracted trace is finite.
% extractRange spans both scoring windows, so a caller can require that an
% event's full baseline and response fall inside its window.

arguments
    eventTable table
    channel (1,1) struct
    params (1,1) struct
    options.baselineWindow_sec (1,2) double {mustBeFinite} = [-2,0]
    options.responseWindow_sec (1,2) double {mustBeFinite} = [0,2]
end

if diff(options.baselineWindow_sec) <= 0 || options.baselineWindow_sec(2) > 0 || ...
        diff(options.responseWindow_sec) <= 0
    error('scoreClampEventRMS:InvalidWindow', ...
        'Use increasing windows, with the baseline ending at or before event delivery.');
end

extractRange = [min(options.baselineWindow_sec(1),options.responseWindow_sec(1)), ...
                max(options.baselineWindow_sec(2),options.responseWindow_sec(2))];
[eventTraces,eventTime] = plotTraces(eventTable.Sample,extractRange,channel.data,params, ...
    signalFs=channel.finalFs,signalSystem=channel.system,eventSystem='ni',plot=false,rmmissing=false);
baselineIdx = eventTime >= options.baselineWindow_sec(1)-1e-9 & eventTime < options.baselineWindow_sec(2)-1e-9;
responseIdx = eventTime >= options.responseWindow_sec(1)-1e-9 & eventTime <= options.responseWindow_sec(2)+1e-9;
if ~any(baselineIdx) || sum(responseIdx) < 2
    error('scoreClampEventRMS:NotEnoughSamples','Not enough DA samples in the scoring windows.');
end

response = eventTraces(:,responseIdx)-mean(eventTraces(:,baselineIdx),2);
eventTable.RMS = sqrt(mean(response.^2,2));
eventTable.Valid = all(isfinite(eventTraces),2) & isfinite(eventTable.RMS);
end
