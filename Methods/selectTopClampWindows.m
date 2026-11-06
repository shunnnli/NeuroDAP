function [topWindows,topTraces] = selectTopClampWindows(windows,timeSeries,params,timeNI,options)
%SELECTTOPCLAMPWINDOWS Keep the best non-overlapping candidate windows per channel.
%
% windows must already be sorted best first and carry Start_sec, End_sec and
% ChannelIndex. Candidates that overlap describe nearly the same stretch of the
% session, so only the best of each overlapping group is kept, and a window is
% kept only when getSampleClampTrace can extract it with finite data.
%
% topWindows keeps every input column and adds Rank (1 = best of that channel),
% ReferenceIdx and TimeRange, so a caller can plot or re-extract the window
% without repeating the lookup. topTraces holds the matching extracted samples.

arguments
    windows table
    timeSeries struct
    params (1,1) struct
    timeNI {mustBeNonempty}
    options.nTop (1,1) double {mustBeInteger,mustBePositive} = 5
    options.requireFiniteClamp (1,1) logical = true
end

required = {'Start_sec','End_sec','ChannelIndex'};
if ~all(ismember(required,windows.Properties.VariableNames))
    error('selectTopClampWindows:MissingColumns', ...
        'windows needs Start_sec, End_sec and ChannelIndex.');
end

windows.Rank = zeros(height(windows),1);
windows.ReferenceIdx = zeros(height(windows),1);
windows.TimeRange = zeros(height(windows),2);
topWindows = windows([],:); topTraces = {};

for ch = reshape(unique(windows.ChannelIndex,'stable'),1,[])
    channelRows = windows(windows.ChannelIndex == ch,:);
    kept = zeros(0,2);
    for i = 1:height(channelRows)
        if size(kept,1) >= options.nTop; break; end
        % Overlaps the window already kept above it, so it adds nothing new
        if any(channelRows.Start_sec(i) < kept(:,2) & channelRows.End_sec(i) > kept(:,1)); continue; end
        [~,referenceIdx] = min(abs(timeNI-(timeNI(1)+channelRows.Start_sec(i))));
        timeRange = [channelRows.Start_sec(i),channelRows.End_sec(i)]-(timeNI(referenceIdx)-timeNI(1));
        try
            sample = getSampleClampTrace(timeSeries,params,referenceIdx,timeRange, ...
                daIdx=ch,removeClampBaseline=false);
        catch ME
            if strcmp(ME.identifier,'getSampleClampTrace:WindowOutsideRecording'); continue; end
            rethrow(ME)
        end
        check = sample.DA;
        if options.requireFiniteClamp; check = [check,sample.blueClamp,sample.redClamp]; end %#ok<AGROW>
        if any(~isfinite(check)); continue; end
        kept(end+1,:) = [channelRows.Start_sec(i),channelRows.End_sec(i)]; %#ok<AGROW>
        keptRow = channelRows(i,:);
        keptRow.Rank = size(kept,1);
        keptRow.ReferenceIdx = referenceIdx;
        keptRow.TimeRange = timeRange;
        topWindows = [topWindows;keptRow]; %#ok<AGROW>
        topTraces{end+1} = sample; %#ok<AGROW>
    end
end
end
