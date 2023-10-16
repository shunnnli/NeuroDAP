function [traces,timestamp] = getTraces(eventInSec,photometry,timeRange,binSize)

% Get photometry traces (zscore) around a specific event
% eventInSec: 1d vector that lists the event time (predivided by nidq.Fs)

timestamp = timeRange(1):binSize:timeRange(2);
traces = zeros(length(eventInSec),length(timestamp));

for i = 1:length(eventInSec)

    if (floor(eventInSec(i)+timeRange(1)) <= 0) || isnan(eventInSec(i))
        traces(i,:) = nan(1,length(timestamp));
        continue
    end

    firstBin = floor((eventInSec(i)+timestamp(1))/binSize);
    lastBin = firstBin + length(timestamp) - 1;
    
    if lastBin > length(photometry)
        trace = nan(1,length(timestamp));
    else
        trace = photometry(firstBin:lastBin);
    end
    traces(i,:) = trace;
end

end