function [rawtraces,timestamp] = getdff(eventInSec,photometry,timeRange,binSize)

% Get photometry traces (dff) around a specific event
% eventInSec: 1d vector that lists the event time (predivided by nidq.Fs)

% timestamp = timeRange(1):binSize:timeRange(2);
% traces = zeros(length(eventInSec),length(timestamp));

% for i = 1:length(eventInSec)
% 
%     if (floor(eventInSec(i)+timeRange(1)) <= 0) || isnan(eventInSec(i))
%         traces(i,:) = nan(1,length(timestamp));
%         continue
%     end
% 
%     firstBin = floor((eventInSec(i)+timestamp(1))/binSize);
%     lastBin = firstBin + length(timestamp) - 1;
%     
%     if lastBin > length(photometry)
%         % baseline = mean(photometry(firstBin:floor((eventInSec(i))/params.downsample_bin)));
%         dff = nan(1,length(timestamp));
%         % dff() = (photometry(firstBin:length(photometry))-baseline)/baseline *100;
%     else 
%         baseline = mean(photometry(firstBin:floor((eventInSec(i))/binSize)));
%         dff = (photometry(firstBin:lastBin)-baseline)/baseline *100;
%     end
%     traces(i,:) = dff;
% end

% Get z-scored/raw/whatever traces
[rawtraces,timestamp] = getTraces(eventInSec,photometry,timeRange,binSize);

baselineBins = 1:abs(timeRange(1)/binSize);
for i = 1:length(eventInSec)
    baseline = mean(rawtraces(i,baselineBins));
    dff = 100*(rawtraces(i,:)-baseline)/baseline;
    rawtraces(i,:) = dff;
end

end