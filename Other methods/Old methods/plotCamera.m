function [traces,timestamp] = plotCamera(timeRange,eventIdx,eye,timeNI,timeCamera,options)

arguments
    timeRange double
    eventIdx double
    eye double
    timeNI double
    timeCamera double
    options.behaviorFs double = 10000
    options.camFs double = 143.5948
end

timestamp = floor(options.camFs*(timeRange(2)-timeRange(1)));
traces = zeros(length(eventIdx),timestamp+1);

% Iterate over all airpuffs
for i = 1:length(eventIdx)
    % Find first & last NI index
    niFirstIdx = eventIdx(i) + floor(timeRange(1)*options.behaviorFs);
    
    % Find closest camera index
    [~, camFirstIdx] = min(abs(timeCamera-timeNI(niFirstIdx)));
    camLastIdx = camFirstIdx + timestamp;

%     disp(['niFirstIdx: ',num2str(timeNI(niFirstIdx))]);
%     disp(['camFirstIdx: ',num2str(timeCamera(camFirstIdx))]);
    
    % Extract eye intensity traces
    if camLastIdx > length(eye)
        trace = nan(1,timestamp+1);
        % trace(camFirstIdx:length(eye)) = eye(camFirstIdx:length(eye),2)';
    else
        trace = eye(camFirstIdx:camLastIdx)';
    end
    traces(i,:) = trace;
end

eventInSec = 

% getTraces
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

end % getEye