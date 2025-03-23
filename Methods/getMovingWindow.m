function windowIndices = getMovingWindow(dataSize,windowSize,options)

arguments
    dataSize double
    windowSize double

    options.reverse logical = 0 % reverse window order (used for fitting last n trials)
end

halfBefore = floor((windowSize-1)/2);
halfAfter  = ceil((windowSize-1)/2);
windowIndices = nan(dataSize,2);

% Get indices for each sliding window position
for i = 1:dataSize
    idx_start = max(1, i - halfBefore);
    idx_end   = min(dataSize, i + halfAfter);
    windowIndices(i,:) = [idx_start, idx_end];
end

% Reverse order if needed
if options.reverse
    windowIndices = flip(windowIndices);
end

end