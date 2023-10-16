function [timeLow,l] = assignTimeStamp(timeLow,timeHigh,idx_low,idx_high,lowFs)

% Fill in time
idx_prev_low = 0;
l = min([length(idx_high),length(idx_low)]);
for i = 1:l
    % Sync the time using the rising edge of the current sync pulse
    timeLow(idx_low(i)) = timeHigh(idx_high(i));
    % Filled in time in between two pulses
    if idx_prev_low < idx_low(i)
        timeLow(idx_prev_low+1 : idx_low(i)) = ...
            timeLow(idx_low(i)) + (1/lowFs)*(idx_prev_low-idx_low(i)+1 :1: 0);
    end
    % Filled in time after the last syncNI pulse
    if i == l
        timeLow(idx_low(i):end) = ...
            timeLow(idx_low(i)) + (1/lowFs)*(0 :1: length(timeLow)-idx_low(i));
    end

    idx_prev_low = idx_low(i);
end

end