function staticWindow = getStaticPeriod(X, Y, options)

arguments
    X double
    Y double
    options.windowDuration double = 30 % in seconds
    options.noMovementThreshold double = 5

    options.Fs double = 20
end

movingWindowInSamples = options.windowDuration * options.Fs;

diffXinWindow = movsum(diff(X),movingWindowInSamples);
diffYinWindow = movsum(diff(Y),movingWindowInSamples);
distanceInWindow = diffXinWindow + diffYinWindow;
still = abs(distanceInWindow) < options.noMovementThreshold;
staticWindow = still;

% dStill = diff([0; still; 0]);   % +1 at block start, -1 at block end
% starts = find(dStill == 1);
% ends   = find(dStill == -1) - 1;
% staticWindow = false(height(data),1); 
% for k = 1:numel(starts)
%     % Duration of this still block:
%     blockDuration = data.time(ends(k)) - data.time(starts(k));
%     if blockDuration >= options.windowDuration/60
%         staticWindow(starts(k):ends(k)) = true;
%     end
% end

end