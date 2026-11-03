function [windowMs, controlWindows] = dmdDetectionWindows(timeMs,windowMs,gapMs)

% Response window and every equally long control window that follows it.
%
% The control windows are exactly as long as the response window. That matters
% once each window is measured around its own median: the de-levelled charge
% density of a drifting trace grows with window length, so a shorter control
% window would make every spot look responsive.
%
% Mirrors detection_windows() in pyNeuroDAP's Shun_dmd_analysis notebook.

arguments
    timeMs double
    windowMs (1,2) double
    gapMs (1,1) double
end

duration = windowMs(2) - windowMs(1);
latest = max(timeMs(:),[],'omitnan');

controlWindows = zeros(0,2);
cursor = windowMs(2) + gapMs;
while cursor + duration <= latest + 1e-9
    controlWindows(end+1,:) = [cursor, cursor+duration]; %#ok<AGROW>
    cursor = cursor + duration;
end

if isempty(controlWindows)
    error('dmdDetectionWindows:noControlWindow',...
          'No %.0f ms control window fits after %.0f ms; the stored window ends at %.0f ms.',...
          duration,windowMs(2),latest);
end

end
