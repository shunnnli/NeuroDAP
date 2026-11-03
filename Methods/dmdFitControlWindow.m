function [ctrlData,ctrlTime] = dmdFitControlWindow(ctrlData,nAnalysisSamples,samplesPerMs)

% The stored control (pre-stimulus) window, trimmed to sit alongside the
% response window, with the time vector that goes with it.
%
% The control and response windows are NOT the same width. The response
% window spans analysisWindowLength ms inclusive of the sample at 0, while
% the control window stops baselineGuard ms short of the pulse so the
% leading edge of the artefact cannot leak into the baseline. At 10 kHz with
% the defaults that is 301 samples against 291. Callers that assumed equal
% widths indexed off the front of the control array.
%
% The time vector is spaced by the real sample interval and ends at 0, so a
% short window reads as a short window rather than being stretched to cover
% the nominal span.

arguments
    ctrlData double
    nAnalysisSamples (1,1) double
    samplesPerMs (1,1) double
end

if isempty(ctrlData)
    ctrlData = [];
    ctrlTime = [];
    return
end

ctrlData = ctrlData(:,dmdTailSamples(size(ctrlData,2),nAnalysisSamples));
ctrlTime = (-(size(ctrlData,2)-1):0) / samplesPerMs;

end
