function idx = dmdTailSamples(n,maxSamples)

% Index of the last maxSamples columns of an n-column array, clamped so it
% never runs off the front.
%
% Written out as n-maxSamples+1 : n, this is a trap: the moment the array is
% shorter than the span asked for, the start index goes zero or negative and
% MATLAB rejects it. That is not a hypothetical -- the stored control window
% is 291 samples at 10 kHz (30 ms less the 1 ms guard before the pulse),
% while several callers ask it for 201 or 301. Clamping here means a short
% window yields a short slice instead of an error.

arguments
    n (1,1) double
    maxSamples (1,1) double
end

idx = max(1,n-maxSamples+1) : n;

end
