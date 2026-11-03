function state = dmdSpotPlotState(isHotspot,nSweepsStored,nSweepsKept)

% How each spot should be drawn, as one code per spot:
%
%    1  a hotspot under the notebook's definition      -> full opacity
%    0  scored, but not a hotspot                      -> transparent
%   -1  stimulated, but no sweep survived quality control -> grey
%
% The third state matters because it is not the same answer as the second.
% A spot whose sweeps all failed the Rs limit has no metrics at all, so it
% is not "not a hotspot" -- it is unknown, and drawing it like a scored
% null spot would claim something the data cannot support.
%
% QC exclusion outranks the verdict: a metric computed from no surviving
% sweeps cannot be believed, whatever it says. A grid position that was
% never stimulated has nothing to draw and is left at 0.

arguments
    isHotspot
    nSweepsStored double
    nSweepsKept double
end

isHotspot = logical(isHotspot(:));
nSweepsStored = nSweepsStored(:);
nSweepsKept = nSweepsKept(:);

state = zeros(numel(isHotspot),1);
state(isHotspot) = 1;
state(nSweepsStored > 0 & nSweepsKept == 0) = -1;

end
