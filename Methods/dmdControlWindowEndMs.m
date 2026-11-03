function endMs = dmdControlWindowEndMs(opts)

% Time of the control window's LAST sample, relative to the stimulus, read
% from the options stored alongside the data.
%
% loadSlicesDMD stops the control window baselineGuard ms short of the pulse
% so the leading edge of the stimulus artifact stays out of the baseline.
% The window therefore ends at -baselineGuard, not at 0, and a plot that
% assumes 0 draws the whole baseline trace that much late -- closing up the
% very interval the guard was added to exclude.
%
% Files written before the guard existed carry no such field. Their control
% window really did run up to the pulse, so 0 is the honest answer for them
% and their figures render as they always have.

arguments
    opts struct
end

if isfield(opts,'baselineGuard')
    endMs = -opts.baselineGuard;
else
    endMs = 0;
end

end
