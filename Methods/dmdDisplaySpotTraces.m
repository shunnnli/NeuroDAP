function displayTraces = dmdDisplaySpotTraces(traces,keptRows)
% Selected sweeps for measured spots; original data for gray QC-only display.
displayTraces = traces;
for s=1:numel(traces)
    if any(keptRows{s})
        displayTraces{s} = traces{s}(keptRows{s},:);
    end
end
end
