function nthSlopes = getNthSlope(DAtrend, n, options)

% Function specific for getting the slope fitting from the last n trials of
% DA recordings
arguments
    DAtrend struct
    n double % find slopes for the last N trials
    
    options.statsType
    options.maxTrial
end

%% Check inputs
if ~isfield(options,'maxTrial')
    % Determine the maximum number of trials across all animals for the given field.
    options.maxTrial = max(arrayfun(@(x) length(x.CueMax), DAtrend));
end

if ~isfield(options,'statsType')
    options.statsType = {'CueMax', 'CueMin', 'CueAvg', 'CueAmp',...
              'CueMax_smoothed','CueMin_smoothed','CueAvg_smoothed','CueAmp_smoothed'};
else
    if isstring(options.statsType)
        options.statsType = {options.statsType};
    end
end

% Preallocate output matrix.
% Columns represent trial numbers starting at 3 (i.e. column 1 is trial 3).
nthSlope = nan(length(DAtrend),1);

%% Loop over trial numbers from 3 to maxTrials.
for s = 1:length(options.statsType)
    st = strcat(options.statsType{s},'_slope');
    for i = 1:length(DAtrend)
        fieldData = DAtrend(i).(st);
        if n > length(fieldData)
            nthSlope(i) = nan;
        else
            nthSlope(i) = fieldData(n);
        end
    end
    nthSlopes.(options.statsType{s}) = nthSlope;
end

end