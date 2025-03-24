function stats = getDAtrend(DAtrend, t1, t2, options)

% Function specific for getting the slope fitting from the last n trials of
% DA recordings
arguments
    DAtrend struct
    t1 {mustBeInteger}
    t2 {mustBeInteger}
    
    options.animalIdx double
    options.mapType string = 'slope'
    options.dataType string = 'raw'
    options.maxTrial
end

%% Check inputs
if ~isfield(options,'maxTrial')
    % Determine the maximum number of trials across all animals for the given field.
    options.maxTrial = max(cell2mat({DAtrend.nTrials}));
end

if sum(contains(options.mapType,["slope","slopes"]))
    mapType = 'slopeMap';
elseif sum(contains(options.mapType,["diff","dif"]))
    mapType = 'diffMap';
elseif sum(contains(options.mapType,["avg","averages","mean"]))
    mapType = 'avgMap';
end

if contains(options.dataType,"raw")
    dataType = 'raw';
elseif contains(options.dataType,"smooth")
    dataType = 'smoothed';
end

if ~isfield(options,'animalIdx')
    options.animalIdx = 1:length(DAtrend);
end

%% Loop over trial numbers from 3 to maxTrials.

DAstatsType = {'max', 'min', 'avg', 'amp'};

% Preallocate output matrix.
% Columns represent trial numbers starting at 3 (i.e. column 1 is trial 3).
cur_stat = nan(length(options.animalIdx),1);

for s = 1:length(DAstatsType)
    % st = DAtrend.(DAstatsType{s},statType);
    for a = 1:length(options.animalIdx)
        animalIdx = options.animalIdx(a);
        fieldData = DAtrend(animalIdx).(DAstatsType{s}).(mapType).(dataType).map;
        if t1 > length(fieldData) || t2 > length(fieldData) || t1 <=0 || t2 <= 0 || t1>t2
            cur_stat(a) = nan;
        else
            cur_stat(a) = fieldData(t1,t2);
        end
    end
    stats.(DAstatsType{s}) = cur_stat;
end

end