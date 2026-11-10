function [weightedQuality,eventWeights,eventRatios] = weightClampQuality(eventSummary,animalList,eventList)
%WEIGHTCLAMPQUALITY Use the same trial-count weights in both conditions.
% For each event, w = min(NClamp,NUnclamp). The combined metric is
% sum(w .* clampMean) / sum(w .* unclampMean), NOT a mean of event ratios.
% Missing event pairs get zero weight. No trials are subsampled.

arguments
    eventSummary table
    animalList (:,1) string
    eventList string = ["Water","Airpuff","Tone"]
end

required = {'Animal','Event','Condition','N','RMSSum','PeakSum'};
if ~isempty(eventSummary) && ~all(ismember(required,eventSummary.Properties.VariableNames))
    error('ClampQuality:EventSummary','Event summary lacks the counts/sums needed for weighting.');
end
nAnimals = numel(animalList);
nEvents = numel(eventList);
values = nan(nAnimals,2,2); % animal x condition x metric
counts = zeros(nAnimals,2);
eventRatios = nan(nAnimals,nEvents,2);
weightCounts = zeros(nAnimals,nEvents,2);
weights = zeros(nAnimals,nEvents);
normalizedWeights = zeros(nAnimals,nEvents);
conditions = ["clamp","unclamp"];

for a = 1:nAnimals
    means = nan(nEvents,2,2); % event x condition x metric
    for e = 1:nEvents
        for c = 1:2
            if isempty(eventSummary); continue; end
            rows = strcmpi(string(eventSummary.Animal),animalList(a)) & ...
                strcmpi(string(eventSummary.Event),eventList(e)) & ...
                strcmpi(string(eventSummary.Condition),conditions(c));
            n = sum(eventSummary.N(rows));
            weightCounts(a,e,c) = n;
            if n <= 0; continue; end
            means(e,c,:) = reshape([sum(eventSummary.RMSSum(rows)), ...
                sum(eventSummary.PeakSum(rows))] / n,1,1,2);
        end
        if all(isfinite(means(e,:,:)),'all')
            weights(a,e) = min(weightCounts(a,e,:));
        end
        for m = 1:2
            if isfinite(means(e,2,m)) && means(e,2,m) > 0
                eventRatios(a,e,m) = means(e,1,m) / means(e,2,m);
            end
        end
    end
    included = weights(a,:) > 0;
    if ~any(included); continue; end
    normalizedWeights(a,:) = weights(a,:) / sum(weights(a,:));
    for c = 1:2
        counts(a,c) = sum(weightCounts(a,included,c));
        for m = 1:2
            values(a,c,m) = sum(normalizedWeights(a,included)' .* means(included,c,m));
        end
    end
end

ratios = nan(nAnimals,2);
for m = 1:2
    valid = isfinite(values(:,2,m)) & values(:,2,m) > 0;
    ratios(valid,m) = values(valid,1,m) ./ values(valid,2,m);
end
weightedQuality = table(animalList,values(:,1,1),values(:,2,1),ratios(:,1), ...
    values(:,1,2),values(:,2,2),ratios(:,2),counts(:,1),counts(:,2), ...
    'VariableNames',{'Animal','ClampRMS','UnclampRMS','RMSRatio','ClampPeak', ...
    'UnclampPeak','PeakRatio','NClamp','NUnclamp'});
weightAnimals = repelem(animalList,nEvents);
eventWeights = table(weightAnimals(:),repmat(eventList(:),nAnimals,1), ...
    reshape(weightCounts(:,:,1)',[],1),reshape(weightCounts(:,:,2)',[],1), ...
    reshape(weights',[],1),reshape(normalizedWeights',[],1), ...
    'VariableNames',{'Animal','Event','NClamp','NUnclamp','Weight','NormalizedWeight'});
end
