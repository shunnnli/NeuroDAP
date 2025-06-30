function results = groupEIbyDA(DAdata, combinedAnimalList, DAanimalList, options)

arguments
    DAdata double
    combinedAnimalList
    DAanimalList

    options.nGroups double = 3
    options.groupLabels double
end

numGroups = options.nGroups;

if ~isfield(options,'groupLabels')
    % Perform k-means clustering with replicates for stability
    [groupIdx, groupCenters] = kmeans(DAdata, numGroups, Replicates=100);
    
    % Sort group centers so that:
    %   Group 1: obvious negative (lowest center)
    %   Group 2: ambiguous (middle center)
    %   Group 3: obvious positive (highest center)
    [~, sortOrder] = sort(groupCenters);
    mapping = zeros(numGroups,1);
    mapping(sortOrder) = 1:numGroups;
    groupLabels = mapping(groupIdx);
    results.groupLabels = groupLabels;
else
    groupLabels = options.groupLabels;
end

% Save grouping results
if numGroups == 2
    downAnimals_kmeans = find(groupLabels == 1);
    upAnimals_kmeans = find(groupLabels == 2);
    results.down.animals = DAanimalList(downAnimals_kmeans);
    results.up.animals = DAanimalList(upAnimals_kmeans);
elseif numGroups == 3
    downAnimals_kmeans = find(groupLabels == 1);
    stableAnimals_kmeans = find(groupLabels == 2);
    upAnimals_kmeans = find(groupLabels == 3);
    results.down.animals = DAanimalList(downAnimals_kmeans);
    results.up.animals = DAanimalList(upAnimals_kmeans);
    results.stable.animals = DAanimalList(stableAnimals_kmeans);
end

% (Analysis) Plot cell EI clustered by kmeans
% Set up getCellIndices function
getCellIndices = @(indices) cell2mat(arrayfun(@(idx) find(strcmpi(combinedAnimalList, DAanimalList{idx})), indices, 'UniformOutput', false));

% Get the animal indices from quant
results.down.cells = getCellIndices(downAnimals_kmeans);
results.up.cells = getCellIndices(upAnimals_kmeans);
if numGroups == 3 
    results.stable.cells = getCellIndices(stableAnimals_kmeans);
end


end