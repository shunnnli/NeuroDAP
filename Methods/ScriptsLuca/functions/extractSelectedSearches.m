function selectedSearchesTable = extractSelectedSearches(experimentSearchFeaturesTable,selectedSearches)

selectedRows = false(height(experimentSearchFeaturesTable), 1);

for iMouse = 1:numel(selectedSearches.mouseNames)
    
    mouseName = selectedSearches.mouseNames{iMouse};
    mouseCellsList = selectedSearches.mouseCells{iMouse};
    mouseEpochsList = selectedSearches.mouseEpochs{iMouse};

    matchMouse = strcmp(experimentSearchFeaturesTable.mouseName, mouseName);
    matchCell = ismember(cell2mat(experimentSearchFeaturesTable.mouseCell), mouseCellsList);
    matchEpoch = false(height(experimentSearchFeaturesTable), 1);
    
    for iEpoch = 1:numel(mouseEpochsList)
        matchEpoch = matchEpoch | ismember(cell2mat(experimentSearchFeaturesTable.mouseEpoch), mouseEpochsList{iEpoch});
    end
    
    selectedRows = selectedRows | (matchMouse & matchCell & matchEpoch);

end

selectedSearchesTable = experimentSearchFeaturesTable(selectedRows, :);

end