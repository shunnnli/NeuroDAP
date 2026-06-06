function searchData = extractSearchData(searchDataDepth)

    searchData = searchDataDepth;
    nSpots = height(searchDataDepth); 
    nPositiveResponsesColumn = zeros(nSpots, 1);
    searchDataDepth.nPositiveResponses = nPositiveResponsesColumn;

    regions = findgroups(searchDataDepth.depth, searchDataDepth.xStart, searchDataDepth.xWidth, searchDataDepth.yStart, searchDataDepth.yHeight);
    searchDataDepth.regions = regions;
    nRegions = numel(unique(regions));

    for iRegion = 1:nRegions

        regionSearchData = searchDataDepth(searchDataDepth.regions == iRegion,:);
        regionSearchDataRows = find(searchDataDepth.regions == iRegion);
        
        if numel(regionSearchDataRows) == 1
            
            if regionSearchData.nPositiveResponses(1) == 0
            
                nPositiveResponsesRegion = regionSearchData.response(1);
                searchData.nPositiveResponses(regionSearchDataRows) = nPositiveResponsesRegion;
                
            end
            
        else

            nPositiveResponsesRegion = regionSearchData.response(end)+regionSearchData.nPositiveResponses(1);
            searchData.nPositiveResponses(regionSearchDataRows) = nPositiveResponsesRegion;

        end
        
    end

    groups = findgroups(searchData.depth, searchData.xStart, searchData.xWidth, searchData.yStart, searchData.yHeight);
    searchData.groups = groups;
    groupsToKeep = unique(searchData.groups(searchData.response == 1));
    rowsToRemove = searchData(ismember(searchData.groups, groupsToKeep) & searchData.response == 0, :);
    rowsToRemove = rowsToRemove(:,:);
    searchData(ismember(searchData, rowsToRemove, 'rows'),:) = [];
    
    [~, uniqueIdx, ~] = unique(searchData, 'rows');
    searchData = searchData(uniqueIdx, :);
    
    searchData = removevars(searchData, {'groups','nPositiveResponses'});
    
end