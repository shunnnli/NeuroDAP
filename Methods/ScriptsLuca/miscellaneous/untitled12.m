matchingGroups = {};

for iGroupsPerProtocol = 1:numel(runGroupsProtocol)

    runGroupsParametersForComparison = removevars(runGroupsParametersTableCell,{'runGroup','sampleSize'});
    parametersToCompare = setdiff(runGroupsParametersForComparison.Properties.VariableNames, {'specificDrugs'});
    referenceGroup = runGroupsProtocol(iGroupsPerProtocol);
    rowReferenceGroup = find(cell2mat(runGroupsParametersTableCell.runGroup) == referenceGroup);
    referenceGroupParameters = runGroupsParametersForComparison(rowReferenceGroup,:);

    for iRowGroup = 1:height(runGroupsParametersForComparison)
        if iRowGroup == rowReferenceGroup; continue; end
        if isequaln(referenceGroupParameters{1, parametersToCompare}, runGroupsParametersForComparison{iRowGroup, parametersToCompare}); 
            matchingGroups{end+1} = runGroupsParametersTableCell.runGroup{iRowGroup};
        end
    end
    
end
