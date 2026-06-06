function selectedRunGroupsAnalysisTable = extractEpochSubset(selectedRunGroupsAnalysisTable, desiredEpochs)

    if isempty(desiredEpochs)

        return;

    elseif size(desiredEpochs,2) == 2

        if desiredEpochs(1) == -1

            selectedEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) <= desiredEpochs(2));

        elseif desiredEpochs(2) == -1

            selectedEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= desiredEpochs(1));

        else
            
            selectedEpochsRows = find(cell2mat(selectedRunGroupsAnalysisTable.epoch) >= desiredEpochs(1) & cell2mat(selectedRunGroupsAnalysisTable.epoch) <= desiredEpochs(2));
        
        end
        
    else 
        
        selectedEpochsRows = find(ismember(cell2mat(selectedRunGroupsAnalysisTable.epoch), desiredEpochs));

    end
        
    selectedRunGroupsAnalysisTable = selectedRunGroupsAnalysisTable(selectedEpochsRows,:);

end