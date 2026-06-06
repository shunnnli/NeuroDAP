function outputTable = extractSelectedMemberFromTable(inputTable, columnName1, selectedMember1, columnName2, selectedMember2)

    if nargin == 3

        selectedMemberRows = find(ismember(cell2mat(inputTable.(columnName1)), selectedMember1));
        outputTable = inputTable(selectedMemberRows,:);
    
    elseif nargin == 5
        
        outputTable = [];
        
        for iMember1 = 1:size(selectedMember1,2)

            selectedMemberRows = find(ismember(cell2mat(inputTable.(columnName1)), iMember1) & ismember(cell2mat(inputTable.(columnName2)), selectedMember2{iMember1}));
            outputSubTable = inputTable(selectedMemberRows,:);
            
            if isempty(outputTable)
                
                outputTable = outputSubTable;
                
            else
                
                outputTable = [outputTable; outputSubTable];
                
            end
            
        end
    
    else
        
        disp('Wrong number of input arguments')
        
    end

end
