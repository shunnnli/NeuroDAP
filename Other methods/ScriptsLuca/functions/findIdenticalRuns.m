function analysisTable = findIdenticalRuns(analysisTable)

    tags = nan(size(analysisTable, 1),1);
    iTag = 1;

    for i = 1:size(analysisTable, 1)
    
        currentStruct = analysisTable.optoParameters{i,1}; %(i)
        identical = false;
    
        for j = 1:i-1
        
            previousStruct = analysisTable.optoParameters{j,1}; %(j)

            if isequaln(currentStruct, previousStruct)
                
                identical = true;
                tags(i) = tags(j);
                break;
            
            end
            
        end
    
        if ~identical
            
            tags(i) = iTag;
            iTag = iTag + 1;
        
        end
        
    end

    analysisTable.runGroup = tags;

end