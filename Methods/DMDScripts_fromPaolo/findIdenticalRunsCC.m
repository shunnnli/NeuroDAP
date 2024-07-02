function analysisTable = findIdenticalRunsCC(analysisTable)

    tags = nan(size(analysisTable, 1),1);
    iTag = 1;

    for i = 1:size(analysisTable, 1)
    
        currentStruct = analysisTable.optoParameters{i,1};
        fieldsToNeglectCurrentStruct = {'delayCurrentPulse', 'delayLightPulse'};

        for iFieldCurrentStruct = 1:numel(fieldsToNeglectCurrentStruct)
            
            if isfield(currentStruct,fieldsToNeglectCurrentStruct{iFieldCurrentStruct})

                currentStruct = rmfield(currentStruct, fieldsToNeglectCurrentStruct{iFieldCurrentStruct});
            
            end
            
        end
        
        identical = false;
    
        for j = 1:i-1
        
            previousStruct = analysisTable.optoParameters{j,1};
            fieldsToNeglectPreviousStruct = {'delayCurrentPulse', 'delayLightPulse'};

            for iFieldPreviousStruct = 1:numel(fieldsToNeglectPreviousStruct)
            
                if isfield(previousStruct,fieldsToNeglectPreviousStruct{iFieldPreviousStruct})
                
                    previousStruct = rmfield(previousStruct, fieldsToNeglectPreviousStruct{iFieldPreviousStruct});

                end
                
            end

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