function phAnalysisStructureRun = phAnalysisExtractData(phAnalysisFile, runEpoch, currentRun)
 
    phAnalysisStructureRun = [];
    epochName = strcat('e', num2str(runEpoch));
    phAnalysisStruct = phAnalysisFile.(epochName);
    phAnalysisFields = fieldnames(phAnalysisStruct);
    
    [~, currentRun, ~] = fileparts(currentRun);
    runIndex = find(strcmp(phAnalysisStruct.acqs, currentRun), 1);
    
    if isempty(runIndex)
        
        disp('Mismatch!')
        
    end
    
    for iField = 1:numel(phAnalysisFields)
        
        fieldName = phAnalysisFields{iField};
        
        if (strcmp(fieldName, 'sagV') || strcmp(fieldName, 'reboundV'))
            
            continue;
            
        end
        
        if isstruct(phAnalysisStruct.(fieldName))
            
            fieldValue = phAnalysisStruct.(fieldName){1,runIndex};
            
        else
            
            fieldValue = phAnalysisStruct.(fieldName)(runIndex);
            
        end
        
        phAnalysisStructureRun = setfield(phAnalysisStructureRun, fieldName, fieldValue);
        
    end

end
