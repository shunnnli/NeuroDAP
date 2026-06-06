function runCyclePosition = findCorrectCyclePosition(excelPath, runEpoch, currentCyclePosition, useCyclePosition)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(excelPath);
    
    runCyclePosition = useCyclePosition;
    
    for row = 1:size(data, 1) % 24 with readcell
        
        epoch = data{row, 3};
        cyclePosition = data{row, 5};

        if epoch == str2double(runEpoch) && cyclePosition == str2double(currentCyclePosition)
            
            runCyclePosition = currentCyclePosition;
            break; 
            
        end
        
    end

end
