function toBeExcluded = findExcludedRuns(excelPath, epochOfInterest, cyclePositionOfInterest)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(excelPath);

    foundString = [];
    
    for row = 1:size(data, 1) % 24 with readcell
        
        epoch = data{row, 3};
        cyclePosition = data{row, 5};

        if epoch == str2double(epochOfInterest) && cyclePosition == str2double(cyclePositionOfInterest)
            
            foundString = data{row, 7};
            break; 
            
        end
    end

    if strcmp(foundString, 'exclude')
        
        toBeExcluded = 1;
        
    else
        
        toBeExcluded = 0;
        
    end

end
