function holdingVoltage = findHoldingVoltage(excelPath, epochOfInterest, cyclePositionOfInterest)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(excelPath);

    holdingVoltage = [];
    
    for row = 1:size(data, 1) %24 with readcell
        
        epoch = data{row, 3};
        cyclePosition = data{row, 5};

        if epoch == str2double(epochOfInterest) && cyclePosition == str2double(cyclePositionOfInterest)
            
            holdingVoltage = data{row, 8};
            
            if ismissing(holdingVoltage)
                
                holdingVoltage = nan;
                
            end
            
            break;
                  
        end
        
    end
    
end
