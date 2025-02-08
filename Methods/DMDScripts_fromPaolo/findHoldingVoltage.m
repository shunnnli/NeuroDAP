function holdingVoltage = findHoldingVoltage(excelPath, epochOfInterest, cyclePositionOfReference)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(excelPath);

    holdingVoltage = [];
    
    for row = 1:size(data, 1)
        
        epoch = data{row, 3};
        cyclePosition = data{row, 5};

        if epoch == str2double(epochOfInterest) && cyclePosition == cyclePositionOfReference
            
            holdingVoltage = data{row, 8};
            
            if ismissing(holdingVoltage)
                
                holdingVoltage = nan;
                
            end
            
            break;
                  
        end
        
    end
    
end
