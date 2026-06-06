function qualityRS0 = findQualityRS0(excelPath, epochOfInterest, cyclePositionOfInterest)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(excelPath);

    qualityRS0 = [];
    
    for row = 1:size(data, 1)
        
        epoch = data{row, 3};
        cyclePosition = data{row, 5};

        if epoch == str2double(epochOfInterest) && cyclePosition == str2double(cyclePositionOfInterest)
            
            qualityRS0 = data{row, 42};
            
            if ismissing(qualityRS0)
                
                qualityRS0 = NaN;
                
            end
            
            break;
                  
        end
        
    end
    
end
