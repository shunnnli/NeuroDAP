function cellTemperature = findCellTemperature(excelPath)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readcell(excelPath);

    foundTemp = data{11,2};
    
    if ismissing(foundTemp)
        
        foundTemp = NaN;
        
    end
    
    cellTemperature = foundTemp;
    
 end