function cellTemperature = findCellTemperatureExcel(data)

    foundTemp = data{11,2};
    
    if ismissing(foundTemp)
        
        foundTemp = NaN;
        
    end
    
    cellTemperature = foundTemp;
    
 end