function cellLocation = findCellLocation(excelPath)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readcell(excelPath);

    foundLocation = data{13,2};
    
    if ismissing(foundLocation)
        
        foundLocation = 'None';
        
    end
    
    upperStr = upper(foundLocation);
    location = split(upperStr, ' ');
    location = erase(location, {',', ';', '-'});
    sortedLocation = sort(location);
    cellLocation = strjoin(sortedLocation, ' ');
    
 end
