function cellType = findCellTypeExcel(data)

    foundCellType = data{14,2};
    
    if ismissing(foundCellType)
        
        foundCellType = 'NONE';
        
    end
    
    upperStr = upper(foundCellType);
    cellType = split(upperStr, ' ');
    cellType = erase(cellType, {',', ';', '-'});
    cellType = strjoin(cellType, ' ');
    
 end