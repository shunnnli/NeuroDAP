function cellFluo = findCellFluoExcel(data)

    foundCellFluo = data{15,2};
    
    if ismissing(foundCellFluo)
        
        foundCellFluo = 'NO';
        
    end
    
    upperStr = upper(foundCellFluo);
    cellFluo = split(upperStr, ' ');
    cellFluo = erase(cellFluo, {',', ';', '-'});
    cellFluo = strjoin(cellFluo, ' ');
    
 end