function cellGFP = findCellGFPExcel(data)

    foundGFP = data{14,2};
    
    if ismissing(foundGFP)
        
        foundGFP = 'NO';
        
    end
    
    upperStr = upper(foundGFP);
    cellGFP = split(upperStr, ' ');
    cellGFP = erase(cellGFP, {',', ';', '-'});
    cellGFP = strjoin(cellGFP, ' ');
    
 end