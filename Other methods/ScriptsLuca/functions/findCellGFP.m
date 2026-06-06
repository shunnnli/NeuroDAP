function cellGFP = findCellGFP(excelPath)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readcell(excelPath);

    foundGFP = data{14,2};
    
    if ismissing(foundGFP)
        
        foundGFP = 'NO';
        
    end
    
    upperStr = upper(foundGFP);
    cellGFP = split(upperStr, ' ');
    cellGFP = erase(cellGFP, {',', ';', '-'});
    cellGFP = strjoin(cellGFP, ' ');
    
 end