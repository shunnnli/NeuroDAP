function cellRFP = findCellRFP(excelPath)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readcell(excelPath);

    foundRFP = data{15,2};
    
    if ismissing(foundRFP)
        
        foundRFP = 'NO';
        
    end
    
    upperStr = upper(foundRFP);
    cellRFP = split(upperStr, ' ');
    cellRFP = erase(cellRFP, {',', ';', '-'});
    cellRFP = strjoin(cellRFP, ' ');
    
 end