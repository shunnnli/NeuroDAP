function generalDrugs = findGeneralDrugs(excelPath)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readcell(excelPath);

    foundDrugs = data{12,2};
    
    if ismissing(foundDrugs)
        
        foundDrugs = 'None';
        
    end
    
    upperStr = upper(foundDrugs);
    drugs = split(upperStr, ' ');
    drugs = erase(drugs, {',', ';', '-'});
    sortedDrugs = sort(drugs);
    generalDrugs = strjoin(sortedDrugs, ' ');
    
 end
