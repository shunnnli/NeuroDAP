function specificDrugs = findSpecificDrugs(excelPath, epochOfInterest, cyclePositionOfInterest)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(excelPath);
    foundDrugs = 'None';
    
    for row = 1:size(data, 1)
        
        epoch = data{row, 3};
        cyclePosition = data{row, 5};

        if epoch == str2double(epochOfInterest) && cyclePosition == str2double(cyclePositionOfInterest)
            
            foundDrugs = data{row, 9};
            foundDrugs = 'None';
            
            if ismissing(foundDrugs)

                foundDrugs = 'None';

            end
            
            break;
                  
        end
        
    end
    
    upperStr = upper(foundDrugs);
    drugs = split(upperStr, ' ');
    drugs = erase(drugs, {',', ';', '-'});
    sortedDrugs = sort(drugs);
    specificDrugs = strjoin(sortedDrugs, ' ');
    
end