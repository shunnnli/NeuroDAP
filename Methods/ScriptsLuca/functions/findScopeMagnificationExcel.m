function scopeMagnification = findScopeMagnificationExcel(data, epochOfInterest, cyclePositionOfInterest)

    foundMagnification = 'None';
    
    for row = 1:size(data, 1)
        
        epoch = data{row, 3};
        cyclePosition = data{row, 5};

        if epoch == str2double(epochOfInterest) && cyclePosition == str2double(cyclePositionOfInterest)
            
            foundMagnification = data{row, 10};
            
            if ismissing(foundMagnification)

                foundMagnification = '40X';

            end
            
            break;
                  
        end
        
    end
    
    upperStr = upper(foundMagnification);
    scopeMagnification = split(upperStr, ' ');
    scopeMagnification = erase(scopeMagnification, {',', ';', '-'});
    scopeMagnification = strjoin(scopeMagnification, ' ');
    
 end