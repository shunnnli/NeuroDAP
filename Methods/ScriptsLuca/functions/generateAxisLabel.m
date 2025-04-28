function label = generateAxisLabel(fieldName)

    if strcmp(fieldName,'height')
        
        label = 'Peak amplitude [pA]';
        
    elseif strcmp(fieldName,'onsetTime')
        
        label = 'Onset time [ms]';
        
    elseif strcmp(fieldName,'riseTime')
        
        label = 'Rise time [ms]';
        
    elseif strcmp(fieldName,'decayTime')
        
        label = 'Decay time [ms]';
        
    elseif strcmp(fieldName,'area')
        
        label = 'Total charge [...C]';
        
    else
        
        label = 'Add label to generateAxisLabel';
        
    end
        
end