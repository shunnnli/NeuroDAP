function experimentAnalysisTable = findIdenticalRunsFastVC(experimentAnalysisTable)

    structArray = vertcat(experimentAnalysisTable.runProtocol{:});
    allParametersTable = struct2table(structArray);
    allParametersTable.activeChannels = cellfun(@(c) strjoin(string(c(:)), ', '),allParametersTable.activeChannels, 'UniformOutput', false);

    for i = 1:height(allParametersTable)
        entry = allParametersTable.options(i); 

        if isstruct(entry)
            entryFields = struct2table(entry);
            allParametersTable.newOptions(i) = strjoin(string(table2cell(entryFields)), ', '); 
        else
            allParametersTable.newOptions(i) = string(entry);
        end
    end

    allParametersTable = removevars(allParametersTable, {'cellCoordinates','delayPulseBlue','delayPulseRed','options'});
    allParametersTable = renamevars(allParametersTable, 'newOptions', 'options');
    
    for i = 1:width(allParametersTable)
        columnName = allParametersTable.Properties.VariableNames{i};
        columnData = allParametersTable.(columnName);
        columnData = string(columnData);
        columnData(ismissing(columnData)) = string("NaN");
        allParametersTable.(columnName) = columnData;
    end

    runGroupArray = findgroups(allParametersTable);
    experimentAnalysisTable.runGroup = runGroupArray;

end