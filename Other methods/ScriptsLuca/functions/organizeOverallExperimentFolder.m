function [experimentAnalysisTable, runGroupsParametersTable, experimentCellsConnectionInfo] = organizeOverallExperimentFolder(matlabDirectory,mainDirectory,overallExperimentName,mouseExperimentNames,recreateExperimentTable)  

    overallExperimentDirectory = [matlabDirectory 'GroupedExperiments' filesep overallExperimentName];
    overallExperimentAnalysisTablePath = [overallExperimentDirectory filesep 'ExperimentAnalysisTable.mat'];
    overallExperimentCellsConnectionInfoPath = [overallExperimentDirectory filesep 'ExperimentCellsConnectionInfo.mat'];
    
    nMouseExperiments = size(mouseExperimentNames,2);

    for iMouse = 1:nMouseExperiments

        mouseDirectory = [mainDirectory mouseExperimentNames{1,iMouse} filesep 'mouseAnalysis'];
        mouseAnalysisTablePath = [mouseDirectory filesep 'MouseAnalysisTable.mat'];
        mouseAnalysisTable = load(mouseAnalysisTablePath);
        mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
        
        cellsConnectionInfoPath = [mouseDirectory filesep 'CellsConnectionInfo.mat'];
        cellsConnectionInfo = load(cellsConnectionInfoPath);
        cellsConnectionInfo = cellsConnectionInfo.cellsConnectionInfo;
        
        saveExperimentTablePath = [overallExperimentDirectory filesep 'ExperimentAnalysisTable.mat'];
        saveExperimentCellsConnectionInfoPath = [overallExperimentDirectory filesep 'ExperimentCellsConnectionInfo.mat'];
              
        if ~isfile(overallExperimentAnalysisTablePath) || ((recreateExperimentTable == 1) && (iMouse == 1))

            repeatedMouseExperiment = repmat({iMouse}, size(mouseAnalysisTable, 1), 1);
            mouseAnalysisTable = addvars(mouseAnalysisTable, repeatedMouseExperiment, 'Before', 1, 'NewVariableNames', 'mouseNumber');

            repeatedMouseExperimentName = repmat({mouseExperimentNames{1,iMouse}}, size(mouseAnalysisTable, 1), 1);
            mouseAnalysisTable = addvars(mouseAnalysisTable, repeatedMouseExperimentName, 'Before', 2, 'NewVariableNames', 'mouseName');
 
            overallCellNames = mouseAnalysisTable.cellName;
            mouseAnalysisTable = addvars(mouseAnalysisTable, overallCellNames, 'Before', 4, 'NewVariableNames', 'overallCellName');
          
            experimentAnalysisTable = mouseAnalysisTable;
            save(saveExperimentTablePath, 'experimentAnalysisTable');
            
            experimentCellsConnectionInfo = [];
            experimentCellsConnectionInfo.blue = [];
            experimentCellsConnectionInfo.blue.connectedCells = {cell2mat(cellsConnectionInfo.blue.connectedCells)};
            experimentCellsConnectionInfo.blue.nonConnectedCells = {cell2mat(cellsConnectionInfo.blue.nonConnectedCells)};
            experimentCellsConnectionInfo.red = [];
            experimentCellsConnectionInfo.red.connectedCells = {cell2mat(cellsConnectionInfo.red.connectedCells)};
            experimentCellsConnectionInfo.red.nonConnectedCells = {cell2mat(cellsConnectionInfo.red.nonConnectedCells)};
            experimentCellsConnectionInfo.doubleConnectedCells = {cellsConnectionInfo.doubleConnectedCells};
            
            experimentCellsConnectionInfo.excitatory = [];
            experimentCellsConnectionInfo.excitatory.connectedCells = [];
            experimentCellsConnectionInfo.excitatory.nonConnectedCells = [];
            experimentCellsConnectionInfo.inhibitory = [];
            experimentCellsConnectionInfo.inhibitory.connectedCells = [];
            experimentCellsConnectionInfo.inhibitory.nonConnectedCells = [];
            
            if strcmp(cellsConnectionInfo.responseType.blue,'excitatory')
                experimentCellsConnectionInfo.excitatory.connectedCells = {cell2mat(cellsConnectionInfo.blue.connectedCells)};
                experimentCellsConnectionInfo.excitatory.nonConnectedCells = {cell2mat(cellsConnectionInfo.blue.nonConnectedCells)};
            end
                
            if strcmp(cellsConnectionInfo.responseType.blue,'inhibitory')
                experimentCellsConnectionInfo.inhibitory.connectedCells = {cell2mat(cellsConnectionInfo.blue.connectedCells)};
                experimentCellsConnectionInfo.inhibitory.nonConnectedCells = {cell2mat(cellsConnectionInfo.blue.nonConnectedCells)};
            end
                
            if strcmp(cellsConnectionInfo.responseType.red,'excitatory')
                experimentCellsConnectionInfo.excitatory.connectedCells = {cell2mat(cellsConnectionInfo.red.connectedCells)};
                experimentCellsConnectionInfo.excitatory.nonConnectedCells = {cell2mat(cellsConnectionInfo.red.nonConnectedCells)};
            end
                
            if strcmp(cellsConnectionInfo.responseType.red,'inhibitory')
                experimentCellsConnectionInfo.inhibitory.connectedCells = {cell2mat(cellsConnectionInfo.red.connectedCells)};
                experimentCellsConnectionInfo.inhibitory.nonConnectedCells = {cell2mat(cellsConnectionInfo.red.nonConnectedCells)};
            end
            
            save(saveExperimentCellsConnectionInfoPath, 'experimentCellsConnectionInfo');

        else 
            
            if iMouse == 1
                experimentAnalysisTable = load(overallExperimentAnalysisTablePath);
                experimentAnalysisTable = experimentAnalysisTable.experimentAnalysisTable;
                experimentCellsConnectionInfo = load(overallExperimentCellsConnectionInfoPath);
                experimentCellsConnectionInfo = experimentCellsConnectionInfo.experimentCellsConnectionInfo;
            end

            if any(strcmp(experimentAnalysisTable.mouseName,mouseExperimentNames{1,iMouse}))
                
                disp([mouseExperimentNames{1,iMouse},' already present in experimentAnalysisTable'])
                continue;

            end
            
            repeatedMouseExperiment = repmat({iMouse}, size(mouseAnalysisTable, 1), 1);
            mouseAnalysisTable = addvars(mouseAnalysisTable, repeatedMouseExperiment, 'Before', 1, 'NewVariableNames', 'mouseNumber');
            
            repeatedMouseExperimentName = repmat({mouseExperimentNames{1,iMouse}}, size(mouseAnalysisTable, 1), 1);
            mouseAnalysisTable = addvars(mouseAnalysisTable, repeatedMouseExperimentName, 'Before', 2, 'NewVariableNames', 'mouseName');
             
            overallCellNames = cell2mat(mouseAnalysisTable.cellName) + max(cell2mat(experimentAnalysisTable.overallCellName)) ;
            mouseAnalysisTable = addvars(mouseAnalysisTable, num2cell(overallCellNames), 'Before', 4, 'NewVariableNames', 'overallCellName');
            
            experimentCellsConnectionInfo.blue.connectedCells{end+1} = cell2mat(cellsConnectionInfo.blue.connectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            experimentCellsConnectionInfo.blue.nonConnectedCells{end+1} = cell2mat(cellsConnectionInfo.blue.nonConnectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            experimentCellsConnectionInfo.red.connectedCells{end+1} = cell2mat(cellsConnectionInfo.red.connectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            experimentCellsConnectionInfo.red.nonConnectedCells{end+1} = cell2mat(cellsConnectionInfo.red.nonConnectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            experimentCellsConnectionInfo.doubleConnectedCells{end+1} = cellsConnectionInfo.doubleConnectedCells + max(cell2mat(experimentAnalysisTable.overallCellName));
            
            if strcmp(cellsConnectionInfo.responseType.blue,'excitatory')
                experimentCellsConnectionInfo.excitatory.connectedCells{end+1} = cell2mat(cellsConnectionInfo.blue.connectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
                experimentCellsConnectionInfo.excitatory.nonConnectedCells{end+1} = cell2mat(cellsConnectionInfo.blue.nonConnectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            end
                
            if strcmp(cellsConnectionInfo.responseType.blue,'inhibitory')
                experimentCellsConnectionInfo.inhibitory.connectedCells{end+1} = cell2mat(cellsConnectionInfo.blue.connectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
                experimentCellsConnectionInfo.inhibitory.nonConnectedCells{end+1} = cell2mat(cellsConnectionInfo.blue.nonConnectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            end
                
            if strcmp(cellsConnectionInfo.responseType.red,'excitatory')
                experimentCellsConnectionInfo.excitatory.connectedCells{end+1} = cell2mat(cellsConnectionInfo.red.connectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
                experimentCellsConnectionInfo.excitatory.nonConnectedCells{end+1} = cell2mat(cellsConnectionInfo.red.nonConnectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            end
                
            if strcmp(cellsConnectionInfo.responseType.red,'inhibitory')
                experimentCellsConnectionInfo.inhibitory.connectedCells{end+1} = cell2mat(cellsConnectionInfo.red.connectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
                experimentCellsConnectionInfo.inhibitory.nonConnectedCells{end+1} = cell2mat(cellsConnectionInfo.red.nonConnectedCells) + max(cell2mat(experimentAnalysisTable.overallCellName));
            end

            experimentAnalysisTable = [experimentAnalysisTable; mouseAnalysisTable];
            
        end

    end
    
    experimentAnalysisTable = removevars(experimentAnalysisTable, 'runGroup');
    experimentAnalysisTable = findIdenticalRunsFastVC(experimentAnalysisTable);
    save(saveExperimentTablePath, 'experimentAnalysisTable');
    save(saveExperimentCellsConnectionInfoPath, 'experimentCellsConnectionInfo');

    runGroups = unique(experimentAnalysisTable.runGroup(:));
    nRunGroups = numel(runGroups);

    columnNames = {'runGroup', 'sampleSize', 'activeChannels', 'holdingVoltage', 'scopeMagnification', 'generalDrugs', 'specificDrugs', 'cellLocation', 'cellTemperature', 'cellType', 'cellFluo', 'cellCoordinates', 'responseType', 'nPulsesBlue', 'isiBlue', 'pulseWidthBlue', 'amplitudeBlue', 'delayPulseBlue', 'functionNameBlue', 'nPulsesRed', 'isiRed', 'pulseWidthRed', 'amplitudeRed', 'delayPulseRed', 'functionNameRed', 'whichPulseFirst', 'pulsesTimeDifference', 'options'};
    dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
    runGroupsParametersTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);

    for iRunGroup = 1:nRunGroups

        rowsRunGroup = find(experimentAnalysisTable.runGroup(:) == iRunGroup);
        sampleSizeRunGroup = size(rowsRunGroup,1);
        runGroupAnalysisTable = experimentAnalysisTable(rowsRunGroup,:);

        runProtocol = runGroupAnalysisTable.runProtocol{1};
        runProtocolCell = struct2cell(runProtocol);

        activeChannels = experimentAnalysisTable.runProtocol{1,1}.activeChannels{1};

        if numel(runGroupAnalysisTable.runProtocol{1,1}.activeChannels) == 2

            secondActiveChannel = runGroupAnalysisTable.runProtocol{1,1}.activeChannels{2};
            activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}]};

        end

        dataCellToAdd = [{iRunGroup}, {sampleSizeRunGroup}, activeChannels, runProtocolCell(2:end)'];

        runGroupsParametersTable = [runGroupsParametersTable; dataCellToAdd];

        saveMouseParametersTablePath = [overallExperimentDirectory filesep 'ExperimentParametersTable.mat'];
        save(saveMouseParametersTablePath, 'runGroupsParametersTable');

    end

end