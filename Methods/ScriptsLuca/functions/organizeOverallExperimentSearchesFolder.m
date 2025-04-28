function experimentSearchFeaturesTable = organizeOverallExperimentSearchesFolder(mainDirectory,dataDirectory,overallExperimentDirectory,mouseNames,recreateExperimentTable)  

overallExperimentTablePath = [overallExperimentDirectory filesep 'ExperimentSearchFeaturesTable.mat'];
nMouseExperiments = size(mouseNames,2);

for iMouse = 1:nMouseExperiments

    mouseDirectory = [dataDirectory mouseNames{1,iMouse}];
    folders = dir(fullfile(mouseDirectory, '*Results*'));
    folderPath = fullfile(mouseDirectory, folders(1).name); 
    searchTableFileName = ['SearchFeaturesTable.mat'];
    searchTablePath = fullfile(folderPath,searchTableFileName);
    
    mouseSearchFeaturesTable = load(searchTablePath);
    mouseSearchFeaturesTable = mouseSearchFeaturesTable.searchFeaturesTable;
    
    if ~isfile(overallExperimentTablePath) || ((recreateExperimentTable == 1) && (iMouse == 1))
        
        repeatedMouseExperiment = repmat({iMouse}, size(mouseSearchFeaturesTable, 1), 1);
        mouseSearchFeaturesTable = addvars(mouseSearchFeaturesTable, repeatedMouseExperiment, 'Before', 1, 'NewVariableNames', 'mouseNumber');

        overallCellNames = mouseSearchFeaturesTable.mouseCell;
        mouseSearchFeaturesTable = addvars(mouseSearchFeaturesTable, overallCellNames, 'Before', 4, 'NewVariableNames', 'overallCellName');

        experimentSearchFeaturesTable = mouseSearchFeaturesTable;
        save(overallExperimentTablePath, 'experimentSearchFeaturesTable');
       
    else 

        if iMouse == 1
            experimentSearchFeaturesTable = load(overallExperimentSearchTablePath);
            experimentSearchFeaturesTable = experimentSearchFeaturesTable.experimentSearchFeaturesTable;
        end

        if any(strcmp(experimentSearchFeaturesTable.mouseName,mouseNames{1,iMouse}))

            disp([mouseNames{1,iMouse},' already present in experimentSearchFeaturesTable'])
            continue;

        end
        
        repeatedMouseExperiment = repmat({iMouse}, size(mouseSearchFeaturesTable, 1), 1);
        mouseSearchFeaturesTable = addvars(mouseSearchFeaturesTable, repeatedMouseExperiment, 'Before', 1, 'NewVariableNames', 'mouseNumber');

        overallCellNames = cell2mat(mouseSearchFeaturesTable.mouseCell) + max(cell2mat(experimentSearchFeaturesTable.overallCellName)) ;
        mouseSearchFeaturesTable = addvars(mouseSearchFeaturesTable, num2cell(overallCellNames), 'Before', 4, 'NewVariableNames', 'overallCellName');
     
        experimentSearchFeaturesTable = [experimentSearchFeaturesTable; mouseSearchFeaturesTable];

    end

end

% Additional features
radiusSomaPixels = 40;
nHotspotsSoma = cellfun(@(c) sum(c < radiusSomaPixels), experimentSearchFeaturesTable.hotspotDistances);
nHotspots = cellfun(@numel,experimentSearchFeaturesTable.hotspotDistances);
hotspotPrcSoma = (nHotspotsSoma./nHotspots)*100;
experimentSearchFeaturesTable = addvars(experimentSearchFeaturesTable, num2cell(hotspotPrcSoma), 'Before', 12, 'NewVariableNames', 'hotspotPrcSoma');
    
save(overallExperimentTablePath, 'experimentSearchFeaturesTable');

end