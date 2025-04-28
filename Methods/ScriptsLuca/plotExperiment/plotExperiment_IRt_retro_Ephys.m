% Clearing workspace
clear all;
close all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Decide saving options
savePlots = 0;

% Indicate experiments to group
overallExperimentName = 'IRt_retro_Ephys';
mouseExperimentNames = {};
mouseExperimentNames{1} = 'IRt_retro_Ephys_1_240401';
mouseExperimentNames{2} = 'IRt_retro_Ephys_2';

% Indicate cells to select
mouseExperimentSelectedCells = {};
mouseExperimentSelectedCells{1} = [1,2];
mouseExperimentSelectedCells{2} = [1,2,3];

% Defining overall experiment directory
overallExperimentDirectory = [matlabDirectory 'GroupedExperiments' filesep overallExperimentName];
overallExperimentAnalysisTablePath = [overallExperimentDirectory filesep 'ExperimentAnalysisTable.mat'];
savePlotPath = overallExperimentDirectory;

nMouseExperiments = size(mouseExperimentNames,2);

for iMouse = 1:nMouseExperiments
    
    mouseDirectory = [mainDirectory mouseExperimentNames{1,iMouse} filesep 'mouseAnalysis'];
    mouseAnalysisTablePath = [mouseDirectory filesep 'MouseAnalysisTable.mat'];
    mouseAnalysisTable = load(mouseAnalysisTablePath);
    mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
    
    saveExperimentTablePath = [overallExperimentDirectory filesep 'ExperimentAnalysisTable.mat'];
    
    if ~isfile(overallExperimentAnalysisTablePath)
        
        mkdir(overallExperimentDirectory);
        addpath(genpath(overallExperimentDirectory));

        repeatedMouseExperiment = repmat({iMouse}, size(mouseAnalysisTable, 1), 1);
        mouseAnalysisTable = addvars(mouseAnalysisTable, repeatedMouseExperiment, 'Before', 1, 'NewVariableNames', 'mouse');
        
        overallCellNames = mouseAnalysisTable.cellName;
        mouseAnalysisTable = addvars(mouseAnalysisTable, overallCellNames, 'Before', 3, 'NewVariableNames', 'overallCellName');
        
        experimentAnalysisTable = mouseAnalysisTable;
        save(saveExperimentTablePath, 'experimentAnalysisTable');
        
    else
             
        experimentAnalysisTable = load(overallExperimentAnalysisTablePath);
        experimentAnalysisTable = experimentAnalysisTable.experimentAnalysisTable;
        
        if any(cell2mat(experimentAnalysisTable.mouse) == iMouse)
        
            continue;
        
        end
        
        repeatedMouseExperiment = repmat({iMouse}, size(mouseAnalysisTable, 1), 1);
        mouseAnalysisTable = addvars(mouseAnalysisTable, repeatedMouseExperiment, 'Before', 1, 'NewVariableNames', 'mouse');
        
        overallCellNames = cell2mat(mouseAnalysisTable.cellName) + cell2mat(experimentAnalysisTable.overallCellName(end)) ;
        mouseAnalysisTable = addvars(mouseAnalysisTable, num2cell(overallCellNames), 'Before', 3, 'NewVariableNames', 'overallCellName');
        
        experimentAnalysisTable = [experimentAnalysisTable; mouseAnalysisTable];
    
    end
    
    experimentAnalysisTable = removevars(experimentAnalysisTable, 'runGroup');
    experimentAnalysisTable = findIdenticalRunsCC(experimentAnalysisTable);
    save(saveExperimentTablePath, 'experimentAnalysisTable');

end

runGroups = unique(experimentAnalysisTable.runGroup(:));
nRunGroups = numel(runGroups);

columnNames = {'runGroup', 'sampleSize', 'activeChannels', 'delayCurrentPulse', 'widthCurrentPulse', 'amplitudeCurrentPulse', 'nLightPulses', 'widthLightPulse', 'amplitudeLightPulse', 'delayLightPulse'};
dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
runGroupsParametersTable = table('Size', [0, numel(columnNames)], 'VariableNames', columnNames, 'VariableTypes', dataTypes);

for iRunGroup = 1:nRunGroups
    
    rowsRunGroup = find(experimentAnalysisTable.runGroup(:) == iRunGroup);
    sampleSizeRunGroup = size(rowsRunGroup,1);
    runGroupAnalysisTable = experimentAnalysisTable(rowsRunGroup,:);
    
    optoParameters = runGroupAnalysisTable.optoParameters{1};
    optoParametersCell = struct2cell(optoParameters);
    
    activeChannels = experimentAnalysisTable.optoParameters{1,1}.activeChannels{1};
    
    if numel(runGroupAnalysisTable.optoParameters{1,1}.activeChannels) == 2
        
        secondActiveChannel = runGroupAnalysisTable.optoParameters{1,1}.activeChannels{2};
        activeChannels = {[activeChannels{1}, ', ', secondActiveChannel{1}]};
    
    end
    
    dataCellToAdd = [{iRunGroup}, {sampleSizeRunGroup}, activeChannels, optoParametersCell(2:end)'];
    
    runGroupsParametersTable = [runGroupsParametersTable; dataCellToAdd];
    
    saveMouseParametersTablePath = [overallExperimentDirectory filesep 'ExperimentParametersTable.mat'];
    save(saveMouseParametersTablePath, 'runGroupsParametersTable');
    
end

%% amplitudeCurrentPulse vs areaVtPulse

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'areaVtPulse';
selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000);
selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'mouse', mouseExperimentNames, 'cellName', mouseExperimentSelectedCells);

generateLinePlot(queryOptoParameter, queryRunFeature, selectedRunGroups, selectedMembersTable, overallExperimentName, savePlotPath, savePlots)

%% process data from pulseAP struct vs amplitudeCurrentPulse

queryOptoParameter =  'amplitudeCurrentPulse';
selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0') & cell2mat(runGroupsParametersTable.widthCurrentPulse) == 6000);
selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'mouse', mouseExperimentNames, 'cellName', mouseExperimentSelectedCells);
selectedRunGroupsAnalysisTable = selectedMembersTable(ismember(selectedMembersTable.runGroup,selectedRunGroups),:);

[featuresStruct, queryOptoParameterValues, queryOptoParameterArray] = extractFeaturesPulseAP(queryOptoParameter, selectedRunGroups, selectedRunGroupsAnalysisTable);

%% nAP

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'nAP';
selectedCellsIDs = unique(cell2mat(selectedRunGroupsAnalysisTable.overallCellName));

generateBarPlotAPpulse(featuresStruct, queryOptoParameter, queryOptoParameterValues, queryRunFeature, selectedRunGroups, selectedCellsIDs, overallExperimentName, savePlotPath, savePlots)

%% timeFirstAP

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'timeFirstAP';
selectedCellsIDs = unique(cell2mat(selectedRunGroupsAnalysisTable.overallCellName));

generateBarPlotAPpulse(featuresStruct, queryOptoParameter, queryOptoParameterValues, queryRunFeature, selectedRunGroups, selectedCellsIDs, overallExperimentName, savePlotPath, savePlots)

%% timeLastAP

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'timeLastAP';
selectedCellsIDs = unique(cell2mat(selectedRunGroupsAnalysisTable.overallCellName));

generateBarPlotAPpulse(featuresStruct, queryOptoParameter, queryOptoParameterValues, queryRunFeature, selectedRunGroups, selectedCellsIDs, overallExperimentName, savePlotPath, savePlots)

%% nAP line

queryOptoParameter =  'amplitudeCurrentPulse';
queryRunFeature = 'nAP';
selectedCellsIDs = unique(cell2mat(selectedRunGroupsAnalysisTable.overallCellName));

generateLinePlotAPpulse(featuresStruct, queryOptoParameter, queryOptoParameterValues, queryRunFeature, selectedRunGroups, selectedCellsIDs, overallExperimentName, savePlotPath, savePlots)
