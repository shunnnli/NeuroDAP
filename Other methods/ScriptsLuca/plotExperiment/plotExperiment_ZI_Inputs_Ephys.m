% Clearing workspace
clear all;
close all;
clc;

% User environment (automatic)
user = getenv('USER'); 

% Indicate main directory and add it to path
matlabDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/'];
mainDirectory = ['/Users/' user '/Desktop/EPFL/Thesis/Matlab/Experiments/']; % including all Matlab toolboxes and necessary paths
addpath(genpath(matlabDirectory));

% Decide saving options
savePlots = 1;

% Indicate experiments to group
overallExperimentName = 'ZI_Inputs_Ephys';
experimentName = overallExperimentName;
mouseExperimentNames = {};
mouseExperimentNames{1} = 'vGAT_ZI_Ephys_11';
mouseExperimentNames{2} = 'SNr_ZI_vGAT_Ephys_50_1';
mouseExperimentNames{3} = 'MLR_ZI_Ephys_1';
mouseExperimentNames{4} = 'MLR_ZI_Ephys_2';
mouseExperimentNames{5} = 'MLR_ZI_Ephys_3';

% Indicate cells to select
mouseExperimentSelectedCells = {};
mouseExperimentSelectedCells{1} = [1,2,3,4,5,6];
mouseExperimentSelectedCells{2} = [1,2,3,4,5,6,7,8];
mouseExperimentSelectedCells{3} = [1,2,3,4,5];
mouseExperimentSelectedCells{4} = [1,3,4,5,6,7];
mouseExperimentSelectedCells{5} = [1,2,3,4,5,6,7,8,9,10,11,12];

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
    experimentAnalysisTable = findIdenticalRunsVC(experimentAnalysisTable);
    save(saveExperimentTablePath, 'experimentAnalysisTable');

end

runGroups = unique(experimentAnalysisTable.runGroup(:));
nRunGroups = numel(runGroups);

columnNames = {'runGroup', 'sampleSize', 'activeChannels', 'holdingVoltage', 'nPulsesBlue', 'pulseWidthBlue', 'amplitudeBlue', 'delayPulseBlue', 'functionNameBlue', 'nPulsesRed', 'pulseWidthRed', 'amplitudeRed', 'delayPulseRed', 'functionNameRed', 'whichPulseFirst', 'pulsesTimeDifference'};
dataTypes = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
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


%% amplitude vs heightPulsePeak

cat1_queryOptoParameter =  'amplitudeRed';
cat2_queryOptoParameter =  'amplitudeRed';
cat3_queryOptoParameter =  'amplitudeRed';
queryRunFeature = 'areaPulse';

plotBaseName = [experimentName '_' queryRunFeature]; 
queryOptoParameterValuesOverall = [5];

cat1_selectedCellsIDs = [9,10,11,12,27,28,29,30,31,32,33,35,36,37];
cat2_selectedCellsIDs = [7,8,13,14,21,22,23];
cat3_selectedCellsIDs = [15,19,34];

cat1_desiredEpochs = [];
cat2_desiredEpochs = [];
cat3_desiredEpochs = [];

groupByCell = 1;

cat1_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
cat2_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));
cat3_selectedRunGroups = find(strcmp(runGroupsParametersTable.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTable.holdingVoltage) == 0 & cell2mat(runGroupsParametersTable.amplitudeRed) == 5 & strcmp(runGroupsParametersTable.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 2, 1)'));

% No need to modify
cat1_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat1_selectedCellsIDs);
cat1_selectedRunGroupsAnalysisTable = cat1_selectedMembersTable(ismember(cat1_selectedMembersTable.runGroup,cat1_selectedRunGroups),:);
cat1_selectedRunGroupsAnalysisTable = extractEpochSubset(cat1_selectedRunGroupsAnalysisTable, cat1_desiredEpochs);

cat2_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat2_selectedCellsIDs);
cat2_selectedRunGroupsAnalysisTable = cat2_selectedMembersTable(ismember(cat2_selectedMembersTable.runGroup,cat2_selectedRunGroups),:);
cat2_selectedRunGroupsAnalysisTable = extractEpochSubset(cat2_selectedRunGroupsAnalysisTable, cat2_desiredEpochs);

cat3_selectedMembersTable = extractSelectedMembersFromTable(experimentAnalysisTable, 'overallCellName', cat3_selectedCellsIDs);
cat3_selectedRunGroupsAnalysisTable = cat3_selectedMembersTable(ismember(cat3_selectedMembersTable.runGroup,cat3_selectedRunGroups),:);
cat3_selectedRunGroupsAnalysisTable = extractEpochSubset(cat3_selectedRunGroupsAnalysisTable, cat3_desiredEpochs);

generateBarPlotVCwith3categories(cat1_queryOptoParameter, cat2_queryOptoParameter, cat3_queryOptoParameter, queryRunFeature, queryOptoParameterValuesOverall, cat1_selectedRunGroupsAnalysisTable, cat2_selectedRunGroupsAnalysisTable, cat3_selectedRunGroupsAnalysisTable, cat1_selectedRunGroups, cat2_selectedRunGroups, cat3_selectedRunGroups, cat1_selectedCellsIDs, cat2_selectedCellsIDs, cat3_selectedCellsIDs, groupByCell, plotBaseName, savePlotPath, savePlots)
