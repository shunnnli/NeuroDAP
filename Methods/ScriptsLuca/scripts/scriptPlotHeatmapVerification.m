disp(['*-*-*-* Running: scriptPlotHeatmapVerification *-*-*-*'])

mouseAnalysisTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseAnalysisTable.mat'];
mouseParametersTablePath = [experimentDirectory filesep 'mouseAnalysis' filesep 'MouseParametersTable.mat'];

mouseAnalysisTable = load(mouseAnalysisTablePath);
mouseAnalysisTable = mouseAnalysisTable.mouseAnalysisTable;
mouseParametersTable = load(mouseParametersTablePath);
runGroupsParametersTable = mouseParametersTable.runGroupsParametersTable;

analyzedCells = unique(cell2mat(mouseAnalysisTable.cellName));

groupByCell = 0;

if ~any(cellfun(@(x) contains(x,'heatmapVerificationBlue'),runGroupsParametersTable.functionNameBlue)) && ~any(cellfun(@(x) contains(x,'heatmapVerificationRed'),runGroupsParametersTable.functionNameRed)); return; end

rowsProtocolBlue = find(cellfun(@(x) contains(x,'heatmapVerificationBlue'),runGroupsParametersTable.functionNameBlue));
rowsProtocolRed = find(cellfun(@(x) contains(x,'heatmapVerificationRed'),runGroupsParametersTable.functionNameRed));

for iCell = 1:numel(analyzedCells)
    
iAnalyzedCell = analyzedCells(iCell);
cellName = ['cell' num2str(iAnalyzedCell)];

plotNames = {};
rawDataPath = [experimentDirectory filesep cellName filesep 'RawData'];

if ~isempty(rowsProtocolBlue) && ~isempty(rowsProtocolRed); nChannels = 2; channels = {'Blue','Red'};
elseif ~isempty(rowsProtocolBlue); nChannels = 1; channels = {'Blue'};
elseif ~isempty(rowsProtocolRed); nChannels = 1; channels = {'Red'};
else; continue; end

for iChannel = 1:nChannels
    
channel = channels{iChannel};
if strcmp(channel,'Blue'); rowsProtocol = rowsProtocolBlue; 
elseif strcmp(channel,'Red'); rowsProtocol = rowsProtocolRed; end

runGroupsProtocol = cell2mat(runGroupsParametersTable.runGroup(rowsProtocol));
rowsAnalysisTableProtocol = find(ismember(mouseAnalysisTable.runGroup,runGroupsProtocol) & cell2mat(mouseAnalysisTable.cellName) == iAnalyzedCell);

if isempty(rowsAnalysisTableProtocol); continue; end

analysisTableProtocol = mouseAnalysisTable(rowsAnalysisTableProtocol,:);
runGroupsProtocolCell = unique(analysisTableProtocol.runGroup);
rowsRunGroupsParametersTableCell = find(ismember(cell2mat(runGroupsParametersTable.runGroup),runGroupsProtocolCell));
runGroupsParametersTableCell = runGroupsParametersTable(rowsRunGroupsParametersTableCell,:);

if strcmp(channel,'Blue'); protocolPatterns = regexp(runGroupsParametersTableCell.functionNameBlue, 'heatmapVerificationBlue_(.*?)_Epoch\d+_Depth\d+', 'tokens');
elseif strcmp(channel,'Red'); protocolPatterns = regexp(runGroupsParametersTableCell.functionNameRed, 'heatmapVerificationRed_(.*?)_Epoch\d+_Depth\d+', 'tokens'); end

protocolPatterns = cellfun(@(x) x{1}, protocolPatterns, 'UniformOutput', true);
uniqueProtocolPatterns = unique(protocolPatterns);
nDifferentProtocolPatterns = numel(unique(protocolPatterns));
protocolEpochs = unique(cell2mat(analysisTableProtocol.epoch));

runGroupsForEachEpoch = {};
for iEpoch = 1:numel(protocolEpochs)
    iEpochProtocol = protocolEpochs(iEpoch);
    rowsEpoch = find(cell2mat(analysisTableProtocol.epoch) == iEpochProtocol);
    runGroupsEpoch = analysisTableProtocol.runGroup(rowsEpoch);
    runGroupsForEachEpoch{end+1} = runGroupsEpoch;
end

protocolReps = cellfun(@(x) mat2str(x), runGroupsForEachEpoch, 'UniformOutput', false);
[uniqueProtocols, protocolIndices, ~] = unique(protocolReps);
nDifferentProtocols = numel(uniqueProtocols);

for iProtocol = 1:nDifferentProtocols

protocolIndex = protocolIndices(iProtocol);
runGroupsProtocol = runGroupsForEachEpoch{protocolIndex};

rowsSingleProtocol = find(ismember(analysisTableProtocol.runGroup,runGroupsProtocol));
analysisTableUniqueProtocol = analysisTableProtocol(rowsSingleProtocol,:);

if strcmp(channel,'Blue'); correspondingSearchName = regexp(analysisTableUniqueProtocol.runProtocol{1,1}.functionNameBlue, 'Epoch\d+_Depth\d+', 'match'); amplitude = analysisTableUniqueProtocol.runProtocol{1}.amplitudeBlue; pulseWidth = analysisTableUniqueProtocol.runProtocol{1}.pulseWidthBlue;
elseif strcmp(channel,'Red'); correspondingSearchName = regexp(analysisTableUniqueProtocol.runProtocol{1,1}.functionNameRed, 'Epoch\d+_Depth\d+', 'match'); amplitude = analysisTableUniqueProtocol.runProtocol{1}.amplitudeRed; pulseWidth = analysisTableUniqueProtocol.runProtocol{1}.pulseWidthRed; end

[~, idxPulsePeak] = max(abs(cell2mat(analysisTableUniqueProtocol.heightPulsePeak)));
maxHeightPulsePeak = max(cell2mat(analysisTableUniqueProtocol.heightPulsePeak));
minHeightPulsePeak = min(cell2mat(analysisTableUniqueProtocol.heightPulsePeak));

[~, idxAreaPulse] = max(abs(cell2mat(analysisTableUniqueProtocol.heightPulsePeak)));
maxAreaPulse = max(cell2mat(analysisTableUniqueProtocol.areaPulse));
minAreaPulse = min(cell2mat(analysisTableUniqueProtocol.areaPulse));

tracesMatrix = cell2mat(analysisTableUniqueProtocol.trace);
tracesMatrix = tracesMatrix(:,4000:7000);
maxTraceValue = max(tracesMatrix(:));
minTraceValue = min(tracesMatrix(:));

initializeFig(0.9,0.9);

nTileRows = 2;
nTileColumns = numel(runGroupsProtocol)*2;
masterLayout = tiledlayout(nTileRows,nTileColumns);
masterLayout.TileSpacing = 'compact';
masterLayout.Padding = 'compact';

colorBlue = [85, 161, 254]./255;
colorRed = [255, 50, 58]./255;
colorGrey = [192, 192, 192]./255;

colorBlueOverlay = [85, 161, 254]./255 * 0.6;
colorRedOverlay = [255, 50, 58]./255 * 0.6;
colorGreyOverlay = [192, 192, 192]./255 * 0.6;

postTextSuptitle = [correspondingSearchName{1}];
    
for iProtocolPattern = 1:numel(runGroupsProtocol)%nDifferentProtocolPatterns

currentRunGroup = runGroupsProtocol(iProtocolPattern);
rowRunGroup = find(cell2mat(runGroupsParametersTableCell.runGroup) == currentRunGroup); % can and below 4

if strcmp(channel,'Blue'); protocolPattern = regexp(runGroupsParametersTableCell.functionNameBlue{rowRunGroup}, 'heatmapVerificationBlue_(.*?)_Epoch\d+_Depth\d+', 'tokens');
elseif strcmp(channel,'Red'); protocolPattern = regexp(runGroupsParametersTableCell.functionNameRed{rowRunGroup}, 'heatmapVerificationRed_(.*?)_Epoch\d+_Depth\d+', 'tokens'); end

protocolPattern = protocolPattern{1}{1};
%protocolPattern = uniqueProtocolPatterns{iProtocolPattern};
protocolRunGroup = runGroupsProtocol(iProtocolPattern);

% if strcmp(channel,'Blue')
%     cat1_rowsSelectedRunGroups = find(strcmp(runGroupsParametersTableCell.activeChannels,'ao0, ao1') & cell2mat(runGroupsParametersTableCell.holdingVoltage) == -70 & strcmp(runGroupsParametersTableCell.functionNameBlue,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
%     cat1_selectedRunGroups = cell2mat(runGroupsParametersTableCell.runGroup(cat1_rowsSelectedRunGroups));
% elseif strcmp(channel,'Red')
%     cat1_rowsSelectedRunGroups = find(strcmp(runGroupsParametersTableCell.activeChannels,'ao0, ao2') & cell2mat(runGroupsParametersTableCell.holdingVoltage) == 0 & strcmp(runGroupsParametersTableCell.functionNameRed,'dmdOutputRect(1, 1, 1, 1, 1, 1)'));
%     cat1_selectedRunGroups = cell2mat(runGroupsParametersTableCell.runGroup(cat1_rowsSelectedRunGroups));
% end

cat1_selectedRunGroupsAnalysisTable = analysisTableUniqueProtocol(ismember(analysisTableUniqueProtocol.runGroup,protocolRunGroup),:);

%% ** Plots **

%% Traces
mainAxes = nexttile(masterLayout,2*(iProtocolPattern-1)+1,[1 2]); axis off;

plotData = generatePlotDataStruct;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

options = generateOptionsStruct(groupByCell);
options.plotTitle = ['Pattern: ', protocolPattern];
options.yLimits = [minTraceValue, maxTraceValue];

if strcmp(channel,'Blue'); mapName = ['heatmapVerificationBlue_',protocolPattern,'_',correspondingSearchName{1}]; options.plotColor.cat1 = colorBlue; options.plotColorOverlay.cat1 = colorBlueOverlay;
elseif strcmp(channel,'Red'); mapName = ['heatmapVerificationRed_',protocolPattern,'_',correspondingSearchName{1}]; options.plotColor.cat1 = colorRed; options.plotColorOverlay.cat1 = colorRedOverlay; end

generatePlotTraces(plotData,options)

hold on
runEpochs = cell2mat(cat1_selectedRunGroupsAnalysisTable.epoch);
fullMapName = [mapName, '_', 'runEpoch', num2str(runEpochs(1))];
fullMapPath = [rawDataPath filesep fullMapName '.mat'];

if exist(fullMapPath,'file') == 2
    
    displayedMap = load(fullMapPath);
    displayedMap = squeeze(displayedMap.displayedPattern(1,:,:));

    mainPosition = mainAxes.Position;
    
    if strcmp(channel,'Blue'); cornerAxes = axes('Position',[mainPosition(1) - 0.008, mainPosition(2) + 0.0017, 0.05, 0.05]);  
    elseif strcmp(channel,'Red'); cornerAxes = axes('Position',[mainPosition(1) - 0.008, mainPosition(2) + 0.35, 0.05, 0.05]); end
    
    imshow(displayedMap, [], 'Parent', cornerAxes);
    colormap([0.9 0.9 0.9; options.plotColor.cat1])
    hold on;
    
    if strcmp(channel,'Blue')
        cellX = 342+round(cat1_selectedRunGroupsAnalysisTable.runProtocol{1,1}.cellCoordinates.blue(1));
        cellY = 304+round(cat1_selectedRunGroupsAnalysisTable.runProtocol{1,1}.cellCoordinates.blue(2));
    elseif strcmp(channel,'Red')
        cellX = 342+round(cat1_selectedRunGroupsAnalysisTable.runProtocol{1,1}.cellCoordinates.red(1));
        cellY = 304+round(cat1_selectedRunGroupsAnalysisTable.runProtocol{1,1}.cellCoordinates.red(2));
    end
    
    scatter(cellX,cellY,30,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    
end

%% Peak current
nexttile(masterLayout,2*(iProtocolPattern-1)+nTileColumns+1,[1 1]); axis off;
cat1_queryRunFeature = 'heightPulsePeak';

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

options = generateOptionsStruct(groupByCell);
options.yLimits = [minHeightPulsePeak, maxHeightPulsePeak];

if strcmp(channel,'Blue'); options.plotColor.cat1 = colorBlue; options.plotColorOverlay.cat1 = colorBlueOverlay;
elseif strcmp(channel,'Red'); options.plotColor.cat1 = colorRed; options.plotColorOverlay.cat1 = colorRedOverlay; end

options.plotTitle = ['Response amplitude'];
options.xLabelText.cat1 = ['Light'];
options.yLabelText = ['Current [pA]'];

generateBarplot(plotData,options)

%% Charge
nexttile(masterLayout,2*(iProtocolPattern-1)+nTileColumns+2,[1 1]); axis off;
cat1_queryRunFeature = 'areaPulse';

plotData = generatePlotDataStruct;
plotData.queryRunFeature.cat1 = cat1_queryRunFeature;
plotData.selectedRunGroupsAnalysisTable.cat1 = cat1_selectedRunGroupsAnalysisTable;

options = generateOptionsStruct(groupByCell);
options.yLimits = [minAreaPulse, maxAreaPulse];

if strcmp(channel,'Blue'); options.plotColor.cat1 = colorBlue; options.plotColorOverlay.cat1 = colorBlueOverlay;
elseif strcmp(channel,'Red'); options.plotColor.cat1 = colorRed; options.plotColorOverlay.cat1 = colorRedOverlay; end

options.plotTitle = ['Response charge'];
options.xLabelText.cat1 = ['Light'];
options.yLabelText = ['Charge [pC]'];

generateBarplot(plotData,options)

end

protocolText = ['Pulse amplitude = ', num2str(amplitude),', Pulse width = ', num2str(pulseWidth)];
plotName = ['cell', num2str(iAnalyzedCell),'_','HeatmapVerification_' experimentName, '_', postTextSuptitle, '_', 'Protocol1'];
plotMainTitle = ['Cell ', num2str(iAnalyzedCell),' - ',experimentName,' - Heatmap Verification - ', postTextSuptitle, ' - ', protocolText];

sgtitle(plotMainTitle,'FontSize',14,'Interpreter', 'none')
protocolNumber = 1;

while any(strcmp(plotNames, plotName))
    
     protocolNumber = protocolNumber + 1;
     plotName = ['cell', num2str(iAnalyzedCell), 'HeatmapVerification_' experimentName, '_', postTextSuptitle, '_', 'Protocol', num2str(protocolNumber)];
  
end

plotNames{end+1} = plotName;
savePlotPath = [experimentDirectory filesep 'mouseAnalysis'];
plotFullPath = fullfile(savePlotPath, plotName);

if savePlots == 1
    exportgraphics(gcf, fullfile(strcat(plotFullPath,'.pdf')), 'ContentType', 'vector');
    saveas(gcf, plotFullPath, 'fig')
end
       
end

end

end

