function [featuresStruct, queryOptoParameterValues, queryOptoParameterArray] = extractFeaturesPulseAP(queryOptoParameter, selectedRunGroups, selectedRunGroupsAnalysisTable)

    queryOptoParameterValues = zeros(numel(selectedRunGroups),1);
    queryOptoParameterArray = [];
    featuresStruct = [];
    featuresStruct.nAPallData = {};
    featuresStruct.timesAPallData = {};
    featuresStruct.timeFirstAPallData = {};
    featuresStruct.timeLastAPallData = {};
    featuresStruct.peakVoltageAPallData = {};
    featuresStruct.thresholdAPallData = {};
    featuresStruct.cellNameData= {};

    groupsCounter = 1;

    for iRunGroup = 1:numel(selectedRunGroups)

        rowsRunGroup = find(selectedRunGroupsAnalysisTable.runGroup(:) == selectedRunGroups(iRunGroup));
        queryOptoParameterValues(groupsCounter) = selectedRunGroupsAnalysisTable.optoParameters{rowsRunGroup(1)}.(queryOptoParameter);

        queryOptoParameterArray = [queryOptoParameterArray; repmat(queryOptoParameterValues(groupsCounter),numel(rowsRunGroup),1)];
        pulseAProwsGroup = selectedRunGroupsAnalysisTable.pulseAP(rowsRunGroup);

        nAProwsGroup = {};
        timesAProwsGroup = {};
        timeFirstAProwsGroup = {};
        timeLastAProwsGroup = {};
        peakVoltageAProwsGroup = {};
        thresholdAProwsGroup = {};

        for iRun = 1:numel(rowsRunGroup)

          if ~isempty(pulseAProwsGroup{iRun})

              peakVoltageAP = pulseAProwsGroup{iRun}.AP_peak_V;
              indicesCorrectPeaks = find(peakVoltageAP>-15);

              nAP = numel(indicesCorrectPeaks);
              timesAP = pulseAProwsGroup{iRun}.AP_peak_time(indicesCorrectPeaks);
              timeFirstAP = timesAP(1);
              timeLastAP = timesAP(end);
              peakVoltageAP = pulseAProwsGroup{iRun}.AP_peak_V(indicesCorrectPeaks);
              thresholdAP = pulseAProwsGroup{iRun}.AP_thresh_V(indicesCorrectPeaks);

          else

              nAP = 0;
              timesAP = nan;
              timeFirstAP = nan;
              timeLastAP = nan;
              peakVoltageAP = nan;
              thresholdAP = nan;

          end

          nAProwsGroup{iRun} = nAP;
          timesAProwsGroup{iRun} = timesAP;
          timeFirstAProwsGroup{iRun} = timeFirstAP;
          timeLastAProwsGroup{iRun} = timeLastAP;
          peakVoltageAProwsGroup{iRun} = peakVoltageAP;
          thresholdAProwsGroup{iRun} = thresholdAP;

        end

        featuresStruct.nAPallData{groupsCounter} = nAProwsGroup;
        featuresStruct.timesAPallData{groupsCounter} = timesAProwsGroup;
        featuresStruct.timeFirstAPallData{groupsCounter} = timeFirstAProwsGroup;
        featuresStruct.timeLastAPallData{groupsCounter} = timeLastAProwsGroup;
        featuresStruct.peakVoltageAPallData{groupsCounter} = peakVoltageAProwsGroup;
        featuresStruct.thresholdAPallData{groupsCounter} = thresholdAProwsGroup;
        featuresStruct.cellNameData{groupsCounter} = cell2mat(selectedRunGroupsAnalysisTable.overallCellName(rowsRunGroup))';

        groupsCounter = groupsCounter + 1;

    end
    
end