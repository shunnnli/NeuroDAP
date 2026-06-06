function [toBeExcluded, holdingVoltage, specificDrugs, qualityRS0, scopeMagnification] = findRunInfoExcel(excelPath, runEpoch, runCyclePosition)

    epochOfInterest = runEpoch;
    cyclePositionOfInterest = runCyclePosition;

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readtable(excelPath);

    toBeExcluded = findExcludedRunsExcel(data, epochOfInterest, cyclePositionOfInterest);
    holdingVoltage = findHoldingVoltageExcel(data, epochOfInterest, cyclePositionOfInterest);
    specificDrugs = findSpecificDrugsExcel(data, epochOfInterest, cyclePositionOfInterest);
    qualityRS0 = findQualityRS0Excel(data, epochOfInterest, cyclePositionOfInterest);
    scopeMagnification = findScopeMagnificationExcel(data, epochOfInterest, cyclePositionOfInterest);

end
