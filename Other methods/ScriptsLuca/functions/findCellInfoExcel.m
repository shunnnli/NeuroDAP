function [generalDrugs, cellLocation, cellTemperature, cellType, cellFluo] = findCellInfoExcel(excelPath)

    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    data = readcell(excelPath);

    generalDrugs = findGeneralDrugsExcel(data);
    cellLocation = findCellLocationExcel(data);
    cellTemperature = findCellTemperatureExcel(data);
    cellType = findCellTypeExcel(data);
    cellFluo = findCellFluoExcel(data);

end