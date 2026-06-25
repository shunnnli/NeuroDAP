function infoTable = readInfoPatchingTable(xlsxPath, maxVariables)
%READINFOPATCHINGTABLE Read the acquisition log from InfoPatching.xlsx.
% The workbook stores its table header on row 25, so the range must be set
% while detecting import options before selecting numbered variables.

arguments
    xlsxPath
    maxVariables (1,1) double {mustBePositive, mustBeInteger}
end

if exist(xlsxPath, 'file') ~= 2
    error("readInfoPatchingTable:FileNotFound", ...
          "InfoPatching.xlsx was not found: %s", char(xlsxPath));
end

try
    opts = detectImportOptions(xlsxPath, "VariableNamesRange", 25);
catch
    opts = detectImportOptions(xlsxPath, "Range", "A25");
end
nVariables = numel(opts.VariableNames);

if nVariables == 0
    error("readInfoPatchingTable:NoVariables", ...
          "No variables were detected in %s using row 25 as the header.", char(xlsxPath));
end

nSelected = min(maxVariables, nVariables);
opts.SelectedVariableNames = opts.VariableNames(1:nSelected);

infoTable = readtable(xlsxPath, opts);

requiredVariables = ["acq_", "epoch", "cyclePos", "holding"];
missingVariables = setdiff(requiredVariables, string(infoTable.Properties.VariableNames));
if ~isempty(missingVariables)
    error("readInfoPatchingTable:MissingVariables", ...
          "InfoPatching.xlsx is missing expected variable(s): %s", ...
          char(strjoin(missingVariables, ", ")));
end

end
