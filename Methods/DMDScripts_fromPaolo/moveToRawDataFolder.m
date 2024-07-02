function moveToRawDataFolder(sourceFolder)

    subfolder = 'RawData';
    matFiles = dir(fullfile(sourceFolder, '*.mat'));

    if ~exist(fullfile(sourceFolder, subfolder), 'dir')

        mkdir(fullfile(sourceFolder, subfolder));

    end

    for i = 1:numel(matFiles)

        sourceFile = fullfile(sourceFolder, matFiles(i).name);
        destinationFile = fullfile(sourceFolder, subfolder, matFiles(i).name);
        movefile(sourceFile, destinationFile);

    end

end
