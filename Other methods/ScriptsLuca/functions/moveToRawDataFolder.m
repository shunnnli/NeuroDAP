function moveToRawDataFolder(sourceFolder)

    subfolder = 'RawData';
    matFiles = dir(fullfile(sourceFolder, '*.mat'));
    figFiles = dir(fullfile(sourceFolder, '*.fig'));    
    allFiles = [matFiles; figFiles];

    if ~exist(fullfile(sourceFolder, subfolder), 'dir')

        mkdir(fullfile(sourceFolder, subfolder));

    end

    parfor i = 1:numel(allFiles)

        sourceFile = fullfile(sourceFolder, allFiles(i).name);
        destinationFile = fullfile(sourceFolder, subfolder, allFiles(i).name);
        movefile(sourceFile, destinationFile);

    end

end
