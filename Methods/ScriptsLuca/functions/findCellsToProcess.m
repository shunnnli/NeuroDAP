function cellsToProcess = findCellsToProcess(experimentPath)

    items = dir(experimentPath);
    cellsToProcess = [];

    for i = 1:length(items)

        if items(i).isdir

            tokens = regexp(items(i).name, '^cell(\d+)$', 'tokens');

            if ~isempty(tokens)
                
                cellPath = fullfile(experimentPath, items(i).name);
                matFiles = dir(fullfile(cellPath, '*.mat'));

                subfolderContents = dir(cellPath);
                isSubfolderPresent = any([subfolderContents.isdir] & ~ismember({subfolderContents.name}, {'.', '..'}));

                if ~isempty(matFiles) || isSubfolderPresent
                    
                    cellsToProcess(end+1) = str2double(tokens{1}{1});
                
                end

            end

        end

    end

end