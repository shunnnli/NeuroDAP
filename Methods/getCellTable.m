function allCells = getCellTable(epochs,options)

% Combine cell-specific statistics from epochs

arguments
    epochs table
    options.save logical = true
    options.saveDataPath
end

%% cells table

% Each row is a cell
% Vhold: every epochs of this cell: nEpochs x 1
% Included: all the included sweeps within each epochs: cell, nEpochs x nSweeps
% Sweep names: matlab names for the sweep: cell, nEpochs x nSweeps
% Raw/processed sweeps: cell, each row contains raw sweeps for this epoch
% Peaks/AUCs: nEpochs x nSweeps
% Rs/Rm/Cm: nEpochs x nSweeps

% If epochs contains only one animal, create cells for that animal
% If epochs contains multiple animals:
    % 1. Create cells for each animal, save them into there corresponding
    % folder
    % 2. Concatenate into allCells and return

%% Create cells table by animal

animalList = unique(epochs{:,"Animal"});
allCells = [];

for i = 1:length(animalList)
    % Find epochs for an animal
    animalEpochs = epochs(epochs.Animal == animalList(i),:);
    cellList = unique(animalEpochs{:,"Cell"});
    dirsplit = split(animalEpochs{1,"Session"},filesep); expName = dirsplit{end};

    % Initialize cells table
    varTypes = {'string','string','string','double',...
                'cell','cell','cell',...
                'cell','cell',...
                'cell','cell',...
                'cell','cell','cell'};
    varNames = {'Session','Animal','Task','Cell',...
                'Vhold','Included','Sweep names',...
                'Raw sweeps','Processed sweeps',...
                'Peaks','AUCs',...
                'Rs','Rm','Cm'};
    cells = table('Size',[length(cellList),length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);

    for c = 1:size(cellList,1)
        cellEpochs = animalEpochs(animalEpochs.Cell == cellList(c),:);

        cells{c,'Session'} = cellEpochs{1,'Session'};
        cells{c,'Animal'} = cellEpochs{1,'Animal'};
        cells{c,'Task'} = cellEpochs{1,'Task'};
        cells{c,'Cell'} = cellList(c);
        cells{c,'Vhold'} = num2cell(cellEpochs.("Vhold epoch mean"),[1 2]); % Save Vhold
        cells{c,'Included'} = {cellEpochs.('Included')};
        cells{c,'Sweep names'} = {cellEpochs.('Sweep names')};
        cells{c,'Raw sweeps'} = {cellEpochs.('Raw sweeps')};
        cells{c,'Processed sweeps'} = {cellEpochs.('Processed sweeps')};
        cells{c,'Peaks'} = {cellEpochs.('Peaks')};
        cells{c,'AUCs'} = {cellEpochs.('AUCs')};
        cells{c,'Rs'} = {cellEpochs.('Rs')};
        cells{c,'Rm'} = {cellEpochs.('Rm')};
        cells{c,'Cm'} = {cellEpochs.('Cm')};
    end

    if options.save
        save(strcat(animalEpochs{1,'Session'},filesep,'cells_',expName),'cells');
        disp(strcat("Created cells.mat: ",expName));
    end
    allCells = [allCells; cells];
end

end
