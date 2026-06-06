[nEpochs, nCyclePositions] = findNumberEpochsCyclePositions('/Users/LucaLiebi/Desktop/Thesis/Matlab/Experiments/vGAT_ZI_Ephys_11_240229/cell1/ProcessedData/Data')

function [nEpochs, nCyclePositions] = findNumberEpochsCyclePositions(path)

    nEpochs = 1;
    
    while true
            
        filename = [path filesep 'Epoch' num2str(nEpochs) '_cyclePosition*.mat'];
        files = dir(filename);
        
        if isempty(files)
            
            nEpochs = nEpochs - 1;
            
            break;
        
        end
        
        nEpochs = nEpochs + 1;
    
    end
    
    nCyclePositions = 1;

    while true

        filename = [path filesep 'Epoch*' '_cyclePosition' num2str(nCyclePositions) '.mat'];
        files = dir(filename);

        if isempty(files)

            nCyclePositions = nCyclePositions - 1;

            break;

        end

        nCyclePositions = nCyclePositions + 1;

    end

    disp(['nEpochs = ', num2str(nEpochs)]);
    disp(['nCyclePositions = ', num2str(nCyclePositions)]);

end