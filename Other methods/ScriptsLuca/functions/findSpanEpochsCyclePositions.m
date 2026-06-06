function  [spanEpochs, spanCyclePositions] = findSpanEpochsCyclePositions(path)

    files = dir(fullfile(path, 'Epoch*_cyclePosition*.mat'));

    smallestEpoch = inf;
    largestEpoch = -inf;
    smallestCyclePosition = inf;
    largestCyclePosition = -inf;

    for i = 1:numel(files)

        fileName = files(i).name;
        numbers = regexp(fileName, 'Epoch(\d+)_cyclePosition(\d+)', 'tokens');
        numbers = str2double([numbers{:}]);

        if ~isempty(numbers)

            smallestEpoch = min(smallestEpoch, numbers(1));
            largestEpoch = max(largestEpoch, numbers(1));
            smallestCyclePosition = min(smallestCyclePosition, numbers(2));
            largestCyclePosition = max(largestCyclePosition, numbers(2));

        end

    end
    
    spanEpochs = [smallestEpoch:largestEpoch];
    spanCyclePositions = [smallestCyclePosition:largestCyclePosition];

end
