function  [spanEpochs, spanSearchDepths, spanSearchRepetitions] = findSpanEpochsRandomSearch(path)

    files = dir(fullfile(path, 'Epoch*_searchDepth*_searchRepetition*.mat'));

    smallestEpoch = inf;
    largestEpoch = -inf;
    smallestDepth = inf;
    largestDepth = -inf;
    smallestRepetition = inf;
    largestRepetition = -inf;

    for i = 1:numel(files)

        fileName = files(i).name;
        numbers = regexp(fileName, 'Epoch(\d+)_searchDepth(\d+)_searchRepetition(\d+)', 'tokens');
        numbers = str2double([numbers{:}]);

        if ~isempty(numbers)

            smallestEpoch = min(smallestEpoch, numbers(1));
            largestEpoch = max(largestEpoch, numbers(1));
            smallestDepth = min(smallestDepth, numbers(2));
            largestDepth = max(largestDepth, numbers(2));
            smallestRepetition = min(smallestRepetition, numbers(3));
            largestRepetition = max(largestRepetition, numbers(3));

        end

    end
    
    spanEpochs = [smallestEpoch:largestEpoch];
    spanSearchDepths = [smallestDepth:largestDepth];
    spanSearchRepetitions = [smallestRepetition:largestRepetition];

end