function cellCoordinates = extractCellCoordinates(tVectorString)

    extractedCoordinates = extractBetween(tVectorString, '[', ']');
    cellCoordinates = str2double(split(extractedCoordinates, ';'));

end