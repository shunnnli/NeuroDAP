function hotspotDistanceStats = computeHotspotDistanceStats(hotspotMapDepth,iSearchDepth,cellX,cellY)

    mapSizeX = 608;
    mapSizeY = 684;

    spotSizeX = mapSizeX/(2^iSearchDepth);
    spotSizeY = mapSizeY/(2^iSearchDepth);
    firstSpotX = spotSizeX/2;
    firstSpotY = spotSizeY/2;
    allSpotsX = firstSpotX:spotSizeX:mapSizeX;
    allSpotsY = firstSpotY:spotSizeY:mapSizeY;

    [X, Y] = meshgrid(allSpotsX, allSpotsY);
    allSpotsXY = round([X(:), Y(:)]);

    responsesSpots = hotspotMapDepth(sub2ind(size(hotspotMapDepth), allSpotsXY(:,2), allSpotsXY(:,1)));
    hotspotsXY = allSpotsXY(logical(responsesSpots),:);
    hotspotsDistance = round(sqrt((hotspotsXY(:,1)-cellX).^2 + (hotspotsXY(:,2)-cellY).^2));
    hotspotNumber = numel(hotspotsDistance);

    hotspotDistanceMean = mean(hotspotsDistance);
    hotspotDistanceStd = std(hotspotsDistance);
    hotspotDistanceMedian = median(hotspotsDistance);
    hotspotDistanceMode = mode(hotspotsDistance);
    hotspotDistancePrctile = prctile(hotspotsDistance,[25,50,75]);

    hotspotDistanceStats = [];
    hotspotDistanceStats.hotspotNumber = hotspotNumber;
    hotspotDistanceStats.hotspotDistances = hotspotsDistance;
    hotspotDistanceStats.hotspotDistanceMean = hotspotDistanceMean;
    hotspotDistanceStats.hotspotDistanceStd = hotspotDistanceStd;
    hotspotDistanceStats.hotspotDistanceMedian = hotspotDistanceMedian;
    hotspotDistanceStats.hotspotDistanceMode = hotspotDistanceMode;
    hotspotDistanceStats.hotspotDistancePrctile = hotspotDistancePrctile;

end
