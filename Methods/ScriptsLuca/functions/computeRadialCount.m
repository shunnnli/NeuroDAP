function [meanRadialResponses, meanRadialResponsesNormalized, radialDistances] = computeRadialCount(isResponseMap,nBins,cellX,cellY)

    dmdWidth = 608;
    dmdHeight = 684;
    dmdX = 1:dmdWidth;
    dmdY = 1:dmdHeight;
    [X,Y] = meshgrid(dmdX,dmdY);

    X = X - cellX;
    Y = Y - cellY;
    spotDistance = sqrt(X.^2+Y.^2);
    maxDistance = max(spotDistance,[],'all');

    binWidth = maxDistance/(nBins);
    radialDistances = (0:binWidth:maxDistance);
    distanceBins = round(spotDistance/binWidth)+1;

    distanceBinsVector = reshape(distanceBins,dmdWidth*dmdHeight,1);
    responseMapVector = reshape(isResponseMap,dmdWidth*dmdHeight,1);

    radialSurfaces = accumarray(distanceBinsVector, 1);
    meanRadialResponses = accumarray(distanceBinsVector,responseMapVector);
    meanRadialResponses = (meanRadialResponses');
    meanRadialResponsesNormalized = (meanRadialResponses./radialSurfaces');

end