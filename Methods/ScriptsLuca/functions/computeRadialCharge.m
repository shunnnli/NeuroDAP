function [meanRadialResponses, meanRadialResponsesNorm, radialDistances] = computeRadialCharge(responseMap,nBins,cellX,cellY,responseType)

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
    responseMapVector = reshape(responseMap,dmdWidth*dmdHeight,1);

    meanRadialResponses = accumarray(distanceBinsVector,responseMapVector,[],@mean);

    if strcmp(responseType,'excitatory')
        meanRadialResponsesNorm = meanRadialResponses'/abs(min(meanRadialResponses));
    elseif strcmp(responseType,'inhibitory')
        meanRadialResponsesNorm = meanRadialResponses'/abs(max(meanRadialResponses));
    end
    
    meanRadialResponses = meanRadialResponses';

end