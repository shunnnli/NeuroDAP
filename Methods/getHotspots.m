function cells = getHotspots(cells,options)

% Calculate a new hotspot map for the given experiment/cell/search 
% using the entered criteria

% If searchIdx is not provided, returns a modified cells table
% else return the hotspotMap for that search

% If baselineType is 'noise', use the noise model of the neuron
% If baselineType is 'online', use the std of the neuron during that sweep

arguments
    cells table

    options.searchIdx double
    options.getMap logical = true % is not

    options.analysisWindowLength double = 50 % in ms
    options.thresholdFactor double = 3
    options.baselineType string = 'noise' % can also be 'online'
end

%% Setup

%% Calcualte hotspots
cellList = unique(cells.Cell);

for cellIdx = 1:size(cells,1)
    % Load cell info
    c = cellList(cellIdx);
    curCell = cells(cells.Cell == c,:);
    searchPerCell = length(curCell.Vhold{1});
    noiseSD = curCell.Stats{1}.noiseSD;
    if isfield(options,'searchIdx'); searchPerCell = options.searchIdx;
    else; searchPerCell = 1:searchPerCell; end

    % Define analysisWindow (0~50)
    maxAnalysisWindowLenth = curCell.Options{1}.spotOptions.timeRange(2);
    if options.analysisWindowLength <= maxAnalysisWindowLenth
        spotAnalysisWindowLength = curCell.Options{1}.spotOptions.analysisWindowLength;
        if spotAnalysisWindowLength == options.analysisWindowLength
            analysisWindow = curCell.Options{1}.spotOptions.analysisWindow;
        else
            outputFs = curCell.Options{1}.outputFs;
            analysisWindowStart = curCell.Options{1}.spotOptions.analysisWindow(1);
            analysisWindowSamples = options.analysisWindowLength * outputFs/1000;
            analysisWindow = analysisWindowStart : analysisWindowStart+analysisWindowSamples;
        end
    else
        warning('getHotspots: analysisWindowLength is longer than timeRange, change to 50ms by default!');
    end

    % Calculate hotspot map for given search
    for searchIdx = searchPerCell
        % Load search info
        curSearch = curCell.Epochs{1}{searchIdx};
        search_depths = curCell.("Response map"){1}.depths{searchIdx};
        search_cmap = curCell.("Response map"){1}.currentMap{searchIdx};
        search_spotLocation = curCell.("Response map"){1}.spotLocation{searchIdx};
        nDepth = length(search_depths);
        search_hotspotMap = zeros(684,608,nDepth);
        
        % Loop through each depth
        for d = 1:nDepth
            curDepth = search_depths(d);
            depthCurrentMap = search_cmap{d};
            depthBaselineSD = curCell.Stats{1}.baseline.std{searchIdx}{d};
            depthHotspotMap = zeros(684,608);
            depthSpotLocation = search_spotLocation{d};
            disp(['Ongoing: calculating hotspot map for search: ', curSearch,' at depth ',num2str(curDepth)]);
            
            % Loop through all spots
            for spot = 1:size(depthCurrentMap,1)
                % Modify hotspot in cells table
                if strcmp(options.baselineType,'noise'); baselineSD = noiseSD;
                else; baselineSD = depthBaselineSD{spot}; end
                response_threshold = baselineSD * options.thresholdFactor;
                spotHotspot = sum(abs(depthCurrentMap{spot}(:,analysisWindow))>=response_threshold,2)>=50;
                cells{cells.Cell == c,'Response map'}{1}.hotspot{searchIdx}{d}{spot} = spotHotspot;

                % Create a new hotspot map
                location = depthSpotLocation(spot,:);
                yRange = location(1):location(2);
                xRange = location(3):location(4);
                % Add to hotspot map
                depthHotspotMap(xRange,yRange) = sum(spotHotspot)>=1;
            end 
            search_hotspotMap(:,:,d) = depthHotspotMap;
        end

        % Save hotspot map
        cells{cells.Cell == c,'Response map'}{1}.hotspotMap{searchIdx} = search_hotspotMap;
    end
end

end