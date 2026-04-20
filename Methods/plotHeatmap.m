function plotHeatmap(data,timestamp,options)

% Plot heatmap of data

arguments
    data double
    timestamp double
    
    options.dataType string = 'trials'
    options.flipYAxis logical = false

    options.plotColorBar logical = true
    options.centerColorMap logical = true

    options.colorlim double
    options.colorBarLabel string
    options.colormap
    
    options.tile
    options.colorBarPosition % left bottom width height
end

%% Check input

% Define colormap
if ~isfield(options,'colormap')
    [~,~,~,~,blueWhiteRed,~,~] = loadColors;

    if options.centerColorMap
        options.colormap = blueWhiteRed;
    else
        % Use only the blue half by default, reversed so:
        % small values -> white, large values -> blue
        nHalf = ceil(size(blueWhiteRed,1)/2);
        blueHalf = blueWhiteRed(1:nHalf,:);
        options.colormap = flipud(blueHalf);
    end
end

% Define color bar label
if ~isfield(options,'colorBarLabel')
    options.colorBarLabel = "Value";
end

% Define color limits
if ~isfield(options,'colorlim')
    if options.centerColorMap
        climit = max(abs(data),[],"all");
        options.colorlim = [-climit, climit];
    else
        options.colorlim = [min(data,[],"all"), max(data,[],"all")];
    end
end

%% Plot heatmap
if strcmpi(options.dataType,'trials')
    h = imagesc(timestamp, 1:size(data,1), data);

    if ~options.flipYAxis
        set(gca,'YDir','normal');
    end

    if ~isempty(options.colorlim)
        clim(gca, options.colorlim);
    end
    colormap(gca, options.colormap);

    if options.plotColorBar
        cb = colorbar;
        cb.Label.String = options.colorBarLabel;
    end

    box off
    return;
    
elseif strcmpi(options.dataType,'heatmap')
    if isfield(options,'tile') 
        h = imagesc(options.tile,data); hold on; box off; 
    else
        h = imagesc(data); hold on; box off; 
    end

    % Plot color bar
    if options.plotColorBar
        if isfield(options,'tile')
            cb = colorbar(options.tile); 
            cb.Label.String = options.colorBarLabel;
            cb.Units = 'normalized';

            options.tile.Units = 'normalized';
            tilePos = options.tile.Position;

            if isfield(options,'colorBarPosition')
                cbWidth  = tilePos(3) * options.colorBarPosition(3);
                cbHeight = tilePos(4) * options.colorBarPosition(4);
                newX = tilePos(1) + tilePos(3) * options.colorBarPosition(1);
                newY = tilePos(2) + tilePos(4) * options.colorBarPosition(2);
                cb.Position = [newX, newY, cbWidth, cbHeight];
            end
        else
            cb = colorbar; 
            cb.Label.String = options.colorBarLabel;
            cb.Location = 'eastoutside';
        end
    end

    clim(options.colorlim);
    colormap(options.colormap);

    % Set NaN to transparent
    alphaMask = ~isnan(data);
    set(h, 'AlphaData', alphaMask);
end

end