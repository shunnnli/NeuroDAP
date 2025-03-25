function plotHeatmap(data,timestamp,options)

% Plot heatmap of data

arguments
    data double
    timestamp double
    
    options.dataType string = 'trials'
    options.flipYAxis logical = false

    options.plotColorBar logical = true
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
    options.colormap = blueWhiteRed;
end

% Define color bar label
if ~isfield(options,'colorBarLabel')
    options.colorBarLabel = "Value";
end

if ~isfield(options,'colorlim')
    climit = max(abs(data),[],"all");
    options.colorlim = [-climit, climit];
end

%% Plot heatmap
if strcmpi(options.dataType,'trials')
    imagesc(timestamp,1:size(data,1),data);
    if ~options.flipYAxis
        set(gca,'YDir','normal');
    end
    colorbar; box off
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
    % Create an alpha mask: 1 for non-NaN, 0 for NaN
    alphaMask = ~isnan(data);
    % Apply the alpha mask to the image
    set(h, 'AlphaData', alphaMask);
end

end