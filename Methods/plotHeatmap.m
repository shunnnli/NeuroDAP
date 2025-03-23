function plotHeatmap(data,timestamp,options)

arguments
    data double
    timestamp double
    
    options.dataType string = 'trials'
    options.flipYAxis logical = false

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

%% Plot heatmap
if strcmpi(options.dataType,'trials')
    imagesc(timestamp,1:size(data,1),data);
    if ~options.flipYAxis
        set(gca,'YDir','normal');
    end
    colorbar; box off
    return;
elseif strcmpi(options.dataType,'heatmap')
    h = imagesc(data); hold on; box off; 
    
    cb = colorbar; 
    cb.Label.String = options.colorBarLabel;

    % Determine cb location
    if isfield(options,'tile')
        cb.Units = 'normalized';
        tilePos = options.tile.Position; 
        if isfield(options,'colorBarPosition')
            cbWidth  = tilePos(3) * options.colorBarPosition(3);
            cbHeight = tilePos(4) * options.colorBarPosition(4);
            newX = tilePos(1) + tilePos(3) * options.colorBarPosition(1);
            newY = tilePos(2) + tilePos(4) * options.colorBarPosition(2);
            cb.Position = [newX, newY, cbWidth, cbHeight];
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