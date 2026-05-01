function plotHeatmap(data,timestamp,options)
% Plot heatmap of data

arguments
    data double
    timestamp double
    options.dataType string = "trials"
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
if ~isfield(options,"colormap")
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
if ~isfield(options,"colorBarLabel")
    options.colorBarLabel = "Value";
end

% Define / sanitize color limits
if ~isfield(options,"colorlim")
    options.colorlim = localAutoColorLimits(data, options.centerColorMap);
else
    options.colorlim = localValidateColorLimits(options.colorlim, data, options.centerColorMap);
end

%% Plot heatmap

if strcmpi(options.dataType,"trials")
    h = imagesc(timestamp, 1:size(data,1), data);

    if ~options.flipYAxis
        set(gca,"YDir","normal");
    end

    if ~isempty(options.colorlim)
        clim(gca, options.colorlim);
    end

    colormap(gca, options.colormap);

    % Make NaN / Inf transparent
    set(h, "AlphaData", isfinite(data));

    if options.plotColorBar
        cb = colorbar;
        cb.Label.String = options.colorBarLabel;
    end

    box off
    return;

elseif strcmpi(options.dataType,"heatmap")

    if isfield(options,"tile")
        h = imagesc(options.tile, data); hold on; box off;
    else
        h = imagesc(data); hold on; box off;
    end

    % Plot color bar
    if options.plotColorBar
        if isfield(options,"tile")
            cb = colorbar(options.tile);
            cb.Label.String = options.colorBarLabel;

            cb.Units = "normalized";
            options.tile.Units = "normalized";
            tilePos = options.tile.Position;

            if isfield(options,"colorBarPosition")
                cbWidth  = tilePos(3) * options.colorBarPosition(3);
                cbHeight = tilePos(4) * options.colorBarPosition(4);
                newX     = tilePos(1) + tilePos(3) * options.colorBarPosition(1);
                newY     = tilePos(2) + tilePos(4) * options.colorBarPosition(2);
                cb.Position = [newX, newY, cbWidth, cbHeight];
            end
        else
            cb = colorbar;
            cb.Label.String = options.colorBarLabel;
            cb.Location = "eastoutside";
        end
    end

    clim(options.colorlim);
    colormap(options.colormap);

    % Set NaN / Inf to transparent
    set(h, "AlphaData", isfinite(data));
end

end

function climits = localAutoColorLimits(data, centerColorMap)
    finiteVals = data(isfinite(data));

    % Fallback for empty / all-NaN / all-Inf
    if isempty(finiteVals)
        if centerColorMap
            climits = [-1, 1];
        else
            climits = [0, 1];
        end
        return;
    end

    if centerColorMap
        maxAbsVal = max(abs(finiteVals), [], "all");

        % Avoid [-0 0]
        if maxAbsVal == 0
            climits = [-1, 1];
        else
            climits = [-maxAbsVal, maxAbsVal];
        end
    else
        lo = min(finiteVals, [], "all");
        hi = max(finiteVals, [], "all");

        % Avoid [x x]
        if hi <= lo
            pad = max([1, abs(lo), abs(hi)]) * 1e-12;
            climits = [lo - pad, hi + pad];
        else
            climits = [lo, hi];
        end
    end
end

function climits = localValidateColorLimits(climitsIn, data, centerColorMap)
    useAuto = isempty(climitsIn) || ...
              ~isnumeric(climitsIn) || ...
              numel(climitsIn) ~= 2 || ...
              any(isnan(climitsIn)) || ...
              climitsIn(2) <= climitsIn(1);

    if useAuto
        climits = localAutoColorLimits(data, centerColorMap);
    else
        climits = reshape(double(climitsIn), 1, 2);
    end
end