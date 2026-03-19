function fig = initializeFig(horizontalPortion, verticalPortion, options)
    %INITIALIZEFIG Initialize a figure at a *reference* pixel size.
    %
    % This function computes a desired figure size in pixels using
    % options.ScreenSize * portions. When the current monitor is smaller than
    % the desired size, many OS/window managers will clamp the *window* size
    % even if the figure is invisible. To make saving reliable, we:
    %   - create the figure invisible and fully OFF-SCREEN when too big
    %   - store the desired pixel size on the figure (appdata: 'TargetPixelSize')
    % so export can use the target size regardless of window clamping.
    
    arguments
        horizontalPortion (1,1) double {mustBeGreaterThan(horizontalPortion,0), mustBeLessThanOrEqual(horizontalPortion,1)} = 0.5
        verticalPortion   (1,1) double {mustBeGreaterThan(verticalPortion,0),   mustBeLessThanOrEqual(verticalPortion,1)}   = 0.5
    
        options.ScreenSize (1,4) double = [0 0 3440 1440]   % reference [left bottom width height] (pixels)
        options.Margin (1,1) double = 0                     % pixels
    
        options.FontSize (1,1) double = 20
        options.AxesLineWidth (1,1) double = 1.5
    
        % 'auto' => Visible='off' if too big for current monitor
        options.Visible {mustBeMember(options.Visible, {'on','off','auto'})} = 'auto'
    end
    
    % --- Desired size from reference ---
    ref = options.ScreenSize;
    if isempty(ref)
        ref = get(0,'ScreenSize');
    end
    
    desiredW = round(ref(3) * horizontalPortion);
    desiredH = round(ref(4) * verticalPortion);
    
    % --- Monitor geometry (for placement/check) ---
    monitors = get(0,'MonitorPositions');
    mon = monitors(1,:); % primary
    
    isTooBig = (desiredW > mon(3)) || (desiredH > mon(4));
    
    % --- Decide visibility ---
    if strcmp(options.Visible,'auto')
        if isTooBig
            finalVisible = 'off';
        else
            finalVisible = 'on';
        end
    else
        finalVisible = options.Visible;
    end
    
    % --- Position ---
    if strcmp(finalVisible,'on')
        monLeft   = mon(1);
        monBottom = mon(2);
        monTop    = mon(2) + mon(4);
    
        left   = monLeft   + round((mon(3) - desiredW)/2);
        bottom = monBottom + round((mon(4) - desiredH)/2);
    
        if desiredW > mon(3)
            left = monLeft + options.Margin;
        end
        if desiredH > mon(4)
            bottom = (monTop - options.Margin) - desiredH;
        end
    
        pos = [left bottom desiredW desiredH];
    else
        % Place fully off-screen (right of all monitors) to avoid negative coords
        maxRight  = max(monitors(:,1) + monitors(:,3));
        minBottom = min(monitors(:,2));
        pos = [maxRight + 50, minBottom + 50, desiredW, desiredH];
    end
    
    % --- Create figure ---
    fig = figure('Color','w', 'Units','pixels', 'Position',pos, ...
        'Visible',finalVisible, 'Resize','off');
    
    % Force undocked/normal window when possible (docked figures often clamp)
    try
        fig.WindowStyle = 'normal';
    catch
    end
    try
        fig.WindowState = 'normal';
    catch
    end
    
    % Store target pixels for saving
    try
        setappdata(fig,'TargetPixelSize',[desiredW desiredH]);
    catch
    end
    
    % Try to enforce size after creation (may still be clamped; saving will not)
    drawnow('nocallbacks');
    try
        fig.Units = 'pixels';
        p = fig.Position;
        p(3:4) = [desiredW desiredH];
        fig.Position = p;
    catch
    end
    
    % Default axes styling
    hold on
    ax = gca;
    set(ax, 'FontSize', options.FontSize, 'LineWidth', options.AxesLineWidth);
    box(ax, 'off');
    
    % For axes layout consistency across versions
    if isprop(ax,'PositionConstraint')
        ax.PositionConstraint = 'innerposition';
    elseif isprop(ax,'ActivePositionProperty')
        ax.ActivePositionProperty = 'position';
    end

end
