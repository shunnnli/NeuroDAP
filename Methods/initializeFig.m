function fig = initializeFig(horizontalPortion, verticalPortion, options)
    %INITIALIZEFIG Create a figure sized as a portion of a *reference* screen.
    %
    % Goal: keep a consistent figure pixel size (based on a reference ScreenSize)
    % even when you run on a smaller monitor.
    %
    % AUTO-VISIBILITY:
    % If the calculated figure size is larger than the current monitor, 
    % the figure is created INVISIBLY ('Visible', 'off') by default. 
    % This prevents MATLAB/OS from resizing it, preserving your label layouts 
    % for export.

    arguments
        horizontalPortion (1,1) double {mustBeGreaterThan(horizontalPortion,0), mustBeLessThanOrEqual(horizontalPortion,1)} = 0.5
        verticalPortion   (1,1) double {mustBeGreaterThan(verticalPortion,0),   mustBeLessThanOrEqual(verticalPortion,1)}   = 0.5
    
        options.ScreenSize (1,4) double = [0 0 3440 1440]   % reference [left bottom width height] (pixels)
        options.Margin (1,1) double = 0                     % pixels; keeps title bar reachable
    
        options.FontSize (1,1) double = 20
        options.AxesLineWidth (1,1) double = 1.5
        
        % 'auto' = turn off visibility if too big for screen
        options.Visible {mustBeMember(options.Visible, {'on', 'off', 'auto'})} = 'auto'
    end
    
    % --- Reference geometry (controls size) ---
    % Use the user-supplied reference screen size for consistent pixel sizing.
    % Fall back to the current screen if none provided.
    ref = options.ScreenSize;
    if isempty(ref)
        try
            ref = get(groot, 'ScreenSize');
        catch
            ref = get(0, 'ScreenSize');
        end
    end
    
    desiredW = round(ref(3) * horizontalPortion);
    desiredH = round(ref(4) * verticalPortion);
    
    % --- Monitor geometry (controls placement/check) ---
    % MonitorPositions is in global desktop coordinates: [left bottom width height]
    try
        monitors = get(groot, 'MonitorPositions');
    catch
        monitors = get(0, 'MonitorPositions');
    end

    % Choose the monitor containing the mouse pointer when possible;
    % otherwise default to the first monitor.
    mon = monitors(1,:);
    try
        ptr = get(groot, 'PointerLocation');
        for k = 1:size(monitors,1)
            m = monitors(k,:);
            if ptr(1) >= m(1) && ptr(1) <= (m(1)+m(3)) && ptr(2) >= m(2) && ptr(2) <= (m(2)+m(4))
                mon = m;
                break;
            end
        end
    catch
        % PointerLocation not available in very old releases; keep default.
    end
    
    monLeft   = mon(1);
    monBottom = mon(2);
    monTop    = mon(2) + mon(4);
    
    % --- Determine Visibility ---
    % Check if desired size exceeds current monitor dimensions
    isTooBig = (desiredW > mon(3)) || (desiredH > mon(4));
    
    if strcmp(options.Visible, 'auto')
        if isTooBig
            finalVisible = 'off';
            fprintf('Info: Figure (%dx%d) exceeds monitor. Created INVISIBLY to preserve layout.\n', desiredW, desiredH);
        else
            finalVisible = 'on';
        end
    else
        finalVisible = options.Visible;
    end

    % --- Calculate Position ---
    if strcmp(finalVisible, 'on')
        % Center on the chosen monitor
        left   = monLeft   + round((mon(3) - desiredW)/2);
        bottom = monBottom + round((mon(4) - desiredH)/2);

        % If bigger than monitor (and visible), anchor top-left so the title bar is reachable.
        if desiredW > mon(3)
            left = monLeft + options.Margin;
        end
        if desiredH > mon(4)
            bottom = (monTop - options.Margin) - desiredH;
        end
        pos = [left, bottom, desiredW, desiredH];
    else
        % IMPORTANT FIX:
        % When invisible (especially if requested size exceeds the current monitor),
        % avoid negative / partially-on-screen coordinates.
        % Some OS/window managers clamp or shrink windows they try to keep "on screen",
        % even if MATLAB sets 'Visible' to 'off'.
        % Putting the window fully off-screen prevents that clamping.
        maxRight  = max(monitors(:,1) + monitors(:,3));
        minBottom = min(monitors(:,2));
        pos = [maxRight + 50, minBottom + 50, desiredW, desiredH];
    end
    
    % --- Create figure ---
    % Keep units in pixels so the requested size is respected.
    % Also disable interactive resizing (helps avoid some window-manager adjustments).
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', pos, ...
                 'Visible', finalVisible, 'Resize', 'off');

    % Force the requested size again after the window is created.
    % This helps on platforms that clamp at creation time.
    drawnow('nocallbacks');
    set(fig, 'Units', 'pixels');
    p = fig.Position;
    p(3:4) = [desiredW desiredH];
    fig.Position = p;
    
    % Make print/saveas sizing independent of current monitor size.
    % Convert pixels -> inches using the detected screen PPI.
    try
        ppi = get(groot, 'ScreenPixelsPerInch');
    catch
        ppi = 96; % reasonable fallback
    end
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 desiredW/ppi desiredH/ppi];
    fig.PaperSize = [desiredW/ppi desiredH/ppi];
    
    % 'manual' tells MATLAB: "Do not look at the screen size when saving.
    % Use the PaperPosition I just gave you."
    fig.PaperPositionMode = 'manual';
    
    % --- CLAMPING (Only if Visible) ---
    % We only clamp to screen edges if the user can actually see the window.
    % If invisible, we strictly maintain the 'pos' calculated above.
    if strcmp(finalVisible, 'on')
        drawnow;
        set(fig, 'Units', 'pixels');
        op = fig.OuterPosition;
        
        % Clamp "top-left reachable":
        if op(1) < monLeft + options.Margin
            op(1) = monLeft + options.Margin;
        end
        if (op(2) + op(4)) > (monTop - options.Margin)
            op(2) = (monTop - options.Margin) - op(4);
        end
        
        fig.OuterPosition = op;
    end
    
    % Default axes styling
    hold on
    ax = gca;
    set(ax, 'FontSize', options.FontSize, 'LineWidth', options.AxesLineWidth);
    % NOTE: PositionConstraint is an AXES property, not a FIGURE property.
    % Use it when available to keep the axes' inner size stable.
    if isprop(ax, 'PositionConstraint')
        ax.PositionConstraint = 'innerposition';
    elseif isprop(ax, 'ActivePositionProperty')
        ax.ActivePositionProperty = 'position';
    end
    box off
end