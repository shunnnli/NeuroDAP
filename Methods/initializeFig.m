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
    ref = options.ScreenSize;
    if isempty(ref)
        ref = get(0, 'ScreenSize');
    end
    
    desiredW = round(ref(3) * horizontalPortion);
    desiredH = round(ref(4) * verticalPortion);
    
    % --- Primary monitor geometry (controls placement/check) ---
    monitors = get(0, 'MonitorPositions'); 
    mon = monitors(1,:);                  
    
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
    % Center on primary monitor
    left   = monLeft   + round((mon(3) - desiredW)/2);
    bottom = monBottom + round((mon(4) - desiredH)/2);
    
    % If bigger than monitor (and we are forcing it visible), anchor top-left 
    % so title bar is reachable.
    if desiredW > mon(3)
        left = monLeft + options.Margin;
    end
    if desiredH > mon(4)
        bottom = (monTop - options.Margin) - desiredH; 
    end
    
    pos = [left, bottom, desiredW, desiredH];
    
    % --- Create figure ---
    fig = figure('Color', 'w', 'Units', 'pixels', 'Position', pos, ...
                 'Visible', finalVisible);
    
    % We use 'points' as a proxy for pixels (1 pt approx 1 pixel at standard 
    % screen DPI) to maintain the font-to-canvas ratio.
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 desiredW desiredH];
    
    % Ensure the paper size is large enough to hold the position
    fig.PaperSize = [desiredW desiredH];
    
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
    set(gca, 'FontSize', options.FontSize, 'LineWidth', options.AxesLineWidth);
    box off
end