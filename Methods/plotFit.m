function plotFit(x, y, options)

arguments
    x double
    y double
    options.color
    options.LineWidth double = 2
    options.plotSEM logical = true
    options.SEMalpha double = 0.2
end

    % Remove NaNs
    validIdx = ~isnan(x) & ~isnan(y);
    x = x(validIdx);
    y = y(validIdx);

    % Need at least 2 points to fit a line
    if numel(x) < 2
        return
    end

    % Fit linear model
    mdl = fitlm(x, y);

    % Grid for smooth fit line
    xFit = linspace(min(x), max(x), 100)';
    [yFit, yCI] = predict(mdl, xFit, 'Alpha', 0.05);

    % Convert CI of fitted mean to SEM of fitted mean
    tcrit = tinv(0.975, mdl.DFE);
    ySEM = (yCI(:,2) - yCI(:,1)) ./ (2 * tcrit);

    % Shaded SEM band
    if options.plotSEM
        fill([xFit; flipud(xFit)], ...
             [yFit - ySEM; flipud(yFit + ySEM)], ...
             options.color, ...
             'FaceAlpha', options.SEMalpha, 'EdgeColor', 'none');
    end

    % Best-fit line
    plot(xFit, yFit, '-', 'Color', options.color, 'LineWidth', options.LineWidth);
end