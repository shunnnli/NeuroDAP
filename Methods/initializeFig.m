function fig = initializeFig(horizontalPortion, verticalPortion, options)
%INITIALIZEFIG Create a figure sized as a portion of a (reference) screen.
%
%   fig = initializeFig(hPortion, vPortion) sizes the figure as a fraction of
%   the *current* monitor's ScreenSize.
%
%   fig = initializeFig(..., 'ScreenSize', [left bottom width height]) sizes
%   the figure as a fraction of a user-provided *reference* ScreenSize.
%
%   Example (use a fixed reference screen size):
%       ref = [1 1 3440 1440];
%       fig = initializeFig(0.6, 0.6, 'ScreenSize', ref);

arguments
    horizontalPortion (1,1) double {mustBeGreaterThan(horizontalPortion,0), mustBeLessThanOrEqual(horizontalPortion,1)} = 0.5
    verticalPortion   (1,1) double {mustBeGreaterThan(verticalPortion,0),   mustBeLessThanOrEqual(verticalPortion,1)}   = 0.5

    options.ScreenSize (1,4) double = [0 0 3440 1440]   % [left bottom width height] (pixels)
    options.FontSize (1,1) double = 20
    options.AxesLineWidth (1,1) double = 1.5
end

% Choose base screen geometry
if isempty(options.ScreenSize)
    baseScreen = get(0, 'ScreenSize');
else
    baseScreen = options.ScreenSize;
end

% Scale width/height by requested portions
pos = baseScreen;
pos(3) = pos(3) * horizontalPortion;
pos(4) = pos(4) * verticalPortion;

% Create figure
fig = figure('Position', pos, 'Color', 'w');

% Default axes styling (applies once axes exist)
hold on
set(gca, 'FontSize', options.FontSize, 'LineWidth', options.AxesLineWidth);
box off

end
