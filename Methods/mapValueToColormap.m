function [mappedIdx,mappedColors] = mapValueToColormap(value,colormap)

% Assume slope is your vector of values to map.
% Determine a symmetric bound around 0.
absBound = max(max(abs(value)), eps);  % eps to avoid division by zero

% Get the number of colors in your custom colormap
nColors = size(colormap, 1);

% Map slope values to indices between 1 and nColors:
%  - When slope equals -absBound, this formula gives index=1.
%  - When slope equals 0, it gives index around mid-color.
%  - When slope equals +absBound, it gives index=nColors.
mappedIdx = round((value + absBound) / (2 * absBound) * (nColors - 1)) + 1;

% Now get the corresponding RGB values for each point.
mappedColors = colormap(mappedIdx, :);


end