function h = plotStackTraces(data, options)
% plotStackTraces  Stack‐plot each row of a trials×time‐points matrix
%
%   h = plotStackTraces(data)
%   h = plotStackTraces(data, Name, Value, …)
%
% Inputs (Name/Value):
%   'x'       – 1×n or n×1 vector of x values (default: 1:n)
%   'color'   – either a single RGB triple (e.g. [0.5 0.5 0.5]) or
%               an m×3 matrix of RGB rows (one color per trace;
%               default: all gray)
%   'spacing' – scalar vertical offset between successive traces
%               (default: 1)
%
% Output:
%   h         – m×1 array of line handles

arguments
    data    (:,:) double
    options.x       double = []
    options.color   = [0.5 0.5 0.5]
    options.spacing (1,1) double = 1
    options.LineWidth double = 2
end

[m,n] = size(data);

% --- x axis ------------------------
if isempty(options.x)
    x = 1:n;
else
    x = options.x;
    assert(numel(x)==n, "Length of x must match #columns of data");
end

% --- colors -----------------------
% expand single RGB to m×3, or check user‐supplied colormap
if isvector(options.color) && numel(options.color)==3
    cmap = repmat(options.color(:)', m, 1);
else
    assert(size(options.color,1)==m && size(options.color,2)==3, ...
        "Color must be 1×3 RGB or m×3 colormap");
    cmap = options.color;
end

% --- spacing ----------------------
dx = options.spacing;

% --- plot -------------------------
hold on; box off;
h = gobjects(m,1);
for k = 1:m
    y = data(k,:) + (k-1)*dx;
    h(k) = plot(x, y, 'Color', cmap(k,:), 'LineWidth', options.LineWidth);
end

end

