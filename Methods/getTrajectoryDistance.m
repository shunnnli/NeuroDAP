function distance = getTrajectoryDistance(X,Y,options)

arguments
    X double
    Y double
    options.filter logical % select frames to calculate distance on
end

% Calculate Euclidean distances between consecutive selected frames.
% A tracking gap (NaN/Inf coordinate) invalidates only the steps that
% touch it; it must not make the entire session distance NaN or connect
% the two valid points on either side of the gap.
x = X(options.filter);
y = Y(options.filter);

if numel(x) < 2
    distance = 0;
    return;
end

validStep = isfinite(x(1:end-1)) & isfinite(y(1:end-1)) & ...
            isfinite(x(2:end))   & isfinite(y(2:end));
d = hypot(diff(x), diff(y));

% Total distance traveled is the sum of valid consecutive steps.
distance = sum(d(validStep));

end
