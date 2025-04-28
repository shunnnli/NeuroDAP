function distance = getTrajectoryDistance(X,Y,options)

arguments
    X double
    Y double
    options.filter logical % select frames to calculate distance on
end

% Calculate Euclidean distances between consecutive selected frames
d = sqrt(diff(X(options.filter)).^2 + diff(Y(options.filter)).^2);

% Total distance traveled is the sum of these distances
distance = sum(d);

end