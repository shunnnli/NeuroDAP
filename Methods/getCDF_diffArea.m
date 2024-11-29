function diffArea = getCDF_diffArea(observed,bootstrapped,options)

arguments
    observed double % a vector saving the observed data
    bootstrapped double % a matrix where each column stores the result of each simulation

    options.simulationInCol logical = true % if false, simulation is stored at each row
    options.abs logical = false % calculate abs diff in area or not
    options.separate logical = true % is false, view bootstrapped matrix as single observation
end

% Input matrices
A = bootstrapped; B = observed;
if ~options.simulationInCol; A = A'; end

% Precompute CDF for vector B
if ~isvector(B); B = reshape(B, [], 1); end
[fB, xB] = ecdf(B);
[xB_unique, idx_unique] = unique(xB); % Ensure unique xB values
fB_unique = fB(idx_unique);           % Adjust fB accordingly

% Initialize an array to store area differences
diffArea = zeros(1, size(A, 2));

if options.separate
    % Loop over each column of A
    for i = 1:size(A, 2)
        % Compute CDF for the current column of A
        [fA, xA] = ecdf(A(:, i));
    
        % Interpolate the CDF of B to match the x-values of the current column
        fB_interp = interp1(xB_unique, fB_unique, xA, 'linear', 'extrap');
    
        % Compute the absolute area difference between the two CDFs
        if options.abs
            diffArea(i) = trapz(xA, abs(fA - fB_interp));
        else
            diffArea(i) = trapz(xA, fA - fB_interp);
        end
    end
else
    A_vector = reshape(A, [], 1);
    [fA, xA] = ecdf(A_vector);
    % Interpolate the CDF of B to match the x-values of the current column
    fB_interp = interp1(xB_unique, fB_unique, xA, 'linear', 'extrap');

    % Compute the absolute area difference between the two CDFs
    if options.abs
        diffArea = trapz(xA, abs(fA - fB_interp));
    else
        diffArea = trapz(xA, fA - fB_interp);
    end
end

end