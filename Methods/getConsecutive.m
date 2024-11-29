function Y = getConsecutive(vector)

if ~isa(vector,'double')
    % Convert vector to double
    vector = double(vector);
end

if size(vector,1) == 1
    d = [true, diff(vector)~=0, true];  % TRUE if values change
elseif size(vector,2) == 1
    d = [true; diff(vector)~=0; true];  % TRUE if values change
end
n = diff(find(d));               % Number of repetitions
Y = repelem(n, n);

end