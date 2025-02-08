function SEM = getSEM(data)

% Data must be a column vector
if size(data,1) == 1
    data = data';
end
SEM = zeros(size(data,2),1);

for i = 1:size(data,2)
    STD = std(data(:,i),"omitnan");
    SEM(i,1) = STD/sqrt(size(data,1));      % Standard Error
end

end % getSEM