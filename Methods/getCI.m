function CI = getCI(data,confidence)

% Data must be a column vector
if size(data,1) == 1
    data = data';
end
CI = zeros(size(data,2),1);

for i = 1:size(data,2)
    STD = std(data(:,i),"omitnan");
    SEM = STD/sqrt(size(data,1));      % Standard Error
    ts = tinv((1-confidence)/2,size(data,1)-1);   % T-Score
    CI(i,1) = ts*SEM; % Total error bar lengths are double the err values
end

end % getCI
