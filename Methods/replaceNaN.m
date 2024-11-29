function result = replaceNaN(array,val)

arguments
    array
    val double
end

result = array;
% Loop through columns
for i = 1:size(result,2)
    if strcmpi(class(result),'table') 
        col = result{:,i};
        if isnumeric(col)
            result{isnan(col),i} = val;
        end
    else
        result(isnan(result(:,i)),i) = val; 
    end
end

end