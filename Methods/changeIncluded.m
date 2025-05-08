function newTable = changeIncluded(combined_cells)

newTable = combined_cells;

for row = 1:size(combined_cells,1)
    included = combined_cells(row,:).Included{1};
    for vhold = 1:length(included)
        if sum(included{vhold}) == 0
            included{vhold} = ones(length(included{vhold}),1);
        end
    end
    newTable(row,:).Included{1} = included;
end

end