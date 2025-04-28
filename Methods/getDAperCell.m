function DA_slope_per_cell = getDAperCell(combined_cells,DAslope,animalList)

[~, cellAnimalIdx] = ismember(combined_cells.Animal, animalList);

DA_slope_per_cell = zeros(size(cellAnimalIdx)); 
DA_slope_per_cell(cellAnimalIdx>0) = DAslope(cellAnimalIdx(cellAnimalIdx>0));

end