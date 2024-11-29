function sortedStruct = sortrows_struct(targetStruct, column)

arguments
    targetStruct struct
    column cell % names of the column
end

T = struct2table(targetStruct);

% % Find column index
% colIdx= [];
% for i = 1:length(column)
%     colIdx = [colIdx find(strcmpi(T.Properties.VariableNames,column{i}))];
% end


sortedT = sortrows(T,column);
sortedStruct = table2struct(sortedT);

end