function combinedTable = concatTables(table1, table2)

arguments
    table1 table
    table2 table
end

% Get variable names from each table
if isempty(table1)
    combinedTable = table2;
    return
else
    table1Vars = table1.Properties.VariableNames;
end
if isempty(table2)
    combinedTable = table1;
    return
else
    table2Vars = table2.Properties.VariableNames;
end

% Union of all variable names
allTableVars = union(table1Vars, table2Vars);

%% Fill missing variables in table1
missingVars_table1 = setdiff(allTableVars, table1Vars);
for m = 1:length(missingVars_table1)
    varName = missingVars_table1{m};
    % Use table2 as reference if it contains the variable.
    if ismember(varName, table2Vars)
        sampleVar = table2.(varName);
        if islogical(sampleVar)
            defaultMissing = false;  % For logical, use false
        elseif isnumeric(sampleVar)
            defaultMissing = NaN;    % For numeric, use NaN
        elseif isdatetime(sampleVar)
            defaultMissing = NaT;    % For datetime, use NaT
        elseif iscellstr(sampleVar) || isstring(sampleVar)
            defaultMissing = missing;
        else
            defaultMissing = missing;
        end
    else
        defaultMissing = missing;
    end
    table1.(varName) = repmat(defaultMissing, height(table1), 1);
end

%% Fill missing variables in table2
missingVars_table2 = setdiff(allTableVars, table2Vars);
for m = 1:length(missingVars_table2)
    varName = missingVars_table2{m};
    % Use table1 as reference if it contains the variable.
    if ismember(varName, table1Vars)
        sampleVar = table1.(varName);
        if islogical(sampleVar)
            defaultMissing = false;
        elseif isnumeric(sampleVar)
            defaultMissing = NaN;
        elseif isdatetime(sampleVar)
            defaultMissing = NaT;
        elseif iscellstr(sampleVar) || isstring(sampleVar)
            defaultMissing = missing;
        else
            defaultMissing = missing;
        end
    else
        defaultMissing = missing;
    end
    table2.(varName) = repmat(defaultMissing, height(table2), 1);
end

%% Reorder columns so both tables have the same order
table1 = table1(:, allTableVars);
table2 = table2(:, allTableVars);

%% Concatenate the tables vertically
combinedTable = [table1; table2];

end