function masterStruct = mergeStructs(structList, options)

% If a cell array of structs are given, where each element is a structure
% with the same set of fields, merge them into one such that each field
% now contains a vector of original values

arguments
    structList cell
    options.summary string = 'None'
end

% Initialize an empty master struct
masterStruct = struct();
% Iterate over the cell array
for i = 1:length(structList)
    currentStruct = structList{i}; % Extract the current struct
    masterStruct = mergeStructFields(masterStruct,currentStruct,...
        length(structList),i,summary=options.summary);
end

end


function masterStruct = mergeStructFields(masterStruct,currentStruct,nStructs,structIdx,options)

    arguments
        masterStruct struct
        currentStruct struct
        nStructs
        structIdx
        options.summary string
    end
    
    fields = fieldnames(currentStruct); % Get the field names
    for j = 1:numel(fields)
        field = fields{j}; % Current field name
        if isstruct(currentStruct.(field))
            % If the field is a struct, call the function recursively
            if ~isfield(masterStruct, field)
                masterStruct.(field) = struct(); % Initialize empty struct if it doesn't exist
            end
            masterStruct.(field) = mergeStructFields(masterStruct.(field), currentStruct.(field),...
                                        nStructs,structIdx,summary=options.summary);
        else
            % If the field is not a struct, append the values
            if ~isfield(masterStruct, field)
                masterStruct.(field) = cell(nStructs,1);
            end
            masterStruct.(field){structIdx} = currentStruct.(field);

            % Convert cell arrays to vectors if possible
            if structIdx == nStructs
                if all(cellfun(@isnumeric, masterStruct.(field)))
                    values = cell2mat(masterStruct.(field));
                    if ~isfield(options,'summary')
                        masterStruct.(field) = values;
                    elseif contains(options.summary,{'average','mean','avg'},IgnoreCase=true)
                        masterStruct.(field) = mean(values);
                    elseif contains(options.summary,{'median','med','mid'},IgnoreCase=true)
                        masterStruct.(field) = median(values);
                    elseif contains(options.summary,{'mode'},IgnoreCase=true)
                        masterStruct.(field) = mode(values);
                    else
                        masterStruct.(field) = values;
                    end
                end
            end
        end
    end
end


%% Past code