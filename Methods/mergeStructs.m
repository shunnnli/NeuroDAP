function masterStruct = mergeStructs(structList, options)

% If a cell array of structs are given, where each element is a structure
% with the same set of fields, merge them into one such that each field
% now contains a vector of original values

arguments
    structList cell
    options.summary string = 'None'
    options.reshape logical = false % normally not used
    options.combine logical = true % combine data from multiple structs
end

% Initialize an empty master struct
masterStruct = struct();
% Iterate over the cell array
for i = 1:length(structList)
    currentStruct = structList{i}; % Extract the current struct
    masterStruct = mergeStructFields(masterStruct,currentStruct,...
                    length(structList),i,...
                    summary=options.summary,...
                    combine=options.combine,...
                    reshape=options.reshape);
end

end


function masterStruct = mergeStructFields(masterStruct,currentStruct,nStructs,structIdx,options)

    arguments
        masterStruct struct
        currentStruct struct
        nStructs
        structIdx
        options.summary string
        options.combine logical = true
        options.reshape logical = false
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
                                        nStructs,structIdx,...
                                        summary=options.summary,...
                                        combine=options.combine,...
                                        reshape=options.reshape);
        else
            % If the field is not a struct, append the values
            if ~isfield(masterStruct, field)
                masterStruct.(field) = cell(nStructs,1);
            end
            masterStruct.(field){structIdx} = currentStruct.(field);

            % Convert cell arrays to vectors if possible
            if structIdx == nStructs
                if options.combine && all(cellfun(@isnumeric, masterStruct.(field)))
                    try
                        vals = masterStruct.(field);
                        % Force each entry to be a row vector (flatten everything)
                        vals = cellfun(@(x) reshape(x, 1, []), vals, 'UniformOutput', false);
                        % Pad to equal length so cell2mat works
                        maxLen = max(cellfun(@numel, vals));
                        vals = cellfun(@(x) [x, nan(1, maxLen - numel(x))], vals, 'UniformOutput', false);
                        values = cell2mat(vals);  % => nStructs x maxLen
                    catch ME
                        disp(getReport(ME));
                        warning(['Field ', field, ' have an error, skipped for now!!!!']);
                        continue
                    end
                    
                    % Normally not needed
                    if options.reshape
                        values = reshape(values,nStructs,[]);
                    end

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