function [] = editStructValue(targetStruct,row,targetField,value)

% UNFINISHED!!!!
% This function edits value a specific subfield of a struct 
arguments
    targetStruct struct % struct to be edited
    row % row number of the struct
    targetField % name of the fields to edit in string or cell array
    value
end

% Check if value and target field are of same length
if length(targetField) ~= length(value) || ~isa(targetField,'cell')
    
    % if targetField is char, find subfields under that field 
    if ischar(targetField)
        subfields = targetStruct.(targetField);
        subfieldsNum = numel(fieldnames(subfields));
        if subfieldsNum ~= length(value)
            error('Number of target fields and values does not match!');
        else
            targets = cell(subfieldsNum,1);
            subfieldNames = fieldnames(subfields);
            for i = 1:subfieldsNum
                targets{i} = subfieldNames{i};
            end
        end
    
    % if targetField is cell (loop through all subfields)
    elseif isa(targetField,'cell')
        subfieldsNum = 0;
        subfields = cell(length(targetField),1);
        for i = 1:length(targetField)
            subfields{i} = targetStruct.(targetField{i});
            subfieldsNum = subfieldsNum + numel(fieldnames(subfields));
        end
        
        % Check if total number of subfields equals value
        if subfieldsNum ~= length(value)
            error('Number of target fields and values does not match!');
        else
            targets = cell(subfieldsNum,1);
            index = 0;
            for i = 1:length(subfields) 
                subfieldNames = fieldnames(subfields{i});
                for j = 1:length(subfieldNames)
                    
                end
            end
            subfieldNames = fieldnames(subfields);
            for i = 1:subfieldsNum
                targets{i} = subfieldNames{i};
            end
        end
    else
        error('targetField should be either cell arrray or a char!');
    end
end

% Get subfield names if nesccessary
if ~isa(targetField,'cell') && subfieldsNum > 1
    
else
    targets = {targetField};
end

% Edit field value
for f = 1:length(targets)
    cur_field = targets{f};
    struct_in_string = getVarName(targetStruct);
    
    % Skip if theres no value or says 'all' in execel sheet
    if ~isnan(value(f)) || ~strcmp(value(f),'all')
        evalin('base',[struct_in_string, '.(', num2str(row),').', cur_field,...
            '=', num2str(value(f))]);
    end
end

end