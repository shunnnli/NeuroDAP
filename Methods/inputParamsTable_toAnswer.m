function Answer = inputParamsTable_toAnswer(T,meta)
%INPUTPARAMSTABLE_TOANSWER Convert an edited parameter table back to a struct array.
%
%   ANSWER = INPUTPARAMSTABLE_TOANSWER(T,META) returns an nSessions-by-1
%   struct array whose field values match what inputsdlg would have returned:
%
%       label -> char
%       check -> logical scalar
%       edit  -> char        (callers run str2double/eval on these)
%       list  -> double, the 1-based index into META.items
%
%   The char and index conventions are load-bearing. Shun_loadSessionData.m:22
%   runs eval() on a char, and :38 indexes a cell array with the list value.
%
%   Pure function, no graphics. See also inputParamsTable, inputParamsTable_build.

arguments
    T table
    meta struct
end

nSessions = height(T);
values = cell(numel(meta),nSessions);

for k = 1:numel(meta)
    col = T.(meta(k).field);
    switch meta(k).kind
        case {'label','edit'}
            values(k,:) = cellstr(string(col))';

        case 'check'
            values(k,:) = num2cell(logical(col))';

        case 'list'
            txt = string(col);
            if any(ismissing(txt))
                error('inputParamsTable:MissingListValue',...
                    'Column ''%s'' has an unset drop-down value.',meta(k).field);
            end
            [ok,loc] = ismember(txt,string(meta(k).items));
            if ~all(ok)
                error('inputParamsTable:UnknownItem',...
                    'Column ''%s'' holds a value outside its item list.',meta(k).field);
            end
            values(k,:) = num2cell(double(loc))';

        otherwise
            error('inputParamsTable:UnknownKind',...
                'Unknown column kind ''%s''.',meta(k).kind);
    end
end

Answer = cell2struct(values,{meta.field},1);
Answer = Answer(:);

end
