function [T,meta] = inputParamsTable_build(Prompt,Formats,DefAns)
%INPUTPARAMSTABLE_BUILD Convert inputsdlg-style dialog inputs into a table.
%
%   [T,META] = INPUTPARAMSTABLE_BUILD(PROMPT,FORMATS,DEFANS) turns the
%   PROMPT/FORMATS/DEFANS triple used by inputsdlg into a MATLAB table with
%   one row per session and one variable per parameter, plus a META struct
%   array describing each column.
%
%   Column types are chosen so that uitable derives the right editor:
%       inactive  -> string    (read-only label, e.g. the session name)
%       check     -> logical   (check box)
%       edit      -> string    (text field)
%       list      -> categorical with categories set to Formats.items (drop-down)
%
%   META has one element per column with fields: field, header, kind, items,
%   editable, width.
%
%   Pure function, no graphics. See also inputParamsTable, inputParamsTable_toAnswer.

arguments
    Prompt cell
    Formats struct
    DefAns struct
end

[headers,fields] = parsePrompt(Prompt);
Formats = Formats(:);

if numel(Formats) ~= numel(fields)
    error('inputParamsTable:SizeMismatch',...
        'Formats has %d entries but Prompt has %d.',numel(Formats),numel(fields));
end

nSessions = numel(DefAns);
if nSessions == 0
    error('inputParamsTable:NoSessions','DefAns is empty; nothing to edit.');
end

meta = struct('field',{},'header',{},'kind',{},'items',{},'editable',{},'width',{});
vars = {};

for k = 1:numel(Formats)
    kind = formatKind(Formats(k),fields{k});
    if strcmp(kind,'none'); continue; end

    if ~isfield(DefAns,fields{k})
        error('inputParamsTable:MissingField',...
            'DefAns has no field ''%s''.',fields{k});
    end
    raw = {DefAns.(fields{k})}';

    items = {};
    switch kind
        case {'label','edit'}
            col = strings(nSessions,1);
            for s = 1:nSessions
                col(s) = toDisplayText(raw{s},fields{k});
            end
            width = textWidth([col; string(headers{k})]);

        case 'check'
            col = false(nSessions,1);
            for s = 1:nSessions
                v = raw{s};
                if ~(islogical(v) || isnumeric(v)) || ~isscalar(v)
                    error('inputParamsTable:BadCheckValue',...
                        'Field ''%s'' is a check box and needs a logical scalar.',fields{k});
                end
                col(s) = logical(v);
            end
            width = max(60,textWidth(string(headers{k})));

        case 'list'
            items = validateItems(Formats(k),fields{k});
            idx = zeros(nSessions,1);
            for s = 1:nSessions
                idx(s) = itemIndex(raw{s},items,fields{k});
            end
            col = categorical(string(items(idx))',string(items));
            width = textWidth([string(items(:)); string(headers{k})]) + 24;
    end

    vars{end+1} = col; %#ok<AGROW>
    meta(end+1) = struct('field',fields{k},'header',headers{k},'kind',kind,...
        'items',{items},'editable',~strcmp(kind,'label'),'width',width); %#ok<AGROW>
end

if isempty(meta)
    error('inputParamsTable:NoColumns','No editable columns were produced.');
end

T = table(vars{:},'VariableNames',{meta.field});

end

%% ------------------------------------------------------------------------

function [headers,fields] = parsePrompt(Prompt)
% inputsdlg takes PROMPT as an N-by-1 (label only) or N-by-2 (label, field
% name) cell. The wrappers pass repmat(names',1,2), so both columns match.

if size(Prompt,2) >= 2
    headers = Prompt(:,1);
    fields  = Prompt(:,2);
else
    headers = Prompt(:);
    fields  = Prompt(:);
end

headers = cellfun(@char,headers,'UniformOutput',false);
fields  = cellfun(@char,fields, 'UniformOutput',false);

if numel(unique(fields)) ~= numel(fields)
    error('inputParamsTable:DuplicateField','Prompt contains duplicate field names.');
end
end

function kind = formatKind(fmt,field)
% Map one inputsdlg Formats entry onto a column kind. Unset struct fields
% come back as [] once the struct array is filled in, hence fld().

if any(strcmp(fld(fmt,'enable',''),{'inactive','off'}))
    kind = 'label'; return
end

type = lower(fld(fmt,'type',''));
switch type
    case 'check', kind = 'check';
    case 'edit',  kind = 'edit';
    case 'list',  kind = 'list';
    case 'none',  kind = 'none';
    case '',      kind = 'edit';  % inputsdlg's own default
    otherwise
        error('inputParamsTable:UnsupportedFormat',...
            ['Field ''%s'' uses Formats.type=''%s'', which inputParamsTable does '...
             'not support. Supported: check, edit, list, none, or enable=''inactive''.'],...
            field,type);
end
end

function v = fld(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function items = validateItems(fmt,field)
items = fld(fmt,'items',{});
if isempty(items)
    error('inputParamsTable:MissingItems',...
        'Field ''%s'' is a list but Formats.items is empty.',field);
end
if isstring(items); items = cellstr(items); end
if ~iscellstr(items) %#ok<ISCLSTR>
    error('inputParamsTable:BadItems','Formats.items for ''%s'' must be text.',field);
end
items = items(:)';
% The list index is recovered on output by matching the cell's text back
% against items, so duplicates would silently return the wrong index.
if numel(unique(items)) ~= numel(items)
    error('inputParamsTable:DuplicateItems',...
        'Formats.items for ''%s'' contains duplicates.',field);
end
end

function i = itemIndex(v,items,field)
% Accept either a numeric index (what the wrappers pass) or the item text.
if isnumeric(v) && isscalar(v)
    i = round(v);
    if i < 1 || i > numel(items)
        error('inputParamsTable:IndexOutOfRange',...
            'Default index %g for ''%s'' is outside 1..%d.',v,field,numel(items));
    end
elseif ischar(v) || isstring(v)
    i = find(strcmp(items,char(v)),1);
    if isempty(i)
        error('inputParamsTable:UnknownItem',...
            'Default ''%s'' for ''%s'' is not in Formats.items.',char(v),field);
    end
else
    error('inputParamsTable:BadListValue',...
        'Default for ''%s'' must be an index or an item name.',field);
end
end

function s = toDisplayText(v,field)
% Callers evaluate edit fields with str2double/eval, so anything that is not
% already text is rendered with num2str rather than left as a number.
if ischar(v)
    s = string(v);
elseif isstring(v) && isscalar(v)
    s = v;
elseif (isnumeric(v) || islogical(v)) && ~isempty(v)
    s = string(num2str(v));
elseif isempty(v)
    s = "";
else
    error('inputParamsTable:BadEditValue',...
        'Default for ''%s'' must be text or numeric.',field);
end
end

function w = textWidth(s)
% Rough pixel width for a column that has to show the longest of s.
n = max([strlength(s); 1]);
w = min(max(double(n)*7.5 + 22, 60), 220);
end
