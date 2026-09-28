function [Answer,Canceled] = inputParamsTable(Prompt,Title,Formats,DefAns)
%INPUTPARAMSTABLE Edit per-session parameters in a scrollable table.
%
%   [ANSWER,CANCELED] = INPUTPARAMSTABLE(PROMPT,TITLE,FORMATS,DEFANS) is a
%   drop-in replacement for the inputsdlg call used by inputAnalysisParams,
%   inputSessionParams and inputSessionParams_singleSlice.
%
%   inputsdlg tiles one column of controls per session, so past roughly eight
%   sessions the figure is wider than the screen and the remaining sessions
%   cannot be reached. This lays the same parameters out as a table with one
%   ROW per session, which scrolls and therefore works at any session count.
%
%   ANSWER is an nSessions-by-1 struct array using the same value conventions
%   as inputsdlg (edit fields as char, drop-downs as a 1-based index), so
%   calling scripts need no changes. CANCELED is 0 on OK and 1 on Cancel or
%   window close; on cancel ANSWER holds the default answers, again matching
%   inputsdlg.
%
%   Select a cell and use the fill buttons to push its value down the whole
%   column, or to just the rows in a multi-row selection.
%
%   See also inputParamsTable_build, inputParamsTable_toAnswer, inputsdlg.

arguments
    Prompt cell
    Title (1,:) char
    Formats struct
    DefAns struct
end

[T,meta] = inputParamsTable_build(Prompt,Formats,DefAns);

Canceled = 1;
Answer = inputParamsTable_toAnswer(T,meta);

nSessions = height(T);
screen = get(groot,'ScreenSize');

fig = uifigure('Name',Title,'Visible','off');

% The figure must be deleted explicitly rather than by an onCleanup object.
% Its callbacks are handles to the nested functions below, which share this
% workspace, so figure -> callback -> workspace -> onCleanup -> figure forms
% a reference cycle that MATLAB never collects. The cleanup would therefore
% never run and the dialog would stay on screen after OK.
try

figW = min(max(sum([meta.width]) + 60, 460), 0.9*screen(3));
figH = min(max(nSessions*22 + 190, 260), 0.85*screen(4));
fig.Position = [0 0 figW figH];
movegui(fig,'center');
try fig.WindowStyle = 'modal'; catch; end %#ok<CTCH>
fig.CloseRequestFcn = @(~,~)onCancel();

gl = uigridlayout(fig,[3 1]);
gl.RowHeight = {'1x',26,30};
gl.ColumnWidth = {'1x'};
gl.RowSpacing = 8;

tbl = uitable(gl,...
    'Data',T,...
    'ColumnName',{meta.header},...
    'ColumnEditable',logical([meta.editable]),...
    'ColumnWidth',num2cell([meta.width]),...
    'RowName',{},...
    'SelectionType','cell',...
    'Multiselect','on');
tbl.SelectionChangedFcn = @(~,~)refreshFillButtons();

fillRow = uigridlayout(gl,[1 3]);
fillRow.ColumnWidth = {215,140,'1x'};
fillRow.Padding = [0 0 0 0];
fillRow.ColumnSpacing = 8;

bFillCol = uibutton(fillRow,'Text','Fill column from selected cell',...
    'Enable','off','ButtonPushedFcn',@(~,~)onFill(false));
bFillSel = uibutton(fillRow,'Text','Fill selected rows',...
    'Enable','off','ButtonPushedFcn',@(~,~)onFill(true));
uilabel(fillRow,'Text',sprintf('%d sessions. Select a cell, then fill.',nSessions));

btnRow = uigridlayout(gl,[1 3]);
btnRow.ColumnWidth = {'1x',90,90};
btnRow.Padding = [0 0 0 0];
btnRow.ColumnSpacing = 8;

uilabel(btnRow,'Text','');
uibutton(btnRow,'Text','Cancel','ButtonPushedFcn',@(~,~)onCancel());
uibutton(btnRow,'Text','OK','ButtonPushedFcn',@(~,~)onOK());

fig.Visible = 'on';
uiwait(fig);

catch ME
    delete(fig);
    rethrow(ME);
end
delete(fig);

%% --- nested callbacks ------------------------------------------------

    function [row,col] = currentCell()
        row = []; col = [];
        sel = tbl.Selection;
        if isempty(sel); return; end
        row = sel(1,1);
        col = sel(1,2);
    end

    function refreshFillButtons()
        [row,col] = currentCell();
        canFill = ~isempty(row) && meta(col).editable;
        multiRow = canFill && numel(unique(tbl.Selection(:,1))) > 1;
        bFillCol.Enable = onOff(canFill);
        bFillSel.Enable = onOff(multiRow);
    end

    function onFill(selectedRowsOnly)
        [row,col] = currentCell();
        if isempty(row) || ~meta(col).editable; return; end
        D = tbl.Data;
        if selectedRowsOnly
            targets = unique(tbl.Selection(:,1));
        else
            targets = (1:height(D))';
        end
        D.(meta(col).field)(targets) = D{row,col};
        tbl.Data = D;
    end

    function onOK()
        try
            Answer = inputParamsTable_toAnswer(tbl.Data,meta);
        catch ME
            uialert(fig,ME.message,'Invalid value');
            return
        end
        Canceled = 0;
        uiresume(fig);
    end

    function onCancel()
        Canceled = 1;
        uiresume(fig);
    end

end

%% ------------------------------------------------------------------------

function s = onOff(tf)
if tf; s = 'on'; else; s = 'off'; end
end
