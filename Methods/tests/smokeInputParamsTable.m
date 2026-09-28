function smokeInputParamsTable()
%SMOKEINPUTPARAMSTABLE Drive the live inputParamsTable dialog end to end.
%
%   Run on demand:  smokeInputParamsTable
%
%   testInputParamsTable covers the pure table<->struct logic. This covers
%   what that cannot: that uitable really does render a categorical column
%   as a drop-down and a logical column as a check box, that the fill
%   buttons work against a live Selection, and that the figure fits on
%   screen at 30 sessions. It opens real windows, so it is deliberately
%   named to fall outside runtests' test-file pattern.
%
%   Timer callbacks swallow errors, so drive() logs to a flushed file and
%   ALWAYS resumes the figure; a watchdog force-closes anything left over.

LOG = fullfile(tempdir,'smokeInputParamsTable.log');
if exist(LOG,'file'); delete(LOG); end
setappdata(groot,'smokeLog',LOG);

say(LOG,'screen = %s',mat2str(get(groot,'ScreenSize')));

sessionList = arrayfun(@(k)sprintf('/data/202401%02d-SL001-g0',k),1:30,'UniformOutput',false)';

% The very first uifigure of a MATLAB session does not expose its children
% to findall until the web-window framework has booted, so burn one dialog
% before asserting anything. Results deliberately ignored.
say(LOG,'--- warm-up (not asserted) ---');
runWith('cancel',@()inputAnalysisParams(sessionList(1)));
% Drop the warm-up's noise so it cannot colour the pass/fail scan below.
delete(LOG);
say(LOG,'screen = %s (warm-up done)',mat2str(get(groot,'ScreenSize')));

say(LOG,'--- 30 sessions, OK path ---');
[A,C] = runWith('ok',@()inputSessionParams(sessionList,paradigm=2));
say(LOG,'Canceled=%d n=%d',C,numel(A));
check(LOG,C==0,'Canceled==0');
check(LOG,numel(A)==30,'30 answers');
check(LOG,all([A.Paradigm]==4),'fill column pushed last item (4) to every row');
check(LOG,ischar(A(1).ReactionTime),'edit field stayed char');
check(LOG,islogical(A(1).redStim),'check field is logical');

say(LOG,'--- 30 sessions, Cancel path ---');
[A2,C2] = runWith('cancel',@()inputSessionParams(sessionList,paradigm=3));
check(LOG,C2==1,'Canceled==1');
check(LOG,all([A2.Paradigm]==3),'cancel returned untouched defaults');

say(LOG,'--- inputSessionParams_singleSlice, 12 sessions ---');
[A4,C4] = runWith('ok',@()inputSessionParams_singleSlice(sessionList(1:12),paradigm=2,animal='SL999'));
check(LOG,C4==0 && numel(A4)==12,'singleSlice returned 12 answers');
check(LOG,all([A4.Paradigm]==5),'singleSlice fill used its own 5-item list');
check(LOG,ischar(A4(1).Animal) && strcmp(A4(1).Animal,'SL999'),'singleSlice Animal is char');
check(LOG,isequal(eval(A4(1).timeRange),[-20 100]),'singleSlice timeRange still evals');

say(LOG,'--- 1 session ---');
[A3,C3] = runWith('ok',@()inputAnalysisParams(sessionList(1)));
check(LOG,C3==0 && isscalar(A3),'single session');
check(LOG,isequal(eval(A3(1).recordLJ),[1 1 0]),'recordLJ still evals');

type(LOG);
txt = fileread(LOG);
if contains(txt,'FAIL') || contains(txt,'ERROR')
    error('smokeInputParamsTable:failed','smoke test failed');
end
fprintf('\nSMOKE OK\n');
end

%% ------------------------------------------------------------------

function [A,C] = runWith(action,fn)
LOG = getappdata(groot,'smokeLog');
t = timer('StartDelay',3,'TimerFcn',@(~,~)drive(action));
w = timer('StartDelay',60,'TimerFcn',@(~,~)watchdog());
start(t); start(w);
try
    [A,C] = fn();
catch ME
    say(LOG,'ERROR calling dialog: %s',ME.message);
    A = []; C = -1;
end
stop(t); delete(t); stop(w); delete(w);
% Regression: the dialog used to stay on screen after OK because an
% onCleanup object sat inside a figure->callback->workspace->figure cycle.
leftover = findall(groot,'Type','figure');
check(LOG,isempty(leftover),sprintf('no figure left open after dialog returned (saw %d)',numel(leftover)));
delete(leftover);
end

function watchdog()
LOG = getappdata(groot,'smokeLog');
f = findall(groot,'Type','figure');
if ~isempty(f)
    say(LOG,'ERROR watchdog fired, force-closing');
    uiresume(f(1)); delete(f);
end
end

function drive(action)
LOG = getappdata(groot,'smokeLog');
f = [];
try
    % On the first use MATLAB's web-window framework brings up an extra
    % figure of its own, so pick the one that actually holds our table
    % rather than whichever findall lists first.
    tbl = [];
    for attempt = 1:60
        figs = findall(groot,'Type','figure');
        for k = 1:numel(figs)
            candidate = findall(figs(k),'Type','uitable');
            if isscalar(candidate); f = figs(k); tbl = candidate; break; end
        end
        if ~isempty(tbl); break; end
        pause(0.5); drawnow limitrate;
    end
    if isempty(tbl)
        % Fall through rather than returning, so the resume block below
        % still runs and the watchdog does not have to fire.
        say(LOG,'ERROR no figure with a uitable (saw %d figures)',numel(figs));
        if ~isempty(figs); f = figs(1); end
    else
        inspect(LOG,f,tbl,action);
    end
catch ME
    say(LOG,'ERROR in drive: %s (%s)',ME.message,ME.identifier);
    if ~isempty(ME.stack)
        say(LOG,'   at %s line %d',ME.stack(1).name,ME.stack(1).line);
    end
end
% Always resume, whatever happened, so the test cannot hang.
if ~isempty(f) && isvalid(f)
    try
        press(f,ternary(strcmp(action,'ok'),'OK','Cancel'));
    catch
        uiresume(f);
    end
end
end

function inspect(LOG,f,tbl,action)
check(LOG,isscalar(tbl),'exactly one uitable');
D = tbl.Data;

scr = get(groot,'ScreenSize');
say(LOG,'figure %gx%g px (screen %gx%g)',f.Position(3),f.Position(4),scr(3),scr(4));
say(LOG,'WindowStyle=%s Visible=%s',f.WindowStyle,f.Visible);
say(LOG,'Data %s %dx%d',class(D),height(D),width(D));
say(LOG,'col classes: %s',strjoin(cellfun(@(v)class(D.(v)),...
    D.Properties.VariableNames,'UniformOutput',false),', '));
say(LOG,'ColumnEditable %s',mat2str(tbl.ColumnEditable));
say(LOG,'SelectionType=%s Multiselect=%s',tbl.SelectionType,tbl.Multiselect);

check(LOG,f.Position(3) <= 0.9*scr(3)+1,'figure within 90% of screen width');
check(LOG,f.Position(4) <= 0.85*scr(4)+1,'figure within 85% of screen height');

bCol = findButton(f,'Fill column from selected cell');
bSel = findButton(f,'Fill selected rows');
check(LOG,strcmp(bCol.Enable,'off'),'fill buttons start disabled');

if ~strcmp(action,'ok'); return; end

if ismember('Paradigm',D.Properties.VariableNames)
    check(LOG,isa(D.Paradigm,'categorical'),'Paradigm is categorical (drop-down)');
    check(LOG,isa(D.redStim,'logical'),'redStim is logical (check box)');

    % Column order and item lists differ between the wrappers, so resolve
    % both from the live table rather than hardcoding them.
    vars = D.Properties.VariableNames;
    cPar = find(strcmp(vars,'Paradigm'));
    cRed = find(strcmp(vars,'redStim'));
    cats = categories(D.Paradigm);
    last = cats{end};
    say(LOG,'Paradigm col %d (%d items, filling "%s"), redStim col %d',...
        cPar,numel(cats),last,cRed);

    D.Paradigm(3) = last; tbl.Data = D;
    tbl.Selection = [3 cPar];
    tbl.SelectionChangedFcn(tbl,[]);
    check(LOG,strcmp(bCol.Enable,'on'),'fill column enables on cell selection');
    check(LOG,strcmp(bSel.Enable,'off'),'fill rows stays off for one row');
    bCol.ButtonPushedFcn(bCol,[]);
    check(LOG,all(tbl.Data.Paradigm==last),'fill column applied to all rows');

    tbl.Data.redStim(:) = false;
    tbl.Data.redStim(1) = true;
    tbl.Selection = [1 cRed; 2 cRed];
    tbl.SelectionChangedFcn(tbl,[]);
    check(LOG,strcmp(bSel.Enable,'on'),'fill rows enables on multi-row selection');
    bSel.ButtonPushedFcn(bSel,[]);
    check(LOG,tbl.Data.redStim(2),'fill rows wrote the selected row');
    check(LOG,~tbl.Data.redStim(5),'fill rows spared unselected rows');
end

tbl.Selection = [1 1];
tbl.SelectionChangedFcn(tbl,[]);
check(LOG,strcmp(bCol.Enable,'off'),'read-only Session column is not fillable');
end

function b = findButton(f,label)
b = findall(f,'Type','uibutton');
b = b(strcmp({b.Text},label));
if ~isscalar(b); error('smoke:button','button "%s" not found',label); end
end

function press(f,label)
b = findButton(f,label);
b.ButtonPushedFcn(b,[]);
end

function say(LOG,fmt,varargin)
fid = fopen(LOG,'a'); fprintf(fid,[fmt '\n'],varargin{:}); fclose(fid);
end

function check(LOG,tf,what)
if tf; say(LOG,'  ok   %s',what); else; say(LOG,'  FAIL %s',what); end
end

function v = ternary(c,a,b)
if c; v = a; else; v = b; end
end
