function [labjack, spikeGLX, livePlot, config] = inputLabjackRecordingConfig(samplerate, configPath)
% inputLabjackRecordingConfig  Select and edit LabJack recording settings.

if nargin < 1 || isempty(samplerate); samplerate = 2000; end
if nargin < 2 || isempty(configPath)
    neuroDAPDir = fileparts(fileparts(mfilename('fullpath')));
    configPath = fullfile(neuroDAPDir,'Labjack','labjack-configs.json');
end

configs = readConfigFiles(configPath);
if isempty(configs); configs = defaultConfigs(); end

channelScanIdx = [1 2 5]; % AIN0, AIN1, AIN10 (PMT)
answer = [];
okPressed = false;

fig = dialog('Name','LabJack recording config','WindowStyle','modal', ...
    'Units','pixels','Position',[100 100 560 365],'Resize','off');
movegui(fig,'center');
set(fig,'CloseRequestFcn',@cancelDialog);

labelX = 105;
labelW = 105;
controlX = 220;
rightEdge = 470;
spikeW = 80;

uicontrol(fig,'Style','text','String','Session name','HorizontalAlignment','right', ...
    'Position',[labelX 325 labelW 20]);
sessionEdit = uicontrol(fig,'Style','edit','HorizontalAlignment','left', ...
    'Position',[controlX 322 rightEdge-controlX 26]);

uicontrol(fig,'Style','text','String','Animal settings','HorizontalAlignment','right', ...
    'Position',[labelX 290 labelW 20]);
animalPopup = uicontrol(fig,'Style','popupmenu','String',{configs.animal}, ...
    'Position',[controlX 287 150 24],'Callback',@selectConfig);
spikeCheck = uicontrol(fig,'Style','checkbox','String','SpikeGLX', ...
    'Position',[rightEdge-spikeW 287 spikeW 24]);

uicontrol(fig,'Style','text','String','Channel','HorizontalAlignment','left', ...
    'Position',[35 250 75 18]);
uicontrol(fig,'Style','text','String','Record','HorizontalAlignment','center', ...
    'Position',[115 250 65 18]);
uicontrol(fig,'Style','text','String','Name','HorizontalAlignment','left', ...
    'Position',[195 250 120 18]);
uicontrol(fig,'Style','text','String','Freq mod','HorizontalAlignment','center', ...
    'Position',[335 250 65 18]);
uicontrol(fig,'Style','text','String','Display','HorizontalAlignment','center', ...
    'Position',[420 250 65 18]);

recordCheck = gobjects(1,3);
nameEdit = gobjects(1,3);
freqCheck = gobjects(1,3);
displayCheck = gobjects(1,3);
for i = 1:3
    y = 250 - i*45;
    uicontrol(fig,'Style','text','String',sprintf('Channel %d',i), ...
        'HorizontalAlignment','left','Position',[35 y+4 75 20]);
    recordCheck(i) = uicontrol(fig,'Style','checkbox','Position',[137 y+4 24 24]);
    nameEdit(i) = uicontrol(fig,'Style','edit','HorizontalAlignment','left', ...
        'Position',[195 y 120 26]);
    freqCheck(i) = uicontrol(fig,'Style','checkbox','Position',[358 y+4 24 24]);
    displayCheck(i) = uicontrol(fig,'Style','checkbox','Position',[442 y+4 24 24]);
end

uicontrol(fig,'Style','pushbutton','String','OK','Position',[315 20 80 30], ...
    'Callback',@okDialog);
uicontrol(fig,'Style','pushbutton','String','Cancel','Position',[410 20 80 30], ...
    'Callback',@cancelDialog);

applyConfig(1);
uiwait(fig);
if ishandle(fig); delete(fig); end

if ~okPressed
    labjack = [];
    spikeGLX = false;
    livePlot = struct('enable',false,'channelIdx',[],'display',false(1,3));
    config = [];
    return
end

config = answer;
spikeGLX = config.spikeGLX;

labjack.name = config.channelNames;
labjack.record = config.record;
labjack.mod = config.freqMod;
labjack.samplerate = samplerate;
labjack.nSignals = sum(labjack.record);
labjack.display = config.display;

livePlot.display = config.display;
livePlot.channelIdx = channelScanIdx(config.display);
livePlot.enable = ~isempty(livePlot.channelIdx);

    function selectConfig(src,~)
        applyConfig(get(src,'Value'));
    end

    function applyConfig(idx)
        cfg = configs(idx);
        set(spikeCheck,'Value',cfg.spikeGLX);
        for c = 1:3
            set(recordCheck(c),'Value',cfg.record(c));
            set(nameEdit(c),'String',cfg.channelNames{c});
            set(freqCheck(c),'Value',cfg.freqMod(c));
            set(displayCheck(c),'Value',cfg.display(c));
        end
    end

    function okDialog(~,~)
        idx = get(animalPopup,'Value');
        cfg = configs(idx);
        cfg.sessionName = strtrim(get(sessionEdit,'String'));
        for c = 1:3
            cfg.record(c) = logical(get(recordCheck(c),'Value'));
            cfg.channelNames{c} = strtrim(get(nameEdit(c),'String'));
            cfg.freqMod(c) = logical(get(freqCheck(c),'Value'));
            cfg.display(c) = logical(get(displayCheck(c),'Value'));
        end
        cfg.spikeGLX = logical(get(spikeCheck,'Value'));

        if isempty(cfg.sessionName)
            errordlg('Please enter a session name.', ...
                'Missing session name','modal');
            return
        end

        if any(cellfun(@isempty,cfg.channelNames))
            errordlg('Please enter a name for all three channels.', ...
                'Missing channel name','modal');
            return
        end

        answer = cfg;
        okPressed = true;
        uiresume(fig);
    end

    function cancelDialog(~,~)
        okPressed = false;
        uiresume(fig);
    end
end

function configs = readConfigFiles(configPath)
items = {};
if isfile(configPath)
    files = dir(configPath);
elseif isfolder(configPath)
    files = dir(fullfile(configPath,'*.json'));
else
    files = [];
end

for f = 1:numel(files)
    filePath = fullfile(files(f).folder,files(f).name);
    raw = jsondecode(fileread(filePath));
    if isfield(raw,'configs'); raw = raw.configs; end
    for i = 1:numel(raw)
        cfg = normalizeConfig(raw(i));
        cfg.sourceFile = filePath;
        items{end+1} = cfg; %#ok<AGROW>
    end
end

configs = [items{:}];
end

function cfg = normalizeConfig(raw)
cfg = blankConfig();
cfg.animal = '';
if isfield(raw,'animal'); cfg.animal = char(raw.animal); end
if isfield(raw,'spikeGLX'); cfg.spikeGLX = logical(raw.spikeGLX); end
if isfield(raw,'record'); cfg.record = logicalVector(raw.record,cfg.record); end
if isfield(raw,'freqMod'); cfg.freqMod = logicalVector(raw.freqMod,cfg.freqMod); end
if isfield(raw,'display'); cfg.display = logicalVector(raw.display,cfg.display); end

if isfield(raw,'channelNames')
    cfg.channelNames = stringCell(raw.channelNames,cfg.channelNames);
end
if isfield(raw,'channels')
    channels = raw.channels;
    for i = 1:min(numel(channels),3)
        if isfield(channels(i),'name'); cfg.channelNames{i} = char(channels(i).name); end
        if isfield(channels(i),'freqMod'); cfg.freqMod(i) = logical(channels(i).freqMod); end
        if isfield(channels(i),'display'); cfg.display(i) = logical(channels(i).display); end
        if isfield(channels(i),'record'); cfg.record(i) = logical(channels(i).record); end
    end
end
if isempty(cfg.animal); cfg.animal = 'Default'; end
end

function cfg = blankConfig()
cfg = struct('sessionName','','animal','Default', ...
    'channelNames',{{'NAc-left','NAc-right','PMT'}}, ...
    'freqMod',false(1,3),'display',false(1,3),'record',[true true false], ...
    'spikeGLX',true,'sourceFile','');
end

function configs = defaultConfigs()
configs = blankConfig();
configs.display = [true true false];
end

function out = logicalVector(value, fallback)
out = fallback;
n = min(numel(value),numel(out));
for i = 1:n
    out(i) = logical(getItem(value,i,out(i)));
end
end

function out = stringCell(value, fallback)
out = fallback;
n = min(numel(value),numel(out));
for i = 1:n
    out{i} = char(getItem(value,i,out{i}));
end
end

function item = getItem(value, idx, fallback)
if iscell(value)
    item = value{idx};
elseif numel(value) >= idx
    item = value(idx);
else
    item = fallback;
end
if isempty(item); item = fallback; end
end
