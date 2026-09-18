function [labjack, spikeGLX, livePlot, config] = inputLabjackRecordingConfig(samplerate, configPath)
% inputLabjackRecordingConfig  Select and edit LabJack recording settings.

if nargin < 1 || isempty(samplerate); samplerate = 2000; end
if nargin < 2 || isempty(configPath)
    neuroDAPDir = fileparts(fileparts(mfilename('fullpath')));
    configPath = fullfile(neuroDAPDir,'Labjack','labjack-configs.json');
end

configs = readConfigFiles(configPath);
if isempty(configs); configs = defaultConfigs(); end

channelScanIdx = [1 2 5 7]; % AIN0, AIN1, AIN10, optional AIN9
channelLabels = {'1: AIN0','2: AIN1','3: AIN10','4: AIN9'};
answer = [];
okPressed = false;

fig = dialog('Name','LabJack recording config','WindowStyle','modal', ...
    'Units','pixels','Position',[100 100 600 460],'Resize','off');
movegui(fig,'center');
set(fig,'CloseRequestFcn',@cancelDialog);

labelX = 75;
labelW = 105;
controlX = 190;
rightEdge = 440;
spikeW = 80;

uicontrol(fig,'Style','text','String','Session name','HorizontalAlignment','right', ...
    'Position',[labelX 420 labelW 20]);
sessionEdit = uicontrol(fig,'Style','edit','HorizontalAlignment','left', ...
    'Position',[controlX 417 rightEdge-controlX 26],'Tag','sessionName');

uicontrol(fig,'Style','text','String','Animal settings','HorizontalAlignment','right', ...
    'Position',[labelX 385 labelW 20]);
animalPopup = uicontrol(fig,'Style','popupmenu','String',{configs.animal}, ...
    'Position',[controlX 382 150 24],'Callback',@selectConfig);
spikeCheck = uicontrol(fig,'Style','checkbox','String','SpikeGLX', ...
    'Position',[rightEdge-spikeW 382 spikeW 24]);

uicontrol(fig,'Style','text','String','Channel','HorizontalAlignment','left', ...
    'Position',[35 345 75 18]);
uicontrol(fig,'Style','text','String','Record','HorizontalAlignment','center', ...
    'Position',[115 345 65 18]);
uicontrol(fig,'Style','text','String','Name','HorizontalAlignment','left', ...
    'Position',[195 345 120 18]);
uicontrol(fig,'Style','text','String','Freq mod','HorizontalAlignment','center', ...
    'Position',[335 345 65 18]);
uicontrol(fig,'Style','text','String','Display','HorizontalAlignment','center', ...
    'Position',[420 345 65 18]);

recordCheck = gobjects(1,4);
nameEdit = gobjects(1,4);
freqCheck = gobjects(1,4);
displayCheck = gobjects(1,4);
for i = 1:4
    y = 345 - i*45;
    uicontrol(fig,'Style','text','String',channelLabels{i}, ...
        'HorizontalAlignment','left','Position',[35 y+4 75 20]);
    recordCheck(i) = uicontrol(fig,'Style','checkbox','Position',[137 y+4 24 24], ...
        'Tag',sprintf('record%d',i));
    nameEdit(i) = uicontrol(fig,'Style','edit','HorizontalAlignment','left', ...
        'Position',[195 y 120 26]);
    freqCheck(i) = uicontrol(fig,'Style','checkbox','Position',[358 y+4 24 24], ...
        'Tag',sprintf('freqMod%d',i));
    displayCheck(i) = uicontrol(fig,'Style','checkbox','Position',[442 y+4 24 24], ...
        'Tag',sprintf('display%d',i));
end

% Each BNC split carries one physical waveform. Both checkboxes remain
% editable, and changing either checkbox updates its shared-DAC partner.
sharedPairs = [1 3; 2 4];
for pairIdx = 1:size(sharedPairs,1)
    callbackPair = sharedPairs(pairIdx,:);
    for sharedChannel = callbackPair
        set(freqCheck(sharedChannel),'Callback',@(src,~) syncSharedMode(src,callbackPair), ...
            'TooltipString','Shared DAC: changing either checkbox updates both.');
    end
end
set(recordCheck(4),'Callback',@selectAIN9Recording);
set(displayCheck(4),'Callback',@selectAIN9Display);
uicontrol(fig,'Style','text','HorizontalAlignment','left', ...
    'Position',[35 65 530 85],'String', ...
    sprintf(['Channels 1 + 3 share DAC0; channels 2 + 4 share DAC1.\n' ...
    'Checking or unchecking either Freq mod box updates both in its pair.\n' ...
    'Recording either input powers BOTH connected LEDs.\n' ...
    'Record selects input data; it cannot switch the split LEDs separately.']));

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
    livePlot = struct('enable',false,'channelIdx',[],'display',false(1,4));
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
        for c = 1:4
            set(recordCheck(c),'Value',cfg.record(c));
            set(nameEdit(c),'String',cfg.channelNames{c});
            set(freqCheck(c),'Value',cfg.freqMod(c));
            set(displayCheck(c),'Value',cfg.display(c));
        end
        % For older presets with different values, either true selects both.
        for p = 1:size(sharedPairs,1)
            pair = sharedPairs(p,:);
            set(freqCheck(pair),'Value',any(cfg.freqMod(pair)));
        end
        selectAIN9Recording();
    end

    function okDialog(~,~)
        idx = get(animalPopup,'Value');
        cfg = configs(idx);
        cfg.sessionName = strtrim(get(sessionEdit,'String'));
        for c = 1:4
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
            errordlg('Please enter a name for all four channels.', ...
                'Missing channel name','modal');
            return
        end

        if cfg.freqMod(2) ~= cfg.freqMod(4)
            errordlg('Channels 2 and 4 share DAC1 and must use the same modulation setting.', ...
                'Conflicting DAC1 settings','modal');
            return
        end
        if cfg.display(4) && ~cfg.record(4)
            errordlg('Select Record for AIN9 before displaying it.', ...
                'AIN9 is not recorded','modal');
            return
        end

        if cfg.freqMod(1) ~= cfg.freqMod(3)
            errordlg('Channels 1 and 3 share DAC0 and must use the same modulation setting.', ...
                'Conflicting DAC0 settings','modal');
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

    function syncSharedMode(src,pair)
        set(freqCheck(pair),'Value',get(src,'Value'));
    end

    function selectAIN9Recording(varargin)
        if ~get(recordCheck(4),'Value'); set(displayCheck(4),'Value',0); end
    end

    function selectAIN9Display(varargin)
        if get(displayCheck(4),'Value'); set(recordCheck(4),'Value',1); end
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
    for i = 1:min(numel(channels),4)
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
    'channelNames',{{'NAc-left','NAc-right','PMT','AIN9'}}, ...
    'freqMod',false(1,4),'display',false(1,4),'record',[true true false false], ...
    'spikeGLX',true,'sourceFile','');
end

function configs = defaultConfigs()
configs = blankConfig();
configs.display = [true true false false];
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
