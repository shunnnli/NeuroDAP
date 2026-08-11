function concatLabjack(sessionpath,options)

% Extract params and concatenate Raw_*.mat files from a given directory
arguments
    sessionpath % path to photometry

    options.plot logical = false % Plot summary figures
    options.save logical = true % save as mat files
    options.record double = [1,1,0] % Recorded channels
    options.rebuildInfo logical = false % rebuild info (old version of the system)

    options.followOriginal logical = true % Follow original recorded information

    options.saveDuplicate logical = true % Save a duplicate into Photometry folder
    options.outputPath string
end

if ~isfield(options,'outputPath'); options.outputPath = sessionpath; end

% Load info.mat
pathPhotometry = strcat(sessionpath,filesep,'Photometry',filesep);
load(strcat(pathPhotometry,'info.mat'));

D = dir(strcat(pathPhotometry,'Raw_*.mat')); 
filename = {D.name}; load(strcat(pathPhotometry,filename{1}));

%% Store mod related signal
if ~exist('labjack','var')
    disp("concatLabjack: 'params' not found, build from current info.");
    options.rebuildInfo = true;

    % Find sample rate
    if exist('samplerate','var'); labjack.samplerate = samplerate;
    else
        warning('Labjack samplerate not found! Use default 2000Hz');
        labjack.samplerate = 2000;
    end

    % Define basic parmas
    labjack.record = [1,0,0];
    labjack.nSignals = sum(labjack.record); % NAc green + NAc red
    labjack.name = {'NAc_green','NAc_red','PMT'};
    labjack.mod = zeros(1,labjack.nSignals);
    labjack.modFreq = zeros(1,labjack.nSignals);
    labjack.LEDpower = zeros(1,labjack.nSignals);

    % Load freq mod settings
    if ~exist('freqMod','var')
        labjack.mod = ones(1,labjack.nSignals);
        disp('     Variable "freqMod" not found, set all to true');
    else
        labjack.mod = ones(1,labjack.nSignals) * freqMod;
    end

    if any(labjack.mod); labjack.modFreq = [171,228,0];
    else; labjack.modFreq = [nan,nan,nan]; end
end

% Add labjack.record if it doesn't exist already (recordings pre 11/12/2023)
if ~isfield(labjack,'record')
    labjack.record = options.record; 
    labjack.nSignals = sum(labjack.record);
end
% If input is different from labjack.record
if sum(labjack.record == options.record) ~= 3 
        disp(['labjack.record: ',num2str(labjack.record)]);
        disp(['options.recordLJ: ',num2str(options.record)]);
    if options.followOriginal
        warning("labjack.record does not agree with recordLJ, reload using labjack.record"); 
    else
        labjack.record = options.record;
        labjack.nSignals = sum(labjack.record);
        labjack.mod(find(~labjack.record)) = [];
        labjack.modFreq(find(~labjack.record)) = [];
        warning("labjack.record does not agree with recordLJ, reload using recordLJ"); 
    end
end
% Replace space in name with underscore
for i = 1:labjack.nSignals
    labjack.name{i} = strrep(labjack.name{i}, ' ', '-');
    labjack.name{i} = strrep(labjack.name{i}, '_', '-');
end

%% Load all data
numChannels = length(temp)/labjack.samplerate;
output = zeros(1,(length(D)*length(temp)));
for i = 1:length(D)
    load(strcat(pathPhotometry,filename{i}));
    output(((i-1)*length(temp)+1):(i*length(temp))) = temp;
end
totalLen = length(output);

% Store sync pulse
sync_labjack = output(mod(1:totalLen,numChannels)==0);  % sync pulse
labjack.sync = sync_labjack;

% Initialize data matrix
labjack.raw = nan(labjack.nSignals,length(sync_labjack));
labjack.modulation = nan(labjack.nSignals,length(sync_labjack));

% Map logical channels to their physical positions in the saved scan.
% Channel 3 is always acquired from AIN10 (scan column 5), regardless of
% its display name. Non-PMT channel 3 shares the DAC0 reference on AIN2.
rawScanIdx = [1,2,5];       % AIN0, AIN1, AIN10
modScanIdx = [3,4,6];       % AIN2, AIN3, AIN11 (PMT galvo copy)
if numel(labjack.name) >= 3 && ...
        ~contains(labjack.name{3},"PMT",IgnoreCase=true)
    modScanIdx(3) = 3;      % non-PMT channel 3 uses the DAC0 reference
end

% Log each recorded signal to its corresponding compacted row.
row = 0;
for i = 1:size(labjack.name,2)
    if labjack.record(i)
        row = row + 1;
        labjack.raw(row,:) = output(mod(1:totalLen,numChannels)==rawScanIdx(i));
        labjack.modulation(row,:) = output(mod(1:totalLen,numChannels)==modScanIdx(i));
    end
end

%% Plot photometry summary plot (skipped by default)
if options.plot
    initializeFig(0.67,0.67); tiledlayout(labjack.nSignals*2 + 1,1);
    nexttile; plot(sync_labjack); title('sync'); box off
    for i = 1:size(labjack.raw,1)
        nexttile; plot(labjack.raw(i,:)); title(strcat(labjack.name{i},'(raw)')); box off
        nexttile; plot(labjack.modulation(i,:)); title(strcat(labjack.name{i},'(mod)')); box off
    end
    saveas(gcf,strcat(options.outputPath,filesep,'Summary_labjack_raw.fig'));
end

%% Save concat files
% raw and modulation are stored with one row per recorded signal. Compact
% the corresponding per-channel metadata in the same way so that row i of
% raw uses mod(i), modFreq(i), and name{i} from the same source channel.
recordedIdx = logical(labjack.record);
if numel(labjack.name) == numel(recordedIdx)
    labjack.name = labjack.name(recordedIdx);
end
if numel(labjack.mod) == numel(recordedIdx)
    labjack.mod = labjack.mod(recordedIdx);
end
if numel(labjack.modFreq) == numel(recordedIdx)
    labjack.modFreq = labjack.modFreq(recordedIdx);
end

% These fields are present in newer acquisition metadata and are also
% defined per source channel.
if isfield(labjack,'LEDpowers') && numel(labjack.LEDpowers) == numel(recordedIdx)
    labjack.LEDpowers = labjack.LEDpowers(recordedIdx);
end
if isfield(labjack,'LEDpowersMin') && numel(labjack.LEDpowersMin) == numel(recordedIdx)
    labjack.LEDpowersMin = labjack.LEDpowersMin(recordedIdx);
end
if isfield(labjack,'LEDpower') && numel(labjack.LEDpower) == numel(recordedIdx)
    labjack.LEDpower = labjack.LEDpower(recordedIdx);
end

labjack.numChannels = numChannels;
labjack.totalLen = totalLen;
labjack.options = options;

if options.save
    save(strcat(options.outputPath,filesep,'data_labjack'),'labjack','sync_labjack','-v7.3');
    if options.saveDuplicate
        save(strcat(options.outputPath,filesep,'Photometry',filesep,'data_labjack'),'labjack','sync_labjack','-v7.3');
    end
end

end
