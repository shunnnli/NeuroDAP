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

% Log signal to corresponding row
row = 0;
for i = 1:size(labjack.name,2)
    if labjack.record(i)
        row = row + 1;
        if contains(labjack.name{i},"PMT",IgnoreCase=true)
            labjack.raw(row,:) = output(mod(1:totalLen,numChannels)==5);
            labjack.modulation(row,:) = output(mod(1:totalLen,numChannels)==6);
        else
            labjack.raw(row,:) = output(mod(1:totalLen,numChannels)==i);
            labjack.modulation(row,:) = output(mod(1:totalLen,numChannels)==i+2);
        end
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
labjack.name(find(~labjack.record)) = [];

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