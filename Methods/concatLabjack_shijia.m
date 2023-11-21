function concatLabjack_shijia(sessionpath,options)

% Extract params and concatenate Raw_*.mat files from a given directory
arguments
    sessionpath % path to photometry

    options.plot logical = false % Plot summary figures
    options.save logical = true % save as mat files
    options.record double = [1,1,0] % Recorded channels
    options.rebuildInfo logical = false % rebuild info (old version of the system)
end

% Load info.mat
pathPhotometry = strcat(sessionpath,filesep,'Photometry',filesep);
% pathPhotometry = strcat(sessionpath,filesep);

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
    labjack.record = [1,1,0];
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

    if any(labjack.mod); labjack.modFreq = [167,223];
    else; labjack.modFreq = [nan,nan]; end
end

% Add labjack.record if it doesn't exist already (recordings pre 11/12/2023)
if ~isfield(labjack,'record')
    labjack.record = options.record; 
    labjack.nSignals = sum(labjack.record);
end
% If input is different from labjack.record
if sum(labjack.record == options.record) ~= 3
    disp(['labjack.record: ',labjack.record]);
    disp(['options.recordLJ: ',options.record]);
    warning("labjack.record does not agree with recordLJ, reload using labjack.record"); 
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

% Store digital pulse
sync_labjack = output(mod(1:totalLen,numChannels)==5);  
labjack.sync = sync_labjack;

cue_labjack = output(mod(1:totalLen,numChannels)==6); 
temp = [false, diff(cue_labjack)];
cue_labjack = (temp==1);
labjack.cue = cue_labjack;

lick_labjack = output(mod(1:totalLen,numChannels)==7);  
temp = [false, diff(lick_labjack)];
lick_labjack = (temp==-1);
labjack.lick = lick_labjack;

solenoid_labjack = output(mod(1:totalLen,numChannels)==0);
temp = [false, diff(solenoid_labjack)];
solenoid_labjack = (temp==1);
labjack.solenoid = solenoid_labjack;

% Initialize data matrix
labjack.raw = nan(labjack.nSignals,length(sync_labjack));
labjack.modulation = nan(labjack.nSignals,length(sync_labjack));

% Log signal to corresponding row
for i = 1:size(labjack.name,2)
    if labjack.record(i)
        if contains(labjack.name{i},"PMT",IgnoreCase=true)
            labjack.raw(3,:) = output(mod(1:totalLen,numChannels)==5);
            labjack.modulation(3,:) = output(mod(1:totalLen,numChannels)==6);
        else
            labjack.raw(i,:) = output(mod(1:totalLen,numChannels)==i);
            labjack.modulation(i,:) = output(mod(1:totalLen,numChannels)==i+2);
        end
    end
end

%% Plot photometry summary plot (skipped)
if options.plot
    initializeFig(0.67,0.67); tiledlayout(labjack.nSignals*2 + 3,1);
    for i = 1:size(labjack.raw,1)
        nexttile; plot(labjack.raw(i,:)); title(strcat(labjack.name{i},'(raw)')); box off
        nexttile; plot(labjack.modulation(i,:)); title(strcat(labjack.name{i},'(mod)')); box off
    end
        nexttile; plot(cue_labjack); title('sync'); box off
        nexttile; plot(lick_labjack); title('sync'); box off
        nexttile; plot(solenoid_labjack); title('sync'); box off
    saveas(gcf,strcat(sessionpath,filesep,'Summary_labjack_raw.fig'));
end

%% Save concat files
labjack.name(find(~labjack.record)) = [];
labjack.numChannels = numChannels;
labjack.totalLen = totalLen;
labjack.options = options;

if options.save
    save(strcat(sessionpath,filesep,'data_labjack'),'labjack','sync_labjack','-v7.3');
end

end