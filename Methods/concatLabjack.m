function concatLabjack(sessionpath,options)

% Extract params and concatenate Raw_*.mat files from a given directory
arguments
    sessionpath % path to photometry

    options.plot logical = false % Plot summary figures
    options.save logical = true % save as mat files
end

% Load info.mat
pathPhotometry = strcat(sessionpath,filesep,'Photometry',filesep);
load(strcat(pathPhotometry,'info.mat'));

D = dir(strcat(pathPhotometry,'Raw_*.mat')); 
filename = {D.name}; load(strcat(pathPhotometry,filename{1}));

% Store mod related signal
if ~exist('labjack','var')
    disp("concatLabjack: 'params' not found, build from current info.");

    % Find sample rate
    if exist('samplerate','var'); labjack.samplerate = samplerate;
    else
        warning('Labjack samplerate not found! Use default 2000Hz');
        labjack.samplerate = 2000;
    end

    % Find number of signals
    labjack.nSignals = 2; % NAc green + NAc red
    labjack.name = {"NAc green","NAc red"};
    labjack.mod = zeros(1,labjack.nSignals);
    labjack.modFreq = zeros(1,labjack.nSignals);
    labjack.LEDpower = zeros(1,labjack.nSignals);

    % Load freq mod settings
    if ~exist('freqMod','var')
        labjack.mod = ones(1,labjack.nSiganls);
        disp('     Variable "freqMod" not found, set all to true');
    else
        labjack.mod = ones(1,labjack.nSiganls) * freqMod;
    end

    if any(labjack.mod); labjack.modFreq = [167,223];
    else; labjack.modFreq = [nan,nan]; end
end

% Load all data
numChannels = length(temp)/labjack.samplerate;
output = zeros(1,(length(D)*length(temp)));
for i = 1:length(D)
    load(strcat(pathPhotometry,filename{i}));
    output(((i-1)*length(temp)+1):(i*length(temp))) = temp;
end
totalLen = length(output);

% Make green and red arrays and demodulate photodiode signals  
sync_labjack = output(mod(1:totalLen,numChannels)==0);  % sync pulse
rawGreen = output(mod(1:totalLen,numChannels)==1);      % photodetector #1 -green NAc
rawGreen2 = output(mod(1:totalLen,numChannels)==2);        % photodetector #2 -green LHb
% rawRed = output(mod(1:totalLen,numChannels)==2);        % photodetector #2 -red
modGreen = output(mod(1:totalLen,numChannels)==3);      % copy of green NAc modulation
modGreen2 = output(mod(1:totalLen,numChannels)==4);        % copy of green LHb modulation
% modRed = output(mod(1:totalLen,numChannels)==4);        % copy of red modulation
rawPMT = output(mod(1:totalLen,numChannels)==5);     % PMT -green
modPMT = output(mod(1:totalLen,numChannels)==6);     % PMT -green galvo (unfinished)

% Initialize data matrix
labjack.raw = nan(labjack.nSignals,length(sync_labjack));
labjack.modulation = nan(labjack.nSignals,length(sync_labjack));
labjack.sync = nan(1,length(sync_labjack));
% Log signal to corresponding row
labjack.raw(1,:) = rawGreen;    labjack.modulation(1,:) = modGreen;
labjack.raw(2,:) = rawGreen2;   labjack.modulation(2,:) = modGreen2;
labjack.raw(3,:) = rawPMT;      labjack.modulation(3,:) = modPMT;
labjack.sync = sync_labjack;

% Plot photometry summary plot (skipped)
if options.plot
    initializeFig(0.67,0.67); tiledlayout(numChannels,1);
    nexttile; plot(sync_labjack); title('sync'); box off
    nexttile; plot(rawGreen); title('green NAc raw'); box off
    nexttile; plot(rawGreen2); title('green LHb raw'); box off
    nexttile; plot(rawPMT); title('green PMT raw'); box off
    nexttile; plot(modGreen); title('green NAc modulation'); box off
    nexttile; plot(modGreen2); title('green LHb modulation'); box off
    nexttile; plot(modPMT); title('green PMT modulation'); box off
    saveas(gcf,strcat(sessionpath,filesep,'photometryLJ_raw.fig'));
end

% Save concat files
if options.save
    save(strcat(sessionpath,filesep,'photometryLJ_raw'),...
        'labjack','numChannels','sync_labjack','-v7.3');
end

disp('Finished: reading photometry data');
end