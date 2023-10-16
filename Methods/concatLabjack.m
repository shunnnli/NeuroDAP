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
numChannels = length(temp)/samplerate;
if ~exist('freqMod','var')
    freqMod = true; 
    disp(['     Variable "freqMod" not found, set to ',num2str(freqMod)]);
end
if ~exist('modFreqGreen','var') && freqMod
    modFreqGreen = 167; 
    disp(['     Variable "modFreqGreen" not found, set to ',num2str(modFreqGreen),' Hz']);
else; modFreqGreen = nan;
end
if ~exist('modFreqRed','var') && freqMod
    modFreqRed = 223; 
    disp(['     Variable "modFreqRed" not found, set to ',num2str(modFreqRed),' Hz']);
else; modFreqRed = nan;
end

% Load all data
output = zeros(1,(length(D)*length(temp)));
for i = 1:length(D)
    load(strcat(pathPhotometry,filename{i}));
    output(((i-1)*length(temp)+1):(i*length(temp))) = temp;
end
totalLen = length(output);

% Make green and red arrays and demodulate photodiode signals  
sync_labjack = output(mod(1:totalLen,numChannels)==0);  % sync pulse
rawGreen = output(mod(1:totalLen,numChannels)==1);      % photodetector #1 -green
rawRed = output(mod(1:totalLen,numChannels)==2);        % photodetector #2 -red
modGreen = output(mod(1:totalLen,numChannels)==3);      % copy of green modulation
modRed = output(mod(1:totalLen,numChannels)==4);        % copy of red modulation
rawGreen2 = output(mod(1:totalLen,numChannels)==5);     % PMT -green
modGreen2 = output(mod(1:totalLen,numChannels)==6);     % PMT -green galvo (unfinished)

% Plot photometry summary plot (skipped)
if options.plot
    initializeFig(0.67,0.67); tiledlayout(numChannels,1);
    nexttile; plot(sync_labjack); title('sync'); box off
    nexttile; plot(rawGreen); title('green raw'); box off
    nexttile; plot(rawRed); title('red raw'); box off
    nexttile; plot(rawGreen2); title('green2 raw'); box off
    nexttile; plot(modGreen); title('green modulation'); box off
    nexttile; plot(modRed); title('red modulation'); box off
    nexttile; plot(modGreen2); title('green2 modulation'); box off
    saveas(gcf,strcat(sessionpath,filesep,'photometryLJ_raw.fig'));
end

% Save concat files
if options.save
    save(strcat(sessionpath,filesep,'photometryLJ_raw'),...
        'samplerate','freqMod','modFreqGreen','modFreqRed','numChannels',...
        'sync_labjack','rawGreen','modGreen','rawRed','modRed','-v7.3');
end

disp('Finished: reading photometry data');
end