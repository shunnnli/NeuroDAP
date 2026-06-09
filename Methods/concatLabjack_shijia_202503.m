%% 2025/03
%% this is the 3rd version for photometry analysis, 2nd version for concatLabjack
%% in the 2nd/last version, options.record double = [1,1], now it's 3 values

function concatLabjack_shijia_202503(sessionpath,options)

% Extract params and concatenate Raw_*.mat files from a given directory
arguments
    sessionpath % path to photometry

    options.plot logical = false % Plot summary figures
    options.save logical = true % save as mat files
    options.record double = [1,1,0] % Recorded channels
    options.rebuildInfo logical = false % rebuild info (old version of the system)
end

disp("*********************Ongoing: concatLabjack: ver 202503*********************");

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
    labjack.name = {'Green','Iso'};
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

% Store sync pulse
sync_labjack = output(mod(1:totalLen,numChannels)==5);  % sync pulse was the 5th channel collected
labjack.sync = sync_labjack;

% Initialize data matrix
labjack.raw = nan(labjack.nSignals,length(sync_labjack));
labjack.modulation = nan(labjack.nSignals,length(sync_labjack));

%% Log signal to corresponding row 
%% Case isosebestic: (Aln0 demodulated with carrier at DAC0, Aln0 demodulated with carrier at DAC1)
row = 0;
for i = 1:length(labjack.name)
    if labjack.record(i)
        row = row + 1;
        % if contains(labjack.name{i}, "Placeholder", 'IgnoreCase', true)
        %     labjack.raw(row,:) = output(mod(1:totalLen, numChannels)==5);
        %     labjack.modulation(row,:) = output(mod(1:totalLen, numChannels)==6);
        % else
            % For the 'Iso' channel, use AIN0 raw data (channel index 1)
            if strcmp(labjack.name{i}, 'Iso')
                labjack.raw(row,:) = output(mod(1:totalLen, numChannels)==1);
            else
                labjack.raw(row,:) = output(mod(1:totalLen, numChannels)==i);
            end
            % Modulation remains the same:
            labjack.modulation(row,:) = output(mod(1:totalLen, numChannels)==i+2);
        % end
    end
end

%% older, working versions
% Case non-isosebestic: (Aln0 demodulated with carrier at DAC0, Aln1 demodulated with carrier at DAC1)
% row = 0;
% for i = 1:size(labjack.name,2)
%     if labjack.record(i)
%         row = row + 1;
%         if contains(labjack.name{i},"Placeholder",IgnoreCase=true)
%             labjack.raw(row,:) = output(mod(1:totalLen,numChannels)==5);
%             labjack.modulation(row,:) = output(mod(1:totalLen,numChannels)==6);
%         else
%             labjack.raw(row,:) = output(mod(1:totalLen,numChannels)==i);
%             labjack.modulation(row,:) = output(mod(1:totalLen,numChannels)==i+2);
%         end
%     end
% end


%% Load all data
% numChannels = length(temp)/labjack.samplerate;
% output = zeros(1,(length(D)*length(temp)));
% for i = 2:length(D)
%     load(strcat(pathPhotometry,filename{i}));
%     output(((i-1)*length(temp)+1):(i*length(temp))) = temp;
% end
% totalLen = length(output);
% 
% % Store digital pulse
% sync_labjack = output(mod(1:totalLen,numChannels)==5);  
% labjack.sync = sync_labjack;
% 
% % clean cue so that cue doesn't have the 30kHz oscillation
% cue_labjack = output(mod(1:totalLen,numChannels)==6);
% cue_labjack(cue_labjack>1000)=0;
% temp = [false, diff(cue_labjack)];
% cue_labjack_allRisingEdges = (temp==1);
% risingEdges = find(cue_labjack_allRisingEdges);
% a = [max(diff(risingEdges)),diff(risingEdges)];
% risingEdges2 = risingEdges(a>(max(diff(risingEdges))/10));
% cue_labjack_use = zeros(size(cue_labjack));
% cue_labjack_use(risingEdges2) = 1;
% labjack.cue = cue_labjack_use;
% 
% % get licks 
% lick_labjack = output(mod(1:totalLen,numChannels)==7);  
% temp = [false, diff(lick_labjack)];
% lick_labjack = (temp==-1); 
% labjack.lick = lick_labjack;
% lick_labjack_tmp = lick_labjack;
% 
% % % clean licks with ITI 90ms - still do it in getTrialTable
% lickFallingEdges = find(lick_labjack);
% previousLickOnset = lickFallingEdges(1);
% lickITICutoffArduino = 90;
% lickFallingEdges_cleaned = [previousLickOnset];
% for p = 1:length(lickFallingEdges)
%     if lickFallingEdges(p)-previousLickOnset>(lickITICutoffArduino/1000*labjack.samplerate)
%         previousLickOnset = lickFallingEdges(p);
%         lickFallingEdges_cleaned = [lickFallingEdges_cleaned previousLickOnset];
%     end
% end
% lick_labjack_zeroed = zeros(size(lick_labjack));
% lick_labjack_zeroed(lickFallingEdges_cleaned)=1;
% labjack.lick = lick_labjack_zeroed;
% % figure;plot(lick_labjack_tmp(25000:30000),'b');hold on;plot(lick_labjack_zeroed(25000:30000),'r');
% 
% % solenoid
% solenoid_labjack = output(mod(1:totalLen,numChannels)==0);
% temp = [false, diff(solenoid_labjack)];
% solenoid_labjack = (temp==1);
% labjack.solenoid = solenoid_labjack;
% 
% % Initialize data matrix
% labjack.raw = nan(labjack.nSignals,length(sync_labjack));
% labjack.modulation = nan(labjack.nSignals,length(sync_labjack));
% 
% % Log signal to corresponding row
% for i = 1:size(labjack.name,2)
%     if labjack.record(i)
% %         if labjack.isosbestic(i)
% %             labjack.raw(i,:) = labjack.raw(i-1,:);
% %             labjack.modulation(i,:) = output(mod(1:totalLen,numChannels)==i+2);
% %         else
%             labjack.raw(i,:) = output(mod(1:totalLen,numChannels)==i);
%             labjack.modulation(i,:) = output(mod(1:totalLen,numChannels)==i+2);
% %         end
%     end
% end

% figure;
% plot(solenoid_labjack); hold on; plot(cue_labjack); hold on;
% plot(labjack.lick);%hold on; plot(lick_labjack_tmp);
% xlim([20000 80000]); title('example cue and solenoid events');
% cd('\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Shijia\BE28 BE31 Task_photometry\BE29 M164-167\20231205_M167L_4W40O40C10D10_g0');
% saveas(gcf, strcat(extractBefore(sessionName,'_g'),'_ExampleTrialEventRisingEdges.tif'));
%% Plot photometry summary plot
if options.plot
    initializeFig(0.67,0.67); tiledlayout(labjack.nSignals*2 + 1,1);
    for i = 1:size(labjack.raw,1)
        nexttile; plot(labjack.raw(i,:)); title(strcat(labjack.name{i},' (raw)')); box off
        nexttile; plot(labjack.modulation(i,:)); title(strcat(labjack.name{i},' (LED driver copy)')); box off
    end
        nexttile; plot(sync_labjack); title('sync'); box off
        % nexttile; plot(lick_labjack); title('sync'); box off
        % nexttile; plot(solenoid_labjack); title('sync'); box off
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
1;
end