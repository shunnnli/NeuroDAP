%{

Setting up stream-in and stream-out together, then reading
stream-in values using .NET.
for info on LabJack methods for Matlab run >> methodsview(LabJack.LJM)

Creates two sine waves on DAC0 and DAC1 for LED modulation.  
DAC0 = 470nm (167Hz) (0.5 offset = about 70uW, used amplitude = 0.275), DAC1 = 565nm (223Hz)
modify frequency and amplitude individually using local variables amplitude, offset and frequency 
negative peak of sine wave should not go below 0.15V as LED begins to turn off/labjack bottoms out

original script from LackJack: example code on MATLAB for LJM (https://labjack.com/support/software/examples/ljm/matlab)

v1.0 Seul Ah Kim and Mike Wallace, Harvard Medical School 
    - added functionality to save stream-in values as .mat files in the
    local folder
    - enabled stream-out of sinusoidal wave on DAC0 and DAC1

v2.0 Shun Li, Harvard Medical School
    - added GUI for starting experiment
    - changed samplerate to 2000 after 10/13/2023
    - cleaned up non-functional live-stream 

v3.0 Shun Li, Harvard Medical School
    - added live plotting function
%}

% concatLabjack('C:\Shun\Recordings\20231107-test-25uW_g0',plot=true);

clear; close all
% addpath(addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Matlab')));
% root_path = 'C:\Shun\Recordings';
addpath(addpath(genpath('E:\Sally\Matlab')));
root_path = 'E:\Sally\Recordings';

%% Set up

% Set frequency modulation and sampling rate
samplerate = 2000;
labjack.samplerate = samplerate;

% Live plot (incoming photometry) settings
enableLivePlot = true;
plotChanIdx = 1;        % 1=AIN0 (photodiode green NAc), 2=AIN1 (photodiode green LHb), etc.
% plotChanIdx = [1 2];  % uncomment to plot both AIN0 and AIN1
plotWindowSec = 10;     % seconds shown in rolling window
plotUpdateSec = 0.1;    % seconds per update (~100 ms)

% Shun's setting
labjack.name = {'dLight','GCaMP8m','PMT'}; 
labjack.record = [1,0,0];
freqMod = true;
LEDpower1 = 0.8; %1.5;%0.5; % power to get 30uW
LEDpower2 = 2.5; % 2.5=30uW
LEDpower1Min = 0.3; %0.3 %0.5 % power to get minimal signal 
LEDpower2Min = 0.2; % power to get minimal signal

% Define modulation params
if freqMod; labjack.mod = [1,1,0]; 
else; labjack.mod = [0,0,0]; end
labjack.modFreq = [200,250,nan]; % NAc green

% Define mod frequency power
labjack.nSignals = sum(labjack.record);
looplength = samplerate*ones(size(labjack.modFreq)) ./ labjack.modFreq; % 200, 250Hz
labjack.LEDpowers = [LEDpower1,LEDpower2];
labjack.LEDpowersMin = [LEDpower1Min,LEDpower2Min];
labjack.Modpowers1 = getModPower(200,2000,LEDpower1,LEDpower1Min);
labjack.Modpowers2 = getModPower(250,2000,LEDpower2,LEDpower1Min);

% Ask for confirmation
names = sprintf(' %s,',labjack.name{:});
quest = {'Please confirm the following labjack settings are correct. If not, click No/Cancel to exit and re-edit.',...
         '',...
         [sprintf('labjack.name: %s', names(1:end-1))],...
         ['labjack.mod: ', num2str(labjack.mod)],...
         ['labjack.record: ', num2str(labjack.record)],...
         ['labjeck.LEDpowers: ', num2str(labjack.LEDpowers)],...
         ['labjack.LEDpowersMin: ', num2str(labjack.LEDpowersMin)]};
answer = questdlg(quest,'Confirm labjack settings');

if ~strcmp(answer,'Yes'); return; end

%% Save metadata
% Choose folder to save photometry data
session = inputdlg({'Enter session name'},"Session name",[1 40]);
session = strcat(session{1},'_g0');
session_path = strcat(root_path,'\',session,'\Photometry\');
script_used = mfilename();

if ~isfile(strcat(session_path,'info.mat'))
    mkdir(session_path);
else
    error("Error: session_path already have recording!");
end

save(strcat(session_path,'info.mat'),'labjack','samplerate','looplength','script_used');

%% Labjack initialization

% Make the LJM .NET assembly visible in MATLAB
ljmAsm = NET.addAssembly('LabJack.LJM');

% Creating an object to nested class LabJack.LJM.CONSTANTS
t = ljmAsm.AssemblyHandle.GetType('LabJack.LJM+CONSTANTS');
LJM_CONSTANTS = System.Activator.CreateInstance(t);
    
% Any device, Any connection, Any identifier
handle = 0;
[ljmError, handle] = LabJack.LJM.OpenS('ANY', 'ANY', 'ANY', handle);
showDeviceInfo(handle);

% Setup stream-out
numAddressesOut = 2; 
aNamesOut = NET.createArray('System.String', numAddressesOut);
aNamesOut(1) = 'DAC0';  
aNamesOut(2) = 'DAC1';

aAddressesOut = NET.createArray('System.Int32', numAddressesOut);
aTypesOut = NET.createArray('System.Int32', numAddressesOut);  % Dummy
LabJack.LJM.NamesToAddresses(numAddressesOut, aNamesOut, ...
    aAddressesOut, aTypesOut);

% Allocate memory for the stream-out buffer DAC0 (OUT0)
LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_TARGET', aAddressesOut(1));
LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_SIZE', 512);
LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_ENABLE', 1);

% Allocate memory for the stream-out buffer DAC1 (OUT1)
LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_TARGET', aAddressesOut(2));
LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_SIZE', 512);
LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_ENABLE', 1);

% Write values to the stream-out buffer
if any(labjack.mod)
    % OUT0
    powers = labjack.Modpowers1;
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_LOOP_SIZE', looplength(1));
    for n = 1:looplength(1)
        LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', powers(n));
    end
    % OUT1 
    powers = labjack.Modpowers2;
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_LOOP_SIZE', looplength(2));
    for n = 1:looplength(2)
        LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', powers(n));
    end

    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_SET_LOOP', 1);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_SET_LOOP', 1);
    [~, value] = LabJack.LJM.eReadName(handle, 'STREAM_OUT2_BUFFER_STATUS', 0);
    disp(['STREAM_OUT2_BUFFER_STATUS = ' num2str(value)])
else
    % OUT0
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_LOOP_SIZE', looplength(1));
    for n = 1:looplength(1)
        LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', LEDpower1);
    end
    % OUT1
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_LOOP_SIZE', looplength(2));
    for n = 1:looplength(2)
        LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', LEDpower2);
    end  

    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_SET_LOOP', 1);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_SET_LOOP', 1);
    [~, value] = LabJack.LJM.eReadName(handle, 'STREAM_OUT2_BUFFER_STATUS', 0);
    disp(['STREAM_OUT2_BUFFER_STATUS = ' num2str(value)])
end

    
% Stream-in  configuration
% Scan list names to stream-in
numAddressesIn = 7; %5;
aScanListNames = NET.createArray('System.String', numAddressesIn);
aScanListNames(1) = 'AIN0';     % photodiode green NAc channel
aScanListNames(2) = 'AIN1';     % photodiode green LHb channel
aScanListNames(3) = 'AIN2';     % copy of DAC0 out to NAc LED
aScanListNames(4) = 'AIN3';     % copy of DAC1 out to LHb LED
aScanListNames(5) = 'AIN10';    % PMT channel
aScanListNames(6) = 'AIN11';    % PMT galvo copy (not working)
aScanListNames(7) = 'DIO0';     % Sync pulse


% Scan list addresses to stream
aScanList = NET.createArray('System.Int32',(numAddressesIn + numAddressesOut));

% Get stream-in addresses
aTypes = NET.createArray('System.Int32', numAddressesIn); % Dummy
LabJack.LJM.NamesToAddresses(numAddressesIn, aScanListNames, ...
    aScanList, aTypes);

% Add the scan list outputs to the end of the scan list.
% STREAM_OUT0 = 4800, STREAM_OUT1 = 4801, ...
aScanList(numAddressesIn+1) = 4800;  % STREAM_OUT0
aScanList(numAddressesIn+2) = 4801;  % STREAM_OUT1

scanRate = double(round(samplerate));  % Scans per second   % SAK previously 2000, now matches samplerate
scansPerRead = max(1, round(scanRate * plotUpdateSec));  % Scans returned by each eStreamRead call

% Stream reads will be stored in aData 
% Needs to be at least numAddressesIn*scansPerRead in size
aData = NET.createArray('System.Double', numAddressesIn*scansPerRead);

%% Run recording

waitdlg = warndlg('Press OK to confirm NI recording started', 'NI started');
waitfor(waitdlg);
waitdlg = warndlg('Press OK to start labjack recording', 'Start session');
waitfor(waitdlg);
% enddlg = questdlg('Press End to end labjack recording','End session','End','Cancel','Cancel');
% endSession = checkEndButton(enddlg);

try
    % When streaming, negative channels and ranges can be configured for
    % individual analog inputs, but the stream has only one settling time
    % and resolution.

    % Ensure triggered stream is disabled.
    LabJack.LJM.eWriteName(handle, 'STREAM_TRIGGER_INDEX', 0);

    % Enabling internally-clocked stream.
    LabJack.LJM.eWriteName(handle, 'STREAM_CLOCK_SOURCE', 0);

    % All negative channels are single-ended, AIN0 and AIN1 ranges are
    % +/-10 V, stream settling is 0 (default) and stream resolution index
    % is 0 (default).
    numFrames = 9; %7;
    aNames = NET.createArray('System.String', numFrames);
    aNames(1) = 'AIN_ALL_NEGATIVE_CH';
    aNames(2) = 'AIN0_RANGE';
    aNames(3) = 'AIN1_RANGE';
    aNames(4) = 'AIN2_RANGE';
    aNames(5) = 'AIN3_RANGE';
    aNames(6) = 'AIN10_RANGE';
    aNames(7) = 'AIN11_RANGE';
    aNames(8) = 'STREAM_SETTLING_US';
    aNames(9) = 'STREAM_RESOLUTION_INDEX';
    aValues = NET.createArray('System.Double', numFrames);
    aValues(1) = LJM_CONSTANTS.GND;
    aValues(2) = 10.0;
    aValues(3) = 10.0;
    aValues(4) = 10.0;
    aValues(5) = 10.0;
    aValues(6) = 10.0;
    aValues(7) = 10.0;
    aValues(8) = 0;
    aValues(9) = 0;

    % Write the analog inputs' negative channels (when applicable), ranges
    % stream settling time and stream resolution configuration.
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, -1);

    % Configure and start stream
    numAddresses = aScanList.Length;
    [~, scanRate] = LabJack.LJM.eStreamStart(handle, scansPerRead, ...
        numAddresses, aScanList, scanRate);
    disp(['Stream started with a scan rate of ' num2str(scanRate) ' Hz.'])

     % ---------------- Live plot + save buffering ----------------
    % Reading smaller chunks so the live plot updates roughly every plotUpdateSec.
    if enableLivePlot
        chanLabels = {'AIN0','AIN1','AIN2','AIN3','AIN10','AIN11','DIO0'};

        % Create a wide + short window so the trace looks "flat" and the title
        % has enough headroom.
        if exist('initializeFig','file') == 2
            liveFig = initializeFig(0.8, 0.4);
            clf(liveFig,'reset');
        else
            liveFig = figure('Units','normalized','Position',[0.01 0.65 0.98 0.22]);
        end

        set(liveFig, 'Name','Photometry (live)', 'NumberTitle','off', ...
            'Color','w', 'MenuBar','none', 'ToolBar','none');

        liveAx = axes('Parent', liveFig, 'Units','normalized');
        hold(liveAx,'on');

        plotLines = gobjects(1, numel(plotChanIdx));
        for k = 1:numel(plotChanIdx)
            plotLines(k) = plot(liveAx, nan, nan);
        end

        lgd = legend(liveAx, chanLabels(plotChanIdx), 'Interpreter','none','Location','northeast');
        lgd.Box = 'off';
        xlabel(liveAx,'Time (s)'); ylabel(liveAx,'Voltage (V)');
        title(liveAx,'Incoming photometry channel (live)');
        % grid(liveAx,'on');

        % Render once so TightInset is up-to-date, then adjust layout so
        % title/labels are not clipped.
        drawnow;
        setLivePlotAxesLayout(liveFig, liveAx);
        set(liveFig, 'ResizeFcn', @(~,~) setLivePlotAxesLayout(liveFig, liveAx));
        drawnow;
    end

    % Plot buffers (rolling window)
    windowSamples = max(1, round(plotWindowSec * scanRate));
    plotTBuf = nan(1, windowSamples);
    plotYBuf = nan(numel(plotChanIdx), windowSamples);
    plotBufPos = 1;
    plotFilled = 0;
    plotSampleCounter = 0;

    % Save buffer: keep writing ~1 second per file (same as original script)
    saveReadsPerFile = max(1, round(samplerate / scansPerRead));  % e.g. 10 reads for 0.1 s reads
    saveBuf = zeros(1, numAddressesIn * scansPerRead * saveReadsPerFile);
    saveBufPos = 1;
    saveRequest = 0;

    % Total acquisition time (minutes)
    recordingMinutes = 120;
    maxRequests = recordingMinutes*60*saveReadsPerFile;  % number of eStreamRead calls
    disp(['Performing ' num2str(maxRequests) ' stream reads (' num2str(recordingMinutes) ' min total).'])
    readRequest = 1;
    curSkippedSamples = 0;
    totalSkippedSamples = 0;

    for readRequest = 1:maxRequests
        [~,devScanBL,ljmScanBL] = LabJack.LJM.eStreamRead(handle,aData,0,0);

        temp = aData.double;  % length = numAddressesIn * scansPerRead

        % ---- Live plot update (roughly plotUpdateSec) ----
        if enableLivePlot && exist('liveFig','var') && ishandle(liveFig)
            % Reshape to [numChannels x scansPerRead] (each column is one scan)
            chunk = reshape(temp, [numAddressesIn, scansPerRead]);
            yChunk = chunk(plotChanIdx, :);
            tChunk = (plotSampleCounter + (0:scansPerRead-1)) ./ scanRate;
            plotSampleCounter = plotSampleCounter + scansPerRead;

            % Circular buffer write
            idxEnd = plotBufPos + scansPerRead - 1;
            if idxEnd <= windowSamples
                plotTBuf(plotBufPos:idxEnd) = tChunk;
                plotYBuf(:, plotBufPos:idxEnd) = yChunk;
            else
                n1 = windowSamples - plotBufPos + 1;
                plotTBuf(plotBufPos:windowSamples) = tChunk(1:n1);
                plotYBuf(:, plotBufPos:windowSamples) = yChunk(:, 1:n1);
                n2 = scansPerRead - n1;
                plotTBuf(1:n2) = tChunk(n1+1:end);
                plotYBuf(:, 1:n2) = yChunk(:, n1+1:end);
            end

            plotBufPos = mod(plotBufPos + scansPerRead - 1, windowSamples) + 1;
            plotFilled = min(windowSamples, plotFilled + scansPerRead);

            if plotFilled < windowSamples
                ord = 1:plotFilled;
            else
                ord = [plotBufPos:windowSamples 1:plotBufPos-1];
            end

            % Push data to the plot
            for k = 1:numel(plotChanIdx)
                set(plotLines(k), 'XData', plotTBuf(ord), 'YData', plotYBuf(k, ord));
            end
            if ~isnan(plotTBuf(ord(end)))
                xlim(liveAx, [max(0, plotTBuf(ord(end)) - plotWindowSec), plotTBuf(ord(end))]);
            end
            drawnow limitrate nocallbacks;
        end

        % ---- Save buffering: write ~1 file per second (original behavior) ----
        nVals = numAddressesIn * scansPerRead;
        saveBuf(saveBufPos:(saveBufPos + nVals - 1)) = temp;
        saveBufPos = saveBufPos + nVals;

        if saveBufPos > numel(saveBuf)
            saveRequest = saveRequest + 1;
            temp = saveBuf;  
            save(strcat(session_path, sprintf('Raw_%d.mat', saveRequest + 1000)), 'temp');
            saveBufPos = 1;
        end

        % Count the skipped samples which are indicated by -9999
        % values. Skipped samples occur after a device's stream buffer
        % overflows and are reported after auto-recover mode ends.
        % When streaming at faster scan rates in MATLAB, try counting
        % the skipped packets outside your eStreamRead loop if you are
        % getting skipped samples/scan.
        curSkippedSamples = sum(double(aData) == -9999.0);
        totalSkippedSamples = totalSkippedSamples + curSkippedSamples;

        % Update endSession and readRequest
        % endSession = checkEndButton(enddlg); disp(endSession);
        % readRequest = readRequest + 1; while ~endSession %
    end

    % Flush any remaining samples in the save buffer (e.g., if you stop early)
    if exist('saveBufPos','var') && saveBufPos > 1
        saveRequest = saveRequest + 1;
        temp = saveBuf(1:saveBufPos-1);  
        save(strcat(session_path, sprintf('Raw_%d.mat', saveRequest + 1000)), 'temp');
    end

    disp(['Stopped: Perfored ' num2str(readRequest) ' stream reads.']);

    LabJack.LJM.eStreamStop(handle);
    LabJack.LJM.Close(handle);
    LabJack.LJM.CloseAll();

catch e
    showErrorMessage(e);
    LabJack.LJM.eStreamStop(handle);
    LabJack.LJM.Close(handle);
    LabJack.LJM.CloseAll();
    close all
end

% %% test getModpower
% 
% powers = getModPower(200,2000,0.2,0.01);
% plot(powers);

%% Calculate freq mod power

function powers = getModPower(freq,samplerate,LEDpower,LEDpowerMin)
    nPointsPerCycle = round(samplerate / freq);
    powers = [];

    for n = 0:nPointsPerCycle-1
        powers = [powers,cos((2*pi/nPointsPerCycle)*n)];
    end

    halfAmp = min(abs(LEDpower-LEDpowerMin), abs(LEDpower-5));
    powers = rescale(powers,LEDpower-halfAmp,LEDpower+halfAmp);

end


%% Set live plot axes

function setLivePlotAxesLayout(fig, ax)
% setLivePlotAxesLayout  Keep a live-plot axes nicely laid out in a short/wide figure.
%
% Usage:
%   setLivePlotAxesLayout(gcf, gca)
%
% This uses TightInset to leave enough room for title/xlabel/ylabel so they
% aren't clipped, and is intended to be called once after creating labels,
% and again from a Figure ResizeFcn.

    if nargin < 2 || ~isgraphics(fig, 'figure') || ~isgraphics(ax, 'axes')
        return;
    end

    % Work in normalized units so the axes scales with the figure.
    ax.Units = 'normalized';

    % TightInset is only accurate after MATLAB has computed text extents.
    % (ResizeFcn calls already happen with callbacks enabled.)
    drawnow;  % do NOT use 'nocallbacks' here

    ti = ax.TightInset;   % [left bottom right top] in normalized units

    % Extra padding beyond TightInset: tweak to taste.
    pad = [0.010 0.020 0.010 0.040]; % [L B R T]  (extra top keeps title visible)

    left   = ti(1) + pad(1);
    bottom = ti(2) + pad(2);
    width  = 1 - (ti(1) + ti(3) + pad(1) + pad(3));
    height = 1 - (ti(2) + ti(4) + pad(2) + pad(4));

    % Failsafe for very small windows
    width  = max(width,  0.10);
    height = max(height, 0.10);

    ax.Position = [left bottom width height];
end
