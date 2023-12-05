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
%}

% concatLabjack('C:\Shun\Recordings\20231107-test-25uW_g0',plot=true);

clear; close all
addpath(addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Matlab')));

%% Set up

% Set frequency modulation and sampling rate
samplerate = 2000;
labjack.samplerate = samplerate;

% Define name
% labjack.name = {'NAc','LHb','PMT'}; labjack.record = [1,1,0];
labjack.name = {'NAc','LHb','PMT'}; labjack.record = [1,0,0];
%labjack.name = {'Lhemi','Rhemi','PMT'}; labjack.record = [1,1,0];
labjack.nSignals = sum(labjack.record);

% Define modulation params
mod = true;
if mod; labjack.mod = [1,1,0]; 
else; labjack.mod = [0,0,0]; end
labjack.modFreq = [200,250,nan]; % NAc green

% Calculate LED power scalar: outputPower = powerAtKnob * (0.5*LEDpowerblue)
LEDpower1 = 0.75; %1.5;%0.5; % power to get 30uW
LEDpower2 = 2.0; %2.5 % power to get 30uW
LEDpower1Min = 0.3; %0.3 %0.5 % power to get minimal signal 
LEDpower2Min = 0.2; % power to get minimal signal
looplength = samplerate*ones(size(labjack.modFreq)) ./ labjack.modFreq; % 200, 250Hz

% Define mod frequency power
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
root_path = 'C:\Shun\Recordings';
session_path = strcat(root_path,'\',session,'\Photometry\');
script_used = mfilename();
mkdir(session_path);
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
scansPerRead = samplerate;  % Scans returned by eStreamRead call

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

    % The number of eStreamRead calls to perform in the stream read loop
    maxRequests = 7200; %1 per sec; 7200s = 120min = 2h
    disp(['Performing ' num2str(maxRequests) ' stream reads.'])

    readRequest = 1;
    curSkippedSamples = 0;
    totalSkippedSamples = 0;

    for readRequest = 1:maxRequests
        [~,devScanBL,ljmScanBL] = LabJack.LJM.eStreamRead(handle,aData,0,0);

        % save each scan data as separate .mat files
        temp = aData.double;
        save(strcat(session_path,sprintf('Raw_%d.mat',readRequest+1000)),'temp');

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
