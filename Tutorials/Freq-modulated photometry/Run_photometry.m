% Run_photometry_streamout_SAK.m
%
% Setting up stream-in and stream-out together, then reading
% stream-in values using .NET.
% for info on LabJack methods for Matlab run >> methodsview(LabJack.LJM)
%
% Creates two sine waves on DAC0 and DAC1 for LED modulation.  
% DAC0 = 470nm (167Hz) (0.5 offset = about 70uW, used amplitude = 0.275), DAC1 = 565nm (223Hz)
% modify frequency and amplitude individually using local variables amplitude, offset and frequency 
% negative peak of sine wave should not go below 0.15V as LED begins to turn off/labjack bottoms out
%
% original script from LackJack: example code on MATLAB for LJM (https://labjack.com/support/software/examples/ljm/matlab)
% 
% v1.0 Seul Ah Kim and Mike Wallace, Harvard Medical School 
%   - added functionality to save stream-in values as .mat files in the
%   local folder
%   - enabled stream-out of sinusoidal wave on DAC0 and DAC1

clc  % Clear the MATLAB command window
clear; close all
addpath(addpath(genpath('\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Shun\Behavior setup\Freq-modulated photometry')));
% Make the LJM .NET assembly visible in MATLAB
ljmAsm = NET.addAssembly('LabJack.LJM');

script_used = mfilename();
samplerate = 2052; % Hz
spectWindow = 200; % window size for frequency calculation in spectrogram (number of samples)
spectOverlap = 180; % overlap between windows in spectrogram; new calculation every 20 samples thus downsampled to 100Hz.
freqRange = [161:5:181;218:5:238]; % frequencies for channel 1 spectrogram (target 171, 228 Hz)

% Choose folder to save data
session = input('Session name: ','s'); session = strcat(session,'_g0');
% root_path = uigetdir('','Choose folder to save photometry data');
root_path = 'C:\Shun\Recordings';
session_path = strcat(root_path,'\',session,'\Photometry\');
mkdir(session_path);

% Calculate LED power scalar: outputPower = powerAtKnob * (0.5*LEDpowerblue)
outputPowerBlue = input('Blue ideal output power (uW): ');
outputPowerGreen = input('Green ideal output power (uW): ');
powerAtKnobBlue = 140; %input('Calibrated blue LED power (uW): ');
powerAtKnobGreen = 123; %input('Calibrated green LED power (uW): ');
LEDpowerblue = 2 * (outputPowerBlue/powerAtKnobBlue);
LEDpowergreen = 2 * (outputPowerGreen/powerAtKnobGreen);
colorarray = [0 1 0;1 0 0];
looplength = [12,9];
% p = gcp;
save(strcat(session_path,'info.mat'),'samplerate','looplength','LEDpowerblue','LEDpowergreen','script_used');

% Creating an object to nested class LabJack.LJM.CONSTANTS
t = ljmAsm.AssemblyHandle.GetType('LabJack.LJM+CONSTANTS');
LJM_CONSTANTS = System.Activator.CreateInstance(t);

handle = 0;
array = zeros(1,6);
array(2) = 8;
values = int32(6);

disp("----- Press ENTER to start recording -----"); pause;


% try
    % Open first found LabJack

    % Any device, Any connection, Any identifier
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
    % OUT0 = 470nm @ 167 Hz
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_LOOP_SIZE', 12);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.6191*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.5687*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.5*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.4313*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.3809*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.3625*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.3809*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.4312*LEDpowerblue); 
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.5000*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.5687*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.6191*LEDpowerblue);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_BUFFER_F32', 0.6375*LEDpowerblue); 
    % OUT1 = 565nm @ 223 Hz  20 uW power 
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_LOOP_SIZE', 9);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (1.2500-0.5)*LEDpowergreen);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (1.1911-0.5)*LEDpowergreen); 
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (1.0422-0.5)*LEDpowergreen); 
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (0.8734-0.5)*LEDpowergreen);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (0.7643-0.5)*LEDpowergreen);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (0.7661-0.5)*LEDpowergreen);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (0.8782-0.5)*LEDpowergreen); 
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (1.0476-0.5)*LEDpowergreen);  
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_BUFFER_F32', (1.1946-0.5)*LEDpowergreen);  
    
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT0_SET_LOOP', 1);
    LabJack.LJM.eWriteName(handle, 'STREAM_OUT1_SET_LOOP', 1);
    [~, value] = LabJack.LJM.eReadName(handle, 'STREAM_OUT2_BUFFER_STATUS', 0);
    disp(['STREAM_OUT2_BUFFER_STATUS = ' num2str(value)])
    
% end

    
    
    % Stream-in  configuration

    % Scan list names to stream-in
    numAddressesIn = 5;
    aScanListNames = NET.createArray('System.String', numAddressesIn);
    aScanListNames(1) = 'AIN0'; % photodiode green channel
    aScanListNames(2) = 'AIN1'; % photodiode red channel
    aScanListNames(3) = 'AIN2'; % copy of DAC0 out to 470 LED
    aScanListNames(4) = 'AIN3'; % copy of DAC1 out to 565 LED
    aScanListNames(5) = 'DIO0'; % Sync pulse

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
        numFrames = 7;
        aNames = NET.createArray('System.String', numFrames);
        aNames(1) = 'AIN_ALL_NEGATIVE_CH';
        aNames(2) = 'AIN0_RANGE';
        aNames(3) = 'AIN1_RANGE';
        aNames(4) = 'AIN2_RANGE';
        aNames(5) = 'AIN3_RANGE';
        aNames(6) = 'STREAM_SETTLING_US';
        aNames(7) = 'STREAM_RESOLUTION_INDEX';
        aValues = NET.createArray('System.Double', numFrames);
        aValues(1) = LJM_CONSTANTS.GND;
        aValues(2) = 10.0;
        aValues(3) = 10.0;
        aValues(4) = 10.0;
        aValues(5) = 10.0;
        aValues(6) = 0;
        aValues(7) = 0;
   %end

    % Write the analog inputs' negative channels (when applicable), ranges
    % stream settling time and stream resolution configuration.
    LabJack.LJM.eWriteNames(handle, numFrames, aNames, aValues, -1);

    % Configure and start stream
    numAddresses = aScanList.Length;
    [~, scanRate] = LabJack.LJM.eStreamStart(handle, scansPerRead, ...
        numAddresses, aScanList, scanRate);
    disp(['Stream started with a scan rate of ' num2str(scanRate) ' Hz.'])

    % The number of eStreamRead calls to perform in the stream read loop
    maxRequests = 7200; %1 per sec

    % Make a cell array out of scan list names. Helps performance when
    % converting and displaying in the loop.
    aScanListNamesML = cell(aScanListNames);

    tic
    disp(['Performing ' num2str(maxRequests) ' stream reads.'])

    % totalScans = 0;
    curSkippedSamples = 0;
    totalSkippedSamples = 0;

    % Initialize data matrix, numAddressesOut not recorded
    DATAM = zeros(maxRequests,samplerate*numAddressesIn); 

    % Added below to support live plotting feature - SAK 11.05.19
    figure('Position',[0 100 600 800]); 
    subplot(7,1,1); hold on; title('Demodulated')
    subplot(7,1,3); hold on; title('Demodulated mBeRFP blue excitation')
    subplot(7,1,4); hold on; title('Sync pulse')
      %subplot(7,1,8);hold on; title('excitation stability')
    drawnow
    speclength = floor((samplerate-spectOverlap)/(spectWindow-spectOverlap));
    rawSig = zeros(2,speclength);
    % ends here

    for i = 1:maxRequests
        [~,devScanBL,ljmScanBL] = LabJack.LJM.eStreamRead(handle,aData,0,0);

        % totalScans = totalScans + scansPerRead;
        % DATAM(i,:)= aData.double;
        temp = aData.double;
        % save each scan data as separate .mat files
        save(strcat(session_path,sprintf('Raw_%d.mat',i+1000)),'temp');

        % added live plotting feature - from Sarah's script. still need
        % to test if this leads to slowed streaming-- check the skipped
        % scan number - SAK 11.05.19
%             if i>2
%             for (j=1:2);
%                 rawSig(j,:)=mean(abs(spectrogram(temp(j:numAddressesIn:end),spectWindow,spectOverlap,freqRange(j,:),samplerate)),1);
%                 subplot(7,1,j);hold on;plot(i*speclength-speclength+1:i*speclength,rawSig(j,:),'Color',colorarray(j,:));axis tight
%             %                     subplot(8,1,j+2);hold on;plot(i*samplerate-samplerate+1:i*samplerate,temp(j:numAddressesIn:end),'Color',colorarray(j,:));axis tight
%                 subplot(7,1,j+5);hold on;plot(i*samplerate-samplerate+1:i*samplerate,temp(j+5:numAddressesIn:end),'Color','b');axis tight
%             %                     hold on;plot(i*samplerate-samplerate+1:i*samplerate,temp(5:numAddressesIn:end),'Color','r');axis tight
%             %                     subplot(8,1,j+4);hold on;plot(i*(samplerate/looplength(j))-(samplerate/looplength(j))+1:i*(samplerate/looplength(j)),temp(j:numAddressesIn*looplength(j):end),'Color',colorarray(j,:));axis tight
%             %                     subplot(7,1,j+7);hold on;plot(i*(samplerate/looplength(j))-(samplerate/looplength(j))+1:i*(samplerate/looplength(j)),temp(j+2:numAddressesIn*looplength(j):end),'Color',colorarray(j,:));axis tight
%                 drawnow  
%             end
%              rawSig(3,:)=mean(abs(spectrogram(temp(2:numAddressesIn:end),spectWindow,spectOverlap,freqRange(1,:),samplerate)),1);
%              subplot(7,1,3);hold on;plot(i*speclength-speclength+1:i*speclength,rawSig(3,:),'Color',colorarray(2,:));axis tight
%              subplot(7,1,4);hold on;plot(i*samplerate-samplerate+1:i*samplerate,temp(8:numAddressesIn:end),'Color','b');axis tight
%              subplot(7,1,5);hold on;plot(i*samplerate-samplerate+1:i*samplerate,temp(5:numAddressesIn:end),'Color','b');axis tight
%             %                            
% 
%             end
%             if mod(i,600)==0;
%             close all
%             figure('Position',[0 100 600 800]); 
%             subplot(7,1,1);hold on; title('demodulated')
%             subplot(7,1,3);hold on; title('demodulated mBeRFP blue excitation')
%             subplot(7,1,4);hold on; title('start')
%             subplot(7,1,5);hold on; title('shock')
%             subplot(7,1,6);hold on; title('sounds')
%             %                 subplot(8,1,8);hold on; title('excitation stability')
%             drawnow
%             end

%             i
        % ends here


        % Count the skipped samples which are indicated by -9999
        % values. Skipped samples occur after a device's stream buffer
        % overflows and are reported after auto-recover mode ends.
        % When streaming at faster scan rates in MATLAB, try counting
        % the skipped packets outside your eStreamRead loop if you are
        % getting skipped samples/scan.
        curSkippedSamples = sum(double(aData) == -9999.0);
        totalSkippedSamples = totalSkippedSamples + curSkippedSamples;

%             disp(['eStreamRead ' num2str(i)])
%             slIndex = 1;
%             for j = 1:scansPerRead*numAddressesIn
%                 fprintf('  %s = %.4f,', ...
%                         char(aScanListNamesML(slIndex)), aData(j))
%                 slIndex = slIndex + 1;
%                 if slIndex > numAddressesIn
%                     slIndex = 1;
%                     fprintf('\n')
%                 end
%             end
%             disp(['  Scans Skipped = ' ...
%                   num2str(curSkippedSamples/numAddressesIn) ...
%                   ', Scan Backlogs: Device = ' num2str(devScanBL) ...
%                   ', LJM = ' num2str(ljmScanBL)])
    end

    disp(['Performing ' num2str(maxRequests) ' stream reads.'])

%         totalScans = 0;
%         curSkippedSamples = 0;
%         totalSkippedSamples = 0;
        
        %timeElapsed = toc;

%         disp(['Total scans = ' num2str(totalScans)])
%         disp(['Skipped Scans = ' ...
%               num2str(totalSkippedSamples/numAddressesIn)])
%         disp(['Time Taken = ' num2str(timeElapsed) ' seconds'])
%         disp(['LJM Scan Rate = ' num2str(scanRate) ' scans/second'])
%         disp(['Timed Scan Rate = ' num2str(totalScans/timeElapsed) ...
%               ' scans/second'])
%         disp(['Sample Rate = ' ...
%               num2str(numAddressesIn*totalScans/timeElapsed) ...
%               ' samples/second'])
%     catch e
%         showErrorMessage(e)
%     end

    disp('Stop Stream')
    LabJack.LJM.eStreamStop(handle);

    % Close handle
    LabJack.LJM.Close(handle);
catch e
    showErrorMessage(e)
    LabJack.LJM.CloseAll();
    close all
end
