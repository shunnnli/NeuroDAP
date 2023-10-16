function processed = demodulatePhotometry(signal,finalFs,params,options)

% Shun Li, 2023/10/15
% Adapted from code written by Bernardo Sabatini

arguments
    signal (1,:) {mustBeNumeric}
    finalFs double
    params struct
    
    options.modFreq = 167; % Hz (green labjack mod frequency)
    options.originalFs = 2000; % Hz
    options.pointsToEstimateCarrier = 1e6; % samples
    options.bandWidth = 3;            % number of frequency steps by Hz.  Eg 1 means analyze center frequency and +/- 1 Hz
    options.detrendWindow = 180;      % seconds
    options.spectralWindow = 966;     % samples (prev 2*9*12)
    
    options.resample = true;          % resample to finalFs when demodulation didn't produce finalFs
    options.plotFFT = false;
end

%% Setup

rawTimeStep = 1/originalFs;
finalTimeStep = 1/finalFs;
options.nSampPerDemodBin = finalTimeStep/rawTimeStep; % =40 if labjackFs=2000Hz; previously finalDownSample

if isnumeric(signal) && ~isfloat(signal); signal = double(signal);
else; error('Data must be numeric'); end

%% Find the carrier frequency with an FFT
endPoint = min(options.pointsToEstimateCarrier, length(signal));
sData_fft = fft(normalize(signal(1:endPoint)));
P2 = abs(sData_fft/endPoint);
P1 = P2(1:endPoint/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% make the frequency bins
fftFreq = originalFs * (0:(endPoint/2))/endPoint;
if options.plotFFT
    figure; plot(fftFreq,P1);
    title('FFT'); 
    set(gca, 'YScale', 'log', 'XScale', 'log');
end

% find the frequency of the peak power
[~, maxFindex] = max(P1);
calculatedModFrequency = fftFreq(maxFindex);
disp(['Discovered modulation frequency is ' num2str(calculatedModFrequency),' Hz']);

if ~isnan(options.modFreq) && (options.modFreq~=calculatedModFrequency)
    warning('Discovered modulation frequency does not equal the user provided frequency!');
    disp(['Using modulation frequency ' num2str(options.modFreq)]);
end

%% do the demodulation for the non-detrended data

% Calculate demodulation params
options.spectralFrequencies = (-options.bandWidth:options.bandWidth)+options.modFreq;
options.spectralWindowOverlap = 2*options.nSampPerBin; % previously: options.spectralWindow-options.nSampPerDemodBin;

% Calculte demodulation window
% ensures that it is an integer multiple of the sampling window
detrendWindowSamples_rawFs = 2*floor(options.detrendWindow*originalFs/2); 
detrendWindowSamples_finalFs = 2*floor(options.detrendWindow*finalFs/2);

% Demod without detrend
[spectVals, ~, dmTimes] = spectrogram(signal, options.spectralWindow, options.spectralWindowOverlap, options.spectralFrequencies, originalFs);
processed.demodData_nodetrend = mean(abs(spectVals),1);    % save raw demodulated (not detrended)
processed.demodTimes = dmTimes;                              % save time points

% Demod with detrend
signal = rollingZ(signal, detrendWindowSamples_rawFs);

[spectVals, ~, ~] = spectrogram(signal, options.spectralWindow, options.spectralWindowOverlap, options.spectralFrequencies, originalFs);
dmData = mean(abs(spectVals),1); % convert spectrogram to power    

dmData = rollingZ(dmData, detrendWindowSamples_finalFs);
dmData(1:(detrendWindowSamples_finalFs/2)) = nan;
dmData((end-detrendWindowSamples_finalFs/2):end) = nan;
processed.demodData = dmData;

% Resample to targetFs if neccessary
if options.resample && (finalFs ~= length(dmData)/params.photometry.totalDuration)
    [p,q] = rat(fianlFs/originalFs);
    % n = 10; beta = 5; % n: length of filter window (default 10); beta: smoothing (default 5)
    processed.demodData = resample(dmData,p,q);
    disp('Finished: resampled demod data to downsampleFs');
end

processed.options = options;

end
