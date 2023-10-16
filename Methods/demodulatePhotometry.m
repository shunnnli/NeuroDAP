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
    
    options.removeTwoEnds = false;    % nan values for the first half of the first spectral window and last half of the last spectural window
    options.resample = true;          % resample to finalFs when demodulation didn't produce finalFs
    options.plotFFT = false;
end

%% Setup
if isnumeric(signal) && ~isfloat(signal); signal = double(signal);
elseif ~isnumeric(signal); error('Demodulation: Data must be numeric'); end

%% Find the carrier frequency with an FFT
endPoint = min(options.pointsToEstimateCarrier, length(signal));
sData_fft = fft(normalize(signal(1:endPoint)));
P2 = abs(sData_fft/endPoint);
P1 = P2(1:endPoint/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% make the frequency bins
fftFreq = options.originalFs * (0:(endPoint/2))/endPoint;
if options.plotFFT
    figure; plot(fftFreq,P1);
    title('FFT'); 
    set(gca, 'YScale', 'log', 'XScale', 'log');
end

% find the frequency of the peak power
[~, maxFindex] = max(P1);
calculatedModFrequency = fftFreq(maxFindex);
disp(['     Demodulation: Discovered modulation frequency is ' num2str(calculatedModFrequency),' Hz']);

if ~isnan(options.modFreq) && (options.modFreq~=calculatedModFrequency)
    disp('     Demodulation: Discovered modulation frequency does not equal the user provided frequency!');
    disp(['     Demodulation: Using modulation frequency ', num2str(options.modFreq), ' Hz']);
end

%% do the demodulation for the non-detrended data

% Calculate demodulation params
options.nSampPerDemodBin = (1/finalFs)*options.originalFs; % =40 if labjackFs=2000Hz; previously finalDownSample
if mod(options.nSampPerDemodBin,1) ~= 0
    warning('Demodulation: nSampPerDemodBin is not integer, use the nearest floor instead!');
    options.nSampPerDemodBin = floor(options.nSampPerDemodBin);
end
options.spectralFrequencies = (-options.bandWidth:options.bandWidth)+options.modFreq;
options.spectralWindow = 2*options.nSampPerDemodBin;
options.spectralWindowOverlap = options.nSampPerDemodBin; % previously: options.spectralWindow-options.nSampPerDemodBin;

% Calculte demodulation window
% ensures that it is an integer multiple of the sampling window
detrendWindowSamples_rawFs = 2*floor(options.detrendWindow*options.originalFs/2); 
detrendWindowSamples_finalFs = 2*floor(options.detrendWindow*finalFs/2);

% Demod without detrend
[spectVals, ~, dmTimes] = spectrogram(signal, options.spectralWindow, options.spectralWindowOverlap, options.spectralFrequencies, options.originalFs);
processed.demodData_nodetrend = mean(abs(spectVals),1);    % save raw demodulated (not detrended)
processed.demodTimes = dmTimes;                              % save time points

% Demod with detrend
signal_detrended = rollingZ(signal, detrendWindowSamples_rawFs);

[spectVals, ~, ~] = spectrogram(signal_detrended, options.spectralWindow, options.spectralWindowOverlap, options.spectralFrequencies, options.originalFs);
dmData = mean(abs(spectVals),1); % convert spectrogram to power    

dmData = rollingZ(dmData, detrendWindowSamples_finalFs);
if options.removeTwoEnds
    dmData(1:(detrendWindowSamples_finalFs/2)) = nan;
    dmData((end-detrendWindowSamples_finalFs/2):end) = nan;
end
processed.demodData = dmData;

% Resample to targetFs if neccessary
demodFs = length(dmData)/params.photometry.totalDuration;
if options.resample && (finalFs ~= demodFs)
    [p,q] = rat(finalFs/demodFs);
    % n = 10; beta = 5; % n: length of filter window (default 10); beta: smoothing (default 5)
    processed.demodData = resample(dmData,p,q);
    disp('     Demodulation: resampled demod data to downsampleFs');
else; options.resample = false; 
end

processed.options = options;

end
