function varargout = plotFFT(signal,options)

arguments
    signal double
    options.Fs double = 2000 % sampling frequency of the signal
    options.timeToEstimateCarrier double = 100 % in seconds
    options.xlim double = [0,inf]
    options.fillmissing logical = true
    options.color double = [0.5 0.5 0.5]
end

options.pointsToEstimateCarrier = options.timeToEstimateCarrier * options.Fs;
endPoint = min(options.pointsToEstimateCarrier, length(signal));

if options.fillmissing
    signal = fillmissing(signal,'next');
end

sData_fft = fft(normalize(signal(1:endPoint)));
P2 = abs(sData_fft/endPoint);
P1 = P2(1:endPoint/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% make the frequency bins
fftFreq = options.Fs * (0:(endPoint/2))/endPoint;

plot(fftFreq,P1,Color=options.color); box off; hold on;
title('FFT'); xlabel('Frequency (Hz)'); ylabel('Power');
xlim(options.xlim);
set(gca, 'YScale', 'log', 'XScale', 'log');

[~, maxFindex] = max(P1);
disp(['Discovered modulation frequency is ' num2str(fftFreq(maxFindex)),' Hz']);


% Optional output
varargout{1} = fftFreq;
varargout{2} = P1;

end