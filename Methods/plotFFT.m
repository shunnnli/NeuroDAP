function varargout = plotFFT(signal,options)

arguments
    signal double
    options.Fs double = 2000 % sampling frequency of the signal
    options.pointsToEstimateCarrier double = 1e6
    options.xlim double = [0,inf]
end

endPoint = min(options.pointsToEstimateCarrier, length(signal));
sData_fft = fft(normalize(signal(1:endPoint)));
P2 = abs(sData_fft/endPoint);
P1 = P2(1:endPoint/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% make the frequency bins
fftFreq = options.Fs * (0:(endPoint/2))/endPoint;

plot(fftFreq,P1); box off
title('FFT'); xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim(options.xlim);
set(gca, 'YScale', 'log', 'XScale', 'log');

[~, maxFindex] = max(P1);
disp(['Discovered modulation frequency is ' num2str(fftFreq(maxFindex)),' Hz']);


% Optional output
varargout{1} = fftFreq;
varargout{2} = P1;

end