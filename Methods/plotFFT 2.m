function varargout = plotFFT(signal,options)

arguments
    signal double                   % can be 1×T or n×T

    options.plot            logical = true
    options.print           logical = true
    
    options.Fs               double = 2000  % sampling frequency
    options.timeToEstimateCarrier double = inf
    options.fillMaxFreq     logical = true % make the max frequency have the power as the bin before
    
    options.xlim             double = [0,inf]
    options.logScale        logical = true
    options.fillmissing     logical = true
    options.color            double = [0.5 0.5 0.5]
    options.plotIndividual  logical = false
end

% compute how many points we'll use for the FFT
options.pointsToEstimateCarrier = options.timeToEstimateCarrier * options.Fs;
L = min(options.pointsToEstimateCarrier, size(signal,2));

if options.fillmissing
    signal = fillmissing(signal,'next');
end

% single‐trial FFT
data_fft = fft(normalize(signal),L,2);
P2 = abs(data_fft/L);
P1 = P2(:,1:L/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
P1(:,end) = P1(:,end-1);

% frequency axis
fftFreq = options.Fs * (0:(L/2))/L;

% plot if asked
if options.plot
    plotSEM(fftFreq,P1,options.color,plotIndividual=options.plotIndividual);
    title('FFT'); xlabel('Frequency (Hz)'); ylabel('Power');
    xlim(options.xlim);
    if options.logScale
        set(gca, 'YScale', 'log');
    end
    box off;
end

% print peak
[~, maxFindex] = max(P1);
if options.print
    disp(['Discovered modulation frequency is ' num2str(fftFreq(maxFindex)),' Hz']);
end

% return varargout
varargout{1} = fftFreq;
varargout{2} = P1;

end
