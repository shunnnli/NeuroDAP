function downsampled = downsamplePhotometry(rawTraces,options)

arguments
    rawTraces double % raw traces
    
    options.targetFs double = 50 % target Fs to downsampled to
    options.originalFs double = 2000 % original sampling freq of the rawTraces
    
    options.movingAvergeFilter logical = false
    options.movingAverageWindowSize double = 2

    options.rollingZ logical = true
    options.rollingWindowTime double = 180

    % Ways to downsample if targetFs and originalFs aren't divided by an
    % integer (nSampPerBin is not an integer)
    % 'closestInteger': select the closestInteger to set nSampPerBin
    % 'resample': use matlab's resample function
    options.dsMethod string = 'resample'
end

% Calculate nSampPerBin
nSampPerBin = (1/targetFs)*originalFs;

% Determine downsample approach based on whether targetFs and
% originalFs can be divided by integer
if mod(nSampPerBin,1) == 0
    % Downsample
    photometry_downsample = zeros(1,floor(length(rawTraces)/nSampPerBin));
    for k = 1:length(photometry_downsample)
        firstBinTime = floor(nSampPerBin*(k-1)+1);
        lastBinTime = floor(nSampPerBin*k);
        photometry_downsample(k) = sum(rawTraces(floor(firstBinTime:lastBinTime)));
    end
    clear nSampPerBin firstBinTime lastBinTime

    if options.movingAvergeFilter
        % Moving average filter
        b = (1/options.movingAverageWindowSize)*ones(1,options.movingAverageWindowSize);
        a = 1;
        photometry_downsample = filter(b,a,photometry_downsample);
    end
    downsampled.dsData = photometry_downsample;

% If sampleFs can't be divided by integer
else
    if strcmp(options.dsMethod,'cloestInteger')
        nSampPerBin = round(nSampPerBin);
        % Downsample
        photometry_downsample = zeros(1,floor(length(rawTraces)/nSampPerBin));
        for k = 1:length(photometry_downsample)
            firstBinTime = floor(nSampPerBin*(k-1)+1);
            lastBinTime = floor(nSampPerBin*k);
            photometry_downsample(k) = sum(rawTraces(floor(firstBinTime:lastBinTime)));
        end
        clear nSampPerBin firstBinTime lastBinTime

        if options.movingAvergeFilter
            % Moving average filter
            b = (1/options.movingAverageWindowSize)*ones(1,options.movingAverageWindowSize);
            a = 1;
            photometry_downsample = filter(b,a,photometry_downsample);
        end
    elseif strcmp(options.dsMethod,'resample')
        % Unfinished, see reason in loadSessions.m
        [p,q] = rat(targetFs/originalFs);
        % n = 10; beta = 5; % n: length of filter window (default 10); beta: smoothing (default 5)
        photometry_downsample = resample(rawTraces,p,q);
    end
    downsampled.dsData = photometry_downsample;
end

% Rolling z score
if options.rollingZ
    % Rolling zscore for demodGreen and demodRed
    options.signalDetrendWindow = floor(options.rollingWindowTime*options.originalFs);
    downsampled.finalData = rollingZ(downsampled.dsData,options.signalDetrendWindow);
end

downsampled.options = options;

end