function photometry_downsample = downsamplePhotometry(rawTraces,targetFs,originalFs,options)

arguments
    rawTraces double % raw traces
    targetFs double % target Fs to downsampled to
    originalFs double % original sampling freq of the rawTraces
    options.movingAvergeFilter logical = false
    options.movingAverageWindowSize double = 2

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
end

end