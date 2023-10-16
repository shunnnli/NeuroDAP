function processed = detrendAndDownsample(rawSignal,options)

arguments
    rawSignal double
    % options.prodecure string = 'detrend-ds-detrend'; % preprocess prodedure
    options.behaviorFs double = 10000
    options.targetFs double = 50
    options.rollingWindowTime double = 60 % in sec
    options.movingAvergeFilter logical = true
    options.movingAverageWindowSize double = 2
    options.dsMethod string = 'resample'
end


% % Rolling z score to detrend
% rollingmean = movmean(rawSignal,rollingSize*options.behaviorFs);
% rollingstd = movstd(rawSignal,rollingSize*options.behaviorFs);
% detrended = (rawSignal - rollingmean)./rollingstd;

% Downsample to 50Hz
ds_photometry = downsamplePhotometry(rawSignal,options.targetFs,options.behaviorFs,...
                                     movingAvergeFilter=options.movingAvergeFilter,...
                                     movingAverageWindowSize=options.movingAverageWindowSize,...
                                     dsMethod=options.dsMethod);

% Rolling z score
rollingSize = options.rollingWindowTime; % rolling window in sec
processed = rollingZ(ds_photometry,rollingSize);
% rollingmean = movmean(ds_photometry,rollingSize*options.targetFs);
% rollingstd = movstd(ds_photometry,rollingSize*options.targetFs);
% processed = (ds_photometry - rollingmean)./rollingstd;

end