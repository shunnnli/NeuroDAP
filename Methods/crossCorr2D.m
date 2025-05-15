function C = crossCorr2D(dataX, dataY, options)
% Computes the trial‐factorized 2D cross‐covariance (or correlation) between
% two 3D arrays of shape (nTrials × T × nNeurons).
%
%   C = compute2dCrossCovFactorized(dataX, dataY, doTimebinZscore)
%
% Inputs:
%   dataX               nTrials×T×nX  array (e.g. DA from region X)
%   dataY               nTrials×T×nY  array (e.g. DA from region Y)
%
% Output:
%   C                   T×T cross‐covariance matrix:
%                       C(t1,t2) = mean_i[ sumX(i,t1)*sumY(i,t2) ]/(nTrials·nX·nY)

arguments
    dataX double
    dataY double

    options.average logical = false % did not care about trial, just use trial averages
    options.zscore logical = false
end

% 0) handle empty inputs
if isempty(dataX) || isempty(dataY)
    C = [];
    return
end

% 1) check dimensions
% Promote 2D → 3D as “1‐neuron” case
if ndims(dataX)==2
    [nTrials,T] = size(dataX);
    dataX = reshape(dataX, [nTrials, T, 1]);
end
if ndims(dataY)==2
    [nTrials2,T2] = size(dataY);
    dataY = reshape(dataY, [nTrials2, T2, 1]);
end

[nTrials, T, nX] = size(dataX);
[nTrials2, T2, nY] = size(dataY);
if ~options.average && (nTrials~=nTrials2 || T~=T2)
    warning('dataX and dataY must have the same (nTrials,T) dimensions. Subsamping to match');
    nX_trials = size(dataX,1);
    nY_trials = size(dataY,1);
    if nX_trials > nY_trials
        idx = randperm(nX_trials, nY_trials);
        dataX = dataX(idx,:,:);
        nTrials = nY_trials;
    elseif nY_trials > nX_trials
        idx = randperm(nY_trials, nX_trials);
        dataY = dataY(idx,:,:);
        nTrials = nX_trials;
    end
end

% trial‐averaged outer-product
if options.average
    % collapse neurons
    sumX = squeeze(sum(dataX,3));  % nTrials×T
    sumY = squeeze(sum(dataY,3));  % nTrials×T

    % mean across all trials
    meanX = mean(sumX,1);          % 1×T
    meanY = mean(sumY,1);          % 1×T

    % zero-mean each waveform
    meanX = meanX - mean(meanX);
    meanY = meanY - mean(meanY);

    % outer product, normalized by #neurons
    C = (meanX' * meanY) / (nX * nY);
    return
end

% 2) compute mean across trials for each (t,neuron)
%    meanX is T×nX, meanY is T×nY
meanX = squeeze(mean(dataX,1));   % collapse trials → [T × nX]
meanY = squeeze(mean(dataY,1));   %              [T × nY]

% 3) subtract those means to center each time‐neuron bin
%    Xp,Yp are still [nTrials × T × nNeurons]
Xp = bsxfun(@minus, dataX, reshape(meanX, [1 T nX]));
Yp = bsxfun(@minus, dataY, reshape(meanY, [1 T nY]));

% 4) optional: z-score each (t,neuron) across trials
if options.zscore
    % MATLAB’s default std uses ddof=0; use normalization by (n−1):
    stdX = squeeze(std(dataX, 0, 1));  % [T × nX]
    stdY = squeeze(std(dataY, 0, 1));  % [T × nY]
    stdX(stdX==0) = eps;
    stdY(stdY==0) = eps;
    Xp = bsxfun(@rdivide, Xp, reshape(stdX, [1 T nX]));
    Yp = bsxfun(@rdivide, Yp, reshape(stdY, [1 T nY]));
end

% 5) sum over neurons to get trial‐by‐time matrices [nTrials × T]
sumX = squeeze(sum(Xp, 3));  % sum over nX
sumY = squeeze(sum(Yp, 3));  % sum over nY

% 6) factorized cross‐covariance: (sumX' * sumY) gives [T × T]
%    then divide by the total number of “samples” (nTrials·nX·nY)
C = (sumX' * sumY) / (nTrials * nX * nY);

end
