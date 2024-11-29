% Caluclate firing rate of given cluster at given time window
% event time: in samples

% ======================================= %
function [spike,spikeRate] = getSpikesSparse(eventIdx,timeRange,binSize,...
                params,timeImec,timeNI)

% The number of spikes for all good units at each time bin (binSize) within
% time before/after the event (timeRange)

% OUTPUT: 3D matrix spikes and spikeRate (x,y,z)
% x: number of neurons
% y: number of events
% z: timepoints

if (1/binSize) > params.ephys.finalFs
    warning('binSize smaller than finalFs binSize, reset to original binSize.');
end

if mod(abs(timeRange(1))/binSize,1) ~= 0
    warning('timeRange not divisible by binSize, consider resetting binSize.');
end

ntimesteps = round((timeRange(2)-timeRange(1))/binSize);
spike = zeros(params.ephys.ap.nGoodClusters,length(eventIdx),ntimesteps);
eventInImec = findCorrespondingTime(eventIdx,timeNI,timeImec);

% Generate spike sparse matrix (downsample to 1/binSize)
params.ephys.finalFs = 500;
spikeIdx = params.ephys.spikeIdx;
spikeDownSample = params.ephys.apFs/(1/binSize); % recording is downsampled by 600x to get a final sample rate of 50 Hz
nDownSamples = floor(params.ephys.ap.totalSampIncluded/spikeDownSample);
spikeTime_downsampled = floor(params.ephys.ap.goodSpikeTimes/spikeDownSample)+1;
% Remove spikes outside of totalSampIncluded
outsideSpikes = find(spikeTime_downsampled > nDownSamples);
if ~isempty(outsideSpikes); spikeIdx(outsideSpikes) = []; spikeTime_downsampled(outsideSpikes) = []; end
nSpikes = length(spikeIdx); val = ones(nSpikes, 1);

% Sparse matrix
spikes = sparse(spikeIdx, spikeTime_downsampled, val, params.ephys.ap.nGoodClusters, nDownSamples, nSpikes);

% Extract event to output spike matrix
for i = 1:length(eventInImec)
    % Find first and last index
    imecFirstIdx = eventInImec(i) - abs(timeRange(1))/binSize;
    imecLastIdx = eventInImec(i) + abs(timeRange(2))/binSize -1;
    disp(imecFirstIdx);
    spike(:,i,:) = spikes(:,imecFirstIdx:imecLastIdx);
end

spikeRate = spike / binSize;

end % getSpikes