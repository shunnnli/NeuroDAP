% Caluclate firing rate of given cluster at given time window
% event time: in samples
% Output: spikeRate -> MxN matrix, average firing rate of M neuron across
% N time points across all events; spikeRate_for_each_cluster -> PxNxM
% matrix, spike rate for M neurons across N time points and P events


% ======================================= %
function [spike,spikeRate,spikeRate_for_each_cluster] = getSpikesByCluster(timeRange,binSize,eventIdx,...
                ap,nidq,timeImec,timeNI)

% The number of spikes for all good units at each time bin (binSize) within
% time before/after the event (timeRange)

spike = zeros(ap.nGoodClusters, round((timeRange(2)-timeRange(1))/binSize));
spike_for_each_event = [];
% spike_for_each_event = zeros(ap.nGoodClusters, round((timeRange(2)-timeRange(1))/binSize),length(eventIdx));

for i = 1:length(eventIdx)
    % Find first & last NI index
    niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidq.Fs);
    if niFirstIdx < 0; niFirstIdx = 1; end
    % Find corresponding imec index
%     imecFirstIdx = round(timeNI(niFirstIdx)*ap.Fs); % make sure timeImec and timeNI is aligned to 0
    [~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
    imecLastIdx = imecFirstIdx + floor(ap.Fs*(timeRange(2)-timeRange(1)));

    % Find spike within timeRange
    spikeTimesinRangeIdx = find...
        (ap.goodSpikeTimes>imecFirstIdx & ap.goodSpikeTimes<imecLastIdx);
    spikeTimesinRange = ap.goodSpikeTimes(spikeTimesinRangeIdx);
    spikeClustersinRange = ap.goodSpikeClusters(spikeTimesinRangeIdx);
    disp(['Found ',num2str(length(spikeTimesinRange)),' spikes for event ', num2str(i)]);

    % Generate spike raster array
    tmp_spike = zeros(ap.nGoodClusters, round((timeRange(2)-timeRange(1))/binSize));
    for s = 1:length(spikeTimesinRange)
        % Find which neuron fires the spike
        cluster_id = spikeClustersinRange(s);
        goodCluster_id = ap.clusterToGoodClusterIndex(cluster_id);
        
        % Find the timebin of the spike
        relativeSpikeSamp = spikeTimesinRange(s)-imecFirstIdx;
        relativeSpikeTime = double(relativeSpikeSamp)/ap.Fs;
        relativeSpikeBin = floor(relativeSpikeTime/binSize)+1;
        
        % Add spike to array
        spike(goodCluster_id,relativeSpikeBin) = 1 + spike(goodCluster_id,relativeSpikeBin);
        tmp_spike(goodCluster_id,relativeSpikeBin) = 1 + tmp_spike(goodCluster_id,relativeSpikeBin);
    end
    spike_for_each_event = cat(3, spike_for_each_event, tmp_spike);
    spikeRate_for_each_event = spike_for_each_event / binSize;
end
spikeRate = spike / (binSize*length(eventIdx));
spikeRate_for_each_cluster = permute(spikeRate_for_each_event,[3 2 1]);

end % getPSTH