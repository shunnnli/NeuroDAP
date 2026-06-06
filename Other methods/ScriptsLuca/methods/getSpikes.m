function [spike,spikeRate,spike_params] = getSpikes(timeRange,binSize,eventIdx,params,options)

% The number of spikes for all good units at each time bin (binSize) within
% time before/after the event (timeRange)

% OUTPUT: 3D matrix spikes and spikeRate (x,y,z)
% x: number of neurons
% y: number of events
% z: timepoints

arguments
    timeRange double
    binSize double
    eventIdx double
    params struct
    options.clusters double = nan % find spikes for selective neurons
    options.maxLatency double = nan % maximum latency of the first spike detected in sec
end

if mod(abs(timeRange(1))/binSize,1) ~= 0
    warning('timeRange not divisible by binSize, consider resetting binSize.');
end

% Initialize params
ap = params.ephys.ap;
nidqFs = params.sync.behaviorFs;
apFs = params.sync.apFs;
timeImec = params.sync.timeImec;
timeNI = params.sync.timeNI;
spike_params.binSize = binSize;
spike_params.timeRange = timeRange;
spike_params.nTrials = length(eventIdx);
spike_params.eventBin = abs(timeRange(1)/binSize) + 1;
spike_params.triggeredSpikeIdx = [];

timesteps = (timeRange(2)-timeRange(1))/binSize;

% Whether to get all spikes or only spikes from selected neurons
if isnan(options.clusters)
    spike_params.clusters = 1:ap.nGoodClusters;
    spike = zeros(ap.nGoodClusters,length(eventIdx),timesteps);
else
    spike_params.clusters = options.clusters;
    spike = zeros(length(options.clusters),length(eventIdx),timesteps);
end
    
% getSpikes
try
    % Get all neuron spikes
    if isnan(options.clusters)
        bar = waitbar(0,'getSpikes...','Name','getSpikes',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(bar,'canceling',0);
        
        for i = 1:length(eventIdx)
        
            % Check for clicked Cancel button
            if getappdata(bar,'canceling')
                break
            end
        
%             % Find first & last NI index
%             niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidqFs);
%             if niFirstIdx < 0; niFirstIdx = 1; end
%             % Find corresponding imec index
%             % imecFirstIdx = round(timeNI(niFirstIdx)*apFs); % make sure timeImec and timeNI is aligned to 0
%             [~, imecFirstIdx] = min(abs(timeImec-timeNI(niFirstIdx)));
%             imecLastIdx = imecFirstIdx + floor(apFs*(timeRange(2)-timeRange(1)));

            % Find imec event sample
%             nSamplesToInclude = apFs*(timeRange(2)-timeRange(1));
            [~,imecEventIdx] = min(abs(timeImec-timeNI(eventIdx(i))));
            imecFirstIdx = imecEventIdx + apFs*timeRange(1);
            imecLastIdx = imecEventIdx + apFs*timeRange(2);
        
            % Find spike within timeRange
            spikeTimesinRangeIdx = find...
                (ap.goodSpikeTimes>imecFirstIdx & ap.goodSpikeTimes<imecLastIdx);
            spikeTimesinRange = ap.goodSpikeTimes(spikeTimesinRangeIdx);
            spikeClustersinRange = ap.goodSpikeClusters(spikeTimesinRangeIdx);
        
            % Generate spike raster array
            for s = 1:length(spikeTimesinRange)
                % Find which neuron fires the spike
                cluster_id = spikeClustersinRange(s);
                goodCluster_id = ap.clusterToGoodClusterIndex(cluster_id);
                
                % Find the timebin of the spike
                relativeSpikeSamp = spikeTimesinRange(s)-imecFirstIdx;
                relativeSpikeTime = double(relativeSpikeSamp)/apFs;
                relativeSpikeBin = floor(relativeSpikeTime/binSize)+1;
                
                % Add spike to array
                spike(goodCluster_id,i,relativeSpikeBin) = 1 + spike(goodCluster_id,i,relativeSpikeBin);

                if ~isnan(options.maxLatency)
                    % Add to event-triggered spike list
                    if imecEventIdx + options.maxLatency*apFs >= spikeTimesinRange(s)
                        spike_params.triggeredSpikeIdx = [spike_params.triggeredSpikeIdx;...
                            goodCluster_id i spikeTimesinRange(s)];
                    end
                end
            end
        
            % Update waitbar and message
            progress = i/length(eventIdx);
            waitbar(progress,bar,['Assigned ',num2str(length(spikeTimesinRange)),' spikes for event ', num2str(i),...
                ' within ', num2str(timeRange(1)),' to ',num2str(timeRange(2)),' sec']);
        end
        
        spikeRate = spike / binSize;
        delete(bar);

    % Get selected spikes
    else
        % XXXX
    end

catch ME
    delete(bar);
    rethrow(ME) % re-issue the error
end

end % getSpikes