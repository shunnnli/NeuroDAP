function [event_wvf,control_wvf] = getSpikeWaveforms(mmfDataX,clusterList,timeRange,params,options)

arguments
    mmfDataX
    clusterList double
    timeRange double
    params struct
    options.event double = zeros(0,0)
    options.maxLatency double = 100000 % maximum latency of the first spike detected in sec
    options.triggeredSpikeIdx double = nan
end


ap = params.ephys.ap;
timeImec = params.sync.timeImec;
timeNI = params.sync.timeNI;
timesteps = floor(ap.Fs*(timeRange(2)-timeRange(1)));

% wvf = number of candidates x number of optostim x timesteps
if isempty(options.event)
    nTrials = 500; nTrialsToAnalyze = nTrials;
    event_wvf = nan(length(clusterList),nTrials,timesteps);
    control_wvf = nan(length(clusterList),nTrials,timesteps);
else
    eventIdx = options.event;
    nTrials = length(eventIdx);
    event_wvf = nan(length(clusterList),nTrials,timesteps);
    control_wvf = nan(length(clusterList),nTrials,timesteps);
end

bar = waitbar(0,'getSpikeWaveforms...','Name','Extracting spike waveforms',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(bar,'canceling',0);

try
    
    for i = 1:length(clusterList)
        % Get best channel of the candidate unit
        cluster_id = ap.goodClusters(clusterList(i));
        best_channel = ap.cluster_info(find(ap.cluster_info(:,1) == cluster_id),6) + 1;
        neuron_spikes_idx = find(ap.goodSpikeClusters == cluster_id);
        neuron_spikes_time = ap.goodSpikeTimes(neuron_spikes_idx);
    
        if ~isempty(options.event)
            if ~isnan(options.triggeredSpikeIdx)
                triggeredSpikeIdx = options.triggeredSpikeIdx(options.triggeredSpikeIdx(:,1) == clusterList(i),3);
                nEventSpikesToPlot = length(triggeredSpikeIdx);
                eventSpikesSamp = triggeredSpikeIdx;
                nTrialsToAnalyze = nTrials + nEventSpikesToPlot;
                for j = 1:length(triggeredSpikeIdx)
                    % Check for clicked Cancel button
                    if getappdata(bar,'canceling'); break; end
                    % Update progress
                    progress = j/nTrialsToAnalyze;
                    waitbar(progress,bar,['Triggered spike ',num2str(j),'/',num2str(length(triggeredSpikeIdx)),...
                        ' for neuron ',num2str(clusterList(i))]);

                    % Find first and last imec index
                    imecFirstIdx = triggeredSpikeIdx(j) + floor(ap.Fs*timeRange(1));
                    imecLastIdx = imecFirstIdx + timesteps - 1;

                    % Read waveform within range
                    wf = mmfDataX(best_channel, imecFirstIdx:imecLastIdx);

                    % Spike waveform for each stim
                    event_wvf(i,j,:) = wf;
                end
            else
                nTrialsToAnalyze = nTrials * 2; nEventSpikesToPlot = nTrials;
                eventSpikesSamp = [];
                % Get opto-triggered waveform
                for j = 1:nTrials
                    % Check for clicked Cancel button
                    if getappdata(bar,'canceling'); break; end
                    % Update waitbar and message
                    progress = j/nTrialsToAnalyze;
                    waitbar(progress,bar,['Event ',num2str(j),'/',num2str(nTrials),...
                        ' for neuron ',num2str(clusterList(i))]);
            
                    % Find Imec stim index
                    niStimIdx = eventIdx(j) + floor(timeRange(1)*params.sync.behaviorFs);
                    [~, imecStimIdx] = min(abs(timeImec-timeNI(niStimIdx)));
            
                    % Find first spike
                    first_spikes_time = neuron_spikes_time(find(neuron_spikes_time > imecStimIdx,1));
                    % disp(1000*(first_spikes_time - imecStimIdx)/ap.Fs);
        
                    % Check spike latency
                    latency = first_spikes_time - imecStimIdx; % in imec samples
                    if latency <= options.maxLatency * ap.Fs
                        % Store plotted spike times
                        eventSpikesSamp = [eventSpikesSamp;first_spikes_time];

                        % Find first and last imec index
                        imecFirstIdx = first_spikes_time + floor(ap.Fs*timeRange(1));
                        imecLastIdx = imecFirstIdx + timesteps - 1;
                
                        % Read waveform within range
                        wf = mmfDataX(best_channel, imecFirstIdx:imecLastIdx);
                    else
                        wf = nan(1,timesteps);
                    end
            
                    % Spike waveform for each stim
                    event_wvf(i,j,:) = wf;
                end
            end
        end
    
        % Randomly select non optotim spikes from this unit
        baseline_spike_time = setdiff(neuron_spikes_time,eventSpikesSamp);
        if length(baseline_spike_time) <= nTrials
            nTrials = length(baseline_spike_time);
            random_spikes_time = baseline_spike_time(randi(nTrials,[nTrials,1]));
        else
            random_spikes_time = baseline_spike_time(randi(nTrials,[nTrials,1]));
        end
        for s = 1:nTrials
            % Check for clicked Cancel button
            if getappdata(bar,'canceling')
                break
            end
            if isempty(options.event); progress = s/nTrials;
            else; progress = (nEventSpikesToPlot+s)/nTrialsToAnalyze; end
            waitbar(progress,bar,['Baseline ',num2str(s),'/',num2str(nTrials),...
                ' for neuron ',num2str(clusterList(i))]);
    
            % Find first and last imec index
            imecFirstIdx = random_spikes_time(s) + floor(ap.Fs*timeRange(1));
            imecLastIdx = imecFirstIdx + timesteps - 1;
            
            % Read waveform within range
            wf = mmfDataX(best_channel, imecFirstIdx:imecLastIdx);
    
            % Spike waveform for each stim
            control_wvf(i,s,:) = wf; 
        end
        delete(bar);
    end

catch ME
    delete(bar);
    rethrow(ME) %re-issue the error
end

% Convert to mV using conversion factor
event_wvf = event_wvf .* 2.3438;
control_wvf = control_wvf .* 2.3438;

end