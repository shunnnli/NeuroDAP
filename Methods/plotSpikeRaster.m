function plotSpikeRaster(spikeRate,clusterList,timeRange,options)

% Draw spikes as rasters around an event

% parse inputs
arguments
    spikeRate double
    clusterList double
    timeRange double
    options.color = [0 0 0]
    options.size = 30
    options.unit string = 's'
    options.originalTimeRange double = timeRange % required if plotting timeRange and getSpikes timeRange is different
end

% Check color is valid
color = options.color;
if ~ischar(color) && sum(color) > 3 && size(color,1) == 1
    color = color ./ 255;
end

% Calculate binSize
binSize = (options.originalTimeRange(2) - options.originalTimeRange(1)) / size(spikeRate,3);
eventIdxInBin = abs(options.originalTimeRange(1) / binSize) + 1;


for i = 1:length(clusterList)
    neuron = clusterList(i);
    psthMatrix = squeeze(spikeRate(neuron,:,:));

    % Initalize progress bar
    textprogressbar(['Plotting spike rasters of neuron ',num2str(neuron),': ']);

    for j = 1:size(psthMatrix,1)
        textprogressbar(j/size(spikeRate,2) * 100);

        spikeTimes = (find(psthMatrix(j,:)) - eventIdxInBin) .* binSize;
        spikeTimes = spikeTimes(spikeTimes>=timeRange(1) & spikeTimes <= timeRange(2));
        if all(size(spikeTimes) == [0 0]); spikeTimes = double.empty(1,0); end

        if strcmp(options.unit,'ms')
            spikeTimes = spikeTimes .* 1000;
            xlabel('Time (ms)'); xlim([timeRange(1)*1000,timeRange(2)*1000]); hold on
        elseif strcmp(options.unit,'s')
            xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]); hold on
        else
            warning('Time unit is neither ms or s, set to ms.');
            spikeTimes = spikeTimes .* 1000;
            xlabel('Time (ms)'); xlim([timeRange(1),timeRange(2)]); hold on
        end
        
        scatter(spikeTimes,j,options.size,'filled','MarkerFaceColor',color); hold on
        % disp(['Finished: plotting raster for trial ', num2str(j), ' of neuron ', num2str(neuron)]);
    end
end

% Set plot related params
ylabel('Trials'); % ylim([0 size(spikeRate,2)]); hold on
box off

textprogressbar('done');
end