function plotSpikes(spikeRate,clusterList,timeRange,color,options)

% Plot spike rate traces around certian event (with SEM)

% parse inputs
arguments
    spikeRate double
    clusterList double
    timeRange double
    color 

    options.average logical = false % Average over clusterList

    options.textOn logical = true
    options.text_source struct = struct([])

    options.unit string = 's' % unit of time (default as ms)

    options.smooth logical = false
    options.smoothWindow double = 10
    options.smoothMethod string = 'movmean'

    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 1
end

% Check color is valid
if ~ischar(color) && any(sum(color)>3) && size(color,1) == 1
    color = color ./ 255;
end

% Set colormap
cmap = colormap(flipud(jet(length(clusterList))));
ncolor = round(linspace(1,size(color,1),length(clusterList)));


% Set axis
if strcmp(options.unit,'ms')
    timeRange = timeRange .* 1000;
    t = linspace(timeRange(1),timeRange(2),size(spikeRate,3));
    xlabel('Time (ms)'); xlim([timeRange(1),timeRange(2)]); hold on
elseif strcmp(options.unit,'s')
    t = linspace(timeRange(1),timeRange(2),size(spikeRate,3));
    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]); hold on
else
    warning('Time unit is neither ms or s, set to ms.');
    timeRange = timeRange .* 1000;
    t = linspace(timeRange(1),timeRange(2),size(spikeRate,3));
    xlabel('Time (ms)'); xlim([timeRange(1),timeRange(2)]); hold on
end
% Set plot related params
set(gca,'ColorOrder',cmap);
ylabel('Spikes/s'); hold on

% Plot traces
if options.average
    % Plot averaged traces
    if options.smooth
        mean_trace = smoothdata(squeeze(mean(spikeRate(clusterList,:,:),2)),2,options.smoothMethod,options.smoothWindow);
    else
        mean_trace = squeeze(mean(spikeRate(clusterList,:,:),2));
    end
    plotSEM(t,mean_trace,color(1,:));
else
    % Plot individual traces
    for i = 1:length(clusterList)
        neuron = clusterList(i);
        if options.smooth
            smoothed_fr = smoothdata(squeeze(spikeRate(neuron,:,:)),2,options.smoothMethod,options.smoothWindow);
        else
            smoothed_fr = squeeze(spikeRate(neuron,:,:));
        end
        plotSEM(t,smoothed_fr,color(ncolor(i),:)); hold on
        % drawnow
        
        % Annotate cluster_id for each traces
        if options.textOn
            mean_trace = mean(smoothed_fr);
            [~,x_max] = max(abs(mean_trace));
            xtxtpos = t(x_max);
            ytxtpos = mean_trace(x_max);
            if isempty(options.text_source)
                text(xtxtpos,ytxtpos,num2str(neuron)); hold on
            else
                text(xtxtpos,ytxtpos,num2str(options.text_source.goodClusters(neuron))); hold on
            end
        end
    end
end

ylim([0 Inf]);
box off

end % plotSpikes