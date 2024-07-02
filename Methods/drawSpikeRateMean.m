function drawSpikeRate(timeRange,spikeRate,clusterList,textOn,ap,colors)

% Set colormap
cmap = colormap(flipud(jet(length(clusterList))));

% Plot individual traces
for i = 1:length(clusterList)
    neuron = clusterList(i);
    if sum(spikeRate(neuron,:)) > 10
        % smoothed_fr = smoothdata(spikeRate(neuron,:),'movmean',10);
        smoothed_fr = spikeRate(neuron,:);
        [~,x_max] = max(abs(smoothed_fr));
        t = linspace(timeRange(1),timeRange(2),size(spikeRate,2));
        if length(clusterList) > length(colors)
            plot(t,smoothed_fr,'-','LineWidth',1.5); hold on
        else
            plot(t,smoothed_fr,'-','LineWidth',1.5,'Color',colors(i)); hold on
        end
        drawnow
        
        if textOn
            xtxtpos = t(x_max);
            ytxtpos = smoothed_fr(x_max);
            text(xtxtpos,ytxtpos,num2str(ap.goodClusters(neuron))); hold on
        end
    end
end

% Set plot related params
set(gca,'ColorOrder',cmap);
xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Spikes/s'); % ylim([0,100]);
box off

end % drawPSTH