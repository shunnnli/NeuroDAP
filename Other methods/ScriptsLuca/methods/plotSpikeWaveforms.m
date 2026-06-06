function [wvf_mean] = plotSpikeWaveforms(wvf,clusterList,timeRange,color,options)

% parse inputs
arguments
    wvf
    clusterList double
    timeRange double
    color
    options.LineWidth double = 1.5
end

t = linspace(timeRange(1),timeRange(2),size(wvf,3)) * 1000;


for j = 1:size(wvf,2)
    plot(t,squeeze(wvf(neuron,j,:)),"Color",color); hold on
end
wvf_mean = plot(t,mean(squeeze(wvf(neuron,:,:))),"Color",color,'LineWidth',options.LineWidth); hold on

xlabel('Time (ms)'); ylabel('Signal'); hold on
box off

end