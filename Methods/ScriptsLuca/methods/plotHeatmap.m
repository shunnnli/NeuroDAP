function plotHeatmap(traces,timestamp,options)

arguments
    traces double
    timestamp double
    options.flipYAxis logical = false
end

imagesc(timestamp,1:size(traces,1),traces);
if ~options.flipYAxis
    set(gca,'YDir','normal');
end
colorbar; box off

end