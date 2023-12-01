function plotHeatmap(traces,t,options)

arguments
    traces double
    t double
    options.colormap double
end

imagesc(t,1:size(traces,1),traces);
set(gca,'YDir','normal');
colorbar; box off

end