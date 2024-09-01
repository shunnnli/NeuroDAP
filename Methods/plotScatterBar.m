function plotScatterBar(data,x,options)

arguments
    data double
    x double % value on the x axis

    options.color 
    options.dotSize double = 70;

    options.MarkerFaceAlpha double = 0.8
    options.XJitter = 'density'
    options.XJitterWidth double = 0.5
    options.LineWidth double = 2
end

box off;

xgroupdata = x * ones(size(data,1),1);
boxchart(xgroupdata,data,'BoxFaceColor',options.color,'WhiskerLineColor',options.color,...
         LineWidth=options.LineWidth); 

hold on;

swarmchart(xgroupdata,data,...
    options.dotSize,options.color,'filled',...
    'MarkerFaceAlpha',options.MarkerFaceAlpha,...
    'XJitter',options.XJitter,'XJitterWidth',options.XJitterWidth); hold on;

end