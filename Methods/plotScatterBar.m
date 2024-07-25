function plotScatterBar(data,x,options)

arguments
    data double
    x double % value on the x axis

    options.color 
    options.dotSize double = 70;

    options.MarkerFaceAlpha double = 0.8
    options.XJitter = 'density'
    options.XJitterWidth double = 0.5

end

xgroupdata = x * ones(size(data,1),1);
boxchart(xgroupdata,data,'BoxFaceColor',options.color); hold on
swarmchart(xgroupdata,data,...
    options.dotSize,options.color,'filled',...
    'MarkerFaceAlpha',options.MarkerFaceAlpha,...
    'XJitter',options.XJitter,'XJitterWidth',options.XJitterWidth); hold on;

end