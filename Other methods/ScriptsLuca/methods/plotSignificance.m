function plotSignificance(pvalue,xloc,yPortion,options)

arguments
    pvalue
    xloc double = [0 1] % must be size [1,2]
    yPortion double = 1 % porportion of the y axis
    
    options.LineWidth double = 2
    options.LineLengthProportion double = 0.75
    options.FontSize double = 12
    options.ylim double = nan
    options.yloc double
end

% determine x
xlineLength = (xloc(2)-xloc(1)) * options.LineLengthProportion;
xlineCenter = mean(xloc);
xlineStart = xlineCenter - xlineLength/2;
xlineEnd = xlineCenter + xlineLength/2;

% determine y
ylimit = ylim; 
if isfield(options,'yloc')
    height = options.yloc;
else
    height = yPortion*(ylimit(2)-ylimit(1)) + ylimit(1);
end

plot([xlineStart xlineEnd],ones(1,2)*height,'Color','k','LineWidth',options.LineWidth);
text(xlineCenter,height,strcat("p = ",num2str(round(pvalue,4))),'FontSize',options.FontSize,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','Bottom');

end