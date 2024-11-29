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

    options.orientation string = 'vertical'
end

% determine x
xlineLength = (xloc(2)-xloc(1)) * options.LineLengthProportion;
xlineCenter = mean(xloc);
xlineStart = xlineCenter - xlineLength/2;
xlineEnd = xlineCenter + xlineLength/2;

% determine y
if strcmpi(options.orientation,'vertical')
    ylimit = ylim; 
    if isfield(options,'yloc')
        height = options.yloc;
    else
        height = yPortion*(ylimit(2)-ylimit(1)) + ylimit(1);
    end
else
    xlimit = xlim; 
    if isfield(options,'yloc')
        height = options.yloc;
    else
        height = yPortion*(xlimit(2)-xlimit(1)) + xlimit(1);
    end
end

if strcmpi(options.orientation,'vertical')
    plot([xlineStart xlineEnd],ones(1,2)*height,'Color','k','LineWidth',options.LineWidth);
    text(xlineCenter,height,strcat("p = ",num2str(round(pvalue,4))),'FontSize',options.FontSize,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','Bottom');
else
    plot(ones(1,2)*height,[xlineStart xlineEnd],'Color','k','LineWidth',options.LineWidth);
    text(height,xlineCenter,strcat("p = ",num2str(round(pvalue,4))),'FontSize',options.FontSize,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','Bottom');
end

end