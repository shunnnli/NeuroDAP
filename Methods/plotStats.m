function plotStats(data1,data2,xloc,options)

arguments
    data1 double
    data2 double
    xloc double = [0 1] % must be size [1,2]
    
    options.testType string = 'kstest'
    options.textType string = '*'
    
    options.LineWidth double = 2 
    options.LineLengthProportion double = 0.75
    options.FontSize double = 12
    options.ylim double = nan
    options.yPortion double = 1 % porportion of the y axis
    options.yloc double
    options.yPortionNoise double = 0.1
end

%% Caculate significance
if strcmpi(options.testType,'kstest')
    [~,pvalue,~] = kstest2(data1,data2);
elseif strcmpi(options.testType,'ranksum')
    pvalue = ranksum(data1,data2);
elseif strcmpi(options.testType,'signrank')
    pvalue = signrank(data1,data2);
elseif contains(options.testType,'ttest')
    [~,pvalue] = ttest2(data1,data2);
end

%% Plot significance bar
% determine x
xlineLength = (xloc(2)-xloc(1)) * options.LineLengthProportion;
xlineCenter = mean(xloc);
xlineStart = xlineCenter - xlineLength/2;
xlineEnd = xlineCenter + xlineLength/2;

% determine y
[ylimit_lower,ylimit_upper] = bounds([data1;data2],'all');
if isfield(options,'yloc')
    height = options.yloc;
else
    yPortion = options.yPortion+rand*options.yPortionNoise;
    if abs(ylimit_upper) >= abs(ylimit_lower)
        height = yPortion*(ylimit_upper-ylimit_lower) + ylimit_lower;
    else
        height = yPortion*(ylimit_lower-ylimit_upper) + ylimit_upper;
    end
end

plot([xlineStart xlineEnd],ones(1,2)*height,'Color','k','LineWidth',options.LineWidth); hold on
if strcmpi(options.textType,'*')
    if pvalue > 0.05
        stars = 'n.s.';
    elseif pvalue <= 1e-4
        stars = '****';
    elseif pvalue <= 1e-3
        stars = '***';
    elseif pvalue <= 1e-2
        stars = '**';
    else
        stars = '*';
    end
    text(xlineCenter,height,stars,'FontSize',options.FontSize,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','Bottom');
else
    text(xlineCenter,height,strcat("p = ",num2str(round(pvalue,4))),'FontSize',options.FontSize,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','Bottom');
end

end