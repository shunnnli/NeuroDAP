function plotEvent(label,duration,options)

% Plot event at PSTH with vertical line, label, and patch indicating
% duration

% parse input
arguments
    label string
    duration double
    options.color = 'r'
    options.x double = 0
    options.unit string = 's'
    options.FaceAlpha double = 0.1
    options.LineWidth double = 1.5
    options.LineStyle string = '-'
    options.shadeOnly logical = false

    options.percentY double = 100
    options.zeroValue double = 0
end

% Check color is valid
if ~ischar(options.color) && sum(options.color) > 3 && size(options.color,1) == 1
    options.color = options.color ./ 255;
end

% Check unit
if strcmpi(options.unit,'ms')
    duration = duration / 1000;
end

% If line do not span the whole y axis (percentY < 100), check whether
% zeroValue is provided
if options.percentY < 100
    if ~isfield(options,'zeroValue')
        warning('plotEvent: value at zero not provided when percentY < 100!');
        options.percentY = 100;
    end
end

% Get current axis limit
xlimit = xlim; left = xlimit(1); right = xlimit(2);
ylimit = ylim; yLength = abs(ylimit(2)-ylimit(1)) * options.percentY / 100;
if options.percentY == 100
    bottom = ylimit(1); top = ylimit(2);
else
    bottom = options.zeroValue - yLength/2; 
    top = options.zeroValue + yLength/2;
end

% Plot events
patch([options.x duration duration options.x],[bottom bottom top top],options.color,...
        'FaceAlpha',options.FaceAlpha,'EdgeColor','none','HandleVisibility','off'); hold on
if ~options.shadeOnly
    if options.percentY == 100
        xline(options.x,options.LineStyle,label,'Color',options.color,'LineWidth',options.LineWidth,'HandleVisibility','off'); hold on
    else
        line([options.x,options.x],[bottom,top],...
                LineStyle=options.LineStyle,Color=options.color,...
                LineWidth=options.LineWidth,HandleVisibility='off'); hold on
        % text(options.x+xlim*0.01,top,label,Color=options.color,Rotation=90); hold on
    end
end

% Adjust x and y axis limits
xlim([left right]);
if options.percentY == 100; ylim([bottom top]); end
box off

end