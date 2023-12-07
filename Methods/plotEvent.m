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
end

% Check color is valid
if ~ischar(options.color) && sum(options.color) > 3 && size(options.color,1) == 1
    options.color = options.color ./ 255;
end

if strcmpi(options.unit,'ms')
    duration = duration / 1000;
end

% Get current axis limit
xlimit = xlim; left = xlimit(1); right = xlimit(2);
ylimit = ylim; bottom = ylimit(1); top = ylimit(2);

% Plot events
patch([options.x duration duration options.x],[bottom bottom top top],options.color,...
        'FaceAlpha',options.FaceAlpha,'EdgeColor','none','HandleVisibility','off'); hold on
xline(options.x,options.LineStyle,label,'Color',options.color,'LineWidth',options.LineWidth,'HandleVisibility','off'); hold on

xlim([left right]);
ylim([bottom top]);
box off

end