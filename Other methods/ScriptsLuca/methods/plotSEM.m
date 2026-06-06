function plotSEM(x,y,color,options)

% parse inputs
arguments
    x double
    y double
    color

    options.smooth double = 0;
    options.smoothMethod string = 'movmean';

    options.meanOnly logical = false
    options.plotIndividual logical = false % plot individual trace in the background
    options.individualColor = 'gray' % Color of individual trace
    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 2
    options.plotStyle string = 'line'
    options.delta double = []
    options.individualAlpha double = 0
end

%% Set up
if isempty(y); return; end

% Check color is valid
if ~isstring(color)
    if any(sum(color)>3) && size(color,1) == 1
        color = color ./ 255;
    end
end

% Check individual color
if ~isstring(options.individualColor) && ~ischar(options.individualColor)
    if any(sum(options.individualColor)>3) && size(options.individualColor,1) == 1
        options.individualColor = options.individualColor ./ 255;
    end
else
    if strcmp(options.individualColor,'same')
        options.individualColor = color;
    elseif strcmp(options.individualColor,'gray')
        options.individualColor = [0.8, 0.8, 0.8];
    end
end

% Smooth data if neccessary
if options.smooth ~= 0
    y = smoothdata(y,2,options.smoothMethod,options.smooth);
end

% N = nnz(~isnan(y(:,1)));   % Number of ‘Experiments’ In Data Set
yMean = mean(y,1,"omitnan");      % Mean Of All Experiments At Each Value Of ‘x’

if isempty(options.delta)
    ySEM = zeros();
    for i=1:size(y,2)
        ySEM(i) = getSEM(y(:,i));
    end
else
    ySEM = options.delta;
end

% Update individual color
if options.individualAlpha ~= 0 
    options.individualColor = [options.individualColor,options.individualAlpha];
end

%% Plot
if strcmp(options.plotStyle,'line')
    if options.plotIndividual
        for i = 1:size(y,1)
            plot(x,y(i,:),'Color',options.individualColor,'LineWidth',options.LineWidth-1,'LineStyle',options.LineStyle,'HandleVisibility','off'); hold on
        end
    end
    plot(x,yMean,'Color',color,'LineWidth',options.LineWidth,'LineStyle',options.LineStyle); hold on
elseif strcmp(options.plotStyle,'stairs')
    stairs(x,yMean,'Color',color,'LineWidth',options.LineWidth,'LineStyle',options.LineStyle); hold on
end

if ~options.meanOnly
    patch('XData',[x, fliplr(x)], 'YData',[yMean-ySEM, fliplr(yMean+ySEM)], ...
        'FaceColor',color,'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
    hold on
end
box off

end % plotSEM