function varargout = plotSEM(x,y,color,options)

arguments
    x double
    y double
    color

    options.smooth double = 0;
    options.smoothMethod string = 'movmean';

    options.plotMean logical = true
    options.plotPatch logical = true % requires plotMean is true
    
    options.opacity double = 1 % 1 is not transparent, 0 is fully transparent

    options.plotIndividual logical = false % plot individual trace in the background
    options.individualColor = 'same' % Color of individual trace
    options.individualAlpha double = 0.3 % 1 is not transparent, 0 is fully transparent

    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 2
    options.plotStyle string = 'line'
    options.delta double = []
    options.label string = ''
    options.plotStyleIndividual string = 'overlay'
end

%% Set up
if isempty(y)
    l = nan; varargout{1} = l;
    return; 
end

% Check color is valid
if ~isstring(color)
    if any(sum(color)>3) && size(color,1) == 1
        color = color ./ 255;
    end
end

% Change color opacity as indiacted
color = 1 - options.opacity*(1-color);

% Check individual color
if ~isstring(options.individualColor) && ~ischar(options.individualColor)
    if any(sum(options.individualColor)>3) && size(options.individualColor,1) == 1
        options.individualColor = options.individualColor ./ 255;
    end
else
    if strcmp(options.individualColor,'same')
        options.individualColor = 1 - options.individualAlpha*(1-color);
    elseif strcmp(options.individualColor,'gray')
        options.individualColor = [0.6, 0.6, 0.6];
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

%% Plot
if strcmp(options.plotStyle,'line')
    if options.plotIndividual
        if strcmp(options.plotStyleIndividual,'overlay')
            for i = 1:size(y,1)
                l = plot(x,y(i,:),'Color',options.individualColor,'LineWidth',max(0.01,options.LineWidth-1),'LineStyle',options.LineStyle,'HandleVisibility','off'); hold on
            end
        elseif strcmp(options.plotStyleIndividual,'stack')
            options.plotMean = false; % By default do not plot mean
            % Get max and min value
        end
    end
    if options.plotMean
        l = plot(x,yMean,'Color',color,'LineWidth',options.LineWidth,'LineStyle',options.LineStyle,'DisplayName',options.label); hold on
    end
elseif strcmp(options.plotStyle,'stairs')
    if options.plotMean
        l = stairs(x,yMean,'Color',color,'LineWidth',options.LineWidth,'LineStyle',options.LineStyle,'DisplayName',options.label); hold on
    end
end

if options.plotMean && options.plotPatch
    patch('XData',[x, fliplr(x)], 'YData',[yMean-ySEM, fliplr(yMean+ySEM)], ...
        'FaceColor',color,'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
    hold on
end
box off

varargout{1} = l;

end % plotSEM