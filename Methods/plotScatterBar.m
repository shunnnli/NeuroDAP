function plotScatterBar(x,data,options)

arguments
    x double % value on the x axis
    data double

    options.color 
    options.dotSize double = 70;

    options.style string = 'box' % can also be bar, which uses a bar
    options.BarFaceOpacity double = 0.2 % 1 is not transparent, 0 is fully transparent
    options.orientation string = 'vertical';
    
    options.connectPairs logical = true
    options.connectColor = [0.75,0.75,0.75]

    options.plotScatter logical = true
    options.MarkerFaceAlpha double = 0.8
    options.XJitter = 'density'
    options.XJitterWidth double = 0.5
    options.LineWidth double = 2
end

%% Check inputs
if isvector(data)
    if isscalar(data)
        if ~isfield(options,'color'); options.color = [.2 .2 .2]; end
    else
        if ~isfield(options,'color')
            options.color = cell(numel(data),1);
            for c = 1:numel(data)
                options.color{c} = rand(1,3);
            end
        end
    end
elseif size(data,2) == 2
    % Check whether x is in correct dimension
    if length(x) ~= 2
        warning('x should be a vector with 2 elements indicating where each column is on the x axis! Reset to x=[1,2]');
        x = [1,2];
    end

    if ~isfield(options,'color')
        warning('No color provided, reset to default red & black');
        options.color = [1 0.196 0.227; .2 .2 .2];
    end

    if ismatrix(options.color) && all(size(options.color) ~= [2,3])
        % Make color the correct dimension
        if all(size(options.color) == [3,2])
            options.color = transpose(options.color);
        else
            warning('Incorrect color dimension, reset to default red & black');
            options.color = [1 0.196 0.227; .2 .2 .2];
        end
    end

    if iscell(options.color) && length(options.color) ~=2
        warning('Incorrect color dimension, reset to default red & black');
        options.color = [1 0.196 0.227; .2 .2 .2];
    end

    % Repackage options.color matrix into cells
    color = cell(2,1);
    for c = 1:size(options.color,1)
        color{c} = options.color(c,:);
    end
    options.color = color;
end

%% Plot data

if isvector(data)
    xgroupdata = x * ones(size(data,1),1);
    if contains(options.style,'box',IgnoreCase=true)
        boxchart(xgroupdata,data,'BoxFaceColor',options.color,'WhiskerLineColor',options.color,...
                 LineWidth=options.LineWidth,...
                 Orientation=options.orientation); 
        hold on;
    elseif contains(options.style,'bar',IgnoreCase=true)
        barFaceColor = 1 - options.BarFaceOpacity*(1-options.color);
        if strcmp(options.orientation,'vertical')
            bar(x,mean(data,'all'),FaceColor=barFaceColor,EdgeColor=options.color,LineWidth=options.LineWidth,HandleVisibility="off");
        else
            barh(x,mean(data,'all'),FaceColor=barFaceColor,EdgeColor=options.color,LineWidth=options.LineWidth,HandleVisibility="off");
        end
        hold on;
    
        SEM = getSEM(data);
        eb = errorbar(x,mean(data,'all'),-SEM,SEM,options.orientation,HandleVisibility="off");
        eb.Color = options.color;                            
        eb.LineWidth = options.LineWidth;  
        hold on;
    end
    
    if options.plotScatter
        if strcmp(options.orientation,'vertical')
            swarmchart(xgroupdata,data,...
                options.dotSize,options.color,'filled',...
                'MarkerFaceAlpha',options.MarkerFaceAlpha,...
                'XJitter',options.XJitter,'XJitterWidth',options.XJitterWidth); 
            hold on;
        else
            swarmchart(data,xgroupdata,...
                    options.dotSize,options.color,'filled',...
                    'MarkerFaceAlpha',options.MarkerFaceAlpha,...
                    'YJitter',options.XJitter,'YJitterWidth',options.XJitterWidth);
            hold on;
        end
    end

elseif size(data,2) == 2
    
    % Plot box or bar for each column
    for col = 1:size(data,2)
        colData = data(:,col);

        xgroupdata = x(col) * ones(size(colData,1),1);
        if contains(options.style,'box',IgnoreCase=true)
            boxchart(xgroupdata,colData,'BoxFaceColor',options.color{col},'WhiskerLineColor',options.color{col},...
                     LineWidth=options.LineWidth,...
                     Orientation=options.orientation,HandleVisibility="off"); 
            hold on;
        elseif contains(options.style,'bar',IgnoreCase=true)
            barFaceColor = 1 - options.BarFaceOpacity*(1-options.color{col});
            if strcmp(options.orientation,'vertical')
                bar(x(col),mean(colData,'all'),FaceColor=barFaceColor,EdgeColor=options.color{col},LineWidth=options.LineWidth,HandleVisibility="off");
            else
                barh(x(col),mean(colData,'all'),FaceColor=barFaceColor,EdgeColor=options.color{col},LineWidth=options.LineWidth,HandleVisibility="off");
            end
            hold on;
        
            SEM = getSEM(colData);
            eb = errorbar(x(col),mean(colData,'all'),-SEM,SEM,options.orientation,HandleVisibility="off");
            eb.Color = options.color{col};                            
            eb.LineWidth = options.LineWidth;  
            hold on;
        end  
    end

    % Pre calculate jitter
    jitterData = cell(1,2);
    for col = 1:2
        colData = data(:,col);
        % Base x position for this group
        xBase = x(col) * ones(size(colData));
        % Compute jitter: here we approximate density jitter with random offsets 
        jitterOffset = (rand(size(colData))-0.5) * options.XJitterWidth;
        jitterData{col} = xBase + jitterOffset;
    end

    % Connect pairs using the computed jittered positions.
    if options.connectPairs
        nPoints = size(data,1);
        for i = 1:nPoints
            if strcmp(options.orientation,'vertical')
                if options.plotScatter
                    x1 = jitterData{1}(i);
                    x2 = jitterData{2}(i);
                else
                    x1 = x(1);
                    x2 = x(2);
                end
                y1 = data(i,1);
                y2 = data(i,2);
                plot([x1, x2], [y1, y2], 'Color', options.connectColor, 'LineWidth', options.LineWidth,HandleVisibility="off");
                hold on;
            else
                if options.plotScatter
                    y1 = jitterData{1}(i);
                    y2 = jitterData{2}(i);
                else
                    y1 = x(1);
                    y2 = x(2);
                end
                x1 = data(i,1);
                x2 = data(i,2);
                plot([x1, x2], [y1, y2], 'Color', options.connectColor, 'LineWidth', options.LineWidth,HandleVisibility="off");
                hold on;
            end
        end
    end

    % Plot scatter points
    if options.plotScatter
        for col = 1:2
            colData = data(:,col);
            if strcmp(options.orientation,'vertical')
                % Use scatter (or you can use plot with marker) to plot at the jittered positions
                scatter(jitterData{col}, colData, options.dotSize, options.color{col}, 'filled',...
                    'MarkerFaceAlpha', options.MarkerFaceAlpha);
                hold on;
            else
                scatter(colData, jitterData{col}, options.dotSize, options.color{col}, 'filled',...
                    'MarkerFaceAlpha', options.MarkerFaceAlpha);
                hold on;
            end
        end
    end

else
    error('Data is neither vector nor matrix with two columns, check again!');
end

box off;
end
