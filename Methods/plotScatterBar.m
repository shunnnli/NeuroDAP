function plotScatterBar(data,x,options)

arguments
    data double
    x double % value on the x axis

    options.color 
    options.dotSize double = 70;

    options.style string = 'box' % can also be bar, which uses a bar
    options.BarFaceOpacity double = 0.2 % 1 is not transparent, 0 is fully transparent
    
    options.connectPairs logical = true
    options.connectColor = [0.75,0.75,0.75]

    options.plotScatter logical = true
    options.MarkerFaceAlpha double = 0.8
    options.XJitter = 'density'
    options.XJitterWidth double = 0.5
    options.LineWidth double = 2
end

%% Check inputs

if size(data,2) == 2
    % Check whether x is in correct dimention
    if length(x) ~= 2
        warning('x should be a vector with 2 elements indicating where each column is on the x axis! Reset to x=[1,2]');
        x = [1,2];
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
                 LineWidth=options.LineWidth); 
        hold on;
    elseif contains(options.style,'bar',IgnoreCase=true)
        barFaceColor = 1 - options.BarFaceOpacity*(1-options.color);
        bar(x,mean(data,'all'),FaceColor=barFaceColor,EdgeColor=options.color,LineWidth=options.LineWidth);
        hold on;
    
        SEM = getSEM(data);
        eb = errorbar(x,mean(data,'all'),-SEM,SEM);
        eb.Color = options.color;                            
        eb.LineWidth = options.LineWidth;  
        hold on;
    end
    
    if options.plotScatter
        swarmchart(xgroupdata,data,...
            options.dotSize,options.color,'filled',...
            'MarkerFaceAlpha',options.MarkerFaceAlpha,...
            'XJitter',options.XJitter,'XJitterWidth',options.XJitterWidth); 
        hold on;
    end

elseif size(data,2) == 2
    
    for col = 1:size(data,2)
        colData = data(:,col);

        xgroupdata = x(col) * ones(size(colData,1),1);
        if contains(options.style,'box',IgnoreCase=true)
            boxchart(xgroupdata,colData,'BoxFaceColor',options.color{col},'WhiskerLineColor',options.color{col},...
                     LineWidth=options.LineWidth); 
            hold on;
        elseif contains(options.style,'bar',IgnoreCase=true)
            barFaceColor = 1 - options.BarFaceOpacity*(1-options.color{col});
            bar(x(col),mean(colData,'all'),FaceColor=barFaceColor,EdgeColor=options.color{col},LineWidth=options.LineWidth);
            hold on;
        
            SEM = getSEM(colData);
            eb = errorbar(x(col),mean(colData,'all'),-SEM,SEM);
            eb.Color = options.color{col};                            
            eb.LineWidth = options.LineWidth;  
            hold on;
        end  
    end

    % Connect pairs
    if options.connectPairs
        plot(x,data,color=options.connectColor,LineWidth=options.LineWidth); hold on;
    end

    % Plot scatter
    if options.plotScatter
        for col = 1:size(data,2)
            colData = data(:,col);
            xgroupdata = x(col) * ones(size(colData,1),1);
            swarmchart(xgroupdata,colData,...
                options.dotSize,options.color{col},'filled',...
                'MarkerFaceAlpha',options.MarkerFaceAlpha,...
                'XJitter',options.XJitter,'XJitterWidth',options.XJitterWidth); 
            hold on;
        end
    end

else
    error('Data is neither vector nor matrix with two columns, check again!');
end

box off;
end