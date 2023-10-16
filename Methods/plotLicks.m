function [lickRate,lickTraces,lickEvents] = plotLicks(eventIdx,timeRange,binSize,color,leftLick,rightLick,params,options)

% Plot lick raster or lick rate traces
arguments
    eventIdx double
    timeRange double
    binSize double
    color % should provide two colors in cell array if side = [1 1]
    leftLick logical
    rightLick logical
    params struct
    options.mode string = 'trace'
    options.side double = [0 1] % only return right lick
    options.leftSolenoid double
    options.rightSolenoid double
    options.airpuff double
end

[lickRate,lickTraces,lickEvents] = getLicks(timeRange,eventIdx,binSize,leftLick,rightLick,...
                                params.sync.behaviorFs,params.sync.timeNI,side=options.side);


%if isempty(lickRate{find(options.side)}); return; end

if strcmp(options.mode,'trace')    
    % draw lick rate traces
    if all(options.side == [1 1])
        t = linspace(timeRange(1),timeRange(2),size(lickRate{1},2));
        plotSEM(t,lickRate{1},color{1});
        plotSEM(t,lickRate{2},color{1});
    else
        t = linspace(timeRange(1),timeRange(2),size(lickRate,2));
        plotSEM(t,lickRate,color);
    end

% elseif strcmp(options.mode,'raster')
    

    
end