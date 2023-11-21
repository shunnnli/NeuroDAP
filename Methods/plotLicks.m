function [lickRate,lickEvents,t] = plotLicks(eventIdx,timeRange,binSize,color,leftLick,rightLick,params,options)

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

    options.plot logical = true
    options.markerSize double = 20
end

if strcmpi(params.session.baselineSystem,'lj')
    Fs = params.sync.labjackFs;
    timeBaseline = params.sync.timePhotometry;
else
    Fs = params.sync.behaviorFs;
    timeBaseline = params.sync.timeNI;
end

% getLicks
[lickRate,~,lickEvents] = getLicks(timeRange,eventIdx,binSize,leftLick,rightLick,...
                                Fs,timeBaseline,side=options.side);

% Create time
if all(options.side == [1 1])
    t = linspace(timeRange(1),timeRange(2),size(lickRate{1},2));
else
    t = linspace(timeRange(1),timeRange(2),size(lickRate,2));
end

%if isempty(lickRate{find(options.side)}); return; end

if options.plot
    if strcmp(options.mode,'trace')    
        % draw lick rate traces
        if all(options.side == [1 1])
            t = linspace(timeRange(1),timeRange(2),size(lickRate{1},2));
            plotSEM(t,lickRate{1},color{1});
            plotSEM(t,lickRate{2},color{2});
        else
            t = linspace(timeRange(1),timeRange(2),size(lickRate,2));
            plotSEM(t,lickRate,color);
        end
    
    elseif strcmp(options.mode,'raster')
        if all(options.side == [1 1])
            for i = 1:size(lickEvents,1)
                scatter(lickEvents{i,1},i,options.markerSize,'filled','MarkerFaceColor',color{1}); hold on
                scatter(eventIdx(i)/params.sync.behaviorFs,i,options.markerSize+10,color{1},'pentagram','filled'); hold on
            end
            for i = 1:size(lickEvents,1)
                scatter(lickEvents{i,2},i,options.markerSize,'filled','MarkerFaceColor',color{2}); hold on
                scatter(eventIdx(i)/params.sync.behaviorFs,i,options.markerSize+10,color{2},'pentagram','filled'); hold on
            end
        else
            for i = 1:size(lickEvents,1)
                scatter(lickEvents{i},i,options.markerSize,'filled','MarkerFaceColor',color); hold on
                scatter(eventIdx(i)/params.sync.behaviorFs,i,options.markerSize+10,color,'pentagram','filled'); hold on
            end
        end
    end
end
    
end