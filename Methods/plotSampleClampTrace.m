function [axDA,axClamp,fig] = plotSampleClampTrace(sample,options)
%PLOTSAMPLECLAMPTRACE Plot one DA (and optionally clamp index) example window.
%
% sample is the struct returned by getSampleClampTrace. Traces are drawn
% against session time, so the interval between deliveries stays visible.
% clampEpochs shades the stretches when the clamp was on, so how much of the
% window was clamped can be read off the DA panel itself.
% Water, airpuff and tone deliveries are marked with plotEvent: the patch
% width is the delivery duration (airpuff defaults to 0.2 s) and the color
% identifies the delivery type. Clamp state appears in parentheses in the label.
%
% Pass ax= to draw DA into existing axes (one tile per window in a summary
% figure); the clamp index needs its own tile, so set plotClampIndex=false.

arguments
    sample (1,1) struct
    options.ax = []                      % target axes for DA; [] opens a figure
    options.plotClampIndex (1,1) logical = true
    options.events table = table()       % Type, IsClamp, Time_sec (+ optional InScore)
    options.clampEpochs double = zeros(0,2) % [start,end] session seconds to shade as clamped
    options.trials table = table()
    options.trialRows double = []
    options.timeNI double = []
    options.xlim double = []             % [] uses the extracted DA window
    options.daColor (1,3) double = [.232 .76 .58]
    options.clampColor (1,3) double = [.232 .76 .58]
    options.unclampColor (1,3) double = [.6 .3 .1]
    options.clampIndexColor (1,3) double = [.3 .3 .3]
    options.waterDuration (1,1) double {mustBeNonnegative} = 0
    options.airpuffDuration (1,1) double {mustBeNonnegative} = 0.2
    options.toneDuration (1,1) double {mustBeNonnegative} = 0
    options.daTitle string = ""          % "" uses the channel name; 2 elements = 2 lines
    options.titleText string = ""        % figure title; skipped when drawing into ax
    options.labelEvents (1,1) logical = true
end

if ~all(isfield(sample,{'DA','sessionTimeDA_sec','DAName','DAUnits'}))
    error('plotSampleClampTrace:InvalidSample','Pass the struct returned by getSampleClampTrace.');
end
if ~isempty(options.trialRows) && (isempty(options.timeNI) || isempty(options.trials))
    error('plotSampleClampTrace:MissingTrials','Trial markers need both trials and timeNI.');
end

% Axes: a new two-panel figure, or the caller's tile for DA alone
fig = []; axClamp = [];
if isempty(options.ax)
    if options.plotClampIndex
        fig = initializeFig(0.65,0.65,ScreenSize=get(groot,'ScreenSize'),Visible='on');
        tiledlayout(2,1); axDA = nexttile;
    else
        fig = initializeFig(0.65,0.4,ScreenSize=get(groot,'ScreenSize'),Visible='on');
        axDA = nexttile(tiledlayout(1,1));
    end
else
    if options.plotClampIndex
        error('plotSampleClampTrace:ClampIndexNeedsLayout', ...
            'Pass plotClampIndex=false when drawing into existing axes.');
    end
    axDA = options.ax;
end

axes(axDA); hold on
plotTraces(sample.DA,sample.sessionTimeDA_sec,color=options.daColor,plotPatch=false);
daTitle = options.daTitle;
if isempty(daTitle) || all(strlength(daTitle) == 0); daTitle = string(sample.DAName); end
title(axDA,daTitle,'Interpreter','none');
xlabel(axDA,'Time (s)'); ylabel(axDA,sample.DAUnits);
xRange = options.xlim;
if isempty(xRange); xRange = [sample.sessionTimeDA_sec(1),sample.sessionTimeDA_sec(end)]; end
xlim(axDA,xRange);

if options.plotClampIndex
    axClamp = nexttile;
    plotTraces(sample.clampIndex,sample.sessionTimeClamp_sec,color=options.clampIndexColor,plotPatch=false);
    if isfield(sample,'removeClampBaseline') && sample.removeClampBaseline
        title(axClamp,'Clamp index (red - blue, baseline subtracted)');
    else
        title(axClamp,'Clamp index (red - blue)');
    end
    xlabel(axClamp,'Time (s)'); ylabel(axClamp,'Clamp index (%)');
    xlim(axClamp,xRange);
end

% Markers go on every panel so the two can be read against each other
[~,~,~,~,~,~,bluePurpleRed] = loadColors;
axList = axDA; if ~isempty(axClamp); axList = [axDA,axClamp]; end
for ax = axList
    axes(ax); hold on %#ok<LAXES>
    for i = 1:numel(options.trialRows)
        cueTime_sec = options.timeNI(options.trials.CueTime(options.trialRows(i)))-options.timeNI(1);
        xline(ax,cueTime_sec,':',sprintf('Trial %g',options.trials.TrialNumber(options.trialRows(i))), ...
            'HandleVisibility','off');
    end
    % Clamp epochs first, so the delivery lines stay on top of the shading
    for i = 1:size(options.clampEpochs,1)
        epoch = [max(options.clampEpochs(i,1),xRange(1)),min(options.clampEpochs(i,2),xRange(2))];
        if diff(epoch) <= 0; continue; end
        plotEvent("",diff(epoch),color=options.clampColor,x=epoch(1), ...
            FaceAlpha=0.12,shadeOnly=true);
    end
    if ~options.labelEvents; continue; end
    for i = 1:height(options.events)
        type = string(options.events.Type(i));
        switch upper(type)
            case "WATER";   duration = options.waterDuration; eventColor = bluePurpleRed(1,:);
            case "AIRPUFF"; duration = options.airpuffDuration; eventColor = [0.2,0.2,0.2];
            case "TONE";    duration = options.toneDuration; eventColor = bluePurpleRed(250,:);
            otherwise;      duration = 0; eventColor = [0.2,0.2,0.2];
        end
        if options.events.IsClamp(i); state = "clamp";
        else; state = "unclamp"; end
        label = sprintf('%s (%s)',type,state);
        % * marks events shown but left out of scoring (partial or missing windows)
        if ismember('InScore',options.events.Properties.VariableNames) && ~options.events.InScore(i)
            label = sprintf('%s *',label);
        end
        plotEvent(label,duration,color=eventColor,x=options.events.Time_sec(i));
    end
    xlim(ax,xRange);
end

if ~isempty(axClamp)
    linkaxes([axDA,axClamp],'x'); xlim(axDA,xRange);
end
if ~isempty(fig) && strlength(options.titleText) > 0
    sgtitle(options.titleText,'Interpreter','none');
end
end
