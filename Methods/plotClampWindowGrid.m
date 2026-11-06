function fig = plotClampWindowGrid(topWindows,topTraces,options)
%PLOTCLAMPWINDOWGRID Screening figure: DA alone for the best candidate windows.
%
% One row per photometry channel and its ranked windows across that row, as
% returned by selectTopClampWindows. Deliveries are labelled with plotEvent and
% clamp epochs are shaded, but the clamp index is not extracted: this figure is
% for picking a window, which the DA / clamp index cell then plots in detail.
%
% A Label column on topWindows adds a second title line per tile (whatever the
% caller ranked on), and ChannelName names the row above its first tile.

arguments
    topWindows table
    topTraces cell
    options.events table = table()
    options.clampEpochs double = zeros(0,2)
    options.channelIdx double = []      % row order; defaults to the order in topWindows
    options.airpuffDuration (1,1) double = 0.2
    options.titleText string = ""
end

if height(topWindows) ~= numel(topTraces)
    error('plotClampWindowGrid:SizeMismatch','Pass one extracted trace per window.');
end
if ~all(ismember({'Start_sec','End_sec','ChannelIndex','Rank'},topWindows.Properties.VariableNames))
    error('plotClampWindowGrid:MissingColumns', ...
        'topWindows needs Start_sec, End_sec, ChannelIndex and Rank.');
end

channels = options.channelIdx;
if isempty(channels); channels = unique(topWindows.ChannelIndex,'stable'); end
nColumns = max(topWindows.Rank);
hasLabel = ismember('Label',topWindows.Properties.VariableNames);
hasName = ismember('ChannelName',topWindows.Properties.VariableNames);

fig = initializeFig(1,0.8,ScreenSize=get(groot,'ScreenSize'),Visible='on');
tiledlayout(numel(channels),nColumns,TileSpacing='compact',Padding='compact');
for r = 1:numel(channels)
    for c = 1:nColumns
        ax = nexttile;
        row = find(topWindows.ChannelIndex == channels(r) & topWindows.Rank == c,1);
        if isempty(row); axis(ax,'off'); continue; end
        windowEvents = options.events;
        if ~isempty(windowEvents)
            windowEvents = windowEvents(windowEvents.Time_sec >= topWindows.Start_sec(row) & ...
                windowEvents.Time_sec <= topWindows.End_sec(row),:);
        end
        tileTitle = string(sprintf('#%d | %.0f-%.0f s',c, ...
            topWindows.Start_sec(row),topWindows.End_sec(row)));
        if hasLabel; tileTitle = [tileTitle;string(topWindows.Label(row))]; end %#ok<AGROW>
        % The channel name goes above the first tile of its row
        if c == 1
            if hasName; channelName = string(topWindows.ChannelName(row));
            else; channelName = string(topTraces{row}.DAName); end
            tileTitle = [channelName;tileTitle]; %#ok<AGROW>
        end
        plotSampleClampTrace(topTraces{row},ax=ax,plotClampIndex=false,events=windowEvents, ...
            clampEpochs=options.clampEpochs,airpuffDuration=options.airpuffDuration, ...
            xlim=[topWindows.Start_sec(row),topWindows.End_sec(row)],daTitle=tileTitle);
    end
end
if strlength(options.titleText) > 0; sgtitle(options.titleText,'Interpreter','none'); end
end
