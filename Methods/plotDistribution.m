function plotDistribution(data,options)

arguments
    data double

    options.plotGroup logical = true % vector, plot this group or not
    options.groupIdx cell

    options.color

    options.nboot double = 5000

    options.masterlayout
    options.tile double = 1
    options.tileSpan double = [1 1];
    options.dotSize double = 100
    options.xlabel string = ''
end

%% Check input
if ~isfield(options,'masterlayout')
    master = tiledlayout(1,1);
else; master = options.masterlayout; 
end

if ~isfield(options,'groupIdx')
    options.groupIdx = {1:size(data,1)};
end

if ~isfield(options,'color')
    options.color = {[.7 .7 .7]};
end

if length(options.plotGroup) ~= size(options.groupIdx)
    error('plotGroup and groupIdx have different size!');
end

%% Plot data
children = tiledlayout(master,4,1);
children.Layout.Tile = options.tile; children.Layout.TileSpan = options.tileSpan;
children.TileSpacing = 'compact'; children.Padding = 'tight'; axis off;
nexttile(children,1);
for group = 1:length(options.plotGroup)
    if options.plotGroup(group)
        groupIdx = options.groupIdx{group};
        groupColor = options.color{group};
        plotScatterBar(data(groupIdx),group,color=groupColor,dotSize=options.dotSize,orientation='horizontal'); hold on;
    end
end
% Plot statistical test
for group = 1:length(options.plotGroup)
    if options.plotGroup(group)
        for nextGroup = group+1:length(options.plotGroup)
            if options.plotGroup(nextGroup)
                [~,p,~] = kstest2(data(options.groupIdx{group}),data(options.groupIdx{nextGroup}));
                plotSignificance(p,[group nextGroup],0.95,orientation='horizontal');
            end
        end
    end
end
h = gca; h.YAxis.Visible = 'off';


nexttile(children,2,[3 1]);
for group = 1:length(options.plotGroup)
    if options.plotGroup(group)
        groupIdx = options.groupIdx{group};
        groupColor = options.color{group};
        bootIdx = groupIdx(round((length(groupIdx)-1).*rand(options.nboot,1) + 1));
        histogram(data(bootIdx),FaceColor=groupColor,EdgeColor=groupColor); hold on;
    end
end
xlabel(options.xlabel); ylabel('Count'); box off;

end