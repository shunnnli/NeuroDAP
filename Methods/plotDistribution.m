function plotDistribution(data,options)

arguments
    data double

    options.plotGroup logical = false % vector, plot this group or not
    options.groupIdx cell

    % Provide bootstrap simulation data
    options.bootstrap logical
    options.bootData double
    options.bootIdx cell

    options.color
    options.opacity double = 0.3

    options.nboot double = 5000
    options.nbins double = 50;

    options.masterlayout
    options.tile double
    options.tileSpan double = [1 1];
    options.dotSize double = 100
    options.LineWidth double = 3
    options.xlabel string

    options.limit double
end

%% Check input
if ~isfield(options,'masterlayout')
    master = tiledlayout(1,1);
else; master = options.masterlayout; 
end

if length(options.plotGroup) ~= size(options.groupIdx)
    error('plotGroup and groupIdx have different size!');
end


if ~isfield(options,'bootstrap')
    if xor(~isfield(options,'bootData'), ~isfield(options,'bootIdx'))
        error('bootIdx or bootData not provided. Need to provide both!');
    elseif isfield(options,'bootIdx') && isfield(options,'bootData')
        options.bootstrap = true;
    else
        options.bootstrap = false;
    end
elseif options.bootstrap
    if ~isfield(options,'bootData') || ~isfield(options,'bootIdx')
        error('bootIdx or bootData not provided. Need to provide both!');
    end
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
        plotScatterBar(group,data(groupIdx),color=groupColor,dotSize=options.dotSize,orientation='horizontal'); hold on;
    end
end

if isfield(options,'limit'); xlim(options.limit); end

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
if isfield(options,'ylim'); ylim(options.ylim); end


nexttile(children,2,[3 1]);
for group = 1:length(options.plotGroup)
    if options.plotGroup(group)
        % Plot bootstrap distribution if provided
        if options.bootstrap
            yyaxis right
            bsIdx = options.bootIdx{group};
            bsColor = 1 - options.opacity*(1-options.color{group});

            % Get average data for each sim (i.e. column)
            sim_avg = mean(options.bootData(bsIdx,:),1);
            histogram(sim_avg,options.nbins,FaceColor=bsColor,EdgeColor=bsColor); hold on;
            observed_avg = mean(data(groupIdx),'all');

            % p value
            left_pval = sum(sim_avg<=observed_avg)/length(sim_avg);
            right_pval = sum(sim_avg>=observed_avg)/length(sim_avg);
            pval = min([left_pval,right_pval]);
            xline(observed_avg,label=['p=',num2str(pval)],...
                    Color=groupColor,LineWidth=options.LineWidth,...
                    LabelOrientation='horizontal');
            xline(mean(data(groupIdx),'all'),Color=groupColor,LineWidth=options.LineWidth);
        else
            % Plot current histogram distribution
            yyaxis left
            groupIdx = options.groupIdx{group};
            groupColor = options.color{group};
            bootIdx = groupIdx(round((size(groupIdx,1)-1).*rand(options.nboot,1) + 1));
            histogram(data(bootIdx),FaceColor=groupColor,EdgeColor=groupColor); hold on;
        end
    end
end

xlabel(options.xlabel);
c = gca; c.YAxis(1).Visible = 'off'; c.YAxis(2).Visible = 'off';
if isfield(options,'limit'); xlim(options.limit); end
box off;

end