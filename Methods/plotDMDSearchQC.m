function fig = plotDMDSearchQC(curCell, searchIdx, options)

% Plot per-sweep DMD search QC metrics from the compact cells_DMD QC cache.

arguments
    curCell table
    searchIdx double

    options.metrics string = ["Rs","Cm","Rm","Verror","Ibaseline","Ibaseline_std"]
    options.save logical = true
    options.savePNG logical = false
    options.savePDF logical = true
    options.saveFIG logical = false
    options.saveDataPath string = 'default'
    options.QCThreshold struct = struct()
end

fig = [];
curSearch = curCell.Epochs{1}{searchIdx};

if strcmp(options.saveDataPath,'default')
    options.saveDataPath = string(curCell.Session);
end

qcTable = localGetSearchQCTable(curCell, searchIdx);
if isempty(qcTable)
    warning('plotDMDSearchQC:MissingQC', ...
            'No compact QC summary found for %s. Re-run loadSlicesDMD to add QC to cells_DMD.', char(curSearch));
    return
end

qcTable.SweepIndex = (1:height(qcTable))';

if isfield(curCell.Options{1},'QCThreshold') && isempty(fieldnames(options.QCThreshold))
    options.QCThreshold = curCell.Options{1}.QCThreshold;
end

metrics = options.metrics(ismember(options.metrics, string(qcTable.Properties.VariableNames)));
if isempty(metrics)
    warning('plotDMDSearchQC:NoMetrics', ...
            'None of the requested QC metrics were found for %s.', char(curSearch));
    return
end

initializeFig(0.85,0.8);
fig = gcf;
tl = tiledlayout(2,ceil(numel(metrics)/2));
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

depths = unique(qcTable.Depth(~isnan(qcTable.Depth)));
depthColors = lines(max(1,numel(depths)));

for m = 1:numel(metrics)
    metric = char(metrics(m));
    nexttile; hold on; box off; grid on;

    for d = 1:numel(depths)
        depthRows = qcTable.Depth == depths(d);
        y = qcTable{depthRows,metric};
        x = qcTable.SweepIndex(depthRows);
        plot(x, y, '-o', ...
             'Color', depthColors(d,:), ...
             'MarkerFaceColor', depthColors(d,:), ...
             'MarkerSize', 4, ...
             'LineWidth', 1.2, ...
             'DisplayName', ['Depth ',num2str(depths(d))]);
    end

    threshold = localGetThreshold(options.QCThreshold, metric);
    if ~isnan(threshold)
        yline(threshold, '--k', 'LineWidth', 1, 'DisplayName', 'QC threshold');
        if strcmpi(metric,'Verror')
            yline(-threshold, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end

    ylabel(localMetricLabel(metric));
    xlabel('Sweep');
    title(localMetricTitle(metric));
    xlim([0.5, height(qcTable)+0.5]);
end

legendTile = nexttile(tl,1);
legend(legendTile, 'Location', 'best');
sgtitle(sprintf('Cell %d %s QC', curCell.Cell, char(curSearch)), 'Interpreter', 'none');

if options.save
    filepath = fullfile(options.saveDataPath, ['cell',num2str(curCell.Cell)], 'Search summary');
    saveFigures(fig, [char(curSearch),'_QC'], filepath, ...
                savePNG=options.savePNG, savePDF=options.savePDF, saveFIG=options.saveFIG);
end

end

function qcTable = localGetSearchQCTable(curCell, searchIdx)
    qcTable = table();
    if ~ismember('QC', curCell.Properties.VariableNames) || isempty(curCell.QC{1})
        return
    end

    cellQC = curCell.QC{1};
    if numel(cellQC) < searchIdx || isempty(cellQC{searchIdx})
        return
    end

    searchQC = cellQC{searchIdx};
    if istable(searchQC)
        qcTable = searchQC;
    elseif iscell(searchQC)
        nonEmpty = cellfun(@(x) istable(x) && height(x) > 0, searchQC);
        if any(nonEmpty)
            qcTable = vertcat(searchQC{nonEmpty});
        end
    end
end

function threshold = localGetThreshold(QCThreshold, metric)
    threshold = nan;
    if isempty(QCThreshold) || ~isstruct(QCThreshold) || ~isfield(QCThreshold, metric)
        return
    end

    value = QCThreshold.(metric);
    if isnumeric(value) && isscalar(value)
        threshold = value;
    end
end

function label = localMetricLabel(metric)
    switch metric
        case {'Rs','Rm'}
            label = [metric, ' (MOhm)'];
        case 'Cm'
            label = 'Cm (pF)';
        case 'tau'
            label = 'Tau (ms)';
        case 'Verror'
            label = 'Verror (mV)';
        case {'Ibaseline','Ibaseline_std'}
            label = [metric, ' (pA)'];
        otherwise
            label = metric;
    end
end

function titleText = localMetricTitle(metric)
    switch metric
        case 'Rs'
            titleText = 'Series resistance';
        case 'Rm'
            titleText = 'Membrane resistance';
        case 'Cm'
            titleText = 'Membrane capacitance';
        case 'tau'
            titleText = 'RC tau';
        case 'Verror'
            titleText = 'Voltage error';
        case 'Ibaseline'
            titleText = 'Baseline current';
        case 'Ibaseline_std'
            titleText = 'Baseline noise';
        otherwise
            titleText = metric;
    end
end
