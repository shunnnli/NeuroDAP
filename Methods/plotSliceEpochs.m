function plotSliceEpochs(epochs,options)

arguments
    epochs

    options.plotEpoch logical = true % plot epoch averaged traces
    options.plotCell logical = true % plot cell averaged traces

    options.colormap double
    
    options.eventSample double = 10000 % in sample
    options.timeRange double = [-20,50] % in ms
    options.outputFs double = 10000

    options.plotAll logical = true % plot all sweeps if true
    options.nArtifactSamples double = 0 % in sample

    options.save logical = true; % save peaks and aucs
    options.resultPath;
end

%% Set up
timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = (plotWindow-1*options.outputFs) ./ (options.outputFs/1000);
analysisWindow = (options.eventSample+options.nArtifactSamples)-timeRangeStartSample : length(plotWindow);

% Fill in resultPath if neccessary
if isempty(options.resultPath)
    options.resultPath = fullfile(epochs{1,'Session'},options.resultPath);
    disp(strcat('Updated resultPath: ', options.resultPath));
end

%% Plot summary raw trace for each cell/epoch

if options.plotCell

    cellList = unique(epochs.Cell);

    initializeFig(1,1); tiledlayout('flow');
    for c = 1:size(cellList,1)
        cellEpochs = epochs{epochs.Cell==cellList(c),["Raw sweeps","Included"]};
        vholdEpochs = epochs{epochs.Cell==cellList(c),"Vhold"};
        haveIncluded = 0;
        for include = 1:size(cellEpochs,1)
            haveIncluded = haveIncluded + sum(cellEpochs{include,2});
        end
        if haveIncluded < 1; continue; end

        nexttile; yMax = nan; yMin = nan; legendList = [];
        nColors = round(linspace(1,size(options.colormap,1),size(cellEpochs,1)));
        for e = 1:size(cellEpochs,1)
            if options.plotAll; included = ones(size(cellEpochs{e,2},1),1);
            else; included = cellEpochs{e,2}; end
            traces = cellEpochs{e,1}(included==1,plotWindow);
            if isempty(traces); continue; end
            plotSEM(timeRangeInms,traces,options.colormap(nColors(e),:),...
                    plotPatch=false,plotIndividual=true);
            xlabel('Time (ms)');
            ylabel('Current (pA)');
            yMinTrace = min(traces(:,analysisWindow),[],"all");
            yMaxTrace = max(traces(:,analysisWindow),[],"all");
    
            if yMaxTrace > yMax || isnan(yMax); yMax = yMaxTrace; end
            if yMinTrace < yMin || isnan(yMin); yMin = yMinTrace; end
    
            vhold_cur = round(vholdEpochs(e));
            legendList = [legendList; strcat(num2str(vhold_cur),"mV (n=",num2str(size(cellEpochs{e,1},1)),")")];
        end
        ylim([yMin-eps,yMax+eps]);
        legend(legendList);
        title(strcat('Cell #',num2str(cellList(c))));
    end
    saveFigures(gcf,'Trace_cell_raw',options.resultPath,saveFIG=true,savePDF=true);
end

if options.plotEpoch
    initializeFig(1,1); tiledlayout('flow');
    for row = 1:size(epochs,1)
        nexttile;
        if options.plotAll; included = ones(size(epochs{row,'Raw sweeps'}{1},1),1);
        else; included = epochs{row,'Included'}{1}; end
        traces = epochs{row,'Raw sweeps'}{1}(included==1,plotWindow);
        if isempty(traces); continue; end
        plotSEM(timeRangeInms,traces,options.colormap(1,:),...
                plotPatch=false,plotIndividual=true);
        xlabel('Time (ms)');
        ylabel('Current (pA)');
        yMin = min(traces(:,analysisWindow),[],"all");
        yMax = max(traces(:,analysisWindow),[],"all");
        ylim([yMin-eps,yMax+eps]);
        title(strcat('Epochs #',num2str(epochs{row,'Epoch'})));
    end
    saveFigures(gcf,'Trace_epoch_raw',options.resultPath,saveFIG=true,savePDF=true);
end


%% Plot summary processed trace for each cell/epoch

if options.plotCell
    initializeFig(1,1); tiledlayout('flow');
    for c = 1:size(cellList,1)
        cellEpochs = epochs{epochs.Cell==cellList(c),["Processed sweeps","Included"]};
        vholdEpochs = epochs{epochs.Cell==cellList(c),"Vhold"};
        haveIncluded = 0;
        for include = 1:size(cellEpochs,1)
            haveIncluded = haveIncluded + sum(cellEpochs{include,2});
        end
        if haveIncluded < 1; continue; end

        nexttile; yMax = nan; yMin = nan; legendList = [];
        nColors = round(linspace(1,size(options.colormap,1),size(cellEpochs,1)));
        for e = 1:size(cellEpochs,1)
            if options.plotAll; included = ones(size(cellEpochs{e,2},1),1);
            else; included = cellEpochs{e,2}; end
            traces = cellEpochs{e,1}(included==1,plotWindow);
            if isempty(traces); continue; end
            plotSEM(timeRangeInms,traces,options.colormap(nColors(e),:),...
                    plotPatch=false,plotIndividual=true);
            xlabel('Time (ms)');
            ylabel('Current (pA)');
            yMinTrace = min(traces(:,analysisWindow),[],"all");
            yMaxTrace = max(traces(:,analysisWindow),[],"all");
    
            if yMaxTrace > yMax || isnan(yMax); yMax = yMaxTrace; end
            if yMinTrace < yMin || isnan(yMin); yMin = yMinTrace; end
    
            vhold_cur = round(vholdEpochs(e));
            legendList = [legendList; strcat(num2str(vhold_cur),"mV (n=",num2str(size(cellEpochs{e,1},1)),")")];
        end
        ylim([yMin-eps,yMax+eps]);
        legend(legendList);
        title(strcat('Cell #',num2str(cellList(c))));
    end
    saveFigures(gcf,'Trace_cell_processed',options.resultPath,saveFIG=true,savePDF=true);
end

if options.plotEpoch
    initializeFig(1,1); tiledlayout('flow');
    for row = 1:size(epochs,1)
        nexttile;
        if options.plotAll; included = ones(size(epochs{row,'Processed sweeps'}{1},1),1);
        else; included = epochs{row,'Included'}{1}; end
        traces = epochs{row,'Processed sweeps'}{1}(included==1,plotWindow);
        if isempty(traces); continue; end
        plotSEM(timeRangeInms,traces,options.colormap(end,:),...
                plotPatch=false,plotIndividual=true);
        xlabel('Time (ms)');
        ylabel('Current (pA)');
        yMin = min(traces(:,analysisWindow),[],"all");
        yMax = max(traces(:,analysisWindow),[],"all");
        ylim([yMin-eps,yMax+eps]);
        title(strcat('Epochs #',num2str(epochs{row,'Epoch'})));
    end
    saveFigures(gcf,'Trace_epoch_processed',options.resultPath,saveFIG=true,savePDF=true);
end

end