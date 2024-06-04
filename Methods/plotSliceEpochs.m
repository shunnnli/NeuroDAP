function plotSliceEpochs(expPaths,options)

arguments
    expPaths % can be cell or string
    % If there's single session: plot everything
    % If there're multiple sessions: combine peaks/aucs and only plot EPSCvsIPSC

    options.groups double = 1; % indicate grouping of sessions if multiple sessions are provided
    options.conditions string = ["Control","Reward pairing","Punishment pairing"];
    options.conditionValues double = [0,1,2];
    options.conditionColors = {[0.6,0.6,0.6];[0.2941,0.3608,0.8];[1,0.1961,0.2275]};
    
    options.eventSample double = 10000 % in sample
    options.timeRange double = [-20,50] % in ms
    options.outputFs double = 10000

    options.plotAll logical = true;
    options.nArtifactSamples double = 0 % in sample
    options.yScalePadding double = 0.1; % in percentage
    options.VholdList double = [-70, 0, 8, 10, 12, 15, 18, 20];

    options.dotSizeScatter double = 100; % dot size for scatter plot
    options.dotSizeBox double = 80; % dot size for box plot
    options.fitScatter logical = true % fit scatter plot

    options.save logical = true; % save peaks and aucs
    options.resultPath string = 'PreQC';
end

%% Load epoch file

% Determine the number of input sessions
if isempty(expPaths); error("ERROR: empty expPaths variable!!");
elseif isstring(expPaths) || ischar(expPaths) || (iscell(expPaths) && isscalar(expPaths))
    expPath = expPaths; clearvars expPaths;
    combine = false;
else
    if length(options.groups) ~= length(expPaths)
        error("ERROR: options.groups has different length from expPaths!!");
    end
    combine = true; 
end

% Load epochs.mat file
if ~combine
    dirsplit = split(expPath,filesep); expName = dirsplit{end};
    if isempty(dir(fullfile(expPath,"epochs_*.mat")))
        error("ERROR: epochs_*.mat file NOT found!");
    else
        load(strcat(expPath,filesep,'epochs_',expName,'.mat'));
    end
    cellList = unique(epochs{:,"Cell"});
    % Whether cell have withVhold
    if contains('Vhold epoch mean',epochs.Properties.VariableNames); withVhold = true; 
    else; withVhold = false; end
else
    epochsList = cell(size(expPaths));
    for exp = 1:length(expPaths)
        expPath = expPaths{exp};
        dirsplit = split(expPath,filesep); expName = dirsplit{end};
        if isempty(dir(fullfile(expPath,"epochs_*.mat")))
            error("ERROR: epochs_*.mat file NOT found!");
        else
            load(strcat(expPath,filesep,'epochs_',expName,'.mat'));
            disp(strcat("Loaded: ", expName));
        end
        epochsList{exp} = epochs;
        withVhold = true;
    end
end

% Set up
timeRangeStartSample = options.eventSample + options.outputFs*options.timeRange(1)/1000;
timeRangeEndSample = options.eventSample + options.outputFs*options.timeRange(2)/1000;
plotWindow = timeRangeStartSample : timeRangeEndSample;
timeRangeInms = (plotWindow-1*options.outputFs) ./ (options.outputFs/1000);
analysisWindow = (options.eventSample+options.nArtifactSamples)-timeRangeStartSample : length(plotWindow);

[~,~,~,~,~,~,bluePurpleRed] = loadColors; 

% Update resultPath
options.resultPath = fullfile(epochs{1,'Session'},options.resultPath);

%% Plot summary raw trace for each cell/epoch

if ~combine
    if withVhold
        initializeFig(1,1); tiledlayout('flow');
        for c = 1:size(cellList,1)
            cellEpochs = epochs{epochs.Cell==cellList(c),["Raw sweeps","Included"]};
            vholdEpochs = epochs{epochs.Cell==cellList(c),"Vhold epoch mean"};
            haveIncluded = 0;
            for include = 1:size(cellEpochs,1)
                haveIncluded = haveIncluded + sum(cellEpochs{include,2});
            end
            if haveIncluded < 1; continue; end
    
            nexttile; yMax = nan; yMin = nan; legendList = [];
            nColors = round(linspace(1,size(bluePurpleRed,1),size(cellEpochs,1)));
            for e = 1:size(cellEpochs,1)
                if options.plotAll; included = 1:size(cellEpochs{1},1);
                else; included = cellEpochs{e,2}; end
                traces = cellEpochs{e,1}(included==1,plotWindow);
                if isempty(traces); continue; end
                plotSEM(timeRangeInms,traces,bluePurpleRed(nColors(e),:),...
                        meanOnly=true,plotIndividual=true);
                xlabel('Time (ms)');
                ylabel('Current (pA)');
                yMinTrace = min(traces(:,analysisWindow),[],"all");
                yMaxTrace = max(traces(:,analysisWindow),[],"all");
        
                if yMaxTrace > yMax || isnan(yMax); yMax = yMaxTrace; end
                if yMinTrace < yMin || isnan(yMin); yMin = yMinTrace; end
        
                vhold_cur = round(vholdEpochs(e));
                legendList = [legendList; strcat(num2str(vhold_cur),"mV (n=",num2str(size(cellEpochs{e,1},1)),")")];
            end
        
            yPad = abs(yMax-yMin)*options.yScalePadding;
            ylim([yMin-yPad,yMax+yPad]);
            legend(legendList);
            title(strcat('Cell #',num2str(cellList(c))));
        end
        saveFigures(gcf,'Trace_cell_raw',options.resultPath,saveFIG=true,savePDF=true);
    end
    
    initializeFig(1,1); tiledlayout('flow');
    for row = 1:size(epochs,1)
        nexttile;
        if options.plotAll; included = 1:size(epochs{row,'Raw sweeps'}{1},1);
        else; included = epochs{row,'Included'}{1}; end
        traces = epochs{row,'Raw sweeps'}{1}(included==1,plotWindow);
        plotSEM(timeRangeInms,traces,bluePurpleRed(1,:),...
                meanOnly=true,plotIndividual=true);
        xlabel('Time (ms)');
        ylabel('Current (pA)');
        yMin = min(traces(:,analysisWindow),[],"all");
        yMax = max(traces(:,analysisWindow),[],"all");
        yPad = abs(yMax-yMin)*options.yScalePadding;
        ylim([yMin-yPad,yMax+yPad]);
        title(strcat('Epochs #',num2str(epochs{row,'Epoch'})));
    end
    saveFigures(gcf,'Trace_epoch_raw',options.resultPath,saveFIG=true,savePDF=true);
end


%% Plot summary processed trace for each cell/epoch
if ~combine
    if withVhold
        initializeFig(1,1); tiledlayout('flow');
        for c = 1:size(cellList,1)
            cellEpochs = epochs{epochs.Cell==cellList(c),["Processed sweeps","Included"]};
            vholdEpochs = epochs{epochs.Cell==cellList(c),"Vhold epoch mean"};
            haveIncluded = 0;
            for include = 1:size(cellEpochs,1)
                haveIncluded = haveIncluded + sum(cellEpochs{include,2});
            end
            if haveIncluded < 1; continue; end
    
            nexttile; yMax = nan; yMin = nan; legendList = [];
            nColors = round(linspace(1,size(bluePurpleRed,1),size(cellEpochs,1)));
            for e = 1:size(cellEpochs,1)
                if options.plotAll; included = 1:size(cellEpochs{1},1);
                else; included = cellEpochs{e,2}; end
                traces = cellEpochs{e,1}(included==1,plotWindow);
                if isempty(traces); continue; end
                plotSEM(timeRangeInms,traces,bluePurpleRed(nColors(e),:),...
                        meanOnly=true,plotIndividual=true);
                xlabel('Time (ms)');
                ylabel('Current (pA)');
                yMinTrace = min(traces(:,analysisWindow),[],"all");
                yMaxTrace = max(traces(:,analysisWindow),[],"all");
        
                if yMaxTrace > yMax || isnan(yMax); yMax = yMaxTrace; end
                if yMinTrace < yMin || isnan(yMin); yMin = yMinTrace; end
        
                vhold_cur = round(vholdEpochs(e));
                legendList = [legendList; strcat(num2str(vhold_cur),"mV (n=",num2str(size(cellEpochs{e,1},1)),")")];
            end
        
            yPad = abs(yMax-yMin)*options.yScalePadding;
            ylim([yMin-yPad,yMax+yPad]);
            legend(legendList);
            title(strcat('Cell #',num2str(cellList(c))));
        end
        saveFigures(gcf,'Trace_cell_processed',options.resultPath,saveFIG=true,savePDF=true);
    end
    
    initializeFig(1,1); tiledlayout('flow');
    for row = 1:size(epochs,1)
        nexttile;
        if options.plotAll; included = 1:size(epochs{row,'Processed sweeps'}{1},1);
        else; included = epochs{row,'Included'}{1}; end
        traces = epochs{row,'Processed sweeps'}{1}(included==1,plotWindow);
        plotSEM(timeRangeInms,traces,bluePurpleRed(end,:),...
                meanOnly=true,plotIndividual=true);
        xlabel('Time (ms)');
        ylabel('Current (pA)');
        yMin = min(traces(:,analysisWindow),[],"all");
        yMax = max(traces(:,analysisWindow),[],"all");
        yPad = abs(yMax-yMin)*options.yScalePadding;
        ylim([yMin-yPad,yMax+yPad]);
        title(strcat('Epochs #',num2str(epochs{row,'Epoch'})));
    end
    saveFigures(gcf,'Trace_epoch_processed',options.resultPath,saveFIG=true,savePDF=true);
end

%% Extract EPSC/IPSC statistics

if ~withVhold; return; end

if ~combine
    peaks = nan(size(cellList,1),2);
    aucs = nan(size(cellList,1),2);
    groups = ones(size(cellList,1),2);
    
    for c = 1:size(cellList,1)
        % Find EPSC and IPSC stats
        % Intuition: if there's multiple epochs for similar range of Vhold, it
        % is because the previous one is slightly wrong (99% is because 0mV is 
        % not the true reversal potential)
        EPSCstats = epochs{epochs.Cell==cellList(c) & epochs.("Vhold epoch mean")<-50, ["Included","Peaks","AUCs"]};
        IPSCstats = epochs{epochs.Cell==cellList(c) & epochs.("Vhold epoch mean")>-10, ["Included","Peaks","AUCs"]};
    
        % Collect EPSC and IPSC stats
        if ~isempty(EPSCstats)
            EPSCstats = EPSCstats(end,:);
            peaks(c,1) = mean(EPSCstats{2}(EPSCstats{1}==1)); 
            aucs(c,1) = mean(EPSCstats{3}(EPSCstats{1}==1)); 
        end
        if ~isempty(IPSCstats)
            IPSCstats = IPSCstats(end,:);
            peaks(c,2) = mean(IPSCstats{2}(IPSCstats{1}==1)); 
            aucs(c,2) = mean(IPSCstats{3}(IPSCstats{1}==1));
        end
    end
    % Save peaks and aucs
    all_peaks = peaks; all_aucs = aucs;
    save(strcat(expPath,filesep,'epochs_',expName),'all_peaks','all_aucs','-append');
else
    all_peaks = []; all_aucs = []; groups = [];

    for exp = 1:length(expPaths)
        epochs = epochsList{exp};
        cellList = unique(epochs{:,"Cell"});
        peaks = nan(size(cellList,1),2);
        aucs = nan(size(cellList,1),2);
        group = ones(size(cellList,1),1) * options.groups(exp);
        
        for c = 1:size(cellList,1)
            % Find EPSC and IPSC stats
            % Intuition: if there's multiple epochs for similar range of Vhold, it
            % is because the previous one is slightly wrong (99% is because 0mV is 
            % not the true reversal potential)
            EPSCstats = epochs{epochs.Cell==cellList(c) & epochs.("Vhold epoch mean")<-50, ["Included","Peaks","AUCs"]};
            IPSCstats = epochs{epochs.Cell==cellList(c) & epochs.("Vhold epoch mean")>-10, ["Included","Peaks","AUCs"]};
        
            % Collect EPSC and IPSC stats
            if ~isempty(EPSCstats)
                EPSCstats = EPSCstats(end,:);
                peaks(c,1) = mean(EPSCstats{2}(EPSCstats{1}==1)); 
                aucs(c,1) = mean(EPSCstats{3}(EPSCstats{1}==1)); 
            end
            if ~isempty(IPSCstats)
                IPSCstats = IPSCstats(end,:);
                peaks(c,2) = mean(IPSCstats{2}(IPSCstats{1}==1)); 
                aucs(c,2) = mean(IPSCstats{3}(IPSCstats{1}==1));
            end
        end
        all_peaks = [all_peaks;peaks];
        all_aucs = [all_aucs;aucs];
        groups = [groups;group];
    end
end

% Remove cells with NaNs
peaks = all_peaks; aucs = all_aucs;
[peaks,rm_row] = rmmissing(peaks); aucs = rmmissing(aucs);
groups = groups(~rm_row);

% Calculate ratios
ratioPeaks = abs(peaks(:,1))./abs(peaks(:,2));
ratioAUCs = abs(aucs(:,1))./abs(aucs(:,2));


%% Plot EPSC/IPSC statistics

groupsList = unique(groups);
initializeFig(0.67,0.5); tiledlayout(1,6);

% Plot scatter plot for EPSC/IPSC peaks
nexttile([1 2]);
for i = 1:length(groupsList)
    groupIdx = (groups==groupsList(i));
    scatter(abs(peaks(groupIdx,1)),abs(peaks(groupIdx,2)),options.dotSizeScatter,options.conditionColors{i},'filled',...
        'MarkerFaceAlpha',1); hold on
    if options.fitScatter
        [p,S] = polyfit(abs(peaks(groupIdx,1)),abs(peaks(groupIdx,2)),1);
        [yfit,~] = polyval(p,abs(peaks(groupIdx,1)),S);
        plot(abs(peaks(groupIdx,1)),yfit,'LineWidth',2,'Color',options.conditionColors{i});
    end
    h = refline(1,0); h.LineWidth = 2; h.Color = [0.3,0.3,0.3]; h.LineStyle = '--';
end
xlabel("EPSC magnitude (pA)"); ylabel("IPSC magnitude (pA)"); 
% xlim([20,1400]); ylim([20,1400]);

% Plot bar plot for EPSC/IPSC peaks
nexttile([1 1]);
ylimMax = 20; 
significanceHeight1 = ylimMax*0.7; 
significanceHeight2 = ylimMax*0.8;
for i = 1:length(groupsList)
    groupIdx = (groups==groupsList(i));
    xgroupdata = ones(sum(groupIdx),1)*groupsList(i);
    boxchart(xgroupdata,ratioPeaks(groupIdx,:),'BoxFaceColor',options.conditionColors{i}); hold on
    swarmchart(xgroupdata,ratioPeaks(groupIdx,:),...
        options.dotSizeBox,options.conditionColors{i},'filled',...
        'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.5);
end
xticks(groupsList); xticklabels(options.conditions(groupsList+1));
ylabel('EPSC/IPSC magnitude'); ylim([0,20]);
% Add significance
if length(groupsList) == 2
    group1Idx = groups==groupsList(1);
    group2Idx = groups==groupsList(2);
    [~,p_peak,~] = kstest2(ratioPeaks(group1Idx,:),ratioPeaks(group2Idx,:));
    text(mean(groupsList)-0.1,15,strcat("p=",num2str(round(p_peak,4))),'FontSize',12);
elseif length(groupsList) == 3
    group1Idx = groups==groupsList(1);
    group2Idx = groups==groupsList(2);
    group3Idx = groups==groupsList(3);
    [~,p_peak1,~] = kstest2(ratioPeaks(group1Idx,:),ratioPeaks(group2Idx,:));
    [~,p_peak2,~] = kstest2(ratioPeaks(group2Idx,:),ratioPeaks(group3Idx,:));
    [~,p_peak3,~] = kstest2(ratioPeaks(group1Idx,:),ratioPeaks(group3Idx,:));

    plotSignificance(p_peak1,groupsList([1 2]),significanceHeight1);
    plotSignificance(p_peak2,groupsList([2 3]),significanceHeight1);
    plotSignificance(p_peak3,groupsList([1 3]),significanceHeight2);
end

% Plot scatter plot for EPSC/IPSC AUCs
nexttile([1 2]);
for i = 1:length(groupsList)
    groupIdx = (groups==groupsList(i));
    scatter(abs(aucs(groupIdx,1)),abs(aucs(groupIdx,2)),options.dotSizeScatter,options.conditionColors{i},'filled',...
        'MarkerFaceAlpha',1); hold on
    if options.fitScatter
        [p,S] = polyfit(abs(aucs(groupIdx,1)),abs(aucs(groupIdx,2)),1);
        [yfit,~] = polyval(p,abs(aucs(groupIdx,1)),S);
        plot(abs(aucs(groupIdx,1)),yfit,'LineWidth',2,'Color',options.conditionColors{i});  
    end
    h = refline(1,0); h.LineWidth = 2; h.Color = [0.3,0.3,0.3]; h.LineStyle = '--';
end
xlabel("EPSC AUC"); ylabel("IPSC AUC"); 
% xlim([20,1400]); ylim([20,1400]);

% Plot bar plot for EPSC/IPSC AUCs
nexttile([1 1]);
ylimMax = 10; 
significanceHeight1 = ylimMax*0.7; 
significanceHeight2 = ylimMax*0.8;
for i = 1:length(groupsList)
    groupIdx = (groups==groupsList(i));
    xgroupdata = ones(sum(groupIdx),1)*groupsList(i);
    boxchart(xgroupdata,ratioAUCs(groupIdx,:),'BoxFaceColor',options.conditionColors{i}); hold on
    swarmchart(xgroupdata,ratioAUCs(groupIdx,:),...
        options.dotSizeBox,options.conditionColors{i},'filled',...
        'MarkerFaceAlpha',0.8,'XJitter','density','XJitterWidth',0.5);
end
if length(groupsList) == 2
    group1Idx = groups==groupsList(1);
    group2Idx = groups==groupsList(2);
    [~,p_auc,~] = kstest2(ratioAUCs(group1Idx,:),ratioAUCs(group2Idx,:));
    text(mean(groupsList)-0.1,15,strcat("p=",num2str(round(p_auc,4))),'FontSize',12);
elseif length(groupsList) == 3
    group1Idx = groups==groupsList(1);
    group2Idx = groups==groupsList(2);
    group3Idx = groups==groupsList(3);
    [~,p_auc1,~] = kstest2(ratioAUCs(group1Idx,:),ratioAUCs(group2Idx,:));
    [~,p_auc2,~] = kstest2(ratioAUCs(group2Idx,:),ratioAUCs(group3Idx,:));
    [~,p_auc3,~] = kstest2(ratioAUCs(group1Idx,:),ratioAUCs(group3Idx,:));
   
    plotSignificance(p_auc1,groupsList([1 2]),significanceHeight1);
    plotSignificance(p_auc2,groupsList([2 3]),significanceHeight1);
    plotSignificance(p_auc3,groupsList([1 3]),significanceHeight2);
end
xticks(groupsList); xticklabels(options.conditions(groupsList+1));
ylabel('EPSC/IPSC AUC'); ylim([0,10]);

% Save
if ~combine; saveFigures(gcf,'Summary_EPSCvsIPSC',options.resultPath,saveFIG=true,savePDF=true);
else
    saveFigures(gcf,'Summary_EPSCvsIPSC',...
        strcat(osPathSwitch(options.resultPath),expName(1:6)),...
        saveFIG=true,savePDF=true); 
    save(strcat(osPathSwitch(options.resultPath),expName(1:6),filesep,'Summary_EPSCvsIPSC'),...
        'ratioPeaks','ratioAUCs','expPaths','groups','-v7.3');
end

end