
function analyzeSlice_OptoPair(expPaths,options)

arguments
    expPaths % can be cell or string
    % If there's single session: plot everything
    % If there're multiple sessions: combine peaks/aucs and only plot EPSCvsIPSC

    options.plotTraces logical = true % just plot traces for epoch/cell
    options.plotStats logical = false % calculate EPSC/IPSC

    options.plotEpoch logical = true
    options.plotCell logical = true

    options.groups double = 1; % indicate grouping of sessions if multiple sessions are provided
    options.conditions string = ["Control","Reward pairing","Punishment pairing"];
    options.conditionValues double = [0,1,2];
    options.conditionColors = {[0.6,0.6,0.6];[0.2941,0.3608,0.8];[1,0.1961,0.2275]};
    
    options.eventSample double = 10000 % in sample
    options.timeRange double = [-20,50] % in ms
    options.outputFs double = 10000

    options.plotAll logical = true
    options.nArtifactSamples double = 0 % in sample
    options.VholdList double = [-70, 0, 8, 10, 12, 15, 18, 20];

    options.dotSizeScatter double = 100; % dot size for scatter plot
    options.dotSizeBox double = 80; % dot size for box plot
    options.fitScatter logical = true % fit scatter plot

    options.save logical = true; % save peaks and aucs
    options.resultPath
end

%% Load epoch file

% Determine the number of input sessions
if isempty(expPaths); error("ERROR: empty expPaths variable!!");
elseif isstring(expPaths) || ischar(expPaths) || (iscell(expPaths) && isscalar(expPaths)) 
    % If only one session is provided
    singleSession = true;
else
    if length(options.groups) ~= length(expPaths)
        error("ERROR: options.groups has different length from expPaths!!");
    end
    singleSession = false; 
end

% Load epochs.mat file
if singleSession
    dirsplit = split(expPaths,filesep); expName = dirsplit{end};
    [allEpochs,allCells] = loadSlices(expPaths,reload=false);
    disp(strcat('Loaded session: ',expName));
else
    allCells = [];
    for exp = 1:length(expPaths)
        dirsplit = split(expPaths{exp},filesep); expName = dirsplit{end};
        [~,sessionCell] = loadSlices(expPaths{exp},reload=false);
        disp(strcat('Loaded session: ',expName, ...
                ' (',num2str(exp),'/',num2str(length(expPaths)),')'));
        allCells = [allCells; sessionCell];
    end
end

% Set up
[~,~,~,~,~,~,bluePurpleRed] = loadColors; 
% Set up resultPath
if ~isfield(options,'resultPath')
    if singleSession; options.resultPath = fullfile(allEpochs{1,'Session'});
    else
        error('Multiple sessions: please specify resultPath!'); 
    end
else
    if singleSession; options.resultPath = fullfile(allEpochs{1,'Session'},options.resultPath);
    elseif ~contains(options.resultPath,filesep)
        error('Multiple sessions: please specify resultPath!');
    end
end

%% Plot response traces for epoch/cell

if singleSession && options.plotTraces
    plotSliceEpochs(allEpochs,colormap=bluePurpleRed,...
            plotCell=options.plotCell,plotEpoch=options.plotEpoch,plotAll=options.plotAll,...
            resultPath=options.resultPath,save=options.save,...
            timeRange=options.timeRange,...
            eventSample=options.eventSample,nArtifactSamples=options.nArtifactSamples,...
            outputFs=options.outputFs);
end

%% ************ UNFINISHED ************** Plot stats

if options.plotStats
    animalList = unique(allCells{:,"Animal"});

    % Initialize stats result matrix
    all_peaks = []; all_aucs = []; groups = [];

    % Loop through animals to get results
    for animal = 1:length(animalList)
        % Find corresponding physiology session
        animalCells = allCells(strcmp(allCells.Animal,animalList(animal)),:);
        cellList = unique(animalCells.Cell);
        peaks = nan(size(cellList,1),2);
        aucs = nan(size(cellList,1),2);
        group = ones(size(cellList,1),1) * options.groups(animal);
        
        for c = 1:size(cellList,1)
            % Find EPSC and IPSC stats
            stats = animalCells{animalCells.Cell==cellList(c),'Stats'}{1};
            peaks(c,1) = stats.EPSC.peakAvg;
            peaks(c,2) = stats.IPSC.peakAvg;
            aucs(c,1) = stats.EPSC.aucAvg;
            aucs(c,2) = stats.IPSC.aucAvg;
        end

        % Save results
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
if singleSession; saveFigures(gcf,'Summary_EPSCvsIPSC',options.resultPath,saveFIG=true,savePDF=true);
else
    saveFigures(gcf,'Summary_EPSCvsIPSC',...
        strcat(osPathSwitch(options.resultPath),expName(1:6)),...
        saveFIG=true,savePDF=true); 
    save(strcat(osPathSwitch(options.resultPath),expName(1:6),filesep,'Summary_EPSCvsIPSC'),...
        'ratioPeaks','ratioAUCs','expPaths','groups','-v7.3');
end

end