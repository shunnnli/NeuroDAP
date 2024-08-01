function analyzeDMDSearch(curCell, searchIdx, options)

% Create spots summary file and plot summary figure

arguments
    curCell table
    searchIdx double

    options.redStim logical = true
    options.depthLineWidth double = [3,2.5,2,1.5,1.1,1,0.5];

    options.save logical = true
    options.savePNG logical = true
    options.savePDF logical = true
    options.saveFIG logical = true

    options.outputFs double = 10000
    options.timeRange double = [-1,50] % in ms
    options.analysisWindowLength double = 50 % in ms after stim onset
    options.controlWindowLength double = 50 % in ms before stim onset
    options.eventSample % in sample
    options.nArtifactSamples double = 0 % in sample
    options.rcCheckRecoveryWindow double = 100 % in ms
    options.peakWindow double = 2 % in ms around the peak to average

    options.feature = 'auc'
    options.thresholdFactor double = 3 % 3*std

    options.resultsPathDate char = 'newest'
    options.saveDataPath string = 'default'
end

%% Setup

[~,~,~,~,blueWhiteRed,~,~] = loadColors;

% Define path
strsplit = split(curCell.Session,'Results');
expPath = strsplit{1};

% Define results path
if strcmp(options.saveDataPath,'default')
    resultsList = sortrows(struct2cell(dir(fullfile(expPath,'Results_*')))',3);
    if strcmp(options.resultsPathDate,'newest')
        resultsPath = fullfile(resultsList{end,2},resultsList{end,1});
    else
        [~,dateIdx] = find(cellfun(@(x) contains(x, options.resultsPathDate), resultsList(:,1)));
        resultsPath = fullfile(resultsList{dateIdx,2},resultsList{dateIdx,1});
    end
    if strcmp(options.saveDataPath,'default')
        options.saveDataPath = resultsPath;
    end
end

% Define colors
if options.redStim; stimColor = blueWhiteRed(end,:);
else; stimColor = blueWhiteRed(1,:); end
blue = blueWhiteRed(1,:);
red = blueWhiteRed(end,:);
purple = [232, 130, 250]./255;

% Define time windows
if options.outputFs ~= curCell.Options{1}.outputFs
    warning('analyzeDMDSearch: Default options.outputFs differs from outputFs extracted from cells_DMD. Using cells_DMD value instead!');
    options.outputFs = curCell.Options{1}.outputFs;
end
if options.timeRange ~= curCell.Options{1}.timeRange
    eventSample = abs(curCell.Options{1}.timeRange(1))*(curCell.Options{1}.outputFs/1000) + 1;
    plotFirstSample = eventSample + options.timeRange(1)*(curCell.Options{1}.outputFs/1000);
    plotLastSample = eventSample + options.timeRange(2)*(curCell.Options{1}.outputFs/1000);
    plotWindowLength = plotLastSample - plotFirstSample + 1;
    plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
    analysisWindow = eventSample : eventSample+options.analysisWindowLength*10;
else
    if isfield(curCell.Options{1},'plotWindowTime')
        plotWindowTime = curCell.Options{1}.plotWindowTime;
        plotWindowLength = length(plotWindowTime);
    else
        if isfield(curCell.Options{1},'plotWindowLength'); plotWindowLength = curCell.Options{1}.plotWindowLength;
        else; plotWindowLength = (options.timeRange(2)-options.timeRange(1))* options.outputFs/1000 + 1;
        end
        plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
    end
    eventSample = find(plotWindowTime==0);
    plotFirstSample = 1; plotLastSample = plotWindowLength;
    analysisWindow = eventSample : eventSample+options.analysisWindowLength*10;
end

%% Load cell info

% Load cell info
curSearch = curCell.Epochs{1}{searchIdx};
search_depths = curCell.("Response map"){1}.depths{searchIdx};
search_rmap = curCell.("Response map"){1}.responseMap{searchIdx};
% search_isResponse = curCell.("Response map"){1}.isResponseMap{searchIdx};
% search_hotspotMap = curCell.("Response map"){1}.hotspotMap{searchIdx};
search_cmap = curCell.("Response map"){1}.currentMap{searchIdx};
search_bmap = curCell.("Response map"){1}.baselineMap{searchIdx};
search_vhold = curCell.Vhold{1}(searchIdx);
search_hotspot = curCell.("Response map"){1}.hotspot{searchIdx};
nDepth = length(search_depths);

% Load cell location (not necessary)
try cellLoc = curCell.Options{1}.cellLocation;
catch
    cellLoc = [nan nan]; 
end

% Load stim duration
if strcmp('Protocol',curCell.Properties.VariableNames)
    stimDuration = curCell.("Protocol"){1}(searchIdx).pulseWidth;
else
    stimDuration = 5;
end

% Define color
if search_vhold < -50
    color = blueWhiteRed(end,:); 
elseif search_vhold > -10
    color = blueWhiteRed(1,:); 
else
    color = purple; 
end

% Load noise model of this neuron
load(fullfile(expPath,['noise_cell',num2str(curCell.Cell),'.mat']),'allNullData');
Ethres = curCell.Stats{1}.Ethres; %-2 * noise_all.sigma;
Ithres = curCell.Stats{1}.Ithres; %2 * noise_all.sigma;

% Load spots.mat for this neuron
% cellResultsPath = fullfile(curCell.Session,['cell',num2str(curCell.Cell)]);
% depthfilename = strcat(curSearch,'_depth*.mat');
% spotsList = dir(fullfile(cellResultsPath,depthfilename));

%% Generate a summary figure for each search

for d = 1:nDepth
    %% Initialization
    curDepth = search_depths(d);
    disp(['Ongoing: calculating statistics for search: ', curSearch,' at depth ',num2str(curDepth)]);

    % Load spots.mat file of this depth
    % depthFilePath = fullfile(spotsList(d).folder,spotsList(d).name);
    % load(depthFilePath,'spotsAtDepth');

    % Load depth specific data
    depthResponseMap = search_rmap(:,:,d);
    % depthHotspotMap = search_isResponse(:,:,d);
    % depthHotspotMap = search_hotspotMap(:,:,d);
    depthCurrentMap = search_cmap{d};
    depthBaselineMap = search_bmap{d};
    spotSequence = repelem(1:height(search_hotspot{d}), cellfun(@numel, search_hotspot{d}))';
    depthHotspot = search_hotspot{d};

    % Determine current plot line width
    if curDepth <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(curDepth);
    else; LineWidth = options.depthLineWidth(end); end

    % Get opto and prestim traces
    optoData = cell2mat(depthCurrentMap); 
    optoData = optoData(:,analysisWindow);
    ctrlData = cell2mat(depthBaselineMap);
    ctrlData = ctrlData(:,size(ctrlData,2)-length(analysisWindow)+1:size(ctrlData,2));
    optoTime = linspace(0,options.analysisWindowLength,size(optoData,2));
    ctrlTime = linspace(-options.controlWindowLength,0,size(ctrlData,2));

    % Calculate max/min current for plotting
    % optoMax = optoData(~isoutlier(optoData)); 
    % optoMin = optoData(~isoutlier(optoData)); 
    optoMax = max(optoData,[],'all'); optoMin = min(optoData,[],'all');

    % Get hotspot vs null spot
    hotspot_spotIdx = cell2mat(cellfun(@(x) sum(x)>=1, search_hotspot{d}, UniformOutput=false));
    hotspot_sweepIdx = cell2mat(cellfun(@(x) repmat(sum(x)>=1,size(x)), search_hotspot{d}, UniformOutput=false));
    hotspotData = optoData(hotspot_sweepIdx,:);
    nullspotData = optoData(~hotspot_sweepIdx,:);

    % Get response statistics for each spot
    spotMax_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.max{searchIdx}{d},UniformOutput=false));
    spotMin_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.min{searchIdx}{d},UniformOutput=false));
    spotMaxTime_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.maxTime{searchIdx}{d},UniformOutput=false));
    spotMinTime_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.minTime{searchIdx}{d},UniformOutput=false));
    spotAUC_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.auc{searchIdx}{d},UniformOutput=false));
    ctrlAUC_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.baseline.auc{searchIdx}{d},UniformOutput=false));
    % ctrlSTD_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.baseline.std{searchIdx}{d},UniformOutput=false));
    ctrlMax_avg = cell2mat(cellfun(@(x) max(x,[],'all'), depthBaselineMap,UniformOutput=false));
    ctrlMin_avg = cell2mat(cellfun(@(x) min(x,[],'all'), depthBaselineMap,UniformOutput=false));
    [~,ctrlMaxTimeIdx_avg] = cellfun(@(x) max(x,[],'all'), depthBaselineMap,UniformOutput=false);
    ctrlMaxTime_avg = cell2mat(ctrlMaxTimeIdx_avg)/options.outputFs*1000;
    [~,ctrlMinTimeIdx_avg] = cellfun(@(x) max(x,[],'all'), depthBaselineMap,UniformOutput=false);
    ctrlMinTime_avg = cell2mat(ctrlMinTimeIdx_avg)/options.outputFs*1000;
    spotAbsAUC = sum(abs(optoData),2)/options.outputFs; 
    ctrlAbsAUC = sum(abs(ctrlData),2)/options.outputFs;
    spotAbsAUC_avg = accumarray(spotSequence(:), spotAbsAUC(:), [], @(x) {x});
    ctrlAbsAUC_avg = accumarray(spotSequence(:), ctrlAbsAUC(:), [], @(x) {x});
    spotAbsAUC_avg = cell2mat(cellfun(@(x) mean(x,"all"), spotAbsAUC_avg,UniformOutput=false));
    ctrlAbsAUC_avg = cell2mat(cellfun(@(x) mean(x,"all"), ctrlAbsAUC_avg,UniformOutput=false));

    % Calculate response rate
    % Prestim success rate
    Iresponse_ctrl = cell2mat(arrayfun(@(x) length(findpeaks(ctrlData(x,:),MinPeakDistance=20,MinPeakProminence=Ithres,MinPeakWidth=options.peakWindow*10)),...
                    1:size(ctrlData,1),UniformOutput=false)');
    Eresponse_ctrl = cell2mat(arrayfun(@(x) length(findpeaks(-ctrlData(x,:),MinPeakDistance=20,MinPeakProminence=-Ethres,MinPeakWidth=options.peakWindow*10)),...
                    1:size(ctrlData,1),UniformOutput=false)');
    ISpotResponse_ctrl = accumarray(spotSequence(:), Iresponse_ctrl(:), [], @(x) {x});
    ESpotResponse_ctrl = accumarray(spotSequence(:), Eresponse_ctrl(:), [], @(x) {x});
    Irate_ctrl = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ISpotResponse_ctrl,UniformOutput=false));
    Erate_ctrl = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ESpotResponse_ctrl,UniformOutput=false));
    % Opto success rate
    Iresponse_opto = cell2mat(arrayfun(@(x) length(findpeaks(optoData(x,:),MinPeakDistance=20,MinPeakProminence=Ithres,MinPeakWidth=options.peakWindow*10)),...
                    1:size(optoData,1),UniformOutput=false)');
    Eresponse_opto = cell2mat(arrayfun(@(x) length(findpeaks(-optoData(x,:),MinPeakDistance=20,MinPeakProminence=-Ethres,MinPeakWidth=options.peakWindow*10)),...
                    1:size(optoData,1),UniformOutput=false)');
    ISpotResponse_opto = accumarray(spotSequence(:), Iresponse_opto(:), [], @(x) {x});
    ESpotResponse_opto = accumarray(spotSequence(:), Eresponse_opto(:), [], @(x) {x});
    Irate_opto = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ISpotResponse_opto,UniformOutput=false));
    Erate_opto = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ESpotResponse_opto,UniformOutput=false));

    %% Plot figure
    disp(['Ongoing: plotting summary for search: ', curSearch,' at depth ',num2str(curDepth)]);
    initializeFig(1,1);
    masterLayout = tiledlayout(4,4);
    masterLayout.TileSpacing = 'compact';
    masterLayout.Padding = 'compact';

    % Plot current trace map
    nexttile(masterLayout,1); axis off;
    depthLayout = tiledlayout(masterLayout,2^curDepth,2^curDepth);
    depthLayout.Layout.Tile = 1;
    depthLayout.Layout.TileSpan = [4 2];
    depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight'; 
    % title(['Depth ', num2str(curDepth),': current responses']);

    for t = 1:4^curDepth
        nexttile(depthLayout,t);
        if isempty(depthCurrentMap{t})
            trace = nan(1,plotWindowLength);
            spotHotspot = false;
        else
            trace = depthCurrentMap{t}(:,plotFirstSample:plotLastSample);
            spotHotspot = sum(depthHotspot{t})>=1;
        end

        if spotHotspot
            plotSEM(plotWindowTime,trace,color,opacity=1,...
                    plotIndividual=true,LineWidth=LineWidth);
        else
            plotSEM(plotWindowTime,trace,color,opacity=0.3,...
                    plotIndividual=true,LineWidth=LineWidth);
        end
        xlim([options.timeRange(1), options.timeRange(2)]); 
        ylim([optoMin, optoMax]);
        plotEvent('',stimDuration,shadeOnly=true,color=stimColor,FaceAlpha=0.3,percentY=30,zeroValue=0);
        % Plot thresholds
        if search_vhold > -10; yline(Ithres,'--',color=blue,Alpha=0.5);
        elseif search_vhold < -50; yline(Ethres,'--',color=red,Alpha=0.5);
        else; yline(Ithres,'--',color=blue,Alpha=0.5); yline(Ethres,'--',color=red,Alpha=0.5);
        end
        % Plot axis at the bottom-left tile
        if t ~= 4^curDepth-2^curDepth+1; axis off
        else
            xlabel('ms'); ylabel('pA'); box off; 
            xticks([0,50]); 
            yticks(sort(unique([round(optoMin)-eps,round(Ethres),round(Ithres),round(optoMax)+eps])));
        end
    end

    % Plot response map
    ax_response = nexttile(masterLayout,3,[2 1]);
    imagesc(depthResponseMap); axis off; hold on;
    cb = colorbar; cb.Label.String = 'Total charge (pC)';
    climit = max(abs(depthResponseMap),[],'all');
    if climit == 0; climit = eps; end 
    clim([-climit, climit]); 
    colormap(ax_response,flip(blueWhiteRed));
    scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    title(['Depth ', num2str(curDepth),': ',curCell.Options{1}.feature]);

    % Plot hotspot map
    % ax_hotspot = nexttile(masterLayout,4,[2 1]); 
    % imagesc(depthHotspotMap); axis off; hold on;
    % clim([0, 1]); colormap(ax_hotspot,[.98,.98,.98; color]);
    % colorbar(Ticks=[0,1],TickLabels={'No response','Responded'});
    % scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    % title(['Depth ', num2str(curDepth),': responded spots']);

    % Plot opto vs prestim trace 
    nexttile(masterLayout,4);
    plotSEM(ctrlTime(end-200:end),ctrlData(:,end-200:end),...
            [.8 .8 .8],plotIndividual=true,...
            label='Baseline');
    plotSEM(optoTime,nullspotData,[.6 .6 .6],plotIndividual=true,label='Nullspot');
    plotSEM(optoTime,hotspotData,color,plotIndividual=true,label='Hotspot');
    xlabel('Time from stim (ms)'); ylabel('pA'), ylim([optoMin, optoMax]);
    legend(Location='best');
    title(['Depth ', num2str(curDepth),': opto vs baseline trace']);

    % Plot noise vs response distribution
    nexttile(masterLayout,8); 
    h_baseline = histogram(allNullData,'Normalization','pdf'); hold on
    h_baseline.FaceColor = [.8,.8,.8]; h_baseline.EdgeColor = [.8,.8,.8];
    if find(~hotspot_spotIdx)
        xline(spotMax_avg(~hotspot_spotIdx),':',color=blueWhiteRed(100,:)); hold on
        xline(spotMin_avg(~hotspot_spotIdx),':',color=blueWhiteRed(400,:)); hold on
    end
    if find(hotspot_spotIdx)
        xline(spotMax_avg(hotspot_spotIdx),'--',color=blue); hold on
        xline(spotMin_avg(hotspot_spotIdx),'--',color=red); hold on
    end
    xline(Ethres,'-',color=red,Label='Exci. threshold',LineWidth=5,LabelHorizontalAlignment='left'); hold on
    xline(Ithres,'-',color=blue,Label='Inhi. threshold',LineWidth=5); hold on; box off;
    xlabel('Current (pA)');
    title(['Depth ', num2str(curDepth),': noise vs response']);

    % Plot statistics
    nexttile(masterLayout,11);
    statLayout1 = tiledlayout(masterLayout,1,5);
    statLayout1.Layout.Tile = 11;
    statLayout1.Layout.TileSpan = [1 2];
    statLayout1.TileSpacing = 'compact'; statLayout1.Padding = 'tight'; axis off;

    % Plot hotspot vs nullspot min response
    nexttile(statLayout1,1,[1 2]);
    if find(hotspot_spotIdx); plotScatterBar(spotMin_avg(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(spotMin_avg(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    % Plot hotspot vs nullspot max response
    if find(hotspot_spotIdx); plotScatterBar(spotMax_avg(hotspot_spotIdx),3,color=blue); end
    if find(~hotspot_spotIdx); plotScatterBar(spotMax_avg(~hotspot_spotIdx),4,color=[.6 .6 .6]); end  
    xticks([1 2 3 4]); 
    xticklabels({'Exci. hotspot','Exci. nullspot','Inhi. hotspot','Inhi. nullspot'});
    ylabel('Current (pA)');
    title('Max response current');

    % Plot ctrl min/max response
    nexttile(statLayout1,3);
    plotScatterBar(ctrlMin_avg,1,color=[.8 .8 .8]);
    plotScatterBar(ctrlMax_avg,2,color=[.8 .8 .8]);
    xticks([1 2]); xticklabels({'Excitatory','Inhibitory'});
    ylabel('Current (pA)');
    title('Max ctrl current');

    % Plot hotspot vs nullspot AUC
    nexttile(statLayout1,4);
    if find(hotspot_spotIdx); plotScatterBar(spotAUC_avg(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(spotAUC_avg(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    plotScatterBar(ctrlAUC_avg,3,color=[.8 .8 .8]);
    xticks([1 2 3]); xticklabels({'Hotspot','Nullspot','Ctrl'});
    ylabel('Net total charge (pC)');
    title('Net total charge');

    % Plot hotspot vs nullspot absolute AUC
    nexttile(statLayout1,5);
    if find(hotspot_spotIdx); plotScatterBar(spotAbsAUC_avg(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(spotAbsAUC_avg(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    plotScatterBar(ctrlAbsAUC_avg,3,color=[.8 .8 .8]);
    xticks([1 2 3]); xticklabels({'Hotspot','Nullspot','Ctrl'});
    ylabel('Absolute total charge (pC)');
    title('Absolute total charge');

    % Other statistics
    nexttile(masterLayout,15);
    statLayout2 = tiledlayout(masterLayout,1,5);
    statLayout2.Layout.Tile = 15;
    statLayout2.Layout.TileSpan = [1 2];
    statLayout2.TileSpacing = 'compact'; statLayout2.Padding = 'tight'; axis off;

    % Plot hotspot vs nullspot minTime
    nexttile(statLayout2,1,[1 2]);
    if find(hotspot_spotIdx); plotScatterBar(spotMinTime_avg(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(spotMinTime_avg(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    % Plot hotspot vs nullspot maxTime
    if find(hotspot_spotIdx); plotScatterBar(spotMaxTime_avg(hotspot_spotIdx),3,color=blue); end
    if find(~hotspot_spotIdx); plotScatterBar(spotMaxTime_avg(~hotspot_spotIdx),4,color=[.6 .6 .6]); end
    xticks([1 2 3 4]); 
    xticklabels({'Exci. hotspot','Exci. nullspot','Inhi. hotspot','Inhi. nullspot'});
    ylabel('Time (ms)');
    title('Time to max response current');

    % Plot control minTime/maxTime
    nexttile(statLayout2,3);
    plotScatterBar(ctrlMinTime_avg,1,color=[.8 .8 .8]);
    plotScatterBar(ctrlMaxTime_avg,2,color=[.8 .8 .8]);
    xticks([1 2 3 4 5 6]); 
    xticklabels({'Excitatory','Inhibitory'});
    ylabel('Time (ms)');
    title('Time to max ctrl current');

    % Plot spot response rate vs noise response rate (excitatory)
    nexttile(statLayout2,4); 
    if find(hotspot_spotIdx); plotScatterBar(Erate_opto(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(Erate_opto(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    plotScatterBar(Erate_ctrl,3,color=[.8 .8 .8]);
    xticks([1 2 3]); xticklabels({'Hotspot opto','Nullspot opto','Ctrl'});
    ylabel('Excitatory response rate');
    title('Excitatory response rate');
    
    % Plot spot response rate vs noise response rate (inhibitory)
    nexttile(statLayout2,5); 
    if find(hotspot_spotIdx); plotScatterBar(Irate_opto(hotspot_spotIdx),1,color=blue); end
    if find(~hotspot_spotIdx); plotScatterBar(Irate_opto(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    plotScatterBar(Irate_ctrl,3,color=[.8 .8 .8]);
    xticks([1 2 3]); xticklabels({'Hotspot opto','Nullspot opto','Ctrl'});
    ylabel('Inhibitory response rate');
    title('Inhibitory response rate');

    %% Save figure
    filename = strcat(curCell.Epochs{1}{searchIdx},'_depth',num2str(curDepth));
    filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Search summary',curCell.Epochs{1}{searchIdx});
    saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath,...
                savePNG=options.savePNG,savePDF=options.savePDF,saveFIG=options.saveFIG);
    disp(['Finished: saving summary figure for search: ', curSearch,' at depth ',num2str(curDepth)]);
end

disp(['Finished: analysis finished for search: ', curSearch]);
close all;
end