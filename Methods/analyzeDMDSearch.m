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
    options.timeRange double = [-20,100] % in ms
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
cellResultsPath = fullfile(curCell.Session,['cell',num2str(curCell.Cell)]);
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
    warning('analyzeSlice_DMD: Default options.outputFs differs from outputFs extracted from cells_DMD. Using cells_DMD value instead!');
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
search_isResponse = curCell.("Response map"){1}.isResponseMap{searchIdx};
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
if search_vhold < -50; color = blueWhiteRed(end,:);
elseif search_vhold > -10; color = blueWhiteRed(1,:); 
else; color = purple; 
end

% Load noise model of this neuron
load(fullfile(expPath,['noise_cell',num2str(curCell.Cell),'.mat']),'allNullData');
Ethres = curCell.Stats{1}.Ethres; %-2 * noise_all.sigma;
Ithres = curCell.Stats{1}.Ithres; %2 * noise_all.sigma;

% Load spots.mat for this neuron
depthfilename = strcat(curSearch,'_depth*.mat');
spotsList = dir(fullfile(cellResultsPath,depthfilename));

%% Generate a summary figure for each search

for d = 3%1:nDepth
    %% Initialization
    curDepth = search_depths(d);
    disp(['Ongoing: calculating statistics for search: ', curSearch,' at depth ',num2str(curDepth)]);

    % Load spots.mat file of this depth
    % depthFilePath = fullfile(spotsList(d).folder,spotsList(d).name);
    % load(depthFilePath,'spotsAtDepth');

    % Load depth specific data
    depthResponseMap = search_rmap(:,:,d);
    isResponseMap_depth = search_isResponse(:,:,d);
    depthCurrentMap = search_cmap{d};
    depthBaselineMap = search_bmap{d};
    spotSequence = repelem(1:height(search_hotspot{d}), cellfun(@numel, search_hotspot{d}))';

    % Determine current plot line width
    if curDepth <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(curDepth);
    else; LineWidth = options.depthLineWidth(end); end

    % Get opto and prestim traces
    optoData = cell2mat(depthCurrentMap); 
    optoData = optoData(:,analysisWindow);
    ctrlData = cell2mat(depthBaselineMap);
    optoTime = linspace(0,options.analysisWindowLength,size(optoData,2));
    ctrlTime = linspace(-options.controlWindowLength,0,size(ctrlData,2));

    % Calculate max/min current for plotting
    optoMax = optoData(~isoutlier(optoData)); 
    optoMin = optoData(~isoutlier(optoData)); 
    optoMax = max(optoData,[],'all'); optoMin = min(optoData,[],'all');

    % Get hotspot vs null spot
    hotspot_spotIdx = cell2mat(cellfun(@(x) sum(x)>=1, search_hotspot{d}, UniformOutput=false));
    hotspot_sweepIdx = cell2mat(cellfun(@(x) repmat(sum(x)>=1,size(x)), search_hotspot{d}, UniformOutput=false));
    hotspotData = optoData(hotspot_sweepIdx,:);
    nullspotData = optoData(~hotspot_sweepIdx,:);

    % Get max/min responses for each spot
    spotMax_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.max{searchIdx}{d},UniformOutput=false));
    spotMin_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.min{searchIdx}{d},UniformOutput=false));
    spotMaxTime_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.maxTime{searchIdx}{d},UniformOutput=false));
    spotMinTime_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.minTime{searchIdx}{d},UniformOutput=false));
    spotAUC_avg = cell2mat(cellfun(@(x) mean(x,"all"), curCell.Stats{1}.auc{searchIdx}{d},UniformOutput=false));
    % ctrlMax_avg = cell2mat(cellfun(@(x) max(x,[],'all'), depthBaselineMap,UniformOutput=false));
    % ctrlMin_avg = cell2mat(cellfun(@(x) min(x,[],'all'), depthBaselineMap,UniformOutput=false));

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
        else
            trace = depthCurrentMap{t}(:,plotFirstSample:plotLastSample);
        end
        plotSEM(plotWindowTime,trace,color,...
                plotIndividual=true,individualColor='same',individualAlpha=0.3,...
                LineWidth=LineWidth);
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
            yticks(sort([ceil(optoMin)-eps,round(Ethres),round(Ithres),floor(optoMax)+eps]));
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

    % Plot isResponse map
    ax_isResponse = nexttile(masterLayout,4,[2 1]); 
    imagesc(isResponseMap_depth); axis off; hold on;
    clim([0, 1]); colormap(ax_isResponse,[.98,.98,.98; color]);
    colorbar(Ticks=[0,1],TickLabels={'No response','Responded'});
    scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    title(['Depth ', num2str(curDepth),': responded spots']);

    % Plot opto vs prestim trace 
    nexttile(masterLayout,11);
    plotSEM(ctrlTime(end-200:end),ctrlData(:,end-200:end),[.8 .8 .8],...
            plotIndividual=true,individualColor='same',individualAlpha=0.3,...
            label='Baseline');
    plotSEM(optoTime,nullspotData,[.6 .6 .6],plotIndividual=true,individualColor='same',individualAlpha=0.3,label='Nullspot');
    plotSEM(optoTime,hotspotData,color,plotIndividual=true,individualColor='same',individualAlpha=0.3,label='Hotspot');
    xlabel('Time from stim (ms)'); ylabel('pA'), ylim([optoMin, optoMax]);
    legend(Location='best');
    title(['Depth ', num2str(curDepth),': opto vs baseline trace']);

    % Plot noise vs response distribution
    nexttile(masterLayout,12); 
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
    nexttile(masterLayout,15);
    statLayout = tiledlayout(masterLayout,1,5);
    statLayout.Layout.Tile = 15;
    statLayout.Layout.TileSpan = [1 2];
    statLayout.TileSpacing = 'compact'; statLayout.Padding = 'tight'; axis off;

    % Plot hotspot vs nullspot minTime
    nexttile(statLayout,1);
    if find(hotspot_spotIdx); plotScatterBar(spotMinTime_avg(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(spotMinTime_avg(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    xticks([1 2]); xticklabels({'Hotspot','Nullspot'});
    ylabel('Time to min (ms)');
    title('Time to minimum');

    % Plot hotspot vs nullspot maxTime
    nexttile(statLayout,2);
    if find(hotspot_spotIdx); plotScatterBar(spotMaxTime_avg(hotspot_spotIdx),1,color=blue); end
    if find(~hotspot_spotIdx); plotScatterBar(spotMaxTime_avg(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    xticks([1 2]); xticklabels({'Hotspot','Nullspot'});
    ylabel('Time to max (ms)');
    title('Time to maximum');

    % Plot hotspot vs nullspot AUC
    nexttile(statLayout,3);
    if find(hotspot_spotIdx); plotScatterBar(spotAUC_avg(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(spotAUC_avg(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    xticks([1 2]); xticklabels({'Hotspot','Nullspot'});
    ylabel('AUC (pC)');
    title('AUC');

    % Plot spot response rate vs noise response rate (excitatory)
    nexttile(statLayout,4); 
    if find(hotspot_spotIdx); plotScatterBar(Erate_opto(hotspot_spotIdx),1,color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(Erate_opto(~hotspot_spotIdx),2,color=[.6 .6 .6]); end
    plotScatterBar(Erate_ctrl,3,color=[.8 .8 .8]);
    xticks([1 2 3]); xticklabels({'Hotspot opto','Nullspot opto','Ctrl'});
    ylabel('Excitatory response rate');
    title('Excitatory response rate');
    
    % Plot spot response rate vs noise response rate (inhibitory)
    nexttile(statLayout,5); 
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

disp(['Finished: analyzeSpot finished for search: ', curSearch]);
close all;
end