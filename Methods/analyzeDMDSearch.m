function analyzeDMDSearch(curCell, searchIdx, options)

% Create spots summary file and plot summary figure

arguments
    curCell table
    searchIdx double

    options.redStim logical = true
    options.depthLineWidth = 'scale' %[3,2.5,2,1.5,1.1,1,0.5];
    options.depthLineWidthAlpha double = 0.3

    options.color
    options.colormap

    options.save logical = true
    options.savePNG logical = true
    options.savePDF logical = true
    options.saveFIG logical = true

    options.outputFs double = 10000
    options.timeRange double = [-10,50] % in ms
    options.analysisWindowLength double = 30 % in ms after stim onset
    options.controlWindowLength double = 30 % in ms before stim onset
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
expPath = osPathSwitch(strsplit{1});

% Define results path
if strcmp(options.saveDataPath,'default')
    resultsList = sortrows(struct2cell(dir(fullfile(expPath,'Results-*')))',3);
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
if ~isfield(options,'colormap'); options.colormap = blueWhiteRed; end
if options.redStim; stimColor = blueWhiteRed(end,:);
else; stimColor = blueWhiteRed(1,:); end
blue = blueWhiteRed(1,:);
red = blueWhiteRed(end,:);
purple = [232, 130, 250]./255;

% Update window lengths
options.analysisWindowLength = curCell.Options{1}.analysisWindowLength;
options.controlWindowLength  = curCell.Options{1}.controlWindowLength;

% Define time windows
if options.outputFs ~= curCell.Options{1}.outputFs
    warning('analyzeDMDSearch: Default options.outputFs differs from outputFs extracted from cells_DMD. Using cells_DMD value instead!');
    options.outputFs = curCell.Options{1}.outputFs;
end

% Sampling helper (samples per ms)
samplesPerMs = options.outputFs/1000;

if options.timeRange ~= curCell.Options{1}.timeRange
    % Use the baseline implied by the *stored* timeRange to locate the event index.
    % Round to ensure integer sample indices.
    eventSample = round(abs(curCell.Options{1}.timeRange(1)) * samplesPerMs) + 1;
    plotFirstSample = round(eventSample + options.timeRange(1) * samplesPerMs);
    plotLastSample  = round(eventSample + options.timeRange(2) * samplesPerMs);
    plotWindowLength = plotLastSample - plotFirstSample + 1;
    plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
    analysisWindow = eventSample : (eventSample + round(options.analysisWindowLength*samplesPerMs));
else
    if isfield(curCell.Options{1},'plotWindowTime')
        plotWindowTime = curCell.Options{1}.plotWindowTime;
        plotWindowLength = length(plotWindowTime);
    else
        if isfield(curCell.Options{1},'plotWindowLength'); plotWindowLength = curCell.Options{1}.plotWindowLength;
        else; plotWindowLength = round((options.timeRange(2)-options.timeRange(1)) * samplesPerMs) + 1;
        end
        plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
    end
    % Avoid exact floating-point comparisons when finding the 0-ms index.
    [~, eventSample] = min(abs(plotWindowTime));
    plotFirstSample = 1; plotLastSample = plotWindowLength;
    analysisWindow = eventSample : (eventSample + round(options.analysisWindowLength*samplesPerMs));
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
search_spotLocation = curCell.("Response map"){1}.spotLocation{searchIdx};
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
if isfield(options,'color')
    color = options.color;
else
    if search_vhold < -50
        color = blueWhiteRed(end,:); 
    elseif search_vhold > -10
        color = blueWhiteRed(1,:); 
    else
        color = purple; 
    end
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

    % ========================= CHANGE (Full grid for analysis) =========================
    % The saved maps (depthCurrentMap/depthBaselineMap/depthHotspot) may only include
    % tiles that were actually sampled (e.g., 4x4 = 16). For analysis/plotting, we
    % expand them to the full 2^curDepth x 2^curDepth grid (length 4^curDepth),
    % leaving unsampled tiles empty.
    [depthCurrentMapFull, ~, depthHotspotFull, ~] = ...
        expandDepthToFullGrid(depthCurrentMap, depthBaselineMap, depthHotspot, search_spotLocation{d}, depthResponseMap, curDepth);

    nCol = 2^curDepth;
    nRow = 2^curDepth;
    % ======================= END CHANGE (Full grid for analysis) =======================


    % Determine current plot line width
    if ~ischar(options.depthLineWidth)
        if curDepth <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(curDepth);
        else; LineWidth = options.depthLineWidth(end); end
    elseif contains(options.depthLineWidth,'scale')
        nTiles = nRow * nCol;   % total subplots on this depth
        nRef  = 4;              % "first depth" reference
        LWref = 3;              % linewidth when nTiles == 4
        alpha = options.depthLineWidthAlpha;            % shrink speed (tune this)
        
        LineWidth = LWref * (nRef / nTiles)^alpha;
        LineWidth = max(0.2, min(3, LineWidth));
    end

    % Get opto and prestim traces
    optoData = cell2mat(depthCurrentMap); 
    ctrlData = cell2mat(depthBaselineMap);

    % Robustly re-align indices if the requested timeRange/analysis window does not fit
    % the actual trace length. This most commonly happens for reconstructed data from .fig
    % files (shorter time range), but can also happen if options.timeRange differs from
    % what was used to generate the stored maps.
    nAvailable = size(optoData, 2);
    needRealign = (plotFirstSample < 1) || (plotLastSample > nAvailable) || (max(analysisWindow) > nAvailable);
    % needRealign = needRealign || (exist('plotWindowLength','var') && plotWindowLength ~= nAvailable);

    if needRealign
        % Infer the event sample using the requested pre-stim baseline implied by timeRange.
        % If the baseline is longer than what we actually have, clamp to the available length.
        baselineSamplesWanted = round(abs(options.timeRange(1)) * samplesPerMs);
        baselineSamples = min(max(baselineSamplesWanted, 0), nAvailable - 1);
        newEventSample = baselineSamples + 1;

        % Optional message for the most common reconstruction case.
        if nAvailable == 601 && options.outputFs == 10000
            disp('[analyzeDMDSearch] Short time-range trace detected (601 samples @ 10 kHz). Re-aligning indices...');
        elseif baselineSamples ~= baselineSamplesWanted
            warning('[analyzeDMDSearch] Requested baseline (%d samples) exceeds available data. Clamping alignment to %d samples.', ...
                baselineSamplesWanted, baselineSamples);
        else
            warning('[analyzeDMDSearch] Time indexing exceeds available data. Re-aligning indices to fit trace length.');
        end

        % 1) Recalculate analysis window (inclusive, in samples)
        analysisLen = round(options.analysisWindowLength * samplesPerMs) + 1;
        analysisWindow = newEventSample : min(nAvailable, newEventSample + analysisLen - 1);

        % 2) Recalculate plotting indices, clamped to available data
        reqFirst = newEventSample + round(options.timeRange(1) * samplesPerMs);
        reqLast  = newEventSample + round(options.timeRange(2) * samplesPerMs);

        plotFirstSample = max(1, reqFirst);
        plotLastSample  = min(nAvailable, reqLast);

        % 3) Update time vector to match clamped samples (ms)
        plotWindowTime = ((plotFirstSample:plotLastSample) - newEventSample) / samplesPerMs;
        plotWindowLength = numel(plotWindowTime);
    end

    % Slice data
    optoData = optoData(:,analysisWindow);
    ctrlData = ctrlData(:,size(ctrlData,2)-length(analysisWindow)+1:size(ctrlData,2));

    % Make control slicing robust (reconstructed data can have shorter baselines)
    % ctrlEnd = size(ctrlData,2);
    % ctrlStart = max(1, ctrlEnd - numel(analysisWindow) + 1);
    % ctrlData = ctrlData(:, ctrlStart:ctrlEnd);

    % Time vectors (ms) based on actual available samples
    optoTime = (0:(size(optoData,2)-1)) / samplesPerMs;
    ctrlTime = (-(size(ctrlData,2)-1):0) / samplesPerMs;

    % Calculate max/min current for plotting
    % optoMax = optoData(~isoutlier(optoData)); 
    % optoMin = optoData(~isoutlier(optoData)); 
    % optoMax = max(optoData,[],'all'); optoMin = min(optoData,[],'all');
    [optoMin, optoMax] = getYlimit(optoData, Ethres=Ethres, Ithres=Ithres, k_sem=5);

    % (Optional) Change hotspot criteria
    % XXXXXXXXX

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
    % Convert peak detection parameters from ms to samples (avoid hard-coded 10 kHz assumptions)
    minPeakDistanceSamples = max(1, round(2 * samplesPerMs));
    peakWidthSamples = max(1, round(options.peakWindow * samplesPerMs));

    Iresponse_ctrl = cell2mat(arrayfun(@(x) length(findpeaks(ctrlData(x,:),MinPeakDistance=minPeakDistanceSamples,MinPeakProminence=Ithres,MinPeakWidth=peakWidthSamples)),...
                    1:size(ctrlData,1),UniformOutput=false)');
    Eresponse_ctrl = cell2mat(arrayfun(@(x) length(findpeaks(-ctrlData(x,:),MinPeakDistance=minPeakDistanceSamples,MinPeakProminence=-Ethres,MinPeakWidth=peakWidthSamples)),...
                    1:size(ctrlData,1),UniformOutput=false)');
    ISpotResponse_ctrl = accumarray(spotSequence(:), Iresponse_ctrl(:), [], @(x) {x});
    ESpotResponse_ctrl = accumarray(spotSequence(:), Eresponse_ctrl(:), [], @(x) {x});
    Irate_ctrl = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ISpotResponse_ctrl,UniformOutput=false));
    Erate_ctrl = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ESpotResponse_ctrl,UniformOutput=false));
    % Opto success rate
    Iresponse_opto = cell2mat(arrayfun(@(x) length(findpeaks(optoData(x,:),MinPeakDistance=minPeakDistanceSamples,MinPeakProminence=Ithres,MinPeakWidth=peakWidthSamples)),...
                    1:size(optoData,1),UniformOutput=false)');
    Eresponse_opto = cell2mat(arrayfun(@(x) length(findpeaks(-optoData(x,:),MinPeakDistance=minPeakDistanceSamples,MinPeakProminence=-Ethres,MinPeakWidth=peakWidthSamples)),...
                    1:size(optoData,1),UniformOutput=false)');
    ISpotResponse_opto = accumarray(spotSequence(:), Iresponse_opto(:), [], @(x) {x});
    ESpotResponse_opto = accumarray(spotSequence(:), Eresponse_opto(:), [], @(x) {x});
    Irate_opto = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ISpotResponse_opto,UniformOutput=false));
    Erate_opto = cell2mat(cellfun(@(x) length(find(x))/size(x,1), ESpotResponse_opto,UniformOutput=false));

    %% Plot figure
    disp(['Ongoing: plotting summary for search: ', curSearch,' at depth ',num2str(curDepth)]);
    close all;
    initializeFig(1,1);
    masterLayout = tiledlayout(4,4);
    masterLayout.TileSpacing = 'compact';
    masterLayout.Padding = 'compact';

    % Plot current trace map
    nexttile(masterLayout,1); axis off;
    depthLayout = tiledlayout(masterLayout,nRow,nCol);
    depthLayout.Layout.Tile = 1;
    depthLayout.Layout.TileSpan = [4 2];
    depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight'; 
    % title(['Depth ', num2str(curDepth),': current responses']);

    % ========================= CHANGE (Full grid for analysis) =========================
    for t = 1:(4^curDepth)
        nexttile(depthLayout,t);
        if isempty(depthCurrentMapFull{t})
            trace = nan(1,plotWindowLength);
            spotHotspot = false;
        else
            trace = depthCurrentMapFull{t}(:,plotFirstSample:plotLastSample);
            spotHotspot = sum(depthHotspotFull{t})>=1;
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
        threshold_lineWidth = max(0.01,LineWidth-1);
        if search_vhold > -10; yline(Ithres,'--',color=blue,Alpha=0.5,LineWidth=threshold_lineWidth);
        elseif search_vhold < -50; yline(Ethres,'--',color=red,Alpha=0.5,LineWidth=threshold_lineWidth);
        else; yline(Ithres,'--',color=blue,Alpha=0.5); yline(Ethres,'--',color=red,Alpha=0.5,LineWidth=threshold_lineWidth);
        end
        % Plot axis at the bottom-left tile
        bottomLeft = (nRow-1)*nCol + 1;
        if t ~= bottomLeft; axis off
        else
            xlabel('ms'); ylabel('pA'); box off; 
            xticks([0,50]); 
            yticks(sort(unique([round(optoMin)-eps,round(Ethres),round(Ithres),round(optoMax)+eps])));
        end
    end
    % ======================= END CHANGE (Full grid for analysis) =======================

    % Plot response map
    ax_response = nexttile(masterLayout,3,[2 1]);
    imagesc(depthResponseMap); axis off; hold on;
    cb = colorbar; cb.Label.String = 'Total charge (pC)';
    applyDMDColormap(ax_response, depthResponseMap, options.colormap, search_vhold);
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
    if find(hotspot_spotIdx); plotScatterBar(1,spotMin_avg(hotspot_spotIdx),color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(2,spotMin_avg(~hotspot_spotIdx),color=[.6 .6 .6]); end
    % Plot hotspot vs nullspot max response
    if find(hotspot_spotIdx); plotScatterBar(3,spotMax_avg(hotspot_spotIdx),color=blue); end
    if find(~hotspot_spotIdx); plotScatterBar(4,spotMax_avg(~hotspot_spotIdx),color=[.6 .6 .6]); end  
    xticks([1 2 3 4]); 
    xticklabels({'Min hotspot','Min nullspot','Max hotspot','Max nullspot'});
    ylabel('Current (pA)');
    title('Max response current');

    % Plot ctrl min/max response
    nexttile(statLayout1,3);
    plotScatterBar(1,ctrlMin_avg,color=[.8 .8 .8]);
    plotScatterBar(2,ctrlMax_avg,color=[.8 .8 .8]);
    xticks([1 2]); xticklabels({'Excitatory','Inhibitory'});
    ylabel('Current (pA)');
    title('Max ctrl current');

    % Plot hotspot vs nullspot AUC
    nexttile(statLayout1,4);
    if find(hotspot_spotIdx); plotScatterBar(1,spotAUC_avg(hotspot_spotIdx),color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(2,spotAUC_avg(~hotspot_spotIdx),color=[.6 .6 .6]); end
    plotScatterBar(3,ctrlAUC_avg,color=[.8 .8 .8]);
    xticks([1 2 3]); xticklabels({'Hotspot','Nullspot','Ctrl'});
    ylabel('Net total charge (pC)');
    title('Net total charge');

    % Plot hotspot vs nullspot absolute AUC
    nexttile(statLayout1,5);
    if find(hotspot_spotIdx); plotScatterBar(1,spotAbsAUC_avg(hotspot_spotIdx),color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(2,spotAbsAUC_avg(~hotspot_spotIdx),color=[.6 .6 .6]); end
    plotScatterBar(3,ctrlAbsAUC_avg,color=[.8 .8 .8]);
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
    if find(hotspot_spotIdx); plotScatterBar(1,spotMinTime_avg(hotspot_spotIdx),color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(2,spotMinTime_avg(~hotspot_spotIdx),color=[.6 .6 .6]); end
    % Plot hotspot vs nullspot maxTime
    if find(hotspot_spotIdx); plotScatterBar(3,spotMaxTime_avg(hotspot_spotIdx),color=blue); end
    if find(~hotspot_spotIdx); plotScatterBar(4,spotMaxTime_avg(~hotspot_spotIdx),color=[.6 .6 .6]); end
    xticks([1 2 3 4]); 
    xticklabels({'Min hotspot','Min nullspot','Max hotspot','Max nullspot'});
    ylabel('Time (ms)');
    title('Time to max response current');

    % Plot control minTime/maxTime
    nexttile(statLayout2,3);
    plotScatterBar(1,ctrlMinTime_avg,color=[.8 .8 .8]);
    plotScatterBar(2,ctrlMaxTime_avg,color=[.8 .8 .8]);
    xticks([1 2 3 4 5 6]); 
    xticklabels({'Excitatory','Inhibitory'});
    ylabel('Time (ms)');
    title('Time to max ctrl current');

    % Plot spot response rate vs noise response rate (excitatory)
    nexttile(statLayout2,4); 
    if find(hotspot_spotIdx); plotScatterBar(1,Erate_opto(hotspot_spotIdx),color=red); end
    if find(~hotspot_spotIdx); plotScatterBar(2,Erate_opto(~hotspot_spotIdx),color=[.6 .6 .6]); end
    plotScatterBar(3,Erate_ctrl,color=[.8 .8 .8]);
    xticks([1 2 3]); xticklabels({'Hotspot opto','Nullspot opto','Ctrl'});
    ylabel('Excitatory response rate');
    title('Excitatory response rate');
    
    % Plot spot response rate vs noise response rate (inhibitory)
    nexttile(statLayout2,5); 
    if find(hotspot_spotIdx); plotScatterBar(1,Irate_opto(hotspot_spotIdx),color=blue); end
    if find(~hotspot_spotIdx); plotScatterBar(2,Irate_opto(~hotspot_spotIdx),color=[.6 .6 .6]); end
    plotScatterBar(3,Irate_ctrl,color=[.8 .8 .8]);
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

%% ========================= CHANGE (Full grid for analysis) =========================
function [curFull, baseFull, hotFull, locFull] = expandDepthToFullGrid(curMap, baseMap, hotMap, spotLoc, respMap, depth)
% Expand per-depth maps to a full 2^depth x 2^depth grid (length 4^depth).
% Unsampled tiles remain empty ([] for cell entries; NaN for locations).
%
% NOTE on coordinate convention in this codebase:
%   location = [colStart colEnd rowStart rowEnd] in 1-based pixel indices, where:
%     - colStart/colEnd index the 2nd dimension (columns) of respMap
%     - rowStart/rowEnd index the 1st dimension (rows) of respMap
%
% This matches how loadSlicesDMD/analyzeDMDSearch index regions.

    nRowsPix = size(respMap, 1);
    nColsPix = size(respMap, 2);

    colStartsFull = localSplitStarts(nColsPix, depth); % 1 x 2^depth
    rowStartsFull = localSplitStarts(nRowsPix, depth); % 1 x 2^depth

    nCol = 2^depth;
    nRow = 2^depth;
    nFull = 4^depth;

    curFull  = cell(nFull, 1);
    baseFull = cell(nFull, 1);
    hotFull  = cell(nFull, 1);
    locFull  = nan(nFull, 4);

    % spotLoc may be stored either as Nx4 numeric or as a cell-wrapped numeric.
    if iscell(spotLoc)
        spotLoc = spotLoc{1};
    end

    nMap = numel(curMap);
    % Ensure spotLoc has at least nMap rows
    if size(spotLoc,1) < nMap
        nMap = size(spotLoc,1);
    end

    for k = 1:nMap
        location = spotLoc(k,:);
        if any(isnan(location))
            continue
        end

        colIdx = find(colStartsFull == location(1), 1);
        rowIdx = find(rowStartsFull == location(3), 1);

        % If an exact match is not found (should be rare), fall back to nearest start.
        if isempty(colIdx)
            [~, colIdx] = min(abs(colStartsFull - location(1)));
        end
        if isempty(rowIdx)
            [~, rowIdx] = min(abs(rowStartsFull - location(3)));
        end

        t = (colIdx - 1) + (rowIdx - 1) * nCol + 1;
        if t < 1 || t > nFull
            continue
        end

        curFull{t}  = curMap{k};
        baseFull{t} = baseMap{k};
        if k <= numel(hotMap)
            hotFull{t}  = hotMap{k};
        end
        locFull(t,:) = location;
    end
end

%% 

function starts = localSplitStarts(totalLen, depth)
% Returns 1-based start indices for the 2^depth segments produced by repeatedly splitting
% each segment into floor(L/2) and (L-floor(L/2)). This matches the typical quadtree
% splitting logic used in DMD search pattern generation.
    lens = totalLen;
    for i = 1:depth
        newLens = zeros(1, numel(lens) * 2);
        for j = 1:numel(lens)
            L = lens(j);
            left = floor(L/2);
            right = L - left;
            newLens(2*j-1) = left;
            newLens(2*j)   = right;
        end
        lens = newLens;
    end
    starts = cumsum([1, lens(1:end-1)]);
end

%% Define optoMax or optoMin

function [optoMin, optoMax] = getYlimit(data, options)

    arguments
        data double
        options.Ethres = nan
        options.Ithres = nan
        options.k_sem double = 3 % room around mean (3*SEM). Try 2–5 depending on how tight you want it.
        options.pad double = 0.1
    end

    % --- Robust y-limits based on mean trace (not outliers) ---
    mu  = mean(data, 1, 'omitnan');
    sig = std(data, 0, 1, 'omitnan');
    sem = sig / sqrt(size(data,1));

    yLo = min(mu - options.k_sem*sem);
    yHi = max(mu + options.k_sem*sem);

    % Make sure thresholds are visible 
    yLo = min([yLo, options.Ethres, options.Ithres], [], 'omitnan');
    yHi = max([yHi, options.Ethres, options.Ithres], [], 'omitnan');
    
    % Add a little extra padding
    pad = options.pad * (yHi - yLo);
    if pad == 0, pad = 1; end
    optoMin = yLo - pad;
    optoMax = yHi + pad;

end
