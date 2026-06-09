function cells = analyzeSlice_DMD_Paolo(expPath,options)

arguments
    expPath string

    options.resultsPathDate char = 'newest'
    options.saveDataPath string = 'default'

    options.redStim logical = true

    options.save logical = true
    options.savePNG logical = true
    options.savePDF logical = true
    options.saveFIG logical = true

    options.plotSearch logical = true
    options.plotPairs logical = false

    options.timeRange double = [-10,50]

    % options.depthLineWidth double = [2,1.5,1,0.5,0.1,0.05,0.01];
end

%% Setup

% Define results path
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

% Decide reload if session has already been loaded
dirsplit = split(expPath,filesep); expName = dirsplit{end};
try load(strcat(resultsPath,filesep,'cells_DMD_',expName,'.mat'));
catch
    error('Error: did not find cells_DMD.mat!');
end

% Remove empty rows of cells
cells = rmmissing(cells,DataVariables='Session');

% % Define time window
% cellsOptions = cells{1,'Options'}{1};
% if options.outputFs ~= cellsOptions.outputFs
%     warning('analyzeSlice_DMD: Default options.outputFs differs from outputFs extracted from cells_DMD. Using cells_DMD value instead!');
%     options.outputFs = cellsOptions.outputFs;
% end
% if options.timeRange ~= cellsOptions.timeRange
%     eventSample = abs(cellsOptions.timeRange(1))*(cellsOptions.outputFs/1000) + 1;
%     plotFirstSample = eventSample + options.timeRange(1)*(cellsOptions.outputFs/1000);
%     plotLastSample = eventSample + options.timeRange(2)*(cellsOptions.outputFs/1000);
%     plotWindowLength = plotLastSample - plotFirstSample + 1;
%     plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
% else
%     if isfield(cellsOptions,'plotWindowTime')
%         plotWindowTime = cellsOptions.plotWindowTime;
%         plotWindowLength = length(plotWindowTime);
%     else
%         if isfield(cellsOptions,'plotWindowLength'); plotWindowLength = cellsOptions.plotWindowLength;
%         else; plotWindowLength = (options.timeRange(2)-options.timeRange(1))* options.outputFs/1000 + 1;
%         end
%         plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
%     end
%     eventSample = find(plotWindowTime==0);
%     plotFirstSample = 1; plotLastSample = plotWindowLength;
% end

%% Plot summary figure for each search

if options.plotSearch
    cellList = unique(cells.Cell);

    for cellIdx = 1:size(cells,1)
        c = cellList(cellIdx);
        disp(['Ongoing: plotting searches for cell',num2str(c)]);
        curCell = cells(cells.Cell == c,:);
        searchPerCell = length(curCell.Vhold{1});

        for searchIdx = 1:searchPerCell
            % Plot a figure for each search depth, organize depth figures of each 
            % search into a single folder
            analyzeDMDSearch_Paolo(curCell,searchIdx,...
                         redStim=options.redStim,...
                         timeRange=options.timeRange,...
                         saveDataPath=options.saveDataPath,...
                         savePNG=false,savePDF=true,saveFIG=false);
        end
    end
end

%% Plot search comparisions

% Plot search pairs with same Vhold and different Vhold separately
if options.plotPairs
    cellList = unique(cells.Cell);

    for cellIdx = 1:size(cells,1)
        c = cellList(cellIdx);
        disp(['Ongoing: plotting search pairs for cell',num2str(c)]);
        curCell = cells(cells.Cell == c,:);
        allPairs = curCell.('Difference map'){1}.diffVhold;

        for pairIdx = 1:length(allPairs)
            % Plot a figure for each search depth, organize depth figures of each 
            % search into a single folder
            analyzeDMDSearchPair(curCell,pairIdx,...
                         redStim=options.redStim,...
                         timeRange=options.timeRange,...
                         saveDataPath=options.saveDataPath,...
                         savePNG=false,savePDF=true,saveFIG=false);
        end
    end
end

close all
end

%% Old code

    % cellList = unique(cells.Cell);
    % for cellIdx = 1:size(cells,1)
    %     c = cellList(cellIdx);
    %     disp(['Ongoing: plotting search pairs for cell',num2str(c)]);
    %     curCell = cells(cells.Cell == c,:);
    %     cellLoc = curCell.Options{1}.cellLocation;
    %     allPairs = curCell.('Difference map'){1}.diffVhold;
    %     diffVholdIdx = find(allPairs);
    %     sameVholdIdx = find(~allPairs);
    % 
    %     % Plot search pairs with different Vhold
    %     for p = 1:length(diffVholdIdx)
    %         diffMap = curCell.('Difference map'){1}.response{diffVholdIdx(p)};
    %         commonSpots = curCell.('Difference map'){1}.commonSpots{diffVholdIdx(p)};
    %         commonDepths = curCell.('Difference map'){1}.commonDepths{diffVholdIdx(p)};
    % 
    %         % Find corresponding current map
    %         pairEpochs = curCell.('Difference map'){1}.pair{diffVholdIdx(p)};
    %         search1_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(1)};
    %         search2_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(2)};
    %         search1_vhold = curCell.Vhold{1}(pairEpochs(1));
    %         search2_vhold = curCell.Vhold{1}(pairEpochs(2));
    % 
    %         % Get stim duration
    %         if strcmp('Protocol',curCell.Properties.VariableNames)
    %             stimDuration1 = curCell.("Protocol"){1}{pairEpochs(1)}{1}.pulseWidth;
    %             stimDuration2 = curCell.("Protocol"){1}{pairEpochs(2)}{1}.pulseWidth;
    %         else
    %             stimDuration1 = 5; stimDuration2 = 5;
    %         end
    % 
    %         % Plot diff map and common map
    %         initializeFig(0.8,1);
    %         masterLayout = tiledlayout(4,size(diffMap,3));
    %         masterLayout.TileSpacing = 'compact';
    %         masterLayout.Padding = 'compact';
    % 
    %         for d = 1:size(diffMap,3)
    %             curDepth = commonDepths(d);
    %             depthResponseMap = diffMap(:,:,d);
    %             depthCommonSpots = commonSpots(:,:,d);
    %             colormap(flip(blueWhiteRed));
    %             depthCurrentMap_search1 = search1_cmap{curDepth};
    %             depthCurrentMap_search2 = search2_cmap{curDepth};
    % 
    %             % Determine current plot line width
    %             if curDepth <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(curDepth);
    %             else; LineWidth = options.depthLineWidth(end); end
    % 
    %             % Plot current response map
    %             nexttile(masterLayout, d); axis off
    %             title(['Depth ', num2str(curDepth),': current responses']); 
    %             depthLayout = tiledlayout(masterLayout,2^curDepth,2^curDepth);
    %             depthLayout.Layout.Tile = d;
    %             depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight';
    % 
    %             maxSpotResponse = cell2mat(cellfun(@(x) max(x,[],'all'),...
    %                 [depthCurrentMap_search1;depthCurrentMap_search2],UniformOutput=false));
    %             minSpotResponse = cell2mat(cellfun(@(x) min(x,[],'all'),...
    %                 [depthCurrentMap_search1;depthCurrentMap_search2],UniformOutput=false));
    %             absMaxCurrent = max([maxSpotResponse -minSpotResponse],[],'all');
    % 
    %             for t = 1:4^curDepth
    %                 nexttile(depthLayout,t);
    %                 if search1_vhold < -50; color1 = blueWhiteRed(end,:);
    %                 elseif search1_vhold > -10; color1 = blueWhiteRed(1,:); 
    %                 else; color1 = purple; 
    %                 end
    %                 if search2_vhold < -50; color2 = blueWhiteRed(end,:);
    %                 elseif search2_vhold > -10; color2 = blueWhiteRed(1,:); 
    %                 else; color2 = purple; 
    %                 end
    %                 if isempty(depthCurrentMap_search1{t})
    %                     trace1 = nan(1,plotWindowLength);
    %                 else
    %                     trace1 = depthCurrentMap_search1{t}(:,plotFirstSample:plotLastSample);
    %                 end
    %                 if isempty(depthCurrentMap_search2{t})
    %                     trace2 = nan(1,plotWindowLength);
    %                 else
    %                     trace2 = depthCurrentMap_search2{t}(:,plotFirstSample:plotLastSample);
    %                 end
    %                 plotSEM(plotWindowTime,trace1,color1,...
    %                     plotIndividual=true,individualColor='same',individualAlpha=0.3,...
    %                     LineWidth=LineWidth);
    %                 plotSEM(plotWindowTime,trace2,color2,...
    %                     plotIndividual=true,individualColor='same',individualAlpha=0.3,...
    %                     LineWidth=LineWidth);
    %                 xlim([options.timeRange(1), options.timeRange(2)]); 
    %                 ylim([-absMaxCurrent-eps, absMaxCurrent+eps]);
    %                 plotEvent('',stimDuration1,shadeOnly=true,color=stimColor,FaceAlpha=0.3,percentY=20,zeroValue=absMaxCurrent/2);
    %                 plotEvent('',stimDuration2,shadeOnly=true,color=stimColor,FaceAlpha=0.3,percentY=20,zeroValue=-absMaxCurrent/2);
    %                 % Plot axis at the bottom-left tile
    %                 if t ~= 4^curDepth-2^curDepth+1; axis off
    %                 else
    %                     xlabel('ms'); ylabel('pA'); box off; 
    %                     xticks([0,50]); yticks([ceil(-absMaxCurrent)-eps,floor(absMaxCurrent)+eps]);
    %                 end
    %             end
    % 
    %             % Plot difference response map
    %             nexttile(masterLayout, d + size(diffMap,3)); 
    %             imagesc(depthResponseMap); axis off; hold on;
    %             cb = colorbar; cb.Label.String = 'Net total charge (pC)';
    %             climit = max(abs(depthResponseMap),[],'all'); 
    %             if climit == 0; climit = eps; end 
    %             clim([-climit, climit]);
    %             scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    %             title(['Depth ', num2str(curDepth),': \Delta',curCell.Options{1}.feature]);
    % 
    %             % Plot difference response map (common spot only)
    %             nexttile(masterLayout, d + 2*size(diffMap,3));
    %             depthCommonResponse = depthResponseMap;
    %             depthCommonResponse(depthCommonResponse & ~depthCommonSpots) = 0;
    %             imagesc(depthCommonResponse); axis off; hold on;
    %             cb = colorbar; cb.Label.String = 'Net total charge (pC)';
    %             climit = max(abs(depthCommonResponse),[],'all'); 
    %             if climit == 0; climit = eps; end 
    %             clim([-climit, climit]);
    %             scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    %             title(['Depth ', num2str(curDepth),': \Delta',curCell.Options{1}.feature, ' (common spot only)']);
    % 
    %             % Plot common spot map
    %             nexttile(masterLayout, d + 3*size(diffMap,3));
    %             depthCommonResponse(depthCommonResponse > 0) = 1;
    %             depthCommonResponse(depthCommonResponse < 0) = -1;
    %             imagesc(depthCommonResponse); axis off; hold on;
    %             clim([-1, 1]); 
    %             colorbar(Ticks=[-1, 0 ,1],TickLabels={'Net excitatory','No response','Net inhibitory'});
    %             scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    %             title(['Depth ', num2str(curDepth),': common responded spots']);
    %         end
    % 
    %         % Save figure
    %         filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Different Vhold pairs');
    %         filename = ['spots_cell',num2str(curCell.Cell),'_NetResponse_epoch',...
    %                     num2str(pairEpochs(1)),'vs',num2str(pairEpochs(2))];
    %         saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath,...
    %                     savePNG=options.savePNG,savePDF=options.savePDF,saveFIG=options.saveFIG);
    %         close all;
    %     end
    % 
    % 
    % 
    %     % Plot search pairs with the same Vhold
    %     for p = 1:length(sameVholdIdx)
    %         diffMap = curCell.('Difference map'){1}.response{sameVholdIdx(p)};
    %         commonSpots = curCell.('Difference map'){1}.commonSpots{sameVholdIdx(p)};
    %         commonDepths = curCell.('Difference map'){1}.commonDepths{diffVholdIdx(p)};
    % 
    %         % Find corresponding current map
    %         pairEpochs = curCell.('Difference map'){1}.pair{sameVholdIdx(p)};
    %         search1_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(1)};
    %         search2_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(2)};
    %         search1_vhold = curCell.Vhold{1}(pairEpochs(1));
    %         search2_vhold = curCell.Vhold{1}(pairEpochs(2));
    % 
    %         % Plot diff map and common map
    %         initializeFig(0.8,1);
    %         masterLayout = tiledlayout(4,size(diffMap,3));
    %         masterLayout.TileSpacing = 'compact';
    %         masterLayout.Padding = 'compact';
    % 
    %         for d = 1:size(diffMap,3)
    %             curDepth = commonDepths(d);
    %             depthResponseMap = diffMap(:,:,d);
    %             depthCommonSpots = commonSpots(:,:,d);
    %             colormap(flip(blueWhiteRed));
    %             depthCurrentMap_search1 = search1_cmap{d};
    %             depthCurrentMap_search2 = search2_cmap{d};
    % 
    %             % Determine current plot line width
    %             if curDepth <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(curDepth);
    %             else; LineWidth = options.depthLineWidth(end); end
    % 
    %             % Plot current response map
    %             nexttile(masterLayout, d); axis off
    %             title(['Depth ', num2str(curDepth),': current responses']); 
    %             depthLayout = tiledlayout(masterLayout,2^d,2^d);
    %             depthLayout.Layout.Tile = d;
    %             depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight';
    % 
    %             maxSpotResponse = cell2mat(cellfun(@(x) max(x,[],'all'),...
    %                 [depthCurrentMap_search1;depthCurrentMap_search2],UniformOutput=false));
    %             minSpotResponse = cell2mat(cellfun(@(x) min(x,[],'all'),...
    %                 [depthCurrentMap_search1;depthCurrentMap_search2],UniformOutput=false));
    %             absMaxCurrent = max([maxSpotResponse -minSpotResponse],[],'all');
    % 
    %             for t = 1:4^curDepth
    %                 nexttile(depthLayout,t);
    %                 if search1_vhold < -50; color1 = blueWhiteRed(end,:);
    %                 else; color1 = blueWhiteRed(1,:); end
    %                 if search2_vhold < -50; color2 = blueWhiteRed(end-150,:);
    %                 else; color2 = blueWhiteRed(150,:); end
    %                 if isempty(depthCurrentMap_search1{t})
    %                     trace1 = nan(1,plotWindowLength);
    %                 else
    %                     trace1 = depthCurrentMap_search1{t}(:,plotFirstSample:plotLastSample);
    %                 end
    %                 if isempty(depthCurrentMap_search2{t})
    %                     trace2 = nan(1,plotWindowLength);
    %                 else
    %                     trace2 = depthCurrentMap_search2{t}(:,plotFirstSample:plotLastSample);
    %                 end
    %                 plotSEM(plotWindowTime,trace1,color1,...
    %                     plotIndividual=true,individualColor='same',individualAlpha=0.3,...
    %                     LineWidth=LineWidth);
    %                 plotSEM(plotWindowTime,trace2,color2,...
    %                     plotIndividual=true,individualColor='same',individualAlpha=0.3,...
    %                     LineWidth=LineWidth);
    %                 xlim([options.timeRange(1), options.timeRange(2)]); 
    %                 ylim([-absMaxCurrent-eps, absMaxCurrent+eps]);
    %                 plotEvent('',stimDuration1,shadeOnly=true,color=stimColor,FaceAlpha=0.3,percentY=20,zeroValue=absMaxCurrent/2);
    %                 plotEvent('',stimDuration2,shadeOnly=true,color=stimColor,FaceAlpha=0.3,percentY=20,zeroValue=-absMaxCurrent/2);
    %                 % Plot axis at the bottom-left tile
    %                 if t ~= 4^curDepth-2^curDepth+1; axis off
    %                 else
    %                     xlabel('ms'); ylabel('pA'); box off; 
    %                     xticks([0,50]); yticks([ceil(-absMaxCurrent)-eps,floor(absMaxCurrent)+eps]);
    %                 end
    %             end
    % 
    %             % Plot difference response map
    %             nexttile(masterLayout, d + size(diffMap,3)); 
    %             imagesc(depthResponseMap); axis off; hold on;
    %             cb = colorbar; cb.Label.String = '\Delta charge (pC)';
    %             climit = max(abs(depthResponseMap),[],'all'); 
    %             if climit == 0; climit = eps; end 
    %             clim([-climit, climit]);
    %             scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    %             title(['Depth ', num2str(curDepth),': \Delta',curCell.Options{1}.feature]);
    % 
    %             % Plot difference response map (common spot only)
    %             nexttile(masterLayout, d + 2*size(diffMap,3));
    %             depthCommonResponse = depthResponseMap;
    %             depthCommonResponse(depthCommonResponse & ~depthCommonSpots) = 0;
    %             imagesc(depthCommonResponse); axis off; hold on; 
    %             cb = colorbar; cb.Label.String = '\Delta charge (pC)';
    %             climit = max(abs(depthCommonResponse),[],'all'); 
    %             if climit == 0; climit = eps; end 
    %             clim([-climit, climit]);
    %             scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    %             title(['Depth ', num2str(curDepth),': \Delta',curCell.Options{1}.feature, ' (common spot only)']);
    % 
    %             % Plot common spot map
    %             nexttile(masterLayout, d + 3*size(diffMap,3));
    %             depthCommonResponse(depthCommonResponse > 0) = 1;
    %             depthCommonResponse(depthCommonResponse < 0) = -1;
    %             imagesc(depthCommonResponse); axis off; hold on;
    %             clim([-1, 1]); 
    %             colorbar(Ticks=[-1, 0 ,1],TickLabels={'Smaller','Same','Larger'});
    %             scatter(cellLoc(1),cellLoc(2),50,'o','filled','MarkerFaceColor','k','MarkerEdgeColor','none');
    %             title(['Depth ', num2str(curDepth),': common responded spots']);
    %         end
    % 
    %         % Save figure
    %         filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Same Vhold pairs');
    %         filename = ['spots_cell',num2str(curCell.Cell),'_DiffResponse_epoch',...
    %                     num2str(pairEpochs(1)),'vs',num2str(pairEpochs(2))];
    %         saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath,...
    %                     savePNG=options.savePNG,savePDF=options.savePDF,saveFIG=options.saveFIG);
    %         close all;
    %     end
    % end