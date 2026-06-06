function cells = analyzeSlice_DMD(expPath,options)

arguments
    expPath string

    options.resultsPathDate char = 'newest'
    options.saveDataPath string = 'default'

    options.save logical = true
    options.savePNG logical = true
    options.savePDF logical = true
    options.saveFIG logical = true

    options.plotQC logical = true
    options.plotSearch logical = true
    options.plotPairs logical = true

    options.timeRange double = [-10,50]
    options.outputFs double = 10000

    options.depthLineWidth double = [2,1.5,1,0.5,0.1,0.05,0.01];
end

%% Setup

[~,~,~,~,blueWhiteRed,~,~] = loadColors;

% Define results path
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

% Decide reload if session has already been loaded
dirsplit = split(expPath,filesep); expName = dirsplit{end};
try load(strcat(resultsPath,filesep,'cells_DMD_',expName,'.mat'));
catch
    error('Error: did not find cells_DMD.mat!');
end

% Define time window
cellsOptions = cells{1,'Options'}{1};
if options.outputFs ~= cellsOptions.outputFs
    warning('analyzeSlice_DMD: Default options.outputFs differs from outputFs extracted from cells_DMD. Using cells_DMD value instead!');
    options.outputFs = cellsOptions.outputFs;
end
if options.timeRange ~= cellsOptions.timeRange
    eventSample = abs(cellsOptions.timeRange(1))*(cellsOptions.outputFs/1000) + 1;
    plotFirstSample = eventSample + options.timeRange(1)*(cellsOptions.outputFs/1000);
    plotLastSample = eventSample + options.timeRange(2)*(cellsOptions.outputFs/1000);
    plotWindowLength = plotLastSample - plotFirstSample + 1;
    plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
else
    if isfield(cellsOptions,'plotWindowTime')
        plotWindowTime = cellsOptions.plotWindowTime;
        plotWindowLength = length(plotWindowTime);
    else
        if isfield(cellsOptions,'plotWindowLength'); plotWindowLength = cellsOptions.plotWindowLength;
        else; plotWindowLength = (options.timeRange(2)-options.timeRange(1))* options.outputFs/1000 + 1;
        end
        plotWindowTime = linspace(options.timeRange(1),options.timeRange(2),plotWindowLength);
    end
    eventSample = find(plotWindowTime==0);
    plotFirstSample = 1; plotLastSample = plotWindowLength;
end

%% Plot quality checks for hotspot detection

if options.plotQC
    analyzeSpots(cells);
end

%% Plot response map for each search

if options.plotSearch
    for c = 1:size(cells,1)
        disp(['Ongoing: plotting searches for cell',num2str(c)]);
        curCell = cells(cells.Cell == c,:);
        searchPerCell = length(curCell.Vhold{1});

        for search = 1:searchPerCell
            search_rmap = curCell.("Response map"){1}.responseMap{search};
            search_isResponse = curCell.("Response map"){1}.isResponseMap{search};
            search_cmap = curCell.("Response map"){1}.currentMap{search};
            search_vhold = curCell.Vhold{1}(search);

            % Initialize figure
            initializeFig(0.8,1); 
            masterLayout = tiledlayout(3,size(search_rmap,3));
            masterLayout.TileSpacing = 'compact';
            masterLayout.Padding = 'compact';

            for d = 1:size(search_rmap,3)
                depthResponseMap = search_rmap(:,:,d);
                isResponseMap_depth = search_isResponse(:,:,d);
                depthCurrentMap = search_cmap{d};

                % Determine current plot line width
                if d <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(d);
                else; LineWidth = options.depthLineWidth(end); end

                % Plot current trace map
                nexttile(masterLayout,d); axis off;
                title(['Depth ', num2str(d),': current responses']); 
                depthLayout = tiledlayout(masterLayout,2^d,2^d);
                depthLayout.Layout.Tile = d;
                depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight';
                maxCurrent = max(abs(depthCurrentMap),[],'all');
                for t = 1:4^d
                    nexttile(depthLayout,t);
                    if search_vhold < -50; color = blueWhiteRed(end,:);
                    else; color = blueWhiteRed(1,:); end
                    plotSEM(plotWindowTime,depthCurrentMap(t,plotFirstSample:plotLastSample),color,LineWidth=LineWidth);
                    xlim([options.timeRange(1), options.timeRange(2)]); ylim([-maxCurrent, maxCurrent]);
                    plotEvent('',5,shadeOnly=true,color=blueWhiteRed(1,:),FaceAlpha=0.3,...
                              percentY=30,zeroValue=depthCurrentMap(t,eventSample));
                    % Plot axis at the bottom-left tile
                    if t ~= 4^d-2^d+1; axis off
                    else
                        xlabel('ms'); ylabel('pA'); box off; 
                        xticks([0,50]); yticks([ceil(-maxCurrent),floor(maxCurrent)]);
                    end
                end

                % Plot response map
                ax_response = nexttile(masterLayout,d + size(search_rmap,3)); 
                imagesc(depthResponseMap); axis off;  
                cb = colorbar; cb.Label.String = 'Total charge (pC)';
                climit = max(abs(depthResponseMap),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': ',curCell.Options{1}.feature]);
        
                % Plot isResponse map
                ax_isResponse = nexttile(masterLayout, d + 2*size(search_rmap,3)); 
                imagesc(isResponseMap_depth); axis off;  clim([0, 1]);
                colorbar(Ticks=[0,1],TickLabels={'No response','Responded'});
                title(['Depth ', num2str(d),': responded spots']);
        
                % Set colormap
                colormap(ax_response,flip(blueWhiteRed));
                if curCell.Vhold{1}(search) < -30; colormap(ax_isResponse,blueWhiteRed(250:500,:)); 
                else; colormap(ax_isResponse,flip(blueWhiteRed(1:250,:))); end
            end

            % Save figure
            filename = curCell.Epochs{1}{search};
            filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Search summary');
            saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath,...
                        savePNG=options.savePNG,savePDF=options.savePDF,saveFIG=options.saveFIG);
            close all;
        end     
    end
end

%% Plot search comparisions

% Plot search pairs with same Vhold and different Vhold separately
if options.plotPairs
    for c = 1:size(cells,1)
        disp(['Ongoing: plotting search pairs for cell',num2str(c)]);
        curCell = cells(cells.Cell == c,:);
        allPairs = curCell.('Difference map'){1}.diffVhold;
        diffVholdIdx = find(allPairs);
        sameVholdIdx = find(~allPairs);
    
    
        % Plot search pairs with different Vhold
        for p = 1:length(diffVholdIdx)
            diffMap = curCell.('Difference map'){1}.response{diffVholdIdx(p)};
            commonSpots = curCell.('Difference map'){1}.commonSpots{diffVholdIdx(p)};

            % Find corresponding current map
            pairEpochs = curCell.('Difference map'){1}.pair{diffVholdIdx(p)};
            search1_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(1)};
            search2_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(2)};
            search1_vhold = curCell.Vhold{1}(pairEpochs(1));
            search2_vhold = curCell.Vhold{1}(pairEpochs(2));
    
            % Plot diff map and common map
            initializeFig(0.8,1);
            masterLayout = tiledlayout(4,size(diffMap,3));
            masterLayout.TileSpacing = 'compact';
            masterLayout.Padding = 'compact';

            for d = 1:size(diffMap,3)
                depthResponseMap = diffMap(:,:,d);
                depthCommonSpots = commonSpots(:,:,d);
                colormap(flip(blueWhiteRed));
                depthCurrentMap_search1 = search1_cmap{d};
                depthCurrentMap_search2 = search2_cmap{d};

                % Determine current plot line width
                if d <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(d);
                else; LineWidth = options.depthLineWidth(end); end

                % Plot current response map
                nexttile(masterLayout, d); axis off
                title(['Depth ', num2str(d),': current responses']); 
                depthLayout = tiledlayout(masterLayout,2^d,2^d);
                depthLayout.Layout.Tile = d;
                depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight';
                maxCurrent = max(abs([depthCurrentMap_search1 depthCurrentMap_search2]),[],'all');
                for t = 1:4^d
                    nexttile(depthLayout,t);
                    if search1_vhold < -50; color1 = blueWhiteRed(end,:);
                    else; color1 = blueWhiteRed(1,:); end
                    if search2_vhold < -50; color2 = blueWhiteRed(end,:);
                    else; color2 = blueWhiteRed(1,:); end
                    plotSEM(plotWindowTime,depthCurrentMap_search1(t,plotFirstSample:plotLastSample),color1,LineWidth=LineWidth);
                    plotSEM(plotWindowTime,depthCurrentMap_search2(t,plotFirstSample:plotLastSample),color2,LineWidth=LineWidth);
                    xlim([options.timeRange(1), options.timeRange(2)]); ylim([-maxCurrent, maxCurrent]);
                    plotEvent('',5,shadeOnly=true,color=blueWhiteRed(1,:),FaceAlpha=0.3,...
                              percentY=30,zeroValue=0);
                    % Plot axis at the bottom-left tile
                    if t ~= 4^d-2^d+1; axis off
                    else
                        xlabel('ms'); ylabel('pA'); box off; 
                        xticks([0,50]); yticks([ceil(-maxCurrent),floor(maxCurrent)]);
                    end
                end
                
                % Plot difference response map
                nexttile(masterLayout, d + size(diffMap,3)); 
                imagesc(depthResponseMap); axis off;  
                cb = colorbar; cb.Label.String = 'Net total charge (pC)';
                climit = max(abs(depthResponseMap),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature]);
                
                % Plot difference response map (common spot only)
                nexttile(masterLayout, d + 2*size(diffMap,3));
                depthCommonResponse = depthResponseMap;
                depthCommonResponse(depthCommonResponse & ~depthCommonSpots) = 0;
                imagesc(depthCommonResponse); axis off;  
                cb = colorbar; cb.Label.String = 'Net total charge (pC)';
                climit = max(abs(depthCommonResponse),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature, ' (common spot only)']);
                
                % Plot common spot map
                nexttile(masterLayout, d + 3*size(diffMap,3));
                depthCommonResponse(depthCommonResponse > 0) = 1;
                depthCommonResponse(depthCommonResponse < 0) = -1;
                imagesc(depthCommonResponse); axis off;  clim([-1, 1]);
                colorbar(Ticks=[-1, 0 ,1],TickLabels={'Net excitatory','No response','Net inhibitory'});
                title(['Depth ', num2str(d),': common responded spots']);
            end
    
            % Save figure
            filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Different Vhold pairs');
            filename = ['spots_cell',num2str(curCell.Cell),'_NetResponse_epoch',...
                        num2str(pairEpochs(1)),'vs',num2str(pairEpochs(2))];
            saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath,...
                        savePNG=options.savePNG,savePDF=options.savePDF,saveFIG=options.saveFIG);
            close all;
        end



        % Plot search pairs with the same Vhold
        for p = 1:length(sameVholdIdx)
            diffMap = curCell.('Difference map'){1}.response{sameVholdIdx(p)};
            commonSpots = curCell.('Difference map'){1}.commonSpots{sameVholdIdx(p)};

            % Find corresponding current map
            pairEpochs = curCell.('Difference map'){1}.pair{sameVholdIdx(p)};
            search1_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(1)};
            search2_cmap = curCell.("Response map"){1}.currentMap{pairEpochs(2)};
            search1_vhold = curCell.Vhold{1}(pairEpochs(1));
            search2_vhold = curCell.Vhold{1}(pairEpochs(2));
    
            % Plot diff map and common map
            initializeFig(0.8,1);
            masterLayout = tiledlayout(4,size(diffMap,3));
            masterLayout.TileSpacing = 'compact';
            masterLayout.Padding = 'compact';

            for d = 1:size(diffMap,3)
                depthResponseMap = diffMap(:,:,d);
                depthCommonSpots = commonSpots(:,:,d);
                colormap(flip(blueWhiteRed));
                depthCurrentMap_search1 = search1_cmap{d};
                depthCurrentMap_search2 = search2_cmap{d};

                % Determine current plot line width
                if d <= length(options.depthLineWidth); LineWidth = options.depthLineWidth(d);
                else; LineWidth = options.depthLineWidth(end); end

                % Plot current response map
                nexttile(masterLayout, d); axis off
                title(['Depth ', num2str(d),': current responses']); 
                depthLayout = tiledlayout(masterLayout,2^d,2^d);
                depthLayout.Layout.Tile = d;
                depthLayout.TileSpacing = 'none'; depthLayout.Padding = 'tight';
                maxCurrent = max(abs([depthCurrentMap_search1 depthCurrentMap_search2]),[],'all');
                for t = 1:4^d
                    nexttile(depthLayout,t);
                    if search1_vhold < -50; color1 = blueWhiteRed(end,:);
                    else; color1 = blueWhiteRed(1,:); end
                    if search2_vhold < -50; color2 = blueWhiteRed(end-150,:);
                    else; color2 = blueWhiteRed(150,:); end
                    plotSEM(plotWindowTime,depthCurrentMap_search1(t,plotFirstSample:plotLastSample),color1,LineWidth=LineWidth);
                    plotSEM(plotWindowTime,depthCurrentMap_search2(t,plotFirstSample:plotLastSample),color2,LineWidth=LineWidth);
                    xlim([options.timeRange(1), options.timeRange(2)]); ylim([-maxCurrent, maxCurrent]);
                    plotEvent('',5,shadeOnly=true,color=blueWhiteRed(1,:),FaceAlpha=0.3,...
                              percentY=30,zeroValue=0);
                    % Plot axis at the bottom-left tile
                    if t ~= 4^d-2^d+1; axis off
                    else
                        xlabel('ms'); ylabel('pA'); box off; 
                        xticks([0,50]); yticks([ceil(-maxCurrent),floor(maxCurrent)]);
                    end
                end
                
                % Plot difference response map
                nexttile(masterLayout, d + size(diffMap,3)); 
                imagesc(depthResponseMap); axis off;  
                cb = colorbar; cb.Label.String = '\Delta charge (pC)';
                climit = max(abs(depthResponseMap),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature]);
                
                % Plot difference response map (common spot only)
                nexttile(masterLayout, d + 2*size(diffMap,3));
                depthCommonResponse = depthResponseMap;
                depthCommonResponse(depthCommonResponse & ~depthCommonSpots) = 0;
                imagesc(depthCommonResponse); axis off;  
                cb = colorbar; cb.Label.String = '\Delta charge (pC)';
                climit = max(abs(depthCommonResponse),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature, ' (common spot only)']);
                
                % Plot common spot map
                nexttile(masterLayout, d + 3*size(diffMap,3));
                depthCommonResponse(depthCommonResponse > 0) = 1;
                depthCommonResponse(depthCommonResponse < 0) = -1;
                imagesc(depthCommonResponse); axis off;  clim([-1, 1]);
                colorbar(Ticks=[-1, 0 ,1],TickLabels={'Smaller','Same','Larger'});
                title(['Depth ', num2str(d),': common responded spots']);
            end
    
            % Save figure
            filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Same Vhold pairs');
            filename = ['spots_cell',num2str(curCell.Cell),'_DiffResponse_epoch',...
                        num2str(pairEpochs(1)),'vs',num2str(pairEpochs(2))];
            saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath,...
                        savePNG=options.savePNG,savePDF=options.savePDF,saveFIG=options.saveFIG);
            close all;
        end
    end
end

close all
end