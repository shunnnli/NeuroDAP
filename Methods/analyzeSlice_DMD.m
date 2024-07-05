function cells = analyzeSlice_DMD(expPath,options)

arguments
    expPath string

    options.save logical = true
    options.saveDataPath string = 'default'

    options.plotSearch logical = true
    options.plotPairs
end

%% Setup

[~,~,~,~,blueWhiteRed,~,~] = loadColors;

if strcmp(options.saveDataPath,'default')
    options.saveDataPath = fullfile(expPath,'Results');
end

% Decide reload if session has already been loaded
dirsplit = split(expPath,filesep); expName = dirsplit{end};
try load(strcat(expPath,filesep,'cells_DMD_',expName,'.mat'));
catch
    error('Error: did not find cells_DMD.mat!');
end

%% Plot response map for each search

if options.plotSearch
    for c = 1:size(cells,1)
        curCell = cells(cells.Cell == c,:);
        searchPerCell = length(curCell.Vhold{1});

        for search = 1:searchPerCell
            search_rmap = curCell.("Response map"){1}.responseMap{search};
            search_isResponse = curCell.("Response map"){1}.isResponseMap{search};

            initializeFig(0.8,1); tiledlayout(2,size(search_rmap,3));
            for d = 1:size(search_rmap,3)
                depthResponseMap = search_rmap(:,:,d);
                isResponseMap_depth = search_isResponse(:,:,d);

                % Plot response map
                ax_response = nexttile(d); 
                imagesc(depthResponseMap); axis off;  
                cb = colorbar; cb.Label.String = 'Total charge (pC)';
                climit = max(abs(depthResponseMap),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': ',curCell.Options{1}.feature]);
        
                % Plot isResponse map
                ax_isResponse = nexttile(d + size(search_rmap,3)); 
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
            saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath);
        end
    end
end

%% Plot search comparisions

% Plot search pairs with same Vhold and different Vhold separately
if options.plotPairs
    for c = 1:size(cells,1)
        curCell = cells(cells.Cell == c,:);
        allPairs = curCell.('Difference map'){1}.Vhold;
        diffVholdIdx = find(allPairs);
        sameVholdIdx = find(~allPairs);
    
        % (Unifinished) Plot search pairs with the same Vhold
        for p = 1:length(sameVholdIdx)
            diffMap = curCell.('Difference map'){1}.response{sameVholdIdx(p)};
            commonSpots = curCell.('Difference map'){1}.commonSpots{sameVholdIdx(p)};
    
            % Plot diff map and common map
            initializeFig(0.8,1); tiledlayout(3,size(diffMap,3));
            for d = 1:size(diffMap,3)
                depthResponseMap = diffMap(:,:,d);
                depthCommonSpots = commonSpots(:,:,d);
                colormap(flip(blueWhiteRed));
                
                % Plot difference response map
                nexttile(d); 
                imagesc(depthResponseMap); axis off;  
                cb = colorbar; cb.Label.String = '\Delta charge (pC)';
                climit = max(abs(depthResponseMap),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature]);
                
                % % Plot difference response map (common spot only)
                nexttile(d + size(diffMap,3));
                depthCommonResponse = depthResponseMap;
                depthCommonResponse(depthCommonResponse & ~depthCommonSpots) = 0;
                imagesc(depthCommonResponse); axis off;  
                cb = colorbar; cb.Label.String = '\Delta charge (pC)';
                climit = max(abs(depthCommonResponse),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature, ' (common spot only)']);
                
                % Plot common spot map
                nexttile(d + 2*size(diffMap,3));
                depthCommonResponse(depthCommonResponse > 0) = 1;
                depthCommonResponse(depthCommonResponse < 0) = -1;
                imagesc(depthCommonResponse); axis off;  clim([-1, 1]);
                colorbar(Ticks=[-1, 0 ,1],TickLabels={'Smaller','Same','Larger'});
                title(['Depth ', num2str(d),': common responded spots']);
            end
    
            % Save figure
            pairEpochs = curCell.('Difference map'){1}.pair{sameVholdIdx(p)};
            filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Same Vhold pairs');
            filename = ['spots_cell',num2str(curCell.Cell),'_DiffResponse_epoch',...
                        num2str(pairEpochs(1)),'vs',num2str(pairEpochs(2))];
            saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath);
        end
    
    
        % Plot search pairs with different Vhold
        for p = 1:length(diffVholdIdx)
            diffMap = curCell.('Difference map'){1}.response{diffVholdIdx(p)};
            commonSpots = curCell.('Difference map'){1}.commonSpots{diffVholdIdx(p)};
    
            % Plot diff map and common map
            initializeFig(0.8,1); tiledlayout(3,size(diffMap,3));
            for d = 1:size(diffMap,3)
                depthResponseMap = diffMap(:,:,d);
                depthCommonSpots = commonSpots(:,:,d);
                colormap(flip(blueWhiteRed));
                
                % Plot difference response map
                nexttile(d); 
                imagesc(depthResponseMap); axis off;  
                cb = colorbar; cb.Label.String = 'Net total charge (pC)';
                climit = max(abs(depthResponseMap),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature]);
                
                % % Plot difference response map (common spot only)
                nexttile(d + size(diffMap,3));
                depthCommonResponse = depthResponseMap;
                depthCommonResponse(depthCommonResponse & ~depthCommonSpots) = 0;
                imagesc(depthCommonResponse); axis off;  
                cb = colorbar; cb.Label.String = 'Net total charge (pC)';
                climit = max(abs(depthCommonResponse),[],'all'); 
                if climit == 0; climit = eps; end 
                clim([-climit, climit]);
                title(['Depth ', num2str(d),': \Delta',curCell.Options{1}.feature, ' (common spot only)']);
                
                % Plot common spot map
                nexttile(d + 2*size(diffMap,3));
                depthCommonResponse(depthCommonResponse > 0) = 1;
                depthCommonResponse(depthCommonResponse < 0) = -1;
                imagesc(depthCommonResponse); axis off;  clim([-1, 1]);
                colorbar(Ticks=[-1, 0 ,1],TickLabels={'Net excitatory','No response','Net inhibitory'});
                title(['Depth ', num2str(d),': common responded spots']);
            end
    
            % Save figure
            pairEpochs = curCell.('Difference map'){1}.pair{diffVholdIdx(p)};
            filepath = fullfile(options.saveDataPath,['cell',num2str(curCell.Cell)],'Different Vhold pairs');
            filename = ['spots_cell',num2str(curCell.Cell),'_NetResponse_epoch',...
                        num2str(pairEpochs(1)),'vs',num2str(pairEpochs(2))];
            saveFigures(gcf,[filename,'_',curCell.Options{1}.feature],filepath);
        end
    end
end

close all
end