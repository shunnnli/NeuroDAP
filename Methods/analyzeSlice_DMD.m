function cells = analyzeSlice_DMD(expPath,options)

arguments
    expPath string
    options.feature = 'auc'

    options.reload logical = false

    options.save logical = true
    options.saveDataPath string = 'default'

    options.plotDepthResponseMap logical = true
end

%% Setup

[~,~,~,~,blueWhiteRed,~,~] = loadColors;
searchList = sortrows(struct2cell(dir(fullfile(expPath,'spots*.mat')))',3);

if strcmp(options.saveDataPath,'default')
    options.saveDataPath = expPath;
end

% Decide reload if session has already been loaded
dirsplit = split(expPath,filesep); expName = dirsplit{end};
if ~isempty(dir(fullfile(expPath,"cells_DMD*.mat")))
    if ~options.reload
        disp('Loading stop: cells_DMD.mat file found.');
        load(strcat(expPath,filesep,'cells_DMD_',expName,'.mat'));
    end
else
    options.reload = true;
end

%% Create cells_DMD.mat

if options.reload
    disp('Ongoing: building cells_DMD.mat');
    %% Initialize cells table
    
    % structure for responseMap
    % For each row, response map is a cell with #searches element
    % within each element, theres a 3-dim matrix (xRange,yRange,depth)
    
    varTypes = {'string','string','string','double','cell','cell',...
                'cell','cell','cell'};
    varNames = {'Session','Animal','Task','Cell','Epochs','Vhold',...
                'Response map','Difference map','Options'};
    cells = table('Size',[1,length(varNames)],...
        'VariableTypes',varTypes,'VariableNames',varNames);
    
    % Initialize other params
    cellResponseMap = {}; isResponseMap_cell = {}; 
    cellVhold = []; cellEpochs = {};
    
    %% Loop through all searches and create cells.mat
    
    for search = 1:length(searchList)
    
        % Find corresponding spots.mat
        load(fullfile(searchList{search,2},searchList{search,1}));
    
        % Initialize params
        depthList = unique(spots.Depth);
        searchResponseMap = zeros(608,684,length(depthList));
        isResponseMap_search = zeros(608,684,length(depthList));
        cellVhold = [cellVhold; spots{1,'Vhold'}];
        cellEpochs{end+1} = {searchList{search,1}(1:end-4)};
        
        % % Plot response map for every depth
        % if options.plotDepthResponseMap
        %     initializeFig(0.75,1); tiledlayout(2,length(depthList));
        % end
    
        % Loop through depth
        for d = 1:length(depthList)
            % Initialization
            depthResponseMap = zeros(608,684);
            isResponseMap_depth = zeros(608,684);
            spotsAtDepth = spots(spots.Depth == d,:);
        
            % Build response map
            for s = 1:size(spotsAtDepth,1)
                % Get response value
                location = spotsAtDepth{s,'Location'}{1};
                xRange = location(1):location(2);
                yRange = location(3):location(4);
                originalValue = mean(depthResponseMap(xRange,yRange),'all');
                newValue = spotsAtDepth{s,'Stats'}{1}.response.(options.feature);
        
                % Add to depthResponseMap
                if originalValue==0; depthResponseMap(xRange,yRange) = newValue;
                else; depthResponseMap(xRange,yRange) = mean([originalValue,newValue]);
                end
    
                % Get isResponse value
                originalValue_isResponse = mode(isResponseMap_depth(xRange,yRange),'all');
                newValue_isResponse = spotsAtDepth{s,'Response'}{1}.isResponse;
        
                % Add to isResponseMap
                if originalValue_isResponse==0; isResponseMap_depth(xRange,yRange) = newValue_isResponse;
                else; isResponseMap_depth(xRange,yRange) = originalValue_isResponse || newValue_isResponse;
                end
            end
        
            % if options.plotDepthResponseMap
            %     % Plot response map
            %     ax_response = nexttile(d); 
            %     imagesc(depthResponseMap); axis off;  
            %     c = colorbar; c.Label.String = 'Total charge (pC)';
            %     climit = max(abs(depthResponseMap),[],'all'); clim([-climit, climit]);
            %     title(['Depth ', num2str(d),': ',options.feature]);
            % 
            %     % Plot isResponse map
            %     ax_isResponse = nexttile(d + length(depthList)); 
            %     imagesc(isResponseMap_depth); axis off;  clim([0, 1]);
            %     title(['Depth ', num2str(d),': responded spots']);
            % 
            %     % Set colormap
            %     colormap(ax_response,flip(blueWhiteRed));
            %     if spots{1,'Vhold'} < -30; colormap(ax_isResponse,blueWhiteRed(250:500,:)); 
            %     else; colormap(ax_isResponse,flip(blueWhiteRed(1:250,:))); end
            % end
    
            % Save depthResponseMap
            searchResponseMap(:,:,d) = depthResponseMap;
            isResponseMap_search(:,:,d) = isResponseMap_depth;
        end
    
        % Save response map
        cellResponseMap{end+1} = searchResponseMap;
        isResponseMap_cell{end+1} = isResponseMap_search;
    
        % if options.plotDepthResponseMap
        %     saveFigures(gcf,[searchList{search,1}(1:end-4),'_',options.feature],expPath);
        % end
    
        % Save to cell
        if search==length(searchList) || spots{1,'Cell'} ~= str2double(searchList{search+1,1}(strfind(searchList{search+1,1},'cell')+4))
            responses.responseMap = cellResponseMap';
            responses.isResponseMap = isResponseMap_cell';
    
            cells{spots{1,'Cell'},'Session'} = spots{1,'Session'};
            cells{spots{1,'Cell'},'Animal'} = spots{1,'Animal'};
            cells{spots{1,'Cell'},'Task'} = spots{1,'Task'};
            cells{spots{1,'Cell'},'Cell'} = spots{1,'Cell'};
            cells{spots{1,'Cell'},'Epochs'} = {cellEpochs'};
            cells{spots{1,'Cell'},'Vhold'} = {cellVhold};
            cells{spots{1,'Cell'},'Response map'} = {responses};
            cells{spots{1,'Cell'},'Options'} = {options};
            disp(['Finished: saving data for cell ',num2str(spots{1,'Cell'})]);
    
            % Update prevCell and reset cellResponseMap
            cellVhold = []; cellEpochs = {};
            cellResponseMap = {}; isResponseMap_cell = {};
        end
    end
    
    %% Cell specific analysis
    
    % The main goal here is to find common spots and calculate there difference
    % in AUC or other feature
    
    for c = 1:size(cells,1)
        cell = cells(cells.Cell == c,:);
        searchPerCell = length(cell.Vhold{1});
        diff_rmap_cell = {}; common_isResponse_cell = {}; 
        diff_pairs = {}; diff_vholds = [];
        
        % Analyze all pairs of searches
        for search1 = 1:searchPerCell
            for search2 = search1+1:searchPerCell
                % Get response map
                search1_rmap = cell.("Response map"){1}.responseMap{search1};
                search2_rmap = cell.("Response map"){1}.responseMap{search2};
                search1_isResponse = cell.("Response map"){1}.isResponseMap{search1};
                search2_isResponse = cell.("Response map"){1}.isResponseMap{search2};
    
                % Check max common depth
                maxCommonDepth = min([size(search1_rmap,3),size(search2_rmap,3)]);
                search1_rmap = search1_rmap(:,:,1:maxCommonDepth);
                search2_rmap = search2_rmap(:,:,1:maxCommonDepth);
                search1_isResponse = search1_isResponse(:,:,1:maxCommonDepth);
                search2_isResponse = search2_isResponse(:,:,1:maxCommonDepth);
    
                % Check whether Vhold are the same
                diffVhold = cell.Vhold{1}(search1) ~= cell.Vhold{1}(search2);
    
                % Initialize difference map
                diff_rmap = zeros(size(search1_rmap));
                common_isResponse = zeros(size(search1_rmap));
    
    
                % Calculate difference map for each depth
                for d = 1:maxCommonDepth
                    if diffVhold
                        % Add two response map together to calculate difference
                        diff_rmap(:,:,d) = search1_rmap(:,:,d) + search2_rmap(:,:,d);
                    else
                        % Subtract two response map to calculate differnece
                        diff_rmap(:,:,d) = search1_rmap(:,:,d) - search2_rmap(:,:,d);
                    end
                    common_isResponse(:,:,d) = search1_isResponse(:,:,d) & search2_isResponse(:,:,d);
                end
    
                % Save difference map
                diff_rmap_cell{end+1} = diff_rmap;
                common_isResponse_cell{end+1} = common_isResponse;
                diff_pairs{end+1} = [search1, search2];
                diff_vholds = [diff_vholds; diffVhold];
            end 
        end
    
        % Save in cells.mat
        diff.response = diff_rmap_cell';
        diff.commonSpots = common_isResponse_cell';
        diff.pair = diff_pairs';
        diff.Vhold = diff_vholds;
        cells{c,'Difference map'} = {diff};
    end
    
    %% Save cells.mat
    
    if options.save
        dirsplit = split(expPath,filesep); expName = dirsplit{end};
        save(strcat(options.saveDataPath,filesep,'cells_DMD_',expName),'cells','-v7.3');
        disp(strcat("Created cells_DMD.mat: ",expName));
    end
end

%% Plot response map for each search

if options.plotDepthResponseMap
    for c = 1:size(cells,1)
        cell = cells(cells.Cell == c,:);
        searchPerCell = length(cell.Vhold{1});

        for search = 1:searchPerCell
            search_rmap = cell.("Response map"){1}.responseMap{search};
            search_isResponse = cell.("Response map"){1}.isResponseMap{search};

            initializeFig(0.75,1); tiledlayout(2,size(search_rmap,3));
            for d = 1:size(search_rmap,3)
                depthResponseMap = search_rmap(:,:,d);
                isResponseMap_depth = search_isResponse(:,:,d);

                % Plot response map
                ax_response = nexttile(d); 
                imagesc(depthResponseMap); axis off;  
                c = colorbar; c.Label.String = 'Total charge (pC)';
                climit = max(abs(depthResponseMap),[],'all'); clim([-climit, climit]);
                title(['Depth ', num2str(d),': ',cell.Options{1}.feature]);
        
                % Plot isResponse map
                ax_isResponse = nexttile(d + size(search_rmap,3)); 
                imagesc(isResponseMap_depth); axis off;  clim([0, 1]);
                colorbar(Ticks=[0,1],TickLabels={'No response','Responded'});
                title(['Depth ', num2str(d),': responded spots']);
        
                % Set colormap
                colormap(ax_response,flip(blueWhiteRed));
                if cell.Vhold{1}(search) < -30; colormap(ax_isResponse,blueWhiteRed(250:500,:)); 
                else; colormap(ax_isResponse,flip(blueWhiteRed(1:250,:))); end
            end

            % Save figure
            filename = cell.Epochs{1}{search};
            saveFigures(gcf,[filename,'_',cell.Options{1}.feature],expPath);
        end
    end
end

%% Plotting code

% initializeFig(0.5,0.6); tiledlayout(2,3);
% d = 4; colormap(flip(blueWhiteRed));
% nexttile; imagesc(search1_rmap(:,:,d)); colorbar;
% climit = max(abs(search1_rmap),[],'all'); clim([-climit, climit]);
% nexttile; imagesc(search2_rmap(:,:,d)); colorbar;
% climit = max(abs(search2_rmap),[],'all'); clim([-climit, climit]);
% nexttile; diffMap = search1_rmap(:,:,d) - search2_rmap(:,:,d);
% imagesc(diffMap); colorbar;
% climit = max(abs(diffMap),[],'all'); clim([-climit, climit]);
% 
% nexttile; imagesc(search1_isResponse(:,:,d)); colorbar;
% nexttile; imagesc(search2_isResponse(:,:,d)); colorbar;
% nexttile; common_isResponse = search1_isResponse(:,:,d) & search2_isResponse(:,:,d);
% imagesc(common_isResponse); colorbar;

close all
end