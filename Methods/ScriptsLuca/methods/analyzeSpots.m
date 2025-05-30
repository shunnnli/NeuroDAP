function analyzeSpots(cells, options)

% Create spots summary file and plot summary figure

arguments
    cells table

    options.outputFs double = 10000
    options.timeRange double = [-20,100] % in ms
    options.analysisWindow double = 50 % in ms after stim onset
    options.controlWindow double = 50 % in ms before stim onset
    options.eventSample % in sample
    options.nArtifactSamples double = 0 % in sample
    options.rcCheckRecoveryWindow double = 100 % in ms
    options.peakWindow double = 1 % in ms around the peak to average

    options.feature = 'auc'
    options.thresholdFactor double = 3 % 3*std

    options.saveDataPath string = 'default'
end

%% Setup

[twoColors,~,~,~,~,~,bluePurpleRed] = loadColors;

% Define results path UNIFINSIHED!!!!!!!!
% resultsList = sortrows(struct2cell(dir(fullfile(expPath,'Results_*')))',3);
% if strcmp(options.resultsPathDate,'newest')
%     resultsPath = fullfile(resultsList{end,2},resultsList{end,1});
% else
%     [~,dateIdx] = find(cellfun(@(x) contains(x, options.resultsPathDate), resultsList(:,1)));
%     resultsPath = fullfile(resultsList{dateIdx,2},resultsList{dateIdx,1});
% end
% if strcmp(options.saveDataPath,'default')
%     options.saveDataPath = resultsPath;
% end

%% Generate summary file/figure for each search
for c = 1:size(cells,1)
    curCell = cells{c,'Cell'};
    
    % Loop through all summary files
    for search = 1:length(cells{curCell,"Epochs"}{1})
        % Get a list of depth for the current cell at current epoch
        curSearch = cells{curCell,'Epochs'}{1}{search};
        curCell = cells{curCell,'Cell'};
        cellResultsPath = fullfile(cells{curCell,'Session'},['cell',num2str(curCell)]);

        depthfilename = strcat(curSearch,'_depth*.mat');
        spotsList = dir(fullfile(cellResultsPath,depthfilename));
        disp(['Ongoing: plotting summary for search: ', curSearch]);

        % Load noise model of this neuron
        % if ~isfield(cells{curCell,'Options'}{1},'Ethres')
        %     load(strcat(cellResultsPath,filesep,'noise_cell',num2str(curCell),'.mat'));
        %     Ethres = 2 * noise_all.sigma;
        %     Ithres = 2 * noise_all.sigma;
        % else
        %     Ethres = cells{curCell,'Options'}{1}.Ethres;
        %     Ithres = cells{curCell,'Options'}{1}.Ithres;
        % end
        load(strcat(cellResultsPath,filesep,'noise_cell',num2str(curCell),'.mat'),...
                    'noise_all','allNullData');
        Ethres = -2 * noise_all.sigma;
        Ithres = 2 * noise_all.sigma;

        % Initialize figures
        initializeFig(0.8,1);
        masterLayout = tiledlayout(1,length(spotsList));
        masterLayout.TileSpacing = 'compact';
        masterLayout.Padding = 'compact';
    
        % Loop through depth
        for d = 1:length(spotsList)
            %% Initialization
            depthFilePath = fullfile(spotsList(d).folder,spotsList(d).name);
            load(depthFilePath,'spotsAtDepth');
            vhold = spotsAtDepth{1,'Vhold'};

            % Load time windows
            analysisWindow = cells{curCell,'Options'}{1}.analysisWindow;
            analysisWindowLength = cells{curCell,'Options'}{1}.analysisWindowLength;
            controlWindowLength = cells{curCell,'Options'}{1}.controlWindowLength;

            % Load color
            if vhold < -50; color = bluePurpleRed(end,:);
            elseif vhold > -10; color = bluePurpleRed(1,:);
            else; color = bluePurpleRed(250,:); 
            end

            % Get max/min responses
            maxResponse = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.max,1:height(spotsAtDepth), UniformOutput=false)');
            minResponse = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Stats'}{1}.response.min,1:height(spotsAtDepth), UniformOutput=false)');

            % Add depth tile
            nexttile(masterLayout,d); axis off;
            title(['Depth ', num2str(d)]); 
            depthLayout = tiledlayout(masterLayout,3,2);
            depthLayout.Layout.Tile = d;
            depthLayout.TileSpacing = 'tight'; depthLayout.Padding = 'tight';

            %% Plot opto vs prestim trace 
            nexttile(depthLayout,1); 
            optoData = cell2mat(arrayfun(@(row) spotsAtDepth{row,'Response'}{1}.processed(analysisWindow), 1:height(spotsAtDepth), UniformOutput=false)');
            ctrlData = cell2mat(arrayfun(@(row) spotsAtDepth{row,'Response'}{1}.control, 1:height(spotsAtDepth), UniformOutput=false)');
            optoTime = linspace(0,analysisWindowLength,size(optoData,2));
            ctrlTime = linspace(-controlWindowLength,0,size(ctrlData,2));

            plotSEM(optoTime,optoData,color,plotIndividual=true,individualColor='same',individualAlpha=0.25);
            plotSEM(ctrlTime,ctrlData,[.3,.3,.3],plotIndividual=true,individualColor='same',individualAlpha=0.25);
            xlabel('ms'); ylabel('pA'), legend({'Opto','Pre-stim'});


            %% Plot hotspot vs nullspot response
            nexttile(depthLayout,2); 
            hotspotIdx = cell2mat(arrayfun(@(x) spotsAtDepth{x,'Response'}{1}.hotspot,1:height(spotsAtDepth), UniformOutput=false)');
            hotspotData = optoData(hotspotIdx,:);
            nullspotData = optoData(~hotspotIdx,:);

            plotSEM(optoTime,hotspotData,color,plotIndividual=true,individualColor='same',individualAlpha=0.25);
            plotSEM(optoTime,nullspotData,[.3,.3,.3],plotIndividual=true,individualColor='same',individualAlpha=0.25);
            xlabel('ms'); ylabel('pA'), legend({'Hotspot','Null spot'});


            %% Bootstrap hotspot max/min response to match #data in allNullData
            if any(hotspotIdx)
                maxHotpot = maxResponse(hotspotIdx); minHotspot = minResponse(hotspotIdx);
                boot_maxHotspot = maxHotpot(randi(length(maxHotpot),[length(allNullData) 1]));
                boot_minHotspot = minHotspot(randi(length(minHotspot),[length(allNullData) 1]));
            else
                boot_maxHotspot = nan; boot_minHotspot = nan;
            end
            if any(~hotspotIdx)
                maxNullspot = maxResponse(~hotspotIdx); minNullspot = minResponse(~hotspotIdx);
                boot_maxNullspot = maxNullspot(randi(length(maxNullspot),[length(allNullData) 1]));
                boot_minNullspot = minNullspot(randi(length(minNullspot),[length(allNullData) 1]));
            else
                boot_maxNullspot = nan; boot_minNullspot = nan;
            end

            %% Plot noise vs hotspot distribution
            nexttile(depthLayout,3); 
            h_baseline = histogram(allNullData,'Normalization','pdf'); hold on
            h_maxHotspot = histogram(boot_maxHotspot,'Normalization','pdf'); hold on
            h_minHotspot = histogram(boot_minHotspot,'Normalization','pdf'); hold on
            h_baseline.FaceColor = [.8,.8,.8]; h_baseline.EdgeColor = [.8,.8,.8];
            h_maxHotspot.FaceColor = bluePurpleRed(1,:); h_maxHotspot.EdgeColor = bluePurpleRed(1,:);
            h_minHotspot.FaceColor = bluePurpleRed(end,:); h_minHotspot.EdgeColor = bluePurpleRed(end,:);
            xline(Ethres,'--',color=[.5,.5,.5],Label='Exci. threshold'); hold on
            xline(Ithres,'--',color=[.5,.5,.5],Label='Inhi. threshold');

            %% Plot noise vs hotspot distribution
            nexttile(depthLayout,4); 
            h_baseline = histogram(allNullData,'Normalization','pdf'); hold on
            h_maxNullspot = histogram(boot_maxNullspot,'Normalization','pdf'); hold on
            h_minNullspot = histogram(boot_minNullspot,'Normalization','pdf'); hold on
            h_baseline.FaceColor = [.8,.8,.8]; h_baseline.EdgeColor = [.8,.8,.8];
            h_baseline.FaceAlpha = 0.25; h_baseline.EdgeAlpha = 0.25;
            h_maxNullspot.FaceColor = bluePurpleRed(1,:); h_maxNullspot.EdgeColor = bluePurpleRed(1,:);
            h_maxNullspot.FaceAlpha = 0.25; h_maxNullspot.EdgeAlpha = 0.25;
            h_minNullspot.FaceColor = bluePurpleRed(end,:); h_minNullspot.EdgeColor = bluePurpleRed(end,:);
            h_minNullspot.FaceAlpha = 0.25; h_minNullspot.EdgeAlpha = 0.25;
            xline(Ethres,'--',color=[.5,.5,.5],Label='Exci. threshold'); hold on
            xline(Ithres,'--',color=[.5,.5,.5],Label='Inhi. threshold');

            %% Calculate response rate

            %% Plot hotspot excitatory response rate vs noise response rate
             

            %% Plot hotspot inhibitory response rate vs noise response rate


            %% Save success rate results if necessary


            % Add data of each spot
            for spot = 1:size(spotsAtDepth,1)
                vhold(spot) = spotsAtDepth{spot,'Vhold'};
                depth(spot) = spotsAtDepth{spot,'Depth'};
                repetition(spot) = spotsAtDepth{spot,'Repetition'};
                location(spot,:) = spotsAtDepth{spot,'Location'}{1};
                hotspot(spot) = spotsAtDepth{spot,'Response'}{1}.hotspot;

                responses(spot,:) = spotsAtDepth{spot,'Response'}{1}.processed(analysisWindow);
                baselines(spot,:) = spotsAtDepth{spot,'Response'}{1}.control;

                maxResponse(spot) = spotsAtDepth{spot,'Stats'}{1}.response.max;
                minResponse(spot) = spotsAtDepth{spot,'Stats'}{1}.response.min;
                maxBaseline(spot) = spotsAtDepth{spot,'Stats'}{1}.baseline.max;
                minBaseline(spot) = spotsAtDepth{spot,'Stats'}{1}.baseline.min;
                EI_response(spot) = spotsAtDepth{spot,'Stats'}{1}.response.EIindex;
                EI_baseline(spot) = spotsAtDepth{spot,'Stats'}{1}.baseline.EIindex;

                % Determine whether theres a success response based on
                % noise model of the neuron
                E_response(spot) = sum(minResponse(spot) <= Ethres);
                E_baseline(spot) = sum(minResponse(spot) <= Ethres);
                I_response(spot) = sum(maxResponse(spot) >= Ithres);
                I_baseline(spot) = sum(maxResponse(spot) >= Ithres);
            end
        end
    end
end

end