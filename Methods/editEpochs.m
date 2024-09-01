function epochs = editEpochs(epochs, options)

arguments
    epochs table
    
    % Can take either a string or a cell
    % If empty, default as editing included, processed sweeps, stats
    options.field

    options.save logical = true
    options.getCellTable logical = false

    % Edit include related params
    options.include
    options.resetOriginal logical = false

    % Edit processed sweeps related params


    % Edit stats related parms
end

%% Setup

if ~isfield(options,field)
    options.field = {'included','processed sweeps','stats'};
else
    if istring(options.field); options.field = {options.field}; end
end

% Include setup
if isfield(options,'include') && isstring(options.include)
    options.include = {options.include};
end

%% Edit selected fields

for f = 1:length(options.field)
    curField = options.field{f};

    % Check field names
    if ~isfield(epochs,curField)
        warning(['Error: ', curField, ' is not in epochs.mat, skip this field!']);
        continue
    end

    for e = 1:size(epochs,1)
        %% Edit included
        if contains(curField,'include',IgnoreCase=true)
            % Remove empty/erraneous sweeps
            if isfield(options,'include')
                % Initialize included if necessary
                if options.resetOriginal; included = ones(length(sweepAcq),1);
                else; included = epochs{e,'Included'}{1}; end 
        
                % Remove sweeps with different cycles
                cycles = epochs{e,'Protocol'}{1}.cycle;
                cyclesCount = tabulate(cycles);
                [~, mostCommonIdx] = max([cyclesCount{:, 2}]);
                if ~all(strcmp(cycles, cycles{mostCommonIdx}))
                    included = strcmp(cycles, cycles{mostCommonIdx});
                end
                
                % Remove sweeps based on provided criteria
                if isfield(options,'include') && ~isempty(options.include)
                    for criterion = options.include
                        included = all([included,eval(criterion{1})],2);
                    end
                    if length(included) ~= length(sweepAcq)
                        warning('included have different size than #sweeps! Reset included to true for all!');
                        included = ones(length(sweepAcq),1);
                    end
                end
                epochs{e,'Included'} = num2cell(included,[1 2]);
            end
    
        %% Edit processed sweeps
        elseif contains(curField,'process',IgnoreCase=true)
           
        %% Edit stats
        elseif contains(curField,'stats',IgnoreCase=true)
            % if ~contains(protocol.cycle,'randomSearch')
            %     % Calculate statistics
            %     % Find auc
            %     stats.response.auc = sum(processed_trace(analysisWindow)) / options.outputFs;
            %     stats.baseline.auc = sum(processed_trace(controlWindow)) / options.outputFs;
            % 
            %     % Find min and max value for stim response
            %     trace = processed_trace(analysisWindow);
            %     [~,maxIdx] = max(trace); [~,minIdx] = min(trace);
            %     % Average around max/min idx to get final value
            %     maxWindowStart = max(1,maxIdx-peakWindowWidth);
            %     maxWindowEnd = min(maxIdx+peakWindowWidth,length(trace));
            %     minWindowStart = max(1,minIdx-peakWindowWidth);
            %     minWindowEnd = min(minIdx+peakWindowWidth,length(trace));
            %     stats.response.max = mean(trace(maxWindowStart:maxWindowEnd));
            %     stats.response.min = mean(trace(minWindowStart:minWindowEnd));
            %     stats.response.maxTime = maxIdx * 1000/options.outputFs;
            %     stats.response.minTime = minIdx * 1000/options.outputFs;
            %     if vhold < -50; stats.response.peak = stats.response.min;
            %     elseif vhold >= 0; stats.response.peak = stats.response.max; 
            %     end
            % 
            %     % Find min and max for control baseline
            %     trace = processed_trace(controlWindow);
            %     [~,maxIdx] = max(trace); [~,minIdx] = min(trace);
            %     % Average around max/min idx to get final value
            %     maxWindowStart = max(1,maxIdx-peakWindowWidth);
            %     maxWindowEnd = min(maxIdx+peakWindowWidth,length(trace));
            %     minWindowStart = max(1,minIdx-peakWindowWidth);
            %     minWindowEnd = min(minIdx+peakWindowWidth,length(trace));
            %     stats.baseline.max = mean(trace(maxWindowStart:maxWindowEnd));
            %     stats.baseline.min = mean(trace(minWindowStart:minWindowEnd));
            %     stats.baseline.maxTime = maxIdx * 1000/options.outputFs;
            %     stats.baseline.minTime = minIdx * 1000/options.outputFs;
            %     if vhold < -50; stats.baseline.peak = stats.baseline.min;
            %     elseif vhold >= 0; stats.baseline.peak = stats.baseline.max; 
            %     end
            % 
            %     % Find E/I index
            %     stats.response.EIindex = abs(stats.response.max)-abs(stats.response.min) / abs(stats.response.max)+abs(stats.response.min);
            %     stats.baseline.EIindex = abs(stats.baseline.max)-abs(stats.baseline.min) / abs(stats.baseline.max)+abs(stats.baseline.min);
            % 
            %     statistics{k} = stats;
            %     epochs{row,'Stats'} = {statistics};
            % end
    
        else
            warning(['Error: ',curField,' not supported to be edited.']);
        end
    end
end

%% Save 
if options.save
    save(strcat(options.saveDataPath,filesep,'epochs_',expName),'epochs','-v7.3');
    disp(strcat("New epochs.mat created & saved: ",expName));
end

%% Get new cell table
if options.getCellTable
    cells = getCellTable(epochs,save=options.save,...
                         timeRange=options.timeRange,...
                         outputFs=options.outputFs,...
                         controlWindowLength=options.controlWindowLength,...
                         nArtifactSamples=options.nArtifactSamples,...
                         peakWindow=options.peakWindow);
end

end