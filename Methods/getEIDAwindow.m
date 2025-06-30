function results = getEIDAwindow(animalEI, DAtrend, nCells, options)

arguments
    animalEI double
    DAtrend struct
    nCells double
    
    options.type string = 'sweep' % or 'moving'
    options.stats string = 'coefficient' % or 'slope'
    options.align string = 'end' % or 'start' (use when type is 'moving')

    options.index double

    options.trialWindowLength double = 10:80
    options.nWindow double = 80
    options.sweepWindow double = 20;
end

%% Set params

% Select index to analyze if provided
if isfield(options,'index')
    animalEI = animalEI(options.index);
    DAtrend = DAtrend(options.index);
    nCells = nCells(options.index);
end

% Initialize params
if contains(options.type, 'sweep')
    trialWindowLength = options.trialWindowLength;
    if contains(options.stats,'slope')
        ending_slopes = zeros(length(trialWindowLength),1);
        ending_pvals = zeros(length(trialWindowLength),1);
        starting_slopes = zeros(length(trialWindowLength),1);
        starting_pvals = zeros(length(trialWindowLength),1);
        ending_vals = zeros(length(trialWindowLength),1);
        ending_vals_pvals = zeros(length(trialWindowLength),1);
        starting_vals = zeros(length(trialWindowLength),1);
        starting_vals_pvals = zeros(length(trialWindowLength),1);
    elseif contains(options.stats,'coeff')
        ending_spearman = zeros(length(trialWindowLength),1);
        ending_spearman_p = zeros(length(trialWindowLength),1);
        starting_spearman = zeros(length(trialWindowLength),1);
        starting_spearman_p = zeros(length(trialWindowLength),1);
        ending_spearman_val = zeros(length(trialWindowLength),1);
        ending_spearman_val_p = zeros(length(trialWindowLength),1);
        starting_spearman_val = zeros(length(trialWindowLength),1);
        starting_spearman_val_p = zeros(length(trialWindowLength),1);
    end
elseif contains(options.type,'mov')
    nWindows  = options.nWindow - 1;
    sweepWindow = options.sweepWindow;
    nWindows  = nWindows-sweepWindow+1;
    sweeping_slopes = zeros(length(nWindows),1);
    sweeping_pvals = zeros(length(nWindows),1);
    sweeping_spearman = zeros(length(nWindows),1);
    sweeping_spearman_p = zeros(length(nWindows),1);
    sweeping_vals = zeros(length(nWindows),1);
    sweeping_vals_p = zeros(length(nWindows),1);
    sweeping_vals_spearman = zeros(length(nWindows),1);
    sweeping_vals_spearman_p = zeros(length(nWindows),1);
end

%% Sweep through windows

if contains(options.type,'sweep')
    for i = 1:length(trialWindowLength)
        % Get slopes
        trialWindow = trialWindowLength(i);
        endingSlope = zeros(length(DAtrend),1); 
        startingSlope = zeros(length(DAtrend),1);
        endingVal = zeros(length(DAtrend),1); 
        startingVal = zeros(length(DAtrend),1);

        for a = 1:length(DAtrend)
            % Get last n trial slope
            startTrial = max(1,DAtrend(a).nTrials-trialWindow);
            endingSlope(a) = DAtrend(a).amp.slopeMap.raw.map(startTrial,end);
            % Get first n trial slope
            endingTrial = min(DAtrend(a).nTrials,trialWindow);
            startingSlope(a) = DAtrend(a).amp.slopeMap.raw.map(1,endingTrial);
    
            % Get last n trial slope
            startTrial = max(1,DAtrend(a).nTrials-trialWindow);
            endingVal(a) = mean(DAtrend(a).amp.raw(startTrial:end));
            % Get first n trial slope
            endingTrial = min(DAtrend(a).nTrials,trialWindow);
            startingVal(a) = mean(DAtrend(a).amp.raw(1:endingTrial));
        end
    
        % Calculate last n trial fitted slopes
        Y = endingSlope;
        if contains(options.stats,'slope')
            [model,ending_pvals(i),~] = fitScatter(animalEI,Y,weights=nCells);
            ending_slopes(i) = model(1); 
        else
            [ending_spearman(i),ending_spearman_p(i)] = corr(animalEI(:), Y(:), 'Type', 'Spearman');
        end
    
        % Calculate first n trial fitted slopes
        Y = startingSlope;
        if contains(options.stats,'slope')
            [model,starting_pvals(i),~] = fitScatter(animalEI,Y,weights=nCells);
            starting_slopes(i) = model(1);
        else
            [starting_spearman(i),starting_spearman_p(i)] = corr(animalEI(:), Y(:), 'Type', 'Spearman');
        end
    
        % Calculate last n trial fitted slopes
        Y = endingVal;
        if contains(options.stats,'slope')
            [model,ending_vals_pvals(i),~] = fitScatter(animalEI,Y,weights=nCells);
            ending_vals(i) = model(1);
        else
            [ending_spearman_val(i),ending_spearman_val_p(i)] = corr(animalEI(:), Y(:), 'Type', 'Spearman');
        end
    
        % Calculate first n trial fitted slopes
        Y = startingVal;
        if contains(options.stats,'slope')
            [model,starting_vals_pvals(i),~] = fitScatter(animalEI,Y,weights=nCells);
            starting_vals(i) = model(1);
        else
            [starting_spearman_val(i),starting_spearman_val_p(i)] = corr(animalEI(:), Y(:), 'Type', 'Spearman');
        end
    end

    % Decide data to return
    if contains(options.stats,'slope')
        results.ending_slopes = ending_slopes;
        results.ending_pvals = ending_pvals;
        results.ending_slopes_maxWindow = getMaxCorrWindow(ending_slopes,ending_pvals,trialWindowLength(1));
        results.starting_slopes = starting_slopes;
        results.starting_pvals = starting_pvals;
        results.starting_slopes_maxWindow = getMaxCorrWindow(starting_slopes,starting_pvals,trialWindowLength(1));
        results.ending_vals = ending_vals;
        results.ending_vals_pvals = ending_vals_pvals;
        results.ending_vals_maxWindow = getMaxCorrWindow(ending_vals,ending_vals_pvals,trialWindowLength(1));
        results.starting_vals = starting_vals;
        results.starting_vals_pvals = starting_vals_pvals;
        results.starting_vals_maxWindow = getMaxCorrWindow(starting_vals,starting_vals_pvals,trialWindowLength(1));
    elseif contains(options.stats,'coeff')
        results.ending_slopes = ending_spearman;
        results.ending_pvals = ending_spearman_p;
        results.ending_slopes_maxWindow = getMaxCorrWindow(ending_spearman,ending_spearman_p,trialWindowLength(1));
        results.starting_slopes = starting_spearman;
        results.starting_pvals = starting_spearman_p;
        results.starting_slopes_maxWindow = getMaxCorrWindow(starting_spearman,starting_spearman_p,trialWindowLength(1));
        results.ending_vals = ending_spearman_val;
        results.ending_vals_pvals = ending_spearman_val_p;
        results.ending_vals_maxWindow = getMaxCorrWindow(ending_spearman_val,ending_spearman_val_p,trialWindowLength(1));
        results.starting_vals = starting_spearman_val;
        results.starting_vals_pvals = starting_spearman_val_p;
        results.starting_vals_maxWindow = getMaxCorrWindow(starting_spearman_val,starting_spearman_val_p,trialWindowLength(1));
    end


%% Moving window
elseif contains(options.type,'mov')
    if strcmpi(options.align,'start')
        % Plot moving 20 trial window (aligned to start)
        nAnimals  = numel(DAtrend);
        movSlopes = nan(nAnimals, nWindows-sweepWindow+1);
        movAvgs   = nan(nAnimals, nWindows-sweepWindow+1);
        % precompute your slope‐kernel once
        xc = (1:sweepWindow)' - mean(1:sweepWindow);
        b  = flipud(xc)/sum(xc.^2);
        for a = 1:nAnimals
            raw = DAtrend(a).amp.smoothed(:);
            % truncate or pad to exactly nWindow samples
            data = raw(1 : min(end,nWindows));
            data(end+1:nWindows,1) = nan;
            
            % compute on the flipped data
            revSlopes = conv(data, b, 'valid'); 
            csum      = cumsum([0; data]);
            revMeans  = (csum(sweepWindow+1:end) - csum(1:end-sweepWindow))/sweepWindow;
            
            % flip the results back
            movSlopes(a,:) = -revSlopes;
            movAvgs(a,:)   = revMeans;
        end
        
        for w = 1:nWindows
            % Calculate fitted slopes with animal EI
            Y = movSlopes(:,w);
            [model,sweeping_pvals(w),~] = fitScatter(X,Y,weights=nCells);
            sweeping_slopes(w) = model(1);
            [X_clean,Y_clean] = removeNaNs(X,Y);
            [sweeping_spearman(w),sweeping_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
            % Calculate average val with animal EI
            Y = movAvgs(:,w);
            [model,sweeping_vals_p(w),~] = fitScatter(X,Y,weights=nCells);
            sweeping_vals(w) = model(1);
            [X_clean,Y_clean] = removeNaNs(X,Y);
            [sweeping_vals_spearman(w),sweeping_vals_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
        end

    elseif strcmpi(options.align,'end')
        % Plot moving 20 trial window (aligned to end)
        nAnimals  = numel(DAtrend);
        movSlopes = nan(nAnimals, nWindows-sweepWindow+1);
        movAvgs   = nan(nAnimals, nWindows-sweepWindow+1);
        % precompute your slope‐kernel once
        xc = (1:sweepWindow)' - mean(1:sweepWindow);
        b  = flipud(xc)/sum(xc.^2);
        for a = 1:nAnimals
            raw = DAtrend(a).amp.smoothed(:);
            L   = min(numel(raw), nWindows);
            data = [nan(nWindows-L,1);raw(1:L)];
            
            % flip the vector end-for-start
            rev = data(end:-1:1);
            
            % compute on the flipped data
            revSlopes = conv(rev, b, 'valid'); 
            csum      = cumsum([0; rev]);
            revMeans  = (csum(sweepWindow+1:end) - csum(1:end-sweepWindow))/sweepWindow;
            
            % flip the results back
            movSlopes(a,:) = -revSlopes;
            movAvgs(a,:)   = revMeans;
        end
        
        for w = 1:nWindows
            % Calculate fitted slopes with animal EI
            Y = movSlopes(:,w);
            [model,sweeping_pvals(w),~] = fitScatter(X,Y,weights=nCells);
            sweeping_slopes(w) = model(1);
            [X_clean,Y_clean] = removeNaNs(X,Y);
            [sweeping_spearman(w),sweeping_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
            % Calculate average val with animal EI
            Y = movAvgs(:,w);
            [model,sweeping_vals_p(w),~] = fitScatter(X,Y,weights=nCells);
            sweeping_vals(w) = model(1);
            [X_clean,Y_clean] = removeNaNs(X,Y);
            [sweeping_vals_spearman(w),sweeping_vals_spearman_p(w)] = corr(X_clean(:), Y_clean(:), 'Type', 'Spearman');
        end
    end

    results.slopes = sweeping_slopes;
    results.slope_pvals = sweeping_pvals;
    results.coeff = sweeping_spearman;
    results.coeff_pvals = sweeping_spearman_p;
    results.vals = sweeping_vals;
    results.vals_pvals = sweeping_vals_p;
    results.vals_coeff = sweeping_vals_spearman;
    results.vals_coeff_pvals = sweeping_vals_spearman_p;
end

end
