function vhold = getSweepVhold(epochPath, sweepName, baselineAvg, options, infoTable)
% GETSWEEPVHOLD Extracts holding voltage with fallback heuristics
% Replaces manual extraction logic in loadSlicesDMD

    vhold = [];
    
    % Parse acquisition number from sweep name (e.g., AD0_15 -> 15)
    acqSplit = split(sweepName, '_'); 
    acqNum = str2double(acqSplit{end});

    % Strategy 1: Excel Lookup
    if strcmpi(options.vholdChannel, 'excel')
        if ~isempty(infoTable)
            try
                % Ensure we find a match
                idx = find(infoTable.acq_ == acqNum, 1);
                if ~isempty(idx)
                    vhold = infoTable.holding(idx);
                end
            catch
                % Keep vhold empty if lookup fails
            end
        end

    % Strategy 2: Load AD channel file (e.g., AD1_15.mat)
    else
        try
            vChanName = [options.vholdChannel, '_', num2str(acqNum)];
            vChanPath = fullfile(epochPath, [vChanName, '.mat']);
            
            if exist(vChanPath, 'file')
                S = load(vChanPath);
                % Handle variable loading without eval
                if isfield(S, vChanName)
                    rawTrace = S.(vChanName);
                    % Assumes extractHoldingVoltage is in your path
                    vhold = extractHoldingVoltage(rawTrace);
                end
            end
        catch
            if ~isempty(infoTable)
                try
                    % Ensure we find a match
                    idx = find(infoTable.acq_ == acqNum, 1);
                    if ~isempty(idx)
                        vhold = infoTable.holding(idx);
                    end
                catch
                    % Keep vhold empty if lookup fails
                end
            end
        end
    end

    % Strategy 3: Heuristic Fallback (based on baseline stats)
    if isempty(vhold) || isnan(vhold)
        if baselineAvg < 0
            vhold = -70;
        else
            vhold = 10; 
        end
    end
end