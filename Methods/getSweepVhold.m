function [vhold, meta] = getSweepVhold(epochPath, sweepName, baselineAvg, vholdChannel, options)
% GETSWEEPVHOLD  Extract holding voltage for a single sweep.
%
%   [vhold, meta] = getSweepVhold(epochPath, sweepName, baselineAvg, vholdChannel)
%   tries, in order:
%     1) Excel lookup (when vholdChannel == 'excel' and options.infoTable provided)
%     2) Load per-sweep AD channel file (e.g. AD2_123.mat) and run extractHoldingVoltage
%     3) Heuristic fallback based on baselineAvg
%
%   meta fields (always present):
%     meta.acqNum   : acquisition number parsed from sweepName (e.g., AD0_123 -> 123)
%     meta.method   : 'excel' | 'adfile' | 'heuristic'
%     meta.trace    : vhold trace vector if loaded from an AD channel file, else []
%     meta.value    : same as vhold (scalar)

    arguments
        epochPath
        sweepName
        baselineAvg
        vholdChannel
        options.infoTable = []
    end

    vhold = [];
    meta = struct('acqNum', NaN, 'method', "", 'trace', [], 'value', NaN);

    infoTable = options.infoTable;

    % Parse acquisition number from sweep name (e.g., AD0_15 -> 15)
    acqSplit = split(string(sweepName), '_');
    acqNum = str2double(acqSplit{end});
    meta.acqNum = acqNum;

    % Strategy 1: Excel Lookup
    if strcmpi(string(vholdChannel), 'excel')
        if ~isempty(infoTable)
            try
                idx = find(infoTable.acq_ == acqNum, 1);
                if ~isempty(idx)
                    vhold = infoTable.holding(idx);
                    meta.method = "excel";
                    meta.value = vhold;
                    return
                end
            catch
                % Keep vhold empty if lookup fails
            end
        end

    % Strategy 2: Load AD channel file (e.g., AD2_15.mat)
    else
        try
            vChanName = string(vholdChannel) + "_" + num2str(acqNum);
            vChanPath = fullfile(epochPath, vChanName + ".mat");

            if exist(vChanPath, 'file')
                S = load(vChanPath);

                if isfield(S, vChanName)
                    rawObj = S.(vChanName);

                    % Extract trace (if available) for downstream bookkeeping
                    if isstruct(rawObj) && isfield(rawObj, 'data')
                        meta.trace = rawObj.data;
                    elseif isnumeric(rawObj)
                        meta.trace = rawObj;
                    else
                        meta.trace = [];
                    end

                    % Assumes extractHoldingVoltage is on your path.
                    % (It can accept either the whole struct or just the data,
                    % depending on your implementation.)
                    try
                        vhold = extractHoldingVoltage(rawObj);
                    catch
                        if ~isempty(meta.trace)
                            vhold = extractHoldingVoltage(meta.trace);
                        else
                            rethrow(lasterror);
                        end
                    end

                    meta.method = "adfile";
                    meta.value = vhold;
                    return
                end
            end
        catch
            % Fall through to Excel (if provided) then heuristic
        end

        % If AD file load failed, try Excel as a fallback (if provided)
        if ~isempty(infoTable)
            try
                idx = find(infoTable.acq_ == acqNum, 1);
                if ~isempty(idx)
                    vhold = infoTable.holding(idx);
                    meta.method = "excel";
                    meta.value = vhold;
                    return
                end
            catch
                % keep empty
            end
        end
    end

    % Strategy 3: Heuristic fallback (based on baseline stats)
    if isempty(vhold) || isnan(vhold)
        if baselineAvg < 0
            vhold = -70;
        else
            vhold = 10;
        end
        meta.method = "heuristic";
        meta.value = vhold;
        meta.trace = [];
    end
end
