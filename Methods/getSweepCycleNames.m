function cycles = getSweepCycleNames(protocolData,nSweeps)

% Return one cycle name per sweep from either supported Protocol format.
% loadSlices stores search protocols as a cell array of scalar structs and
% merged full-field protocols as a scalar struct whose cycle field is a
% cell array. This helper keeps downstream code independent of that storage
% detail.

if nargin < 2
    nSweeps = inferSweepCount(protocolData);
end

cycles = strings(nSweeps,1);

if iscell(protocolData)
    nProtocol = min(nSweeps,numel(protocolData));
    for sweep = 1:nProtocol
        protocol = protocolData{sweep};
        if isstruct(protocol) && isfield(protocol,'cycle') && ~isempty(protocol.cycle)
            cycleValue = string(protocol.cycle);
            cycles(sweep) = cycleValue(1);
        end
    end
elseif isstruct(protocolData) && ~isempty(protocolData)
    if ~isscalar(protocolData)
        nProtocol = min(nSweeps,numel(protocolData));
        for sweep = 1:nProtocol
            if isfield(protocolData(sweep),'cycle') && ~isempty(protocolData(sweep).cycle)
                cycleValue = string(protocolData(sweep).cycle);
                cycles(sweep) = cycleValue(1);
            end
        end
    elseif isfield(protocolData,'cycle')
        cycleValues = string(protocolData.cycle);
        cycleValues = cycleValues(:);
        nProtocol = min(nSweeps,numel(cycleValues));
        cycles(1:nProtocol) = cycleValues(1:nProtocol);
    end
end

end


function nSweeps = inferSweepCount(protocolData)

if iscell(protocolData) || (isstruct(protocolData) && ~isscalar(protocolData))
    nSweeps = numel(protocolData);
elseif isstruct(protocolData) && isscalar(protocolData) && isfield(protocolData,'cycle')
    cycleValue = protocolData.cycle;
    if ischar(cycleValue) || (isstring(cycleValue) && isscalar(cycleValue))
        nSweeps = 1;
    else
        nSweeps = numel(cycleValue);
    end
else
    nSweeps = 0;
end

end
