function [EPSC_traces, IPSC_traces] = getResponseTraces(combined_cells,options)

arguments
    combined_cells table

    options.plot logical = false
    options.plotIndividual logical = true
    options.ylim double
    
    options.timeRange double = [-10,50]
    options.Fs double = 10000;

    options.animalRange string = 'All'
end

timeRange = options.timeRange; 
time2sample = options.Fs/1000; 
timeRangeSamples = timeRange * time2sample;
timeWindowLength = timeRangeSamples(2) - timeRangeSamples(1) + 1;
timeRangeInms = linspace(-10,50,timeWindowLength);

% Define cell to plot
cellList = 1:size(combined_cells,1);
if ~strcmpi(options.animalRange,'All')
    cellList = cellList(contains(combined_cells.Animal,options.animalRange));
end


EPSC_traces = []; IPSC_traces = [];

if options.plot
    initializeFig(1,1); tiledlayout('flow');
    [~,~,~,~,~,~,bluePurpleRed] = loadColors;
end

for c = 1:length(cellList)
    cellIdx = cellList(c);
    curCell = combined_cells(cellIdx,:); 

    % Get rows
    EPSCrows = curCell.Vhold{1} == -70;
    IPSCrows = curCell.Vhold{1} >= 0;

    % Get included
    if sum(EPSCrows)>0; EPSCincluded = curCell.Included{1}{EPSCrows};
    else; EPSCincluded = []; end
    if sum(IPSCrows)>0; IPSCincluded = curCell.Included{1}{IPSCrows};
    else; IPSCincluded = []; end
    if sum(EPSCincluded) == 0; EPSCincluded = ones(length(EPSCincluded),1); end
    if sum(IPSCincluded) == 0; IPSCincluded = ones(length(IPSCincluded),1); end
    EPSCincluded = logical(EPSCincluded);
    IPSCincluded = logical(IPSCincluded);

    % Concatenate traces
    if sum(EPSCrows) ~= 0
        EPSCIdx = find(EPSCrows); curTrace_EPSC = [];
        for epsc = EPSCIdx'
            stimOnset = curCell.Protocol{1}{epsc}.stimOnset(1);
            timeRangeWindow = stimOnset+timeRangeSamples(1) : stimOnset+timeRangeSamples(2);
            curTrace_EPSC = [curTrace_EPSC;curCell.('Processed sweeps'){1}{epsc}(:,timeRangeWindow)];
            curTrace_EPSC(~EPSCincluded,:) = [];
        end
        EPSC_traces = [EPSC_traces;mean(curTrace_EPSC,1)];
    else
        EPSC_traces = [EPSC_traces;nan(1,size(EPSC_traces,2))];
    end

    if sum(IPSCrows) ~= 0
        IPSCIdx = find(IPSCrows); curTrace_IPSC = [];
        for ipsc = IPSCIdx'
            stimOnset = curCell.Protocol{1}{ipsc}.stimOnset(1);
            timeRangeWindow = stimOnset+timeRangeSamples(1) : stimOnset+timeRangeSamples(2);
            curTrace_IPSC = [curTrace_IPSC;curCell.('Processed sweeps'){1}{ipsc}(:,timeRangeWindow)];
            curTrace_IPSC(~IPSCincluded,:) = [];
        end
        IPSC_traces = [IPSC_traces;mean(curTrace_IPSC,1)];
    else
        IPSC_traces = [IPSC_traces;nan(1,size(IPSC_traces,2))];
    end

    if options.plot
        nexttile;
        plotSEM(timeRangeInms,curTrace_EPSC,bluePurpleRed(end,:),plotIndividual=options.plotIndividual,label='EPSC');
        plotSEM(timeRangeInms,curTrace_IPSC,bluePurpleRed(1,:),plotIndividual=options.plotIndividual,label='IPSC');
        xlabel('Time (ms)'); ylabel('Current (pA)'); legend; 
        if isfield(options,'ylim'); ylim(options.ylim); end
        title(strcat(curCell.Task,": Cell ",num2str(curCell.Cell)," (",curCell.Animal,")"));
    end
end

end