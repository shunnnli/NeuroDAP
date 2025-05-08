function [EPSC_traces, IPSC_traces] = getResponseTraces(combined_cells,options)

arguments
    combined_cells table

    options.plot logical = false
    options.plotIndividual logical = true
    options.color
    options.ylim double
    options.newFig logical = true % initialize a new figure or not
    options.overlay logical = false % plot all cells overlaying in same panel or not
    options.plotLegend logical = true
    options.plotNormalized logical = false
    
    options.timeRange double = [-10,50]
    options.Fs double = 10000;

    options.animalRange string = 'All'
    options.cellList double
end

timeRange = options.timeRange; 
time2sample = options.Fs/1000; 
timeRangeSamples = timeRange * time2sample;
timeWindowLength = timeRangeSamples(2) - timeRangeSamples(1) + 1;
timeRangeInms = linspace(-10,50,timeWindowLength);


% Define cell to plot
if ~isfield(options,'cellList')
    cellList = 1:size(combined_cells,1);
    if ~strcmpi(options.animalRange,'All')
        cellList = cellList(contains(combined_cells.Animal,options.animalRange));
    end
else
    cellList = options.cellList;
end

EPSC_traces = []; IPSC_traces = [];

if options.newFig
    initializeFig(1,1); tiledlayout('flow');
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
        if ~isfield(options,'color')
            % [~,~,~,~,~,~,bluePurpleRed] = loadColors;
            EPSCcolor = [255 157 33]/255;
            IPSCcolor = [71 144 253]/255;
        else
            if ~iscell(options.color)
                EPSCcolor = options.color;
                IPSCcolor = options.color;
            elseif numel(options.color) == 1
                EPSCcolor = options.color{1};
                IPSCcolor = options.color{1};
            elseif numel(options.color) == 2
                EPSCcolor = options.color{1};
                IPSCcolor = options.color{2};
            end
        end

        if ~options.overlay; nexttile; end
        if options.plotNormalized
            curEPSC_peaks = curCell.Stats{1}.summary.EPSC.peakAvg;
            curIPSC_peaks = curCell.Stats{1}.summary.IPSC.peakAvg;
            totalCurrents = abs(curEPSC_peaks)+abs(curIPSC_peaks);
            curTrace_EPSC = curTrace_EPSC / totalCurrents;
            curTrace_IPSC = curTrace_IPSC / totalCurrents;
        end
        plotSEM(timeRangeInms,curTrace_EPSC,EPSCcolor,plotIndividual=options.plotIndividual,label='EPSC');
        plotSEM(timeRangeInms,curTrace_IPSC,IPSCcolor,plotIndividual=options.plotIndividual,label='IPSC');
         
        if options.plotLegend; legend; end
        if isfield(options,'ylim'); ylim(options.ylim); end
        xlabel('Time (ms)'); ylabel('Current (pA)');
        title(strcat(curCell.Task,": Cell ",num2str(curCell.Cell)," (",curCell.Animal,")"));
    end
end

end