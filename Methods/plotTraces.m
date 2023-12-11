function [traces,timestamp] = plotTraces(varargin,options)

%% Notes
% plotTraces.m, Shun Li
% 2023/09/02: change to generalize to event or signal streams of choice
% 2023/09/05: fixed minor error to make it applicable for camera and
% downsampled LJ photometry

%% Function
arguments (Repeating)
    varargin
end
arguments
    % eventIdx double
    % timeRange double
    % signal double
    % color
    % params struct
    
    % Common options
    options.signalFs double = 50
    options.baseline double = 0 % subtracted by x secs before event time

    % Extract options
    options.extract logical = true    
    options.eventIdx double
    options.timeRange double
    options.signal doube
    options.params struct
    options.eventSystem string = 'ni'
    options.signalSystem string = 'lj'
    options.rmmissing logical = true % remove nan rows

    % Plot options
    options.plot logical = true % whether or not to plot traces
    options.traces double
    options.timestamp double

    options.color = [0.75, 0.75, 0.75]
    options.LineStyle (1,1) string = "-"
    options.LineWidth (1,1) {mustBeNumeric} = 2

    options.plotShuffled logical = false
    options.shuffledColor double = [.75, .75, .75]
    
    options.plotIndividual logical = false
    options.individualColor double = [.75, .75, .75]

    options.plotStyle string = 'line'

    options.smooth double = 0; % 0: no smoothing, else is samples of smooth data
    options.smoothMethod string = 'movmean';
end

%% Check inputs
if nargin == 4
    if ~options.extract
        disp('plotTraces: extract is false but 4 inputs are provided, changed to true!'); 
        options.extract = true;
    end
    options.eventIdx = varargin{1};
    options.timeRange = varargin{2};
    options.signal = varargin{3};
    options.params = varargin{4};
    if options.plot
        disp(['plotTraces: ',num2str(nargin),' inputs are detected, execute extract and plot']);
    else
        disp(['plotTraces: ',num2str(nargin),' inputs are detected, execute extract']);
    end
elseif nargin == 5
    if ~options.extract
        disp('plotTraces: extract is false but 4 inputs are provided, changed to true!'); 
        options.extract = true;
    end
    options.eventIdx = varargin{1};
    options.timeRange = varargin{2};
    options.signal = varargin{3};
    options.params = varargin{5};
    if isfield(options,'color'); options.color = varargin{4}; end
    if options.plot
        disp(['plotTraces: ',num2str(nargin),' inputs are detected, execute extract and plot']);
    else
        disp(['plotTraces: ',num2str(nargin),' inputs are detected, execute extract']);
    end
elseif nargin == 2
    if ~options.plot
        disp('plotTraces: plot is false but 2 inputs are provided, changed to true!'); 
        options.plot = true;
    end
    if options.extract
        disp('plotTraces: extract is true but 2 inputs are provided, changed to false!');
        options.extract = false;
    end
    options.traces = varargin{1};
    options.timestamp = varargin{2};
    disp(['plotTraces: ',num2str(nargin),' inputs are detected, execute plot']);
end

%% Check extract & plot
if ~options.extract && ~options.plot
    if isfield(options,'traces'); options.plot = true; end
    if isfield(options,'params'); options.extract = true; end
    warning('plotTraces: extract and plot can not be all false, changed to the following:');
    disp(['options.extract: ',num2str(options.extract)]);
    disp(['options.plot: ',num2str(options.plot)]);
elseif options.plot && ~options.extract
    traces = options.traces;
    timestamp = options.timestamp;
    if ~isfield(options,'traces'); error('plotTraces: options.traces not provided!'); end
    if ~isfield(options,'timestamp'); error('plotTraces: options.timestamp not provided!'); end
elseif options.extract
    if ~isfield(options,'eventIdx'); error('plotTraces: options.eventIdx not provided!'); end
    if ~isfield(options,'timeRange'); error('plotTraces: options.timeRange not provided!'); end
    if ~isfield(options,'signal'); error('plotTraces: options.signal not provided!'); end
    if ~isfield(options,'params'); error('plotTraces: options.params not provided!'); end
end


%% Get traces if needed
if options.extract
    [traces,timestamp] = getTraces(options.eventIdx,options.timeRange,...
                                        options.signal,options.params,...
                                        signalSystem=options.signalSystem,...
                                        signalFs=options.signalFs,...
                                        eventSystem=options.eventSystem,...
                                        rmmissing=options.rmmissing,...
                                        baseline=options.baseline);
end

%% Subtract baseline if needed

if ~options.extract && any(options.baseline ~= 0)
    % Check input
    if options.baseline(1) < timestamp(1)
        options.baseline(1) = timestamp(1);
        warning('plotTraces: baseline is longer than input timeRange, changed to the first timestamp instead');
        disp(['baseline: ',num2str(options.baseline)]);
    end

    % Convert to bin
    [~,eventBin] = min(abs(timestamp-0));
    firstBin_baseline = floor(eventBin + options.baseline(1)*options.signalFs);
    lastBin_baseline = floor(eventBin + options.baseline(2)*options.signalFs);
    baseline = mean(traces(:,firstBin_baseline:lastBin_baseline),2);
    traces = traces - baseline;
end

%% Plot traces if needed

if options.plot
    % 4.1 Smooth data if neccessary
    if options.smooth ~= 0
        if options.smooth == 1; smoothWindow = 0.1*signalFs;
        else; smoothWindow = options.smooth;
        end
    else
        smoothWindow = 0;
    end

    % 4.2 Plot shuffled data
    if options.plotShuffled
        shuffled = shuffleTraces(traces);
        plotSEM(timestamp,shuffled,options.shuffledColor,smooth=smoothWindow,smoothMethod=options.smoothMethod,...
        LineStyle=options.LineStyle,LineWidth=options.LineWidth,plotStyle=options.plotStyle); hold on
    end

    % 4.3 Plot SEM
    plotSEM(timestamp,traces,options.color,smooth=smoothWindow,smoothMethod=options.smoothMethod,...
        LineStyle=options.LineStyle,LineWidth=options.LineWidth,...
        plotIndividual=options.plotIndividual,individualColor=options.individualColor,...
        plotStyle=options.plotStyle);
    xlabel('Time (s)'); ylabel('z-score'); hold on
end

end

