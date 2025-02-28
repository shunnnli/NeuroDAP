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
    % Common options
    options.signalFs double = 50
    options.baseline double = [0,0] % subtracted by x secs before event time
    options.dff logical = false

    % Extract options
    options.extract logical = true    
    options.eventIdx double
    options.timeRange double
    options.signal doube
    options.params struct
    options.sameSystem logical = false
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
    options.individualColor = 'same'
    options.individualAlpha double = 0.3 %1 is not transparent, 0 is fully transparent

    options.plotStyle string = 'line'
    options.plotPatch logical = true

    options.smooth double = 0 % 0: no smoothing, else is samples of smooth data
    options.smoothMethod string = 'movmean';

    options.print logical = false % print progress message
    options.ylabel string = 'z-score'
    options.xlabel string = 'Time (s)'
end

%% Parse required inputs

% Determine number of inputs (option inputs started as a string)
nInputs = find(cellfun(@isstring,varargin),1) - 1;

if isempty(nInputs); nInputs = nargin;
else; error('Input format error: check whether options is correct'); end

if nInputs == 4
    if ~options.extract
        disp('plotTraces: extract is false but 4 inputs are provided, changed to true!'); 
        options.extract = true;
    end
    options.eventIdx = varargin{1};
    % if size(options.eventIdx,1) == 1; options.eventIdx = options.eventIdx'; end
    options.timeRange = varargin{2};
    options.signal = varargin{3};
    options.params = varargin{4};
    if options.print
        if options.plot
            disp(['plotTraces: ',num2str(nInputs),' inputs are detected, execute extract and plot']);
        else
            disp(['plotTraces: ',num2str(nInputs),' inputs are detected, execute extract']);
        end
    end
elseif nInputs == 5
    if ~options.extract
        if options.print
            disp('plotTraces: extract is false but 4 inputs are provided, changed to true!'); 
        end
        options.extract = true;
    end
    options.eventIdx = varargin{1};
    % if size(options.eventIdx,1) == 1; options.eventIdx = options.eventIdx'; end
    options.timeRange = varargin{2};
    options.signal = varargin{3};
    options.params = varargin{5};
    if isfield(options,'color'); options.color = varargin{4}; end
    if options.print
        if options.plot
            disp(['plotTraces: ',num2str(nInputs),' inputs are detected, execute extract and plot']);
        else
            disp(['plotTraces: ',num2str(nInputs),' inputs are detected, execute extract']);
        end
    end
elseif nInputs <= 2
    if ~options.plot
        if options.print
            disp('plotTraces: plot is false but 2 inputs are provided, changed to true!'); 
        end
        options.plot = true;
    end
    if options.extract
        if options.print
            disp('plotTraces: extract is true but 2 inputs are provided, changed to false!');
        end
        options.extract = false;
    end
    options.traces = varargin{1};
    if nInputs == 2
        options.timestamp = varargin{2};
    else
        options.timestamp = 1:size(options.traces,2);
    end
    if options.print
        disp(['plotTraces: ',num2str(nInputs),' inputs are detected, execute plot']);
    end
elseif nInputs == 3
    if ~options.extract
        if options.print
            disp('plotTraces: extract is false but 3 inputs are provided, changed to true!'); 
        end
        options.extract = true;
    end
    options.eventIdx = varargin{1};
    % if size(options.eventIdx,1) == 1; options.eventIdx = options.eventIdx'; end
    options.timeRange = varargin{2};
    options.signal = varargin{3};
    options.params = struct([]);
    if options.print
        if options.plot
            disp(['plotTraces: ',num2str(nInputs),' inputs are detected, execute extract and plot']);
        else
            disp(['plotTraces: ',num2str(nInputs),' inputs are detected, execute extract']);
        end
    end
else
    error(['Input format error: ', num2str(nInputs),' inputs are detected, should be either 2, 4, 5']);
end

%% Parse options inputs

% optionIdx = find(cellfun(@isstring,varargin));
% 
% for i = 1:length(optionIdx)
%     if ~isfield(options,varargin{optionIdx(i)})
%         error(['option input error: can not find corresponding options for ',convertStringsToChars(varargin{optionIdx(i)})]);
%     end
%     options.(varargin{optionIdx(i)}) = varargin{optionIdx(i)+1};
% end

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
                                        options.signal,params=options.params,...
                                        signalSystem=options.signalSystem,...
                                        signalFs=options.signalFs,...
                                        eventSystem=options.eventSystem,...
                                        sameSystem=options.sameSystem,...
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
    lastBin_baseline = floor(eventBin + min(options.baseline(2),0)*options.signalFs);
    baseline = mean(traces(:,firstBin_baseline:lastBin_baseline),2);

    % Subtract
    traces = traces - baseline;
end

%% Turn into df/f (if input signal is raw)
if options.dff
    % Check input
    if options.baseline(1) < timestamp(1)
        options.baseline(1) = timestamp(1);
        warning('plotTraces: baseline is longer than input timeRange, changed to the first timestamp instead');
        disp(['baseline: ',num2str(options.baseline)]);
    elseif options.baseline(1) == 0
        options.baseline(1) = timestamp(1);
    end

    % Convert to bin
    [~,eventBin] = min(abs(timestamp-0));
    firstBin_baseline = floor(eventBin + options.baseline(1)*options.signalFs);
    lastBin_baseline = floor(eventBin + min(options.baseline(2),0)*options.signalFs);
    baseline = mean(traces(:,firstBin_baseline:lastBin_baseline),2);

    % Convert to dff
    traces = (traces-baseline)./baseline;
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
        % shuffled = shuffleTraces(traces);
        shuffledEvents = randi([1, length(timestamp)], length(options.eventIdx), 1);
        [shuffled,~] = getTraces(shuffledEvents,options.timeRange,...
                                options.signal,params=options.params,...
                                signalSystem=options.signalSystem,...
                                signalFs=options.signalFs,...
                                eventSystem=options.eventSystem,...
                                sameSystem=options.sameSystem,...
                                rmmissing=options.rmmissing,...
                                baseline=options.baseline);
        plotSEM(timestamp,shuffled,options.shuffledColor,smooth=smoothWindow,smoothMethod=options.smoothMethod,...
        LineStyle=options.LineStyle,LineWidth=options.LineWidth,...
        plotStyle=options.plotStyle,plotPatch=options.plotPatch); hold on
    end

    % 4.3 Plot SEM
    plotSEM(timestamp,traces,options.color,smooth=smoothWindow,smoothMethod=options.smoothMethod,...
        LineStyle=options.LineStyle,LineWidth=options.LineWidth,...
        plotIndividual=options.plotIndividual,individualColor=options.individualColor,...
        plotStyle=options.plotStyle,plotPatch=options.plotPatch);
    
    xlabel(options.xlabel); 
    ylabel(options.ylabel); 
    
    hold on
end

end

