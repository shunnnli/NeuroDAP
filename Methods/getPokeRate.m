function poke_rate = getPokeRate(data,options)

arguments
    data table

    options.method string = 'moving window'
    options.binWindow double = 10 % in min
    options.smoothWindow double = 3 % in bins
end

% Extract poking data from table
event = data.Event;
data.timestamp = data.MM_DD_YYYYHh_mm_ss;
timeSinceStart = minutes(data.timestamp - data.timestamp(1));
lastevent = max(timeSinceStart);

% extract poke events and isolate rows, also correct time format
poke_data = data(ismember(data.Event, ["Right", "Left"]), :);
poke_times = minutes(poke_data.timestamp - poke_data.timestamp(1));

% calculate poke rate using bin
binSize = options.binWindow;
edges = 0:binSize:lastevent; % form bins
poke_totals = histcounts(timeSinceStart,edges);
poke_rate = poke_totals/ binSize;

% Calculate poke rate
if sum(contains(options.method,["moving window","moving","mov"]))
    poke_rate = movsum(poke_rate, options.smoothWindow);
end

end