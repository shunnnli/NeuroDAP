function shuffled = shuffleTraces(traces,options)

arguments
    traces double
    options.nShuffle double = nan
end

if isnan(options.nShuffle); options.nShuffle = size(traces,1); end
shuffled = zeros(options.nShuffle,size(traces,2));
% Shuffle data
for i = 1:size(shuffled,1)
    shuffled(i,:) = traces(randi(size(traces,1)),randperm(size(traces,2)));
end

end