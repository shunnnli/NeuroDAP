function maps = getStatsMap(data,options)

% Given a vectorized data, compute the pairwise stats
arguments
    data double
    
    options.stat string = 'diff'
    options.pval logical = false
    options.nboot double = 1000
    options.reverse logical = false
end

%% Check input

if contains(options.stat,'diff',IgnoreCase=true)
    stat = 'diff';
elseif contains(options.stat,'slope',IgnoreCase=true)
    stat = 'slope';
end

% Reverse if neccessary
if options.reverse
    data = flip(data);
end

%% Generate map

map_stat = nan(length(data),length(data));
map_pval = nan(length(data),length(data));

for i = 1:length(data)
    for j = i:length(data)
        if strcmpi(stat,'slope') && (j-i>=2)
            Y = data(i:j);
            X = (1:length(i:j))';
            fit_obs = polyfit(X, Y, 1); obs = fit_obs(1);
            map_stat(i,j) = obs;
            
            if options.pval && j == length(data)
                fit_boot = bootstrp(options.nboot, @(n) polyfit(X, Y(n), 1), 1:length(Y));
                map_pval(i,j) = sum(abs(fit_boot(:,1)) >= abs(obs))/options.nboot;
            end

        elseif strcmpi(stat,'diff')
            obs = data(j) - data(i);
            map_stat(i,j) = obs;
        end
    end
end

% Reverse if neccessary
if options.reverse
    map_stat = -map_stat;
end

% Store results
maps.map = map_stat;
if options.pval; maps.pval = map_pval; end

end