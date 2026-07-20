function index = getClampIdx(blue,red,options)

    arguments
        blue double
        red double
        options.type string = 'diff'

        options.removeBaseline logical = true
        options.timestamp double
        options.baselineWindow double % in samples
    end
    
    if options.removeBaseline
        if ~isfield(options,'timestamp') && ~isfield(options,'baselineWindow')
            error('Provide either timestamp or baselineWindow in format of [1 25]!!');
        elseif isfield(options,'timestamp') && ~isfield(options,'baselineWindow')
            options.baselineWindow = find(options.timestamp<0 & options.timestamp >= -0.2);
        end

        % Remove baseline for each row
        final_red = red - median(red(:,options.baselineWindow),2);
        final_blue = blue - median(blue(:,options.baselineWindow),2);
    else
        final_red = red;
        final_blue = blue;
    end

    if strcmpi(options.type,'diff')
        index = final_red - final_blue;
    elseif strcmpi(options.type,'norm')
        index = (final_red - final_blue) ./ (final_red + final_blue);
    end

end