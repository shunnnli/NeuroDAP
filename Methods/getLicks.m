function [lickRate,lickTraces,lickEvents] = getLicks(timeRange,eventIdx,binSize,...
                leftLick,rightLick,Fs,timeNI,options)

% Get licks around event

arguments
    timeRange double
    eventIdx double
    binSize double
    leftLick 
    rightLick 
    Fs double
    timeNI double
    options.side = [0,1]
    options.inputLickIdx logical = false % If true, then don't need to find the lick once again
    options.getRate logical = true
    options.behaviorFs double = 10000 % sample freq for behavior (NI)
end

% Initialize lickTraces (number of licks in bin) & lickRate 
lickTraces = cell(1,2);
lickTraces_left = zeros(length(eventIdx),round((timeRange(2)-timeRange(1))/binSize));
lickTraces_right = zeros(length(eventIdx),round((timeRange(2)-timeRange(1))/binSize));
lickRate = lickTraces;
% Initialize lickEvents (relative time of every lick for every trial)
% Left side: first column, right side: second column
lickEvents = cell(length(eventIdx),2);

% Return if there's no licks
if isempty(find(leftLick,1)) & isempty(find(rightLick, 1)); return; end

% Find licks
if isempty(rightLick)
    % Find licking events
    if options.inputLickIdx; leftLickOnIdx = leftLick;
    else; leftLickOnIdx = find(leftLick==1); 
    end

    for i = 1:length(eventIdx)
        % Find first & last NI index
        niFirstIdx = eventIdx(i) + floor(timeRange(1)*Fs);
        niLastIdx = eventIdx(i) + floor(timeRange(2)*Fs);

        % Find licks within timeRange
        leftLickTimesinRange = leftLickOnIdx(leftLickOnIdx>niFirstIdx & leftLickOnIdx<niLastIdx);
        relativeLeftLickTime = timeNI(leftLickTimesinRange)-timeNI(eventIdx(i));

        if options.getRate
            % Calculate lick rate
            relativeLeftLickBin = floor((relativeLeftLickTime - timeRange(1))/binSize)+1;
            lickTraces_left(i,relativeLeftLickBin) = 1 + lickTraces_left(i,relativeLeftLickBin);
        end
        
        % lickEvents
        lickEvents{i,1} = relativeLeftLickTime;

    end

elseif isempty(leftLick)
    % Find licking events
    if options.inputLickIdx; rightLickOnIdx = rightLick;
    else; rightLickOnIdx = find(rightLick==1); 
    end

    for i = 1:length(eventIdx)
        % Find first & last NI index
        niFirstIdx = eventIdx(i) + floor(timeRange(1)*Fs);
        niLastIdx = eventIdx(i) + floor(timeRange(2)*Fs);

        % Find licks within timeRange
        rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
        relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i)); % Can through out error is timeNI(0)

        if options.getRate
            % Calculate lick rate
            relativeRightLickBin = floor((relativeRightLickTime - timeRange(1))/binSize)+1;
            lickTraces_right(i,relativeRightLickBin) = 1 + lickTraces_right(i,relativeRightLickBin);
        end
        
        % lickEvents
        lickEvents{i,2} = relativeRightLickTime;
    end

else
    % Find licking events
    if options.inputLickIdx
        rightLickOnIdx = rightLick;
        leftLickOnIdx = leftLick;
    else
        rightLickOnIdx = find(rightLick==1);
        leftLickOnIdx = find(leftLick==1);
    end

    for i = 1:length(eventIdx)
        % Find first & last NI index
        niFirstIdx = eventIdx(i) + floor(timeRange(1)*Fs);
        niLastIdx = eventIdx(i) + floor(timeRange(2)*Fs);

        % Find licks within timeRange
        leftLickTimesinRange = leftLickOnIdx(leftLickOnIdx>niFirstIdx & leftLickOnIdx<niLastIdx);
        rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
        relativeLeftLickTime = timeNI(leftLickTimesinRange)-timeNI(eventIdx(i));
        relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i));

        if options.getRate
            % Calculate lick rate
            relativeLeftLickBin = floor((relativeLeftLickTime - timeRange(1))/binSize)+1;
            lickTraces_left(i,relativeLeftLickBin) = 1 + lickTraces_left(i,relativeLeftLickBin);
            relativeRightLickBin = floor((relativeRightLickTime - timeRange(1))/binSize)+1;
            lickTraces_right(i,relativeRightLickBin) = 1 + lickTraces_right(i,relativeRightLickBin);
        end
        
        % lickEvents
        lickEvents{i,1} = relativeLeftLickTime;
        lickEvents{i,2} = relativeRightLickTime;
    end
end

if all(options.side == [1 1])
    lickTraces{1} = lickTraces_left;
    lickTraces{2} = lickTraces_right;
    lickRate{1} = lickTraces{1} / binSize;
    lickRate{2} = lickTraces{2} / binSize;
elseif all(options.side == [1 0])
    lickTraces = lickTraces_left;
    lickRate = lickTraces / binSize;
    lickEvents = lickEvents(:,1);
elseif all(options.side == [0 1])
    lickTraces = lickTraces_right;
    lickRate = lickTraces / binSize;
    lickEvents = lickEvents(:,2);
else 
    warning('side format is wrong. Should be a 1x2 matrix (e.g. [0 1])');
end

end % getLicks