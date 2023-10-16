function [] = drawTrials_optoPair(timeRange,eventIdx,...
                leftLick,rightLick,nidq,timeNI,varargin)

if ~isempty(varargin)
    leftSolenoid = find(varargin{1});
    rightSolenoid = find(varargin{2});
    airpuff = find(varargin{3});
end


% Find licking events
leftLickOnIdx = find(leftLick==1);
rightLickOnIdx = find(rightLick==1);

for i = 1:length(eventIdx)
    % Find first & last NI index
    niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidq.Fs);
    niLastIdx = eventIdx(i) + floor(timeRange(2)*nidq.Fs);

    % Find licks within timeRange
    leftLickTimesinRange = leftLickOnIdx(leftLickOnIdx>niFirstIdx & leftLickOnIdx<niLastIdx);
    rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
    relativeLeftLickTime = timeNI(leftLickTimesinRange)-timeNI(eventIdx(i));
    relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i));

    % Plot lick raster plot
    scatter(relativeLeftLickTime,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
    scatter(relativeRightLickTime,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on

    if ~isempty(varargin)
        % Find reward/punishment within timeRange
        leftSolenoidTimesinRange = leftSolenoid(leftSolenoid>niFirstIdx & leftSolenoid<niLastIdx);
        relativeLeftSolenoidTime = timeNI(leftSolenoidTimesinRange)-timeNI(eventIdx(i));
        rightSolenoidTimesinRange = rightSolenoid(rightSolenoid>niFirstIdx & rightSolenoid<niLastIdx);
        relativeRightSolenoidTime = timeNI(rightSolenoidTimesinRange)-timeNI(eventIdx(i));
        punishTimesinRange = airpuff(airpuff>niFirstIdx & airpuff<niLastIdx);
        relativePunishTime = timeNI(punishTimesinRange)-timeNI(eventIdx(i));

        if sum(size(relativePunishTime)==0) == 2; relativePunishTime = zeros(1,0); end % if matrix A is of size 0x0, make it 1x0

        scatter(relativeLeftSolenoidTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',1.5); hold on
        scatter(relativeRightSolenoidTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',1.5); hold on
        scatter(relativePunishTime,i,[],[0.7647,0.3333,0.7373],"x",'LineWidth',2); hold on
    end
end

xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Trial'); ylim([0,length(eventIdx)]);

end % drawTrials