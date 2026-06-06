function drawLicks(timeRange,eventIdx,leftLick,rightLick,nidqFs,timeNI,varargin)

if ~isempty(varargin)
    leftSolenoid = find(varargin{1});
    rightSolenoid = find(varargin{2});
    airpuff = find(varargin{3});
end

if isempty(rightLick)
    % Find licking events
    leftLickOnIdx = find(leftLick==1);

    for i = 1:length(eventIdx)
        % Find first & last NI index
        niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidqFs);
        niLastIdx = eventIdx(i) + floor(timeRange(2)*nidqFs);

        % Find licks within timeRange
        leftLickTimesinRange = leftLickOnIdx(leftLickOnIdx>niFirstIdx & leftLickOnIdx<niLastIdx);
        relativeLeftLickTime = timeNI(leftLickTimesinRange)-timeNI(eventIdx(i));
        % Plot lick raster plot
        scatter(relativeLeftLickTime,i,'filled','MarkerFaceColor','#6DBAA1'); hold on
        % drawnow
        
        if ~isempty(varargin)
            % Find reward/punishment within timeRange
            leftSolenoidTimesinRange = leftSolenoid(leftSolenoid>niFirstIdx & leftSolenoid<niLastIdx);
            relativeLeftSolenoidTime = timeNI(leftSolenoidTimesinRange)-timeNI(eventIdx(i));
            punishTimesinRange = airpuff(airpuff>niFirstIdx & airpuff<niLastIdx);
            relativePunishTime = timeNI(punishTimesinRange)-timeNI(eventIdx(i));

            scatter(relativeLeftSolenoidTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',1.5); hold on
            scatter(relativePunishTime,i,[],[0.7647,0.3333,0.7373],"x",'LineWidth',2); hold on
        end
    end

    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
    ylabel('Trial'); ylim([0,length(eventIdx)]);
    %xline(0,'-',eventName,'Color','r','LineWidth',2,'HandleVisibility','off');
    %box off
    
elseif isempty(leftLick)
    % Find licking events
    rightLickOnIdx = find(rightLick==1);
    
    for i = 1:length(eventIdx)
        % Find first & last NI index
        niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidqFs);
        niLastIdx = eventIdx(i) + floor(timeRange(2)*nidqFs);

        % Find licks within timeRange
        rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
        relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i));
        % Plot lick raster plot
        scatter(relativeRightLickTime,i,'filled','MarkerFaceColor','#6DBAA1'); hold on
        
        if ~isempty(varargin)
            % Find reward/punishment within timeRange
            rightSolenoidTimesinRange = rightSolenoid(rightSolenoid>niFirstIdx & rightSolenoid<niLastIdx);
            relativeRightSolenoidTime = timeNI(rightSolenoidTimesinRange)-timeNI(eventIdx(i));
            punishTimesinRange = airpuff(airpuff>niFirstIdx & airpuff<niLastIdx);
            relativePunishTime = timeNI(punishTimesinRange)-timeNI(eventIdx(i));
            %disp(relativeRightSolenoidTime);

            scatter(relativeRightSolenoidTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',1.5); hold on
            scatter(relativePunishTime,i,[],[0.7647,0.3333,0.7373],"x",'LineWidth',2); hold on
        end
    end

    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
    ylabel('Trial'); ylim([0,length(eventIdx)]);
else
    % Find licking events
    leftLickOnIdx = find(leftLick==1);
    rightLickOnIdx = find(rightLick==1);

    for i = 1:length(eventIdx)
        % Find first & last NI index
        niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidqFs);
        niLastIdx = eventIdx(i) + floor(timeRange(2)*nidqFs);

        % Find licks within timeRange
        leftLickTimesinRange = leftLickOnIdx(leftLickOnIdx>niFirstIdx & leftLickOnIdx<niLastIdx);
        rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
        relativeLeftLickTime = timeNI(leftLickTimesinRange)-timeNI(eventIdx(i));
        relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i));
        % Plot lick raster plot
        scatter(relativeLeftLickTime,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
        scatter(relativeRightLickTime,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on
        % drawnow
        
        if ~isempty(varargin)
            disp(i)
            % Find reward/punishment within timeRange
            leftSolenoidTimesinRange = leftSolenoid(leftSolenoid>niFirstIdx & leftSolenoid<niLastIdx);
            relativeLeftSolenoidTime = timeNI(leftSolenoidTimesinRange)-timeNI(eventIdx(i));
            rightSolenoidTimesinRange = rightSolenoid(rightSolenoid>niFirstIdx & rightSolenoid<niLastIdx);
            relativeRightSolenoidTime = timeNI(rightSolenoidTimesinRange)-timeNI(eventIdx(i));
            punishTimesinRange = airpuff(airpuff>niFirstIdx & airpuff<niLastIdx);
            relativePunishTime = timeNI(punishTimesinRange)-timeNI(eventIdx(i));

            scatter(relativeLeftSolenoidTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',1.5); hold on
            scatter(relativeRightSolenoidTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',1.5); hold on
            scatter(relativePunishTime,i,[],[0.7647,0.3333,0.7373],"x",'LineWidth',2); hold on
        end
    end

    xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
    ylabel('Trial'); ylim([0,length(eventIdx)]);
end



end % drawLicks