function [] = drawTrials_riskTaking(timeRange,eventIdx,trials,...
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
    % Find cooresponding row in trials table
    [~, row] = min(abs(trials{:,'CueTime'}-eventIdx(i)));
    outcome = trials{row,"Outcome"};
    outcomeSamp = trials{row,"OutcomeTime"};
    outcomeTime = outcomeSamp/nidq.Fs;

    % Find first & last NI index
    niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidq.Fs);
    niLastIdx = eventIdx(i) + floor(timeRange(2)*nidq.Fs);

    % Find licks within timeRange
    leftLickTimesinRange = leftLickOnIdx(leftLickOnIdx>niFirstIdx & leftLickOnIdx<niLastIdx);
    rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
    relativeLeftLickTime = timeNI(leftLickTimesinRange)-timeNI(eventIdx(i));
    relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i));

    choice_l_licks = relativeLeftLickTime(relativeLeftLickTime <= outcomeTime);
    choice_r_licks = relativeRightLickTime(relativeRightLickTime <= outcomeTime);
    if sum(size(choice_l_licks)==0) == 2; choice_l_licks = zeros(1,0); end % if matrix A is of size 0x0, make it 1x0
    if sum(size(choice_r_licks)==0) == 2; choice_r_licks = zeros(1,0); end % if matrix A is of size 0x0, make it 1x0
    scatter(choice_l_licks,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
    scatter(choice_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on

    outcome_l_licks = setdiff(relativeLeftLickTime,choice_l_licks);
    outcome_r_licks = setdiff(relativeRightLickTime,choice_r_licks);
    if contains(outcome,'F')
        if strcmp(outcome,'LF')
            line([outcomeTime timeRange(2)],[row row],"LineStyle","--",'Color','#FFC25C','LineWidth',1); hold on
        elseif strcmp(outcome,'RF')
            line([outcomeTime timeRange(2)],[row row],"LineStyle","--",'Color','#6DBAA1','LineWidth',1); hold on
        end
        % Plot lick raster plot
        scatter(outcome_l_licks,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
        scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on
    elseif contains(outcome,'R') && ~contains(outcome,'F')
        % Plot lick raster plot
        scatter(outcome_l_licks,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
        scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on
    else
        % Plot lick raster plot
        scatter(outcome_l_licks,i,20,'filled','MarkerFaceColor','#FFC25C',...
           'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); hold on
        scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1',...
           'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); hold on
        
        if strcmp(outcome,'NO')
            scatter(outcomeSamp/nidq.Fs,i,[],[0.0784,0.2471,0.3216],"o",'LineWidth',1.5); hold on
        end
    end

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