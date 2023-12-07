function [] = drawTrials_flexibleLearning(timeRange,eventIdx,trials,...
                ~,rightLick,nidq,timeNI,varargin)

if ~isempty(varargin)
    % leftSolenoid = find(varargin{1});
    rightSolenoid = find(varargin{2});
    airpuff = find(varargin{3});
end

% Find licking events
rightLickOnIdx = find(rightLick==1);

for i = 1:length(eventIdx)
    % Find cooresponding row in trials table
    [~, row] = min(abs(trials{:,'CueTime'}-eventIdx(i)));
    TrialType = trials{row,"TrialType"};
    outcome = trials{row,"Outcome"};
    outcomeSamp = trials{row,"OutcomeTime"};
    outcomeTime = outcomeSamp/nidq.Fs;

    % Find first & last NI index
    niFirstIdx = eventIdx(i) + floor(timeRange(1)*nidq.Fs);
    niLastIdx = eventIdx(i) + floor(timeRange(2)*nidq.Fs);
    %disp(['trial: ', num2str(i)])
    %disp(['niFirst: ', num2str(niFirstIdx)])
    %disp(['niLast: ', num2str(niLastIdx)])
    

    % Find licks within timeRange
    rightLickTimesinRange = rightLickOnIdx(rightLickOnIdx>niFirstIdx & rightLickOnIdx<niLastIdx);
    relativeRightLickTime = timeNI(rightLickTimesinRange)-timeNI(eventIdx(i));

    choice_r_licks = relativeRightLickTime(relativeRightLickTime <= outcomeTime);
    if sum(size(choice_r_licks)==0) == 2; choice_r_licks = zeros(1,0); end % if matrix A is of size 0x0, make it 1x0
    if TrialType == 0
        scatter(choice_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on
    else
        scatter(choice_r_licks,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
    end

    outcome_r_licks = setdiff(relativeRightLickTime,choice_r_licks);
    if contains(outcome,'H')
        if TrialType == 0
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on
        else
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
        end
    elseif contains(outcome,'M')
        if TrialType == 0
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1',...
               'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); hold on
        else
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#FFC25C',...
               'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); hold on
        end
    elseif contains(outcome,'FA')
        if TrialType == 0
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1',...
               'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); hold on
        else
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#FFC25C',...
               'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5); hold on
        end
    elseif contains(outcome,'CR')
        if TrialType == 0
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#6DBAA1'); hold on
        else
            scatter(outcome_r_licks,i,20,'filled','MarkerFaceColor','#FFC25C'); hold on
        end
    end

    if ~isempty(varargin)
        % Find reward/punishment within timeRange
        rightSolenoidTimesinRange = rightSolenoid(rightSolenoid>niFirstIdx & rightSolenoid<niLastIdx);
        relativeRightSolenoidTime = timeNI(rightSolenoidTimesinRange)-timeNI(eventIdx(i));
        %disp(relativeRightSolenoidTime);
        %disp('----')
        punishTimesinRange = airpuff(airpuff>niFirstIdx & airpuff<niLastIdx);
        relativePunishTime = timeNI(punishTimesinRange)-timeNI(eventIdx(i));
        if sum(size(relativePunishTime)==0) == 2; relativePunishTime = zeros(1,0); end % if matrix A is of size 0x0, make it 1x0

        if trials{row,"RewardSize"} ~= 0
            scatter(relativeRightSolenoidTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',trials{row,"RewardSize"}); hold on
        end
        %scatter(outcomeTime,i,[],[0.2588,0.5294,0.9608],"+",'LineWidth',trials{row,"RewardSize"}); hold on
        
        if trials{row,"PunishSize"} ~= 0
            scatter(relativePunishTime,i,[],[0.7647,0.3333,0.7373],"x",'LineWidth',2); hold on
        end
    end
end

xlabel('Time (s)'); xlim([timeRange(1),timeRange(2)]);
ylabel('Trial'); ylim([0,length(eventIdx)]);

end % drawTrials