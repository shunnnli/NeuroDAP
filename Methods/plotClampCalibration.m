function fig = plotClampCalibration(calibration)
%PLOTCLAMPCALIBRATION Mean +/- SEM by power and late-on dose-response curves.
fig = initializeFig(0.67,0.8);
tiledlayout(2,2,'TileSpacing','compact');
channels = ["red","blue"];
labels = ["Red / excitation","Blue / inhibition"];
baseColors = [0.85,0.2,0.15;0.15,0.4,0.85];
for c = 1:2
    selected = arrayfun(@(group) group.channel == channels(c),calibration.groups);
    groups = calibration.groups(selected);
    traceAx = nexttile(c);
    hold(traceAx,'on');
    responseAx = nexttile(c+2);
    hold(responseAx,'on');
    set([traceAx,responseAx],'FontSize',14,'LineWidth',1.2);
    nLevels = numel(groups);
    for g = 1:nLevels
        group = groups(g);
        color = 1-(0.4+0.6*g/max(1,nLevels))*(1-baseColors(c,:));
        if group.nTrials == 0; continue; end
        t = group.time_sec;
        avg = group.mean_dff;
        sem = group.sem_dff;
        fill(traceAx,[t fliplr(t)],[avg-sem fliplr(avg+sem)],color, ...
            'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
        plot(traceAx,t,avg,'Color',color,'LineWidth',2, ...
            'DisplayName',sprintf('%g%% (n=%d)',group.power_pct,group.nTrials));
        duration = median(group.duration_sec);
        xline(traceAx,duration,':','Color',color,'HandleVisibility','off');
        trials = group.lateOnMedian_dff;
        scatter(responseAx,repmat(group.power_pct,size(trials)),trials,24,color, ...
            'filled','MarkerFaceAlpha',0.3,'HandleVisibility','off');
    end
    xline(traceAx,0,'k--','HandleVisibility','off');
    yline(traceAx,0,':','HandleVisibility','off');
    title(traceAx,labels(c));
    xlabel(traceAx,'Time from laser onset (s)');
    ylabel(traceAx,'Photometry (\DeltaF/F_0)');
    if any([groups.nTrials] > 0)
        legend(traceAx,'Location','best');
    else
        text(traceAx,0.5,0.5,'No valid calibration trials','Units','normalized', ...
            'HorizontalAlignment','center');
    end
    if ~isempty(groups)
        powers = [groups.power_pct];
        responses = [groups.medianResponse_dff];
        plot(responseAx,powers,responses,'o-','Color',baseColors(c,:), ...
            'LineWidth',2,'DisplayName','Median across trials');
        fit = calibration.fits(c);
        finite = isfinite(responses);
        if isfinite(fit.slope)
            x = [min(powers(finite)),max(powers(finite))];
            plot(responseAx,x,fit.slope*x+fit.intercept,'--','Color',baseColors(c,:), ...
                'DisplayName',sprintf('Fit: slope %.3g, R^2 %.3f',fit.slope,fit.r2));
        end
        if any(finite); legend(responseAx,'Location','best'); end
    end
    title(responseAx,labels(c)+": power response");
    xlabel(responseAx,'Laser command (% of DAC range)');
    ylabel(responseAx,'Late-on median (\DeltaF/F_0)');
    grid(responseAx,'on');
    box(traceAx,'off'); box(responseAx,'off');
end
label = 'Calibration';
if isfield(calibration,'name'); label = sprintf('Calibration: %s (LabJack)',calibration.name); end
sgtitle(sprintf('%s: mean +/- SEM; F0 from %.3g s before onset', ...
    label,calibration.options.preTime),'Interpreter','none');
end
