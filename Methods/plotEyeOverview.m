function plotEyeOverview(clampIdx,unclampIdx,trials,eyeArea_detrend,pupilArea_detrend,params,options)

arguments
    clampIdx double
    unclampIdx double
    trials table
    eyeArea_detrend double
    pupilArea_detrend double
    params struct

    options.cameraTimeRange double = [3,5]
    options.groupSize double = 20
    options.smooth double = 15
    options.eventDuration double = 0.25
    options.sessionpath string = ''
    options.save logical = true
    options.figureName string = 'Behavior_EyeOverview'
end

[~,~,~,~,~,~,bluePurpleRed] = loadColors;
cameraTimeRange = options.cameraTimeRange;

initializeFig(0.67,0.67); tiledlayout(2,4);

%% Eye trace heatmap
nexttile([2 2]);
[allEyeTraces,t_cam] = plotTraces(trials{1:end-1,"CueTime"},cameraTimeRange,...
    eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
plotHeatmap(allEyeTraces,t_cam,centerColorMap=false);
set(gca,'YDir','normal');
colorbar; box off
plotEvent('Cue',options.eventDuration);
xlabel('Time (s)'); ylabel('Trials');

%% Extract traces by trial type
[clampEyeArea,t_cam] = plotTraces(clampIdx,cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
[unclampEyeArea,~] = plotTraces(unclampIdx,cameraTimeRange,eyeArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
[clampPupilArea,~] = plotTraces(clampIdx,cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');
[unclampPupilArea,~] = plotTraces(unclampIdx,cameraTimeRange,pupilArea_detrend,[0,0,0],params,plot=false,signalSystem='camera');

%% Eye and pupil area across session
signalPairs = {
    {clampEyeArea, unclampEyeArea}, 'Eye area (a.u.)', [3, 7];
    {clampPupilArea, unclampPupilArea}, 'Pupil area (a.u.)', [4, 8];
};
labels = {'Clamp','Unclamp'};

for sp = 1:size(signalPairs,1)
    traces = signalPairs{sp,1};
    yLabel = signalPairs{sp,2};
    tilePositions = signalPairs{sp,3};

    for event = 1:length(traces)
        nexttile(tilePositions(event));
        legendList = plotGroupTraces(traces{event},t_cam,bluePurpleRed,...
            groupSize=options.groupSize);
        plotEvent(labels{event},options.eventDuration);
        xlabel('Time (s)'); ylabel(yLabel);
        legend(legendList);
    end
end

%% Save
if options.save && strlength(options.sessionpath) > 0
    saveFigures(gcf,options.figureName,options.sessionpath,savePDF=true);
end

end