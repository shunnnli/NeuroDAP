function [] = drawCSUS(trialtype)
if trialtype == 0
    ylimit = ylim;
    patch([0 1 1 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
        'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
    xline(0,'-','Reward cue','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
    xline(2,'-','Reward','Color','r','LineWidth',1.5,'HandleVisibility','off');
    box off
elseif trialtype == 1
    ylimit = ylim;
    patch([0 1 1 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
        'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
    xline(0,'-','Airpuff cue','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
    xline(2,'-','Airpuff','Color','r','LineWidth',1.5,'HandleVisibility','off');
    box off
else % both reward and airpuff
    ylimit = ylim;
    patch([0 1 1 0],[ylimit(1) ylimit(1) ylimit(2) ylimit(2)],[.5 .5 .5],...
        'FaceAlpha',0.1,'EdgeColor','none','HandleVisibility','off');
    xline(0,'-','CS','Color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off');
    xline(2,'-','US','Color','r','LineWidth',1.5,'HandleVisibility','off');
    box off
end

end