function plotCI(x,y,color)

% N = nnz(~isnan(y(:,1)));   % Number of ‘Experiments’ In Data Set
yMean = mean(y,1,"omitnan");      % Mean Of All Experiments At Each Value Of ‘x’
yCI95 = zeros();
for i=1:size(y,2)
    yCI95(i) = getCI(y(:,i),0.95);
end

finalCI = yCI95;

plot(x,yMean,'Color',color,'LineWidth',2); hold on
if size(y,1) >= 5
    patch('XData',[x, fliplr(x)], 'YData',[yMean-finalCI, fliplr(yMean+finalCI)], ...
        'FaceColor',color,'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
end
hold on
box off

end % plotCI