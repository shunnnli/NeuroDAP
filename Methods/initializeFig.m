function fig = initializeFig(horizontalPortion,verticalPortion,options)

% Initialize figure of certain size
% Orientation: vertical or horizontal
% portion: portion of the screensize of the other direction
% E.g.: vertical, 1/2 means the horizontal width is half of the screen, the
% vertical height is the full size of the screen

arguments
    horizontalPortion double = 0.5
    verticalPortion double = 0.5
    
    options.FontSize double = 20
    options.AxesLineWidth double = 1.5
end

narginchk(0,2);
if nargin < 1
    horizontalPortion = 0.5;
    verticalPortion = 0.5;
end

screenSize = get(0,'Screensize');

if horizontalPortion > 1 || horizontalPortion <= 0 || verticalPortion > 1 || verticalPortion <= 0
    error('Portion should between 0 and 1!');
end

screenSize(3) = screenSize(3) * horizontalPortion;
screenSize(4) = screenSize(4) * verticalPortion;

fig = figure('Position',screenSize);
set(gcf,'Color','w'); box off;
% set(gca,'FontSize',options.FontSize);
% set(gca,'LineWidth',options.AxesLineWidth);


hold on
end