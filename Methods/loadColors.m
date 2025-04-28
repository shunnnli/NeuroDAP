function varargout = loadColors(options)

arguments
    options.type string
end

% RGB settings: 
red = [255, 50, 58]; % FF323A

% blue = [75, 92, 204]; % #4B5CCC
blue = [85, 161, 254]; % #55A1FE
lightBlue = [7 136 225]; % 0788E1

% white = [];
lightPurple = [241 160 255]; % #F1A0FF
darkPurple = [37 12 98]; % #250C62

% Color palettes
twoColors = string(['#0788E1';'#EC4943']);
colors = string(['#4B5CCC';'#C355BC';'#FF6198';...
                 '#FF8C72';'#FFC25C';'#F9F871']);
blueRedYellow = string(['#4b5ccc'; '#6a5bcc'; '#8359ca'; '#9958c7';...
    '#ac56c3'; '#cf55b7'; '#ea57a8'; '#ff6098'; ...
    '#ff7d7b'; '#ffa564'; '#ffcf5c'; '#f9f871']);

% Color gradients
blueGreenYellow = getColormap(blue,red,500,'midCol',[27 227 61]);
blueWhiteRed = getColormap(blue,red,500,'midCol',[255 255 255]);
% blueDarkPurpleRed = getColormap(blue,red,500,'midCol',darkPurple);
bluePurpleRed = getColormap(blue,red,500,'midCol',lightPurple);
purpleWhiteRed = getColormap([213, 90, 233],[255, 50, 58],500,'midCol',[255, 255, 255]);
pinkPurpleCyan = getColormap([75, 217, 230],[225 83 161],500,'midCol',[170, 139, 252]);

varargout{1} = twoColors;
varargout{2} = colors;
varargout{3} = blueRedYellow;
varargout{4} = blueGreenYellow;
varargout{5} = blueWhiteRed;
varargout{6} = pinkPurpleCyan;
varargout{7} = bluePurpleRed;
varargout{8} = purpleWhiteRed;

end

