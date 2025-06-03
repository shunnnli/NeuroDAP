% shijia's own colors, 09/2024
function varargout = loadColorsShijia()
% RGB settings: 
black = [0 0 0]; %#000000
red = [255, 50, 58]; % FF323A
blue = [75, 92, 204]; % 4B5CCC
% lightBlue = [7 136 225]; % 0788E1
orange = [252, 137, 35]; % fc8923
green = [58, 161, 85]; % 3aa155
% lightPurple = [241 160 255]; % #F1A0FF
darkPurple = [37 12 98]; % #250C62
darkGray = [84 76 74]; % #544c4a
lightGray = [153 157 160]; % #999DA0

% Color gradients
blueWhiteRed = getColormap(blue,red,500,'midCol',[255 255 255]);
blueDarkPurpleRed = getColormap(blue,red,500,'midCol',darkPurple);
purpleWhiteRed = getColormap([213 90 233],[255, 50, 58],500,'midCol',[255 255 255]);
purpleWhiteOrange = getColormap(darkPurple,orange,500,'midCol',[255 255 255]);
purpleWhiteGreen = getColormap(darkPurple,green,500,'midCol',[255 255 255]);
blackDarkgrayLightgray = getColormap(black,darkGray,500,'midCol',lightGray);

% Color palettes
fourColors = string([ ...
    '#4B5CCC'; % blue, rewarded
    '#FF6198'; % pink, omission
    '#C967ED'; % purple, light in ACC
    '#82ED82'; % green, light in OFC
    ]); 

colors = string([ ...
    '#FF6198'; % pink, omission 
    '#fc8923'; % orange
    '#FFC25C'; % yellow
    '#3aa155'; % green
    '#12c7c1'; % cyan
    '#4B5CCC'; % blue, rewarded
    '#c541cc'; % purple
    '#999DA0'; % gray
    '#d12156']); % red

varargout{1} = fourColors;
varargout{2} = colors;
varargout{3} = blackDarkgrayLightgray;
varargout{4} = purpleWhiteOrange;
varargout{5} = blueWhiteRed;
varargout{6} = blueDarkPurpleRed;
varargout{7} = purpleWhiteGreen;
varargout{8} = purpleWhiteRed;

end

% colors = string(['#4B5CCC';'#C355BC';'#FF6198';...
%                  '#FF8C72';'#FFC25C';'#F9F871']);

% blueRedYellow = string(['#4b5ccc'; '#6a5bcc'; '#8359ca'; '#9958c7';...
%     '#ac56c3'; '#cf55b7'; '#ea57a8'; '#ff6098'; ...
%     '#ff7d7b'; '#ffa564'; '#ffcf5c'; '#f9f871']);

% this is how i get the purple and green colors in fourColors palette
    % % Define colors
    % lighterGray = [0.75, 0.75, 0.75]; % Lightened gray for no light condition
    % darkGray = [0.5, 0.5, 0.5]; % Slightly darker gray for contrast
    % 
    % % Highly contrasting purple
    % highContrastPurple = [0.7, 0.15, 0.9]; % Boost blue, reduce green for maximum distinction
    % % highContrastGreen = [0.1, 0.8, 0.2]; % Boost green, reduce red and blue for maximum distinction
    % highContrastGreen = [0.3, 0.9, 0.3]; % A vibrant and bright green (balance R, high G, and moderate B)
    % 
    % if strcmp(animalPrefix,'EF')
    %     highContrastColor = highContrastGreen;
    % else
    %     highContrastColor = highContrastPurple;
    % end
    % 
    % pastelWeight = 0.3; % Keep pastel effect for vibrancy
    % 
    % % Generate colors
    % pastelColorL = pastelWeight * [1, 1, 1] + (1 - pastelWeight) * highContrastColor; % Vibrant lighter purple
    % darkColorL = highContrastColor * 0.85; % Darker, richer purple for contrast
    % pastelColorRO = pastelWeight * [1, 1, 1] + (1 - pastelWeight) * lighterGray; % Softer light gray
    % darkColorRO = darkGray; % Slightly darker gray for differentiation
