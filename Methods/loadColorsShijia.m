% shijia's own colors, 09/2024
function varargout = loadColorsShijia()

% RGB settings:
black = [0 0 0]; %#000000
red = [255, 50, 58]; % FF323A
blue = [75, 92, 204]; % 4B5CCC
orange = [252, 137, 35]; % fc8923
green = [58, 161, 85]; % 3aa155
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

% Color palettes (HEX)
mainColors = string([ ...
    '#4B5CCC'; % (1) blue, rewarded
    '#FF6198'; % (2) pink, omission
    '#C967ED'; % (3) purple, light in ACC
    '#82ED82'; % (4) green, light in OFC
    '#CC1A99'; % (5) red-purple, light in PL
    '#008B8B'; % (6) teal, light in SC and InsCtx
    '#8C564B'; % (7) brown, light in IRt
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

% ============================================================
% NEW: block-task 4-color set (small/big x rew/omit)
%   smallR  = lightBLUE
%   smallO  = lightPINK
%   bigR    = darkBLUE
%   bigO    = darkPINK
% ============================================================
amount = 0.2;  % same as your snippet

% mainColors(1) = blue, mainColors(2) = pink
[lightBLUE, darkBLUE] = adjustHexBrightness(mainColors(1), amount); % returns RGB in [0,1]
[lightPINK, darkPINK] = adjustHexBrightness(mainColors(2), amount);

% convenient 4x3 RGB stack in YOUR requested order:
blockColors = [
    lightBLUE;   % smallR
    lightPINK;   % smallO
    darkBLUE;    % bigR
    darkPINK     % bigO
];


% ---- outputs (keep existing numbering; add new ones at the end) ----
varargout{1}  = mainColors;
varargout{2}  = colors;
varargout{3}  = blackDarkgrayLightgray;
varargout{4}  = purpleWhiteOrange;
varargout{5}  = blueWhiteRed;
varargout{6}  = blueDarkPurpleRed;
varargout{7}  = purpleWhiteGreen;
varargout{8}  = purpleWhiteRed;

% NEW outputs:
varargout{9}  = blockColors;   % 4x3 RGB, order: [smallR; smallO; bigR; bigO]
end