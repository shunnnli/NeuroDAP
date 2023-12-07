function yFilt = fftFilter(y,Fs,CutOff,FilterType)
% 
% INPUT
%   y: data to be filtered (matrix M samples x N channels)
%   Fs: Sampling rate
%   CutOff: Cutoff frequency
%   FilterType: FilterType: low-pass = 1, high-pass = 2. 
% OUTPUT
%   yFilt: Filtered data

%   Created: SRO 2/19/10
%   Modified: SRO 6/16/10 (fixed edge artifacts)


% Do NOT pad beginning and end vector to avoid edge artifacts

% Get length of original data
ly = size(y,1);

% Pad with nPad points
nPad = ceil(ly/1000);

% Temp matrix for storing padded data
temp = zeros(ly+2*nPad,size(y,2))+1;

% Start pad
p1 = y(1:nPad,:);

% End pad
p2 = y(end-nPad+1:end,:);

% Flip up to down
p1 = flipud(p1);
p2 = flipud(p2);

% Fill in temp matrix
temp(1:nPad,:) = p1;
temp(nPad+1+ly:end,:) = p2;
temp(nPad+1:nPad+ly,:) = y;

y = temp;
n = length(y);                                  % number of points
t = n/Fs;                                       % duration in seconds
freqPos = (0:n/2-1)/t;                          % frequencies
freqNeg = sort(freqPos(2:end),'descend');
freqAll = [freqPos n/2/t freqNeg];              % vector with all frequencies
CutOffInd = (freqAll < CutOff)';                % get indices of frequencies > 500 Hz
CutOffInd = repmat(CutOffInd,1,size(y,2));      % replicate vector to match size of matrix y
Fcoef = fft(y);                                 % FFT

switch FilterType
    case 1    % low-pass; set fourier coefficients above CutOff to zero
        Fcoef(~CutOffInd) = 0;
    case 2    % high-pass; set fourier coeffients below CutOff to zero
        Fcoef(CutOffInd) = 0;
end

yFilt = ifft(Fcoef);                            % inverse FFT

% Get unpadded data
yFilt = real(yFilt(nPad+1:nPad+ly,:));

