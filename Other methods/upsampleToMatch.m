function [scalar,upsampled] = upsampleToMatch(lowFs,highFs,lowData)

% Upsample to match sampling rate
scalar = round(highFs/lowFs);
upsampled = upsample(lowData,scalar);
% Fill upsampled steps with 1 if during sync pulse
for i = 1:length(upsampled)-scalar
    if upsampled(i) == 1 && upsampled(i+scalar) == 1
        for j = 1:scalar-1
            upsampled(i+j) = 1;
        end
    end
end

end