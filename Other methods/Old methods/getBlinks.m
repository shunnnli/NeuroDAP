function blink_occurance = getBlinks(eyeTraces,blink_thresh)

blink_occurance = zeros(size(eyeTraces,1),size(eyeTraces,2));
for i = 1:size(eyeTraces,1)
    blink_occurance(i,:) = eyeTraces(i,:) > blink_thresh;
end

end