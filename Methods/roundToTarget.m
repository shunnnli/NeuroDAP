function rounded = roundToTarget(list,target)
% Round individual element in list to nearest element in target

rounded = zeros(size(list));
for i = 1:length(list)
    [~, idx] = min(abs(list(i)*ones(size(target)) - target));
    rounded(i) = target(idx);
end

end